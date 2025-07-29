import os
from flask import Flask, jsonify, request
import firebase_admin
from firebase_admin import credentials, firestore, storage
from flask_cors import CORS
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
import re
import uuid
from google.oauth2 import service_account
from googleapiclient.discovery import build
from googleapiclient.http import MediaFileUpload
import tempfile
import base64
import datetime
from collections import Counter
import tempfile
import sys
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "*"}})

# 🔍 DEBUG: Logging environment variable presence
print("🔥 Checking for FIREBASE_KEY in environment variables...")
if "FIREBASE_KEY" not in os.environ:
    print("❌ FIREBASE_KEY is missing from environment", file=sys.stderr)
    sys.exit(3)

try:
    firebase_key_json = os.environ["FIREBASE_KEY"]
    print("✅ FIREBASE_KEY loaded from environment")

    # Write to temporary file
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".json") as temp_key_file:
        temp_key_file.write(firebase_key_json)
        firebase_key_path = temp_key_file.name
    print(f"📁 Firebase key written to temp file at: {firebase_key_path}")

    # Initialize Firebase
    cred = credentials.Certificate(firebase_key_path)
    firebase_admin.initialize_app(cred, {
        'storageBucket': 'polaris-test-3b8f8.firebasestorage.app'
    })
    db = firestore.client()
    print("✅ Firebase Admin initialized successfully")

except Exception as e:
    print("❌ Failed to initialize Firebase Admin SDK:", str(e), file=sys.stderr)
    sys.exit(3)


# Only treat FN - FNG and FNG - FN as equivalent, not all transitions
FN_FNG_EQUIV = {('FN', 'FNG'), ('FNG', 'FN')}

def canonicalize_transition(transition: str) -> str:
    """
    HELPER FUNCTION: Standardizes phase transition names
    Takes a phase transition like "FN - FNG" and makes sure it's written consistently
    This helps avoid duplicates when searching for similar compounds
    """
    parts = [p.strip().upper() for p in transition.split('-')]
    if len(parts) == 2:
        pair = (parts[0], parts[1])
        if pair in FN_FNG_EQUIV:
            return 'FN - FNG'
        return f'{parts[0]} - {parts[1]}'
    return transition.strip().upper()

def extract_transitions_with_temps(phase_map_str: str):
    """
    HELPER FUNCTION: Extracts phase transitions and temperatures from text
    Takes a string like "CR-120-N;N-140-I" and converts it to a list
    of standardized transitions with temperatures
    """
    if not isinstance(phase_map_str, str):
        return []

    transitions = []
    segments = phase_map_str.split(";")

    valid_phases = {
        "CR", "N", "I", "FN", "NF", "SMA", "SMX", "FNG", "NX", "GLASS", "DECOMP", "SMZA"
    }

    for segment in segments:
        # Remove any parenthetical notes (e.g., (no crystals observed))
        segment = re.sub(r"\([^)]*\)", "", segment)
        parts = re.split(r'\s*-\s*', segment.strip())
        parsed = []

        for part in parts:
            clean = part.strip().upper().rstrip("?")
            if clean in valid_phases:
                parsed.append(("phase", clean))
            elif re.match(r"[~>]?[-+]?[0-9]+(?:\.[0-9]+)?", clean):
                temp_match = re.search(r"[-+]?[0-9]+(?:\.[0-9]+)?", clean)
                if temp_match:
                    parsed.append(("temp", temp_match.group()))

        i = 0
        while i + 2 < len(parsed):
            if parsed[i][0] == "phase" and parsed[i+1][0] == "temp" and parsed[i+2][0] == "phase":
                phase1 = parsed[i][1]
                temp = parsed[i+1][1]
                phase2 = parsed[i+2][1]
                if phase1 != phase2:
                    canonical = canonicalize_transition(f"{phase1} - {phase2}")
                    transitions.append(f"{canonical} @ {temp}")
                i += 2
            else:
                i += 1

    return transitions


@app.route('/compounds')
def get_compounds():
    """
    API ENDPOINT: Get all compounds from the database
    Returns a list of all compounds stored in Firestore
    Used by the frontend to display the compounds table
    """
    compounds_ref = db.collection('compounds')
    docs = compounds_ref.stream()
    compounds = [doc.to_dict() for doc in docs]
    return jsonify(compounds)


@app.route("/compounds/<compound_id>")
def get_single_compound(compound_id):
    """
    API ENDPOINT: Get details for one specific compound
    Takes a compound ID and returns all the data for that compound
    Used when viewing compound details or editing a compound
    """
    try:
        doc = db.collection("compounds").document(compound_id).get()
        if not doc.exists:
            return jsonify({"error": "Compound not found"}), 404
        return jsonify(doc.to_dict()), 200
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/search-substructure", methods=["POST"])
def search_substructure():
    """
    API ENDPOINT: Find compounds that contain a specific molecular structure
    Takes a SMILES string (chemical structure code) and finds all compounds
    that contain that structure as a part of their molecule
    """
    data = request.get_json()
    print("Received JSON:", data)
    query_smiles = data.get("query")

    if not query_smiles:
        return jsonify({"error": "No SMILES provided"}), 400

    try:
        query_mol = Chem.MolFromSmiles(query_smiles)
        if query_mol is None:
            return jsonify([])

        matched_compounds = []
        compounds_ref = db.collection("compounds")
        docs = compounds_ref.stream()

        for doc in docs:
            compound = doc.to_dict()
            smiles = compound.get("smiles", "")
            mol = Chem.MolFromSmiles(smiles)
            if mol and mol.HasSubstructMatch(query_mol):
                matched_compounds.append({ "id": doc.id, **compound })

        return jsonify(matched_compounds)

    except Exception as e:
        return jsonify({"error": str(e)}), 500
    
@app.route("/search-substructure-formulations", methods=["POST"])
def search_substructure_formulations():
    """
    API ENDPOINT: Find formulations that contain compounds with a specific structure
    First finds compounds with the given structure, then finds all formulations
    that use those compounds as ingredients
    """
    data = request.get_json()
    query_smiles = data.get("query")

    if not query_smiles:
        return jsonify({"error": "No SMILES provided"}), 400

    try:
        query_mol = Chem.MolFromSmiles(query_smiles)
        if query_mol is None:
            return jsonify([])

        # Step 1: Find matching compounds
        matched_ids = set()
        compounds_ref = db.collection("compounds")
        for doc in compounds_ref.stream():
            compound = doc.to_dict()
            smiles = compound.get("smiles", "")
            mol = Chem.MolFromSmiles(smiles)
            if mol and mol.HasSubstructMatch(query_mol):
                matched_ids.add(compound.get("id"))

        if not matched_ids:
            return jsonify([])

        # Step 2: Filter formulations by matching compoundId in components
        matching_formulations = []
        formulations_ref = db.collection("formulations")
        for doc in formulations_ref.stream():
            formulation = doc.to_dict()
            components = formulation.get("components", [])
            for comp in components:
                if comp.get("compoundId") in matched_ids:
                    matching_formulations.append({ "id": doc.id, **formulation })
                    break  # One match is enough

        return jsonify(matching_formulations)

    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/add-compound", methods=["POST", "OPTIONS"])
def add_compound():
    """
    API ENDPOINT: Add a new compound to the database
    Takes compound data from the frontend, generates molecular images,
    calculates molecular weight, and saves everything to Firestore
    """
    print("🧪 Received method:", request.method)
    if request.method == "OPTIONS":
        return '', 200

    data = request.get_json()

    if "original_id" in data:
        return jsonify({"error": "Cannot add a compound that originates from a lot."}), 400

    if "id" not in data:
        return jsonify({"error": "Compound must have an 'id' field"}), 400

    data["createdAt"] = firestore.SERVER_TIMESTAMP

    phase_map_str = data.get("phase map", "")
    transitions = extract_transitions_with_temps(phase_map_str)
    data["parsed_phase_transitions"] = transitions

    smiles = data.get("smiles")
    if smiles and not data.get("imageUrl"):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                data["MW"] = round(Descriptors.MolWt(mol), 2)  # Optional: round to 2 decimal places
            if mol:
                rdDepictor.Compute2DCoords(mol)

                # Use even larger canvas and adjust draw options for bigger molecule and slightly smaller font
                width, height = 2400, 1200  # Keep large canvas
                drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
                options = drawer.drawOptions()
                num_atoms = mol.GetNumAtoms()
                # Dynamically set both minFontSize and maxFontSize for atom labels
                if num_atoms <= 20:
                    options.minFontSize = 120
                    options.maxFontSize = 120
                elif num_atoms <= 40:
                    options.minFontSize = 90
                    options.maxFontSize = 90
                else:
                    options.minFontSize = 58
                    options.maxFontSize = 58
                options.bondLineWidth = 4
                options.padding = 0.01    # Minimal padding to fill image
                options.fixedScale = True # Force molecule to fill canvas

                # options.useBWAtomPalette()

                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()

                # Save image using WriteDrawingText (RDKit >=2022)
                with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
                    drawer.WriteDrawingText(tmp.name)
                    tmp.flush()

                    bucket = storage.bucket()
                    filename = f"compound_images/{data['id']}.png"
                    blob = bucket.blob(filename)
                    blob.upload_from_filename(tmp.name, content_type="image/png")
                    blob.make_public()
                    data["imageUrl"] = blob.public_url

                os.unlink(tmp.name)
                # Generate true black and white image for printing
                drawer_bw = rdMolDraw2D.MolDraw2DCairo(width, height)
                options_bw = drawer_bw.drawOptions()
                options_bw.minFontSize = options.minFontSize
                options_bw.maxFontSize = options.maxFontSize
                options_bw.bondLineWidth = options.bondLineWidth
                options_bw.padding = options.padding
                options_bw.fixedScale = options.fixedScale
                options_bw.useBWAtomPalette()  # Force black atom labels and bonds
                drawer_bw.DrawMolecule(mol)
                drawer_bw.FinishDrawing()
                with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp_bw:
                    drawer_bw.WriteDrawingText(tmp_bw.name)
                    tmp_bw.flush()
                    filename_bw = f"compound_images/{data['id']}_bw.png"
                    blob_bw = bucket.blob(filename_bw)
                    blob_bw.upload_from_filename(tmp_bw.name, content_type="image/png")
                    blob_bw.make_public()
                    data["bwImageUrl"] = blob_bw.public_url
                os.unlink(tmp_bw.name)

        except Exception as e:
            print(f"Failed to generate/upload image for {data['id']}: {e}")

    db.collection("compounds").document(data["id"]).set(data, merge=True)
    return jsonify({"success": True}), 200


@app.route("/compute-mw", methods=["POST"])
def compute_mw():
    """
    API ENDPOINT: Calculate molecular weight from SMILES
    Takes a SMILES string (chemical structure code) and returns
    the calculated molecular weight of that compound
    """
    data = request.get_json()
    smiles = data.get("smiles", "")
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mw = round(Descriptors.MolWt(mol), 2)
            return jsonify({"MW": mw})
        else:
            return jsonify({"error": "Invalid SMILES"}), 400
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/update-compound", methods=["POST"])
def update_compound():
    """
    API ENDPOINT: Update an existing compound's information
    Takes compound data with changes and saves the updated
    information back to the database
    """
    data = request.get_json()
    if "id" not in data:
        return jsonify({"error": "Compound must have an 'id' field"}), 400

    compound_id = data["id"]
    try:
        db.collection("compounds").document(compound_id).set(data, merge=True)
        return jsonify({"success": True}), 200
    except Exception as e:
        return jsonify({"error": str(e)}), 500
    

@app.route("/delete-compound/<string:compound_id>", methods=["DELETE"])
def delete_compound(compound_id):
    """
    API ENDPOINT: Delete a compound from the database
    Permanently removes a compound and all its data from Firestore
    Used when a compound is no longer needed
    """
    try:
        db.collection("compounds").document(compound_id).delete()
        return jsonify({"success": True}), 200
    except Exception as e:
        return jsonify({"error": str(e)}), 500
    

@app.route("/similar-compounds", methods=["POST"])
def similar_compounds():
    """
    API ENDPOINT: Find compounds with similar molecular structures
    Uses chemical fingerprinting to calculate similarity scores
    and returns compounds that are structurally similar to the input
    """
    data = request.get_json()
    smiles = data.get("smiles")
    threshold = float(data.get("threshold", 0.7))  # Default threshold for similarity

    if not smiles:
        return jsonify({"error": "No SMILES provided"}), 400

    query_mol = Chem.MolFromSmiles(smiles)
    if query_mol is None:
        return jsonify([])

    query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=2048)

    similar = []
    for doc in db.collection("compounds").stream():
        compound = doc.to_dict()
        compound_smiles = compound.get("smiles", "")
        mol = Chem.MolFromSmiles(compound_smiles)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            similarity = DataStructs.TanimotoSimilarity(query_fp, fp)
            if similarity >= threshold and compound["id"] != data.get("id"):
                compound["similarity"] = similarity
                similar.append(compound)

    similar.sort(key=lambda x: x["similarity"], reverse=True)
    return jsonify(similar)
    

@app.route("/filter-phase-map", methods=["POST"])
def filter_phase_map():
    """
    API ENDPOINT: Filter compounds by phase transitions and temperatures
    Searches for compounds that have specific phase transitions
    (like crystal to liquid) at specific temperatures
    """
    data = request.get_json()
    selected_transition = data.get("transition")
    target_temp = data.get("temperature")

    print("Received filter:", selected_transition, target_temp)

    results = []
    canonical_selected = canonicalize_transition(selected_transition)
    for doc in db.collection("compounds").stream():
        compound = doc.to_dict()
        parsed_transitions = compound.get("parsed_phase_transitions", [])

        # Use canonical form for matching
        if any(item.strip().upper().startswith(canonical_selected) for item in parsed_transitions):
            if target_temp is None or check_temp_near(parsed_transitions, selected_transition, target_temp):
                results.append(compound)

    return jsonify(results)


def normalize_transitions(phase_map_str):
    """
    HELPER FUNCTION: Clean up and standardize phase transition text
    Takes raw phase map text and converts it to a standard format
    for easier searching and comparison
    """
    if not isinstance(phase_map_str, str):
        return []
    
    transitions = set()
    entries = phase_map_str.split(";")
    for entry in entries:
        parts = re.split(r'\s*-\s*', entry.strip())
        phases = [
            p.strip().upper().rstrip("?") for p in parts
            if re.match(r'^[A-Z]{1,4}$|^RT$|^GLASS$', p.strip(), re.IGNORECASE)
        ]
        for i in range(len(phases) - 1):
            transitions.add(f"{phases[i]} - {phases[i+1]}")
    return transitions


def check_temp_near(transition_list, transition, input_temp_str, default_tolerance=10):
    """
    HELPER FUNCTION: Check if a compound has a transition near a target temperature
    Looks through a compound's phase transitions to see if any match
    the requested transition type and temperature range
    """
    try:
        input_temp_str = str(input_temp_str).strip()
        if "-" in input_temp_str and not input_temp_str.endswith("-"):
            low, high = map(float, input_temp_str.split("-"))
            temp_range = lambda t: low <= t <= high
        elif input_temp_str.endswith("+"):
            low = float(input_temp_str[:-1])
            temp_range = lambda t: t >= low
        elif input_temp_str.endswith("-"):
            high = float(input_temp_str[:-1])
            temp_range = lambda t: t <= high
        else:
            target = float(input_temp_str)
            temp_range = lambda t: abs(t - target) <= default_tolerance
    except Exception as e:
        print("Temp parse error:", e)
        return False

    for item in transition_list:
        if item.strip().upper().startswith(transition.strip().upper()):
            match = re.search(r"@\s*([-+]?[0-9]+(?:\.[0-9]+)?)", item)
            if match:
                try:
                    temp_val = float(match.group(1))
                    if temp_range(temp_val):
                        return True
                except ValueError:
                    continue
    return False


@app.route("/lots", methods=["GET"])
def list_lots():
    """
    API ENDPOINT: Get list of all lot names
    Returns just the names of all lots in the database
    Used to populate dropdown menus and lot selection lists
    """
    lots_ref = db.collection("lots")
    docs = lots_ref.stream()
    lots = [doc.id for doc in docs]
    return jsonify(lots)


@app.route("/lot/<lot_id>", methods=["GET"])
def get_lot_compounds(lot_id):
    """
    API ENDPOINT: Get all compounds in a specific lot
    Takes a lot name and returns all the compound samples
    that belong to that manufacturing lot
    """
    lot_ref = db.collection("lots").document(lot_id).collection("compounds")
    docs = lot_ref.stream()
    compounds = [doc.to_dict() for doc in docs]
    return jsonify(compounds)


@app.route("/create-lot", methods=["POST"])
def create_lot():
    """
    API ENDPOINT: Create a new manufacturing lot
    Takes a compound and creates a new lot for testing samples
    Copies the compound data but clears testing-specific fields
    """
    data = request.get_json()
    compound_ids = data.get("compoundIds")
    lot_name = data.get("lotName")

    if not compound_ids or len(compound_ids) != 1:
        return jsonify({"error": "Only one compound can be added per lot"}), 400

    if not lot_name or lot_name.strip() == "":
        return jsonify({"error": "Lot name is required"}), 400

    compound_id = compound_ids[0]
    lot_name = lot_name.strip()

    # Check if a lot with the same name already exists
    if db.collection("lots").document(lot_name).get().exists:
        return jsonify({"error": "A lot with that name already exists"}), 400

    compound_doc = db.collection("compounds").document(compound_id).get()
    if not compound_doc.exists:
        return jsonify({"error": "Compound does not exist"}), 400

    compound = compound_doc.to_dict()

    # Clear out specific fields for retesting
    excluded_fields = {
        "J/g DSC melt (total)", "kJ/mol DSC melt (total)",
        "Refractive index (ne/no)", "Notes",
        "Lambda Max (DCM/AcCN)", "Lambda Max (neat film)",
        "phase map", "r33", "attachments"
    }
    for field in excluded_fields:
        compound[field] = ""

    new_id = str(uuid.uuid4())
    compound["id"] = new_id
    compound["original_id"] = compound_id

    db.collection("lots").document(lot_name).set({
        "compound_id": compound_id,
        "created_at": firestore.SERVER_TIMESTAMP
    })

    db.collection("lots").document(lot_name).collection("compounds").document(new_id).set(compound)

    return jsonify({"success": True, "lotName": lot_name, "newCompoundId": new_id}), 200


@app.route("/lots-for-compound/<compound_id>")
def lots_for_compound(compound_id):
    """
    API ENDPOINT: Find all lots that contain a specific compound
    Takes a compound ID and returns all the manufacturing lots
    that were made from that base compound
    """
    compound_id = compound_id.lower()  # normalize
    matched_lots = []
    for doc in db.collection("lots").stream():
        lot_data = doc.to_dict()
        if lot_data.get("compound_id", "").lower() == compound_id:
            matched_lots.append(doc.id)
    return jsonify(matched_lots)


@app.route("/update-lot-compound", methods=["POST"])
def update_lot_compound():
    """
    API ENDPOINT: Update test results for a compound in a lot
    Saves new testing data (like phase maps, optical properties)
    for a specific compound sample in a manufacturing lot
    """
    data = request.get_json()
    lot_id = data.get("lotId")
    compound_id = data.get("id")

    if not lot_id or not compound_id:
        return jsonify({"error": "Missing lotId or compound id"}), 400

    try:
        db.collection("lots").document(lot_id).collection("compounds").document(compound_id).set(data, merge=True)
        return jsonify({"success": True}), 200
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/delete-lot-compound/<lot_id>/<compound_id>", methods=["DELETE"])
def delete_lot_compound(lot_id, compound_id):
    """
    API ENDPOINT: Remove a compound sample from a lot
    Deletes a specific compound from a manufacturing lot
    If it's the last compound in the lot, deletes the entire lot
    """
    try:
        # Delete the compound from the lot
        compound_ref = db.collection("lots").document(lot_id).collection("compounds").document(compound_id)
        compound_ref.delete()

        # Check if the lot is now empty
        remaining = list(db.collection("lots").document(lot_id).collection("compounds").stream())
        if len(remaining) == 0:
            # Delete the parent lot document as well
            db.collection("lots").document(lot_id).delete()

        return jsonify({"success": True}), 200

    except Exception as e:
        return jsonify({"error": str(e)}), 500
    

@app.route("/batch-lots-for-compounds", methods=["POST"])
def batch_lots_for_compounds():
    """
    API ENDPOINT: Get lot information for multiple compounds at once
    Takes a list of compound IDs and returns which lots contain
    each compound - more efficient than making individual requests
    """
    data = request.get_json() or {}
    compound_ids = data.get("compound_ids")

    if not isinstance(compound_ids, list):
        return jsonify({"error": "compound_ids must be a list"}), 400

    compound_ids = [cid.lower() for cid in compound_ids if isinstance(cid, str)]

    lot_map = {}
    for doc in db.collection("lots").stream():
        lot = doc.to_dict()
        compound_id = lot.get("compound_id", "").lower()
        if compound_id in compound_ids:
            lot_map.setdefault(compound_id, []).append(doc.id)

    return jsonify(lot_map)


@app.route("/upload-image-to-firebase", methods=["POST"])
def upload_image_to_firebase():
    """
    API ENDPOINT: Upload images to cloud storage
    Takes image files (photos, spectra, etc.) and stores them
    in Firebase Storage, returning a public URL for the image
    """
    if 'file' not in request.files:
        return jsonify({"error": "No file part"}), 400

    uploaded_file = request.files['file']
    note = request.form.get("note", "")
    filename = uploaded_file.filename

    if filename == "":
        return jsonify({"error": "No selected file"}), 400

    try:
        bucket = storage.bucket()
        blob = bucket.blob(f"data_images/{filename}")
        blob.upload_from_file(uploaded_file, content_type=uploaded_file.content_type)
        blob.make_public()  # optional: allows public access

        return jsonify({
            "fileUrl": blob.public_url,
            "note": note,
        }), 200

    except Exception as e:
        return jsonify({"error": str(e)}), 500
    
@app.route("/regenerate-images", methods=["POST"])
def regenerate_images():
    """
    API ENDPOINT: Regenerates structure images for a compound after SMILES change
    Replaces imageUrl and bwImageUrl in Firestore and Storage
    """
    try:
        data = request.get_json()
        compound_id = data.get("id")
        smiles = data.get("smiles")

        if not compound_id or not smiles:
            return jsonify({"error": "Missing compound ID or SMILES"}), 400

        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return jsonify({"error": "Invalid SMILES"}), 400

        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D
        import tempfile
        import os

        rdDepictor.Compute2DCoords(mol)
        width, height = 2400, 1200
        bucket = storage.bucket()

        # Regenerate color image
        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        options = drawer.drawOptions()
        options.minFontSize = 90
        options.maxFontSize = 90
        options.bondLineWidth = 4
        options.padding = 0.01
        options.fixedScale = True
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            drawer.WriteDrawingText(tmp.name)
            blob = bucket.blob(f"compound_images/{compound_id}.png")
            blob.upload_from_filename(tmp.name, content_type="image/png")
            blob.make_public()
            imageUrl = blob.public_url
            os.unlink(tmp.name)

        # Regenerate BW image
        drawer_bw = rdMolDraw2D.MolDraw2DCairo(width, height)
        opts_bw = drawer_bw.drawOptions()
        opts_bw.minFontSize = 90
        opts_bw.maxFontSize = 90
        opts_bw.bondLineWidth = 4
        opts_bw.padding = 0.01
        opts_bw.fixedScale = True
        opts_bw.useBWAtomPalette()
        drawer_bw.DrawMolecule(mol)
        drawer_bw.FinishDrawing()
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp_bw:
            drawer_bw.WriteDrawingText(tmp_bw.name)
            blob_bw = bucket.blob(f"compound_images/{compound_id}_bw.png")
            blob_bw.upload_from_filename(tmp_bw.name, content_type="image/png")
            blob_bw.make_public()
            bwImageUrl = blob_bw.public_url
            os.unlink(tmp_bw.name)

        # Update Firestore
        db.collection("compounds").document(compound_id).update({
            "imageUrl": imageUrl,
            "bwImageUrl": bwImageUrl
        })

        return jsonify({"imageUrl": imageUrl, "bwImageUrl": bwImageUrl}), 200

    except Exception as e:
        print("Error regenerating images:", e)
        return jsonify({"error": str(e)}), 500

    

@app.route("/create-formulation", methods=["POST"])
def create_formulation():
    """
    API ENDPOINT: Create a new formulation recipe
    Takes a list of compounds with percentages and creates
    a formulation recipe that can be used for manufacturing
    """
    try:
        data = request.get_json()

        if not data.get("components"):
            return jsonify({"error": "No components provided"}), 400

        formulation_ref = db.collection("formulations").document()
        data["createdAt"] = firestore.SERVER_TIMESTAMP
        formulation_ref.set(data)

        return jsonify({"success": True, "id": formulation_ref.id}), 200
    except Exception as e:
        return jsonify({"error": str(e)}), 500
    



@app.route("/formulations", methods=["GET"])
def list_formulations():
    """
    API ENDPOINT: Get all formulation recipes
    Returns a list of all formulations stored in the database
    Used to display the formulations table on the frontend
    """
    try:
        docs = db.collection("formulations").stream()
        formulations = [{ "id": doc.id, **doc.to_dict() } for doc in docs]
        return jsonify(formulations), 200
    except Exception as e:
        return jsonify({"error": str(e)}), 500
    

@app.route("/delete-formulation/<formulation_id>", methods=["DELETE"])
def delete_formulation(formulation_id):
    """
    API ENDPOINT: Delete a formulation recipe
    Permanently removes a formulation and all its data from the database
    Used when a formulation is no longer needed
    """
    try:
        db.collection("formulations").document(formulation_id).delete()
        return jsonify({"success": True}), 200
    except Exception as e:
        return jsonify({"error": str(e)}), 500
    

@app.route("/update-formulation/<formulation_id>", methods=["POST"])
def update_formulation(formulation_id):
    """
    API ENDPOINT: Update an existing formulation recipe
    Takes formulation data with changes and saves the updated
    recipe information back to the database
    """
    try:
        data = request.get_json()
        db.collection("formulations").document(formulation_id).update(data)
        return jsonify({"success": True}), 200
    except Exception as e:
        return jsonify({"error": str(e)}), 500
    

# Dashboard logic
@app.route('/compound-analytics', methods=['GET'])
def compound_analytics():
    """
    API ENDPOINT: Generate analytics dashboard data
    Calculates statistics about compounds (total count, average molecular weight,
    creation timeline) for displaying charts and graphs on the dashboard
    """
    try:
        compounds_ref = db.collection('compounds')
        docs = compounds_ref.stream()
        compounds = [doc.to_dict() for doc in docs]

        if not compounds:
            return jsonify({"error": "No compounds found"}), 404

        total_compounds = len(compounds)
        most_viewed = max(compounds, key=lambda c: c.get("views", 0))

        # Only include compounds that have MW
        mw_values = [c.get("MW") for c in compounds if "MW" in c and isinstance(c.get("MW"), (int, float))]

        avg_mw = round(sum(mw_values) / len(mw_values), 2) if mw_values else "N/A"

        # Build MW histogram (simple bins)
        bins = ["0-500", "501-700", "701-900", "901-1100", "1101+"]
        bin_counts = {bin: 0 for bin in bins}

        for c in compounds:
            mw = c.get("MW")
            if mw is None:
                continue
            if mw <= 500:
                bin_counts["0-500"] += 1
            elif mw <= 700:
                bin_counts["501-700"] += 1
            elif mw <= 900:
                bin_counts["701-900"] += 1
            elif mw <= 1100:
                bin_counts["901-1100"] += 1
            else:
                bin_counts["1101+"] += 1

        mw_distribution = [{"range": bin, "count": count} for bin, count in bin_counts.items()]

        # ✅ Build Creation Timeline
        creation_data = {}
        for c in compounds:
            created_at = c.get("createdAt")
            if created_at and isinstance(created_at, dict) and "seconds" in created_at:
                timestamp = datetime.datetime.fromtimestamp(created_at["seconds"])
                month_str = timestamp.strftime("%Y-%m")
                creation_data[month_str] = creation_data.get(month_str, 0) + 1

        timeline = [{"month": k, "count": v} for k, v in sorted(creation_data.items())]

        return jsonify({
            "totalCompounds": total_compounds,
            "mostViewedCompound": most_viewed.get("id", "N/A"),
            "averageMW": avg_mw,
            "mwDistribution": mw_distribution,
            "creationData": timeline
        })
    except Exception as e:
        return jsonify({"error": str(e)}), 500



    
@app.route('/save-starred', methods=['POST'])
def save_starred():
    """
    API ENDPOINT: Save user's favorite compounds
    Stores which compounds a user has marked as favorites
    so they can quickly find them later
    """
    data = request.get_json()
    email = data.get('email')
    starred = data.get('starred')
    if not email or not isinstance(starred, list):
        return jsonify({'error': 'Missing email or starred list'}), 400
    try:
        # Save in Firestore under collection 'user_starred', document id is email
        db.collection('user_starred').document(email).set({
            'starred': starred,
            'updated_at': firestore.SERVER_TIMESTAMP
        }, merge=True)
        return jsonify({'success': True}), 200
    except Exception as e:
        print('Failed to save starred compounds:', e)
        return jsonify({'error': str(e)}), 500
    

    # Get starred compounds for a user
@app.route('/get-starred', methods=['GET'])
def get_starred():
    """
    API ENDPOINT: Get user's favorite compounds
    Returns the list of compounds that a user has
    previously marked as favorites
    """
    email = request.args.get('email')
    if not email:
        return jsonify({'error': 'Missing email'}), 400
    try:
        doc = db.collection('user_starred').document(email).get()
        if not doc.exists:
            return jsonify({'starred': []}), 200
        data = doc.to_dict()
        return jsonify({'starred': data.get('starred', [])}), 200
    except Exception as e:
        print('Failed to fetch starred compounds:', e)
        return jsonify({'error': str(e)}), 500
    
    
@app.route('/save-starred-formulations', methods=['POST'])
def save_starred_formulations():
    """
    API ENDPOINT: Save user's favorite formulations
    Stores which formulation recipes a user has marked as favorites
    so they can quickly find them later
    """
    data = request.get_json()
    email = data.get('email')
    starred = data.get('starred')
    if not email or not isinstance(starred, list):
        return jsonify({'error': 'Missing email or starred list'}), 400
    try:
        db.collection('user_starred_formulations').document(email).set({
            'starred': starred,
            'updated_at': firestore.SERVER_TIMESTAMP
        }, merge=True)
        return jsonify({'success': True}), 200
    except Exception as e:
        print('Failed to save starred formulations:', e)
        return jsonify({'error': str(e)}), 500


@app.route('/get-starred-formulations', methods=['GET'])
def get_starred_formulations():
    """
    API ENDPOINT: Get user's favorite formulations
    Returns the list of formulation recipes that a user has
    previously marked as favorites
    """
    email = request.args.get('email')
    if not email:
        return jsonify({'error': 'Missing email'}), 400
    try:
        doc = db.collection('user_starred_formulations').document(email).get()
        if not doc.exists:
            return jsonify({'starred': []}), 200
        data = doc.to_dict()
        return jsonify({'starred': data.get('starred', [])}), 200
    except Exception as e:
        print('Failed to fetch starred formulations:', e)
        return jsonify({'error': str(e)}), 500


@app.route("/plans", methods=["GET"])
def list_plans():
    """
    API ENDPOINT: Get all synthesis plans
    Returns a list of all planned compound syntheses
    Used to display the synthesis planning table
    """
    try:
        docs = db.collection("plans").stream()
        plans = [{"id": doc.id, **doc.to_dict()} for doc in docs]
        return jsonify(plans), 200
    except Exception as e:
        return jsonify({"error": str(e)}), 500
    

@app.route("/create-plan", methods=["POST"])
def create_plan():
    """
    API ENDPOINT: Create a new synthesis plan
    Takes planned compound data and creates a synthesis plan entry
    Generates molecular images and calculates molecular weight if SMILES provided
    """
    try:
        data = request.get_json()
        # Ensure completed is always set to False by default
        if "completed" not in data:
            data["completed"] = False
        # Generate and upload structure image if smiles is present and no imageUrl
        smiles = data.get("smiles")
        if smiles and not data.get("imageUrl"):
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    data["MW"] = round(Descriptors.MolWt(mol), 2)
                if mol:
                    rdDepictor.Compute2DCoords(mol)

                    # Use even larger canvas and adjust draw options for bigger molecule and slightly smaller font
                    width, height = 2400, 1200  # Keep large canvas
                    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
                    options = drawer.drawOptions()
                    num_atoms = mol.GetNumAtoms()
                    # Dynamically set both minFontSize and maxFontSize for atom labels
                    if num_atoms <= 20:
                        options.minFontSize = 120
                        options.maxFontSize = 120
                    elif num_atoms <= 40:
                        options.minFontSize = 90
                        options.maxFontSize = 90
                    else:
                        options.minFontSize = 58
                        options.maxFontSize = 58
                    options.bondLineWidth = 4
                    options.padding = 0.01    # Minimal padding to fill image
                    options.fixedScale = True # Force molecule to fill canvas

                    # options.useBWAtomPalette()

                    drawer.DrawMolecule(mol)
                    drawer.FinishDrawing()
                    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
                        drawer.WriteDrawingText(tmp.name)
                        tmp.flush()

                        bucket = storage.bucket()
                        filename = f"plan_images/{data.get('id', 'plan')}.png"
                        blob = bucket.blob(filename)
                        blob.upload_from_filename(tmp.name, content_type="image/png")
                        blob.make_public()
                        data["imageUrl"] = blob.public_url
                    os.unlink(tmp.name)
            except Exception as e:
                print(f"Failed to generate/upload image for plan {data.get('id')}: {e}")
        # Use the provided id as the Firestore document ID
        plan_id = data.get("id")
        if not plan_id:
            plan_ref = db.collection("plans").document()  # fallback to auto-id if not provided
        else:
            plan_ref = db.collection("plans").document(str(plan_id))
        data["createdAt"] = firestore.SERVER_TIMESTAMP
        plan_ref.set(data)
        return jsonify({"success": True}), 200
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/delete-plan/<plan_id>", methods=["DELETE"])
def delete_plan(plan_id):
    """
    API ENDPOINT: Delete a synthesis plan
    Permanently removes a planned synthesis from the database
    Used when a synthesis plan is no longer needed
    """
    try:
        db.collection("plans").document(plan_id).delete()
        return jsonify({"success": True}), 200
    except Exception as e:
        return jsonify({"error": str(e)}), 500
    

@app.route("/update-plan/<plan_id>", methods=["POST"])
def update_plan(plan_id):
    """
    API ENDPOINT: Update an existing synthesis plan
    Takes plan data with changes (like marking as completed)
    and saves the updated information back to the database
    """
    try:
        data = request.get_json()
        if not data:
            return jsonify({"error": "No data provided"}), 400

        doc_ref = db.collection("plans").document(plan_id)
        if not doc_ref.get().exists:
            return jsonify({"error": "Plan not found"}), 404

        doc_ref.update(data)  # Use update instead of set
        return jsonify({"success": True}), 200
    except Exception as e:
        return jsonify({"error": str(e)}), 500


if __name__ == '__main__':
    app.run(debug=True)


#   Copyright © 2025 Polaris Electro Optics
#   This code is the property of Polaris Electro Optics and may not be reused, modified, or distributed without explicit permission.

