import os
from flask import Flask, jsonify, request
import firebase_admin
from firebase_admin import credentials, firestore, storage
from flask_cors import CORS
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem, DataStructs
import re
import uuid
from google.oauth2 import service_account
from googleapiclient.discovery import build
from googleapiclient.http import MediaFileUpload
import tempfile
import base64

app = Flask(__name__)
# CORS(app, origins="*", allow_headers="*", supports_credentials=True, methods=["GET", "POST", "OPTIONS", "DELETE"])
CORS(app, resources={r"/*": {"origins": "http://localhost:3000"}})  # Match frontend port



# Initialize Firebase Admin SDK
cred = credentials.Certificate("firebase-key.json")
firebase_admin.initialize_app(cred, {
    'storageBucket': 'polaris-test-3b8f8.firebasestorage.app'
})
db = firestore.client()




# Only treat FN - FNG and FNG - FN as equivalent, not all transitions
FN_FNG_EQUIV = {('FN', 'FNG'), ('FNG', 'FN')}

def canonicalize_transition(transition: str) -> str:
    """
    Only treat 'FN - FNG' and 'FNG - FN' as the same. All other transitions keep their order.
    """
    parts = [p.strip().upper() for p in transition.split('-')]
    if len(parts) == 2:
        pair = (parts[0], parts[1])
        if pair in FN_FNG_EQUIV:
            return 'FN - FNG'
        return f'{parts[0]} - {parts[1]}'
    return transition.strip().upper()

def extract_transitions_with_temps(phase_map_str: str):
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
    compounds_ref = db.collection('compounds')
    docs = compounds_ref.stream()
    compounds = [doc.to_dict() for doc in docs]
    return jsonify(compounds)


@app.route("/compounds/<compound_id>")
def get_single_compound(compound_id):
    try:
        doc = db.collection("compounds").document(compound_id).get()
        if not doc.exists:
            return jsonify({"error": "Compound not found"}), 404
        return jsonify(doc.to_dict()), 200
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/search-substructure", methods=["POST"])
def search_substructure():
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
    
@app.route("/add-compound", methods=["POST", "OPTIONS"])
def add_compound():
    print("🧪 Received method:", request.method)
    if request.method == "OPTIONS":
        return '', 200

    data = request.get_json()

    # 🚫 Prevent saving lot-based compounds directly
    if "original_id" in data:
        return jsonify({"error": "Cannot add a compound that originates from a lot."}), 400

    if "id" not in data:
        return jsonify({"error": "Compound must have an 'id' field"}), 400

    phase_map_str = data.get("phase map", "")
    transitions = extract_transitions_with_temps(phase_map_str)
    data["parsed_phase_transitions"] = transitions

    # Generate and upload structure image if smiles is present and no imageUrl
    smiles = data.get("smiles")
    if smiles and not data.get("imageUrl"):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # Generate a high-resolution image
                img = Draw.MolToImage(mol, size=(1200, 600), dpi=300)
                with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
                    img.save(tmp.name, dpi=(300, 300))
                    tmp.flush()
                    # Upload to Firebase Storage
                    bucket = storage.bucket()
                    filename = f"compound_images/{data['id']}.png"
                    blob = bucket.blob(filename)
                    blob.upload_from_filename(tmp.name, content_type="image/png")
                    blob.make_public()
                    data["imageUrl"] = blob.public_url
                os.unlink(tmp.name)
        except Exception as e:
            print(f"Failed to generate/upload image for {data['id']}: {e}")

    db.collection("compounds").document(data["id"]).set(data, merge=True)
    return jsonify({"success": True}), 200



@app.route("/update-compound", methods=["POST"])
def update_compound():
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
    try:
        db.collection("compounds").document(compound_id).delete()
        return jsonify({"success": True}), 200
    except Exception as e:
        return jsonify({"error": str(e)}), 500
    
@app.route("/similar-compounds", methods=["POST"])
def similar_compounds():
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
    lots_ref = db.collection("lots")
    docs = lots_ref.stream()
    lots = [doc.id for doc in docs]
    return jsonify(lots)


@app.route("/lot/<lot_id>", methods=["GET"])
def get_lot_compounds(lot_id):
    lot_ref = db.collection("lots").document(lot_id).collection("compounds")
    docs = lot_ref.stream()
    compounds = [doc.to_dict() for doc in docs]
    return jsonify(compounds)


@app.route("/create-lot", methods=["POST"])
def create_lot():
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
    compound_id = compound_id.lower()  # normalize
    matched_lots = []
    for doc in db.collection("lots").stream():
        lot_data = doc.to_dict()
        if lot_data.get("compound_id", "").lower() == compound_id:
            matched_lots.append(doc.id)
    return jsonify(matched_lots)





@app.route("/update-lot-compound", methods=["POST"])
def update_lot_compound():
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



    


@app.route("/create-formulation", methods=["POST"])
def create_formulation():
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
    try:
        docs = db.collection("formulations").stream()
        formulations = [{ "id": doc.id, **doc.to_dict() } for doc in docs]
        return jsonify(formulations), 200
    except Exception as e:
        return jsonify({"error": str(e)}), 500
    

@app.route("/delete-formulation/<formulation_id>", methods=["DELETE"])
def delete_formulation(formulation_id):
    try:
        db.collection("formulations").document(formulation_id).delete()
        return jsonify({"success": True}), 200
    except Exception as e:
        return jsonify({"error": str(e)}), 500
    
@app.route("/update-formulation/<formulation_id>", methods=["POST"])
def update_formulation(formulation_id):
    try:
        data = request.get_json()
        db.collection("formulations").document(formulation_id).update(data)
        return jsonify({"success": True}), 200
    except Exception as e:
        return jsonify({"error": str(e)}), 500











if __name__ == '__main__':
    app.run(debug=True)
