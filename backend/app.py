import os
from flask import Flask, jsonify, request
import firebase_admin
from firebase_admin import credentials, firestore
from flask_cors import CORS
from rdkit import Chem
import re

app = Flask(__name__)
# CORS(app, origins="*", allow_headers="*", supports_credentials=True, methods=["GET", "POST", "OPTIONS", "DELETE"])
CORS(app, resources={r"/*": {"origins": "http://localhost:3000"}})  # Match frontend port



# Initialize Firebase Admin SDK
cred = credentials.Certificate("firebase-key.json")
firebase_admin.initialize_app(cred)
db = firestore.client()

@app.route('/compounds')
def get_compounds():
    compounds_ref = db.collection('compounds')
    docs = compounds_ref.stream()
    compounds = [doc.to_dict() for doc in docs]
    return jsonify(compounds)


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
        return '', 200  # Preflight success response

    data = request.get_json()
    if "id" not in data:
        return jsonify({"error": "Compound must have an 'id' field"}), 400

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
    
@app.route("/filter-phase-map", methods=["POST"])
def filter_phase_map():
    data = request.get_json()
    selected_transition = data.get("transition")
    target_temp = data.get("temperature")

    print("Received filter:", selected_transition, target_temp)

    results = []
    for doc in db.collection("compounds").stream():
        compound = doc.to_dict()
        parsed_transitions = compound.get("parsed_phase_transitions", [])

        if any(item.strip().upper().startswith(selected_transition.strip().upper()) for item in parsed_transitions):
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






if __name__ == '__main__':
    app.run(debug=True)



