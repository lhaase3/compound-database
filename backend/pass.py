import os
from flask import Flask, jsonify, request
import firebase_admin
from firebase_admin import credentials, firestore
from flask_cors import CORS
from rdkit import Chem
import csv
import re


app = Flask(__name__)
CORS(app, origins="*", allow_headers="*", supports_credentials=True, methods=["GET", "POST", "OPTIONS", "DELETE"])


# Initialize Firebase Admin SDK
cred = credentials.Certificate("firebase-key.json")
firebase_admin.initialize_app(cred)
db = firestore.client()

# === Parsing Function ===
def extract_transitions_with_temps(phase_map_str: str):
    if not isinstance(phase_map_str, str):
        return []

    results = []
    entries = phase_map_str.split(";")
    for entry in entries:
        parts = re.split(r'\s*-\s*', entry.strip())
        i = 0
        while i < len(parts) - 1:
            phase1 = parts[i].strip().upper().rstrip("?")
            i += 1

            temp = None
            if i < len(parts):
                temp_match = re.search(r"[-+]?[0-9]+(?:\.[0-9]+)?", parts[i])
                if temp_match:
                    temp = temp_match.group()
                    i += 1

            if i < len(parts):
                phase2 = parts[i].strip().upper().rstrip("?")
                if re.match(r'^[A-Z]{1,4}$|^RT$', phase1) and re.match(r'^[A-Z]{1,4}$|^RT$', phase2):
                    label = f"{phase1} - {phase2}"
                    if temp:
                        label += f" @ {temp}"
                    results.append(label)
            i += 1
    return results

# === Main update function ===
def update_firebase_with_transitions(csv_path):
    with open(csv_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            cid = row.get("id")
            phase_map = row.get("phase map", "")
            if cid:
                transitions = extract_transitions_with_temps(phase_map)
                doc_ref = db.collection("compounds").document(cid)
                doc_ref.update({"parsed_phase_transitions": transitions})
                print(f"✅ Updated {cid}")

# === Run if this file is executed ===
if __name__ == "__main__":
    update_firebase_with_transitions("phase_map - Sheet1.csv")  # Update with your actual filename
