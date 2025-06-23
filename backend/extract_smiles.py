import firebase_admin
from firebase_admin import credentials, firestore
import csv

# Initialize Firebase
cred = credentials.Certificate("firebase-key.json")
firebase_admin.initialize_app(cred)
db = firestore.client()

output_file = "compounds_missing_images.csv"

compounds_ref = db.collection("compounds")
missing = []

for doc in compounds_ref.stream():
    data = doc.to_dict()
    if not data.get("imageUrl") and data.get("smiles"):
        missing.append({"id": data.get("id", doc.id), "smiles": data["smiles"]})

with open(output_file, "w", newline="") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=["id", "smiles"])
    writer.writeheader()
    writer.writerows(missing)

print(f"Exported {len(missing)} compounds to {output_file}")