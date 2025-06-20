import firebase_admin
from firebase_admin import credentials, firestore

# 🔐 Path to your service account key
cred = credentials.Certificate("firebase-key.json")
firebase_admin.initialize_app(cred)

db = firestore.client()

field_to_remove = "r33 (neat, assuming ne=1.7, no=1.5)"

# Go through each lot
lots_ref = db.collection("lots")
lots = lots_ref.stream()

for lot in lots:
    lot_id = lot.id
    compounds_ref = lots_ref.document(lot_id).collection("compounds")
    compounds = compounds_ref.stream()

    for compound in compounds:
        compound_id = compound.id
        compound_data = compound.to_dict()

        if field_to_remove in compound_data:
            print(f"Deleting '{field_to_remove}' from {compound_id} in lot {lot_id}")
            compounds_ref.document(compound_id).update({
                field_to_remove: firestore.DELETE_FIELD
            })

print("Cleanup complete.")
