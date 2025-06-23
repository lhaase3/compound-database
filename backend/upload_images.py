import os
import firebase_admin
from firebase_admin import credentials, firestore, storage

# --- CONFIG ---
IMAGES_DIR = "./compound_images"  # Folder with images named by compound id (e.g., PEO-0100.png)
FIREBASE_KEY = "path/to/firebase-key.json"
BUCKET_NAME = "polaris-test-3b8f8.firebasestorage.app"  # Update to your bucket

# --- INIT FIREBASE ---
cred = credentials.Certificate("firebase-key.json")
firebase_admin.initialize_app(cred, {'storageBucket': BUCKET_NAME})
db = firestore.client()
bucket = storage.bucket()

# --- MAIN SCRIPT ---
for filename in os.listdir(IMAGES_DIR):
    if not (filename.endswith(".png") or filename.endswith(".jpg") or filename.endswith(".jpeg")):
        continue
    compound_id = os.path.splitext(filename)[0]
    local_path = os.path.join(IMAGES_DIR, filename)
    blob = bucket.blob(f"compound_images/{filename}")
    blob.upload_from_filename(local_path)
    blob.make_public()
    image_url = blob.public_url

    # Update Firestore
    doc_ref = db.collection("compounds").document(compound_id)
    doc_ref.set({"imageUrl": image_url}, merge=True)
    print(f"Updated {compound_id} with image {image_url}")

print("Done.")