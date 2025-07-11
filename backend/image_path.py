import firebase_admin
from firebase_admin import credentials, firestore

# Initialize Firebase Admin SDK
cred = credentials.Certificate("firebase-key.json")
firebase_admin.initialize_app(cred)
db = firestore.client()

BUCKET = "polaris-test-3b8f8.firebasestorage.app"

def make_firebase_url(compound_id, bw=False):
    filename = f"{compound_id}_bw.png" if bw else f"{compound_id}.png"
    # URL encode the path
    object_path = f"compound_images%2F{filename}"
    return f"https://firebasestorage.googleapis.com/v0/b/{BUCKET}/o/{object_path}?alt=media"

def update_compound_image_urls():
    compounds_ref = db.collection('compounds')
    docs = compounds_ref.stream()
    for doc in docs:
        compound_id = doc.id
        image_url = make_firebase_url(compound_id)
        bw_image_url = make_firebase_url(compound_id, bw=True)
        update_data = {
            "imageUrl": image_url,
            "bwImageUrl": bw_image_url
        }
        db.collection('compounds').document(compound_id).update(update_data)
        print(f"Updated {compound_id}: imageUrl -> {image_url}, bwImageUrl -> {bw_image_url}")

if __name__ == "__main__":
    update_compound_image_urls()