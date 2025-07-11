import firebase_admin
from firebase_admin import credentials, firestore, storage
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
import tempfile
import os

# Initialize Firebase Admin SDK
cred = credentials.Certificate("firebase-key.json")
firebase_admin.initialize_app(cred, {
    'storageBucket': 'polaris-test-3b8f8.firebasestorage.app'
})
db = firestore.client()
bucket = storage.bucket()

def regenerate_and_upload_images():
    compounds_ref = db.collection('compounds')
    docs = compounds_ref.stream()
    for doc in docs:
        compound = doc.to_dict()
        compound_id = doc.id
        smiles = compound.get('smiles')
        if not smiles:
            print(f"Compound {compound_id} missing SMILES, skipping.")
            continue

        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            print(f"Compound {compound_id} invalid SMILES, skipping.")
            continue
        rdDepictor.Compute2DCoords(mol)

        width, height = 2400, 1200
        # Color image
        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        options = drawer.drawOptions()
        num_atoms = mol.GetNumAtoms()
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
        options.padding = 0.01
        options.fixedScale = True
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            drawer.WriteDrawingText(tmp.name)
            tmp.flush()
            blob = bucket.blob(f"compound_images/{compound_id}.png")
            blob.upload_from_filename(tmp.name, content_type="image/png")
            blob.make_public()
        os.unlink(tmp.name)

        # Black-and-white image
        drawer_bw = rdMolDraw2D.MolDraw2DCairo(width, height)
        options_bw = drawer_bw.drawOptions()
        options_bw.minFontSize = options.minFontSize
        options_bw.maxFontSize = options.maxFontSize
        options_bw.bondLineWidth = options.bondLineWidth
        options_bw.padding = options.padding
        options_bw.fixedScale = options.fixedScale
        options_bw.useBWAtomPalette()
        drawer_bw.DrawMolecule(mol)
        drawer_bw.FinishDrawing()
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp_bw:
            drawer_bw.WriteDrawingText(tmp_bw.name)
            tmp_bw.flush()
            blob_bw = bucket.blob(f"compound_images/{compound_id}_bw.png")
            blob_bw.upload_from_filename(tmp_bw.name, content_type="image/png")
            blob_bw.make_public()
        os.unlink(tmp_bw.name)

        print(f"Uploaded new images for compound {compound_id}")

if __name__ == "__main__":
    regenerate_and_upload_images()