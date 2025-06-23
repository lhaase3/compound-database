from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import os

def smiles_to_png(smiles, out_path, mol_size=(400, 200)):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol, size=mol_size)
        img.save(out_path)
        return True
    return False

# Example usage with the CSV from above
df = pd.read_csv("compounds_missing_images.csv")
os.makedirs("compound_images", exist_ok=True)

for _, row in df.iterrows():
    out_file = os.path.join("compound_images", f"{row['id']}.png")
    if smiles_to_png(row["smiles"], out_file):
        print(f"Saved {out_file}")
    else:
        print(f"Failed to render {row['id']}")