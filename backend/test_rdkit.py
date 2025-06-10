# # from rdkit import Chem
# # from rdkit.Chem import Draw

# # mol = Chem.MolFromSmiles("CC(C)OC1=CC(=CC=C1C(=O)OC2=CC=C(C(=C2)OC)C(=O)OC3=C(C(=C(C(=C3))I=CC4=CC=C(C=C4)[N+](=O)[O-])OC)OC)N=NC5=C(C=C(C=C5OC)N=NC6=C(C=C(C=C6)N(C)C)OC)F")
# # img = Draw.MolToImage(mol)
# # img.show()


# from rdkit import Chem
# from rdkit.Chem import Draw

# smiles = "CCCOc1cc(OC)ccc1C(=O)Oc2ccc(C(=O)Oc3ccc(C(=O)Oc4ccc(C(=O)Oc5ccc(C(=O)Oc6ccc(cc6)[N+](=O)[O-])c(OCCC)c5)c(OCCC)c4)c(OCCC)c3)c(OCCC)c2"
# mol = Chem.MolFromSmiles(smiles)

# # Optional: Compute 2D coordinates
# from rdkit.Chem import AllChem
# AllChem.Compute2DCoords(mol)

# # Create a high-resolution image
# img = Draw.MolToImage(mol, size=(1600, 600))  # Increase resolution here
# img.save("high_quality.png")



import csv
import os
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

# CONFIG
input_csv = "transport_res.csv"
output_folder = "output_images"
image_size = (600, 300)

os.makedirs(output_folder, exist_ok=True)

with open(input_csv, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        filename = row["filename"]
        smiles = row["smiles"]

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"[!] Could not parse SMILES for {filename}")
            continue

        AllChem.Compute2DCoords(mol)
        img = Draw.MolToImage(mol, size=image_size)
        img.save(os.path.join(output_folder, filename))
        print(f"[✓] Saved {filename}")



