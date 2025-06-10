# import os
# import subprocess
# import pandas as pd

# image_folder = "compound_test"
# output_csv = "compound_test_output.csv"

# results = []

# for filename in os.listdir(image_folder):
#     if filename.lower().endswith((".png", ".jpg", ".jpeg")):
#         image_path = os.path.join(image_folder, filename)
#         try:
#             result = subprocess.run(["osra", image_path], capture_output=True, text=True)
#             smiles = result.stdout.strip()
#             print(f"[✔] {filename} → {smiles}")
#             results.append({"filename": filename, "smiles": smiles})
#         except Exception as e:
#             print(f"[✘] Failed on {filename}: {e}")
#             results.append({"filename": filename, "smiles": "ERROR"})

# df = pd.DataFrame(results)
# df.to_csv(output_csv, index=False)
# print(f"\n✅ SMILES saved to: {output_csv}")



import os
import pandas as pd
from rdkit import Chem

# === CONFIGURATION ===
mol_folder = "mol_files" 
output_csv = "smiles_from_mol.csv"

results = []

# Process each .mol file
for filename in os.listdir(mol_folder):
    if filename.endswith(".mol"):
        filepath = os.path.join(mol_folder, filename)
        try:
            mol = Chem.MolFromMolFile(filepath)
            if mol:
                smiles = Chem.MolToSmiles(mol)
                print(f" {filename} → {smiles}")
                results.append({"filename": filename, "smiles": smiles})
            else:
                print(f" Invalid molecule: {filename}")
                results.append({"filename": filename, "smiles": "ERROR"})
        except Exception as e:
            print(f" Error loading {filename}: {e}")
            results.append({"filename": filename, "smiles": "ERROR"})

# Save results to CSV
df = pd.DataFrame(results)
df.to_csv(output_csv, index=False)
print(f"\n SMILES saved to: {output_csv}")

