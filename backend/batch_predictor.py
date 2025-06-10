import os
import csv
from DECIMER import predict_SMILES

# === CONFIG ===
image_folder = "sing_test"
output_csv = "decimer_216_res.csv"

if not os.path.exists(image_folder):
    print(f"❌ Folder not found: {image_folder}")
    exit(1)

results = []
for filename in os.listdir(image_folder):
    if filename.lower().endswith((".png", ".jpg", ".jpeg")):
        filepath = os.path.join(image_folder, filename)
        print(f"🔍 Predicting: {filename}")

        try:
            smiles = predict_SMILES(filepath)
        except Exception as e:
            print(f"✘ Error with {filename}: {e}")
            smiles = "ERROR"

        results.append((filename, smiles))

with open(output_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["filename", "smiles"])
    writer.writerows(results)

print(f"\n✅ All done. Results saved to: {output_csv}")
