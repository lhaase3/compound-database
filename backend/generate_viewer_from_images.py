import csv
import os
import shutil
import base64

# === CONFIG ===
input_csv = "transport_res.csv"
image_folder = "output_images"
html_file = "chemical_structures.html"

# === Generate HTML ===
html = ["<html><head><title>Chemical Structures</title></head><body>"]
html.append("<h1>Chemical Structures</h1><table border='1' cellspacing='5' cellpadding='10'>")
html.append("<tr><th>Image</th><th>Filename</th><th>SMILES</th></tr>")

with open(input_csv, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        filename = row["filename"]
        smiles = row["smiles"]
        image_path = os.path.join(image_folder, filename)

        if os.path.exists(image_path):
            with open(image_path, "rb") as image_file:
                encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
                img_tag = f"<img src='data:image/png;base64,{encoded_string}' width='350'>"
            html.append(f"<tr><td>{img_tag}</td><td>{filename}</td><td>{smiles}</td></tr>")
        else:
            print(f"[!] Missing image for {filename}, skipping.")

html.append("</table></body></html>")

# Save HTML inside the image folder
html_path = os.path.join(image_folder, html_file)
with open(html_path, "w", encoding="utf-8") as f:
    f.write("\n".join(html))

# Create ZIP of the output_images folder
shutil.make_archive("chemical_package_test", 'zip', root_dir='.', base_dir=image_folder)
print("[✓] chemical_package_test.zip created and ready to send.")
