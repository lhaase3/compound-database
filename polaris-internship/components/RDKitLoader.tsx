"use client";

import Script from "next/script";

export default function RDKitLoader() {
  return (
    <Script
      src="https://unpkg.com/@rdkit/rdkit@latest/Code/MinimalLib/dist/RDKit_minimal.js"
      strategy="beforeInteractive"
      onLoad={() => console.log("✅ RDKit loaded")}
    />
  );
}
