"use client";

import { useEffect, useRef, useState } from "react";

declare global {
  interface Window {
    RDKit: any;
    initRDKitModule: () => Promise<any>;
  }
}

type Props = {
  smiles: string;
};

export default function SMILESRenderer({ smiles }: Props) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [rdkitLoaded, setRdkitLoaded] = useState(false);

  useEffect(() => {
    if (!window.RDKit) {
      window.initRDKitModule()
        .then((RDKit) => {
          window.RDKit = RDKit;
          setRdkitLoaded(true);
        })
        .catch((err) => {
          console.error("Failed to load RDKit:", err);
        });
    } else {
      setRdkitLoaded(true);
    }
  }, []);

  useEffect(() => {
    if (!rdkitLoaded || !canvasRef.current) return;

    try {
      const mol = window.RDKit.get_mol(smiles);
      const svg = mol.get_svg(); // no sizing support in JS API
      const blob = new Blob([svg], { type: "image/svg+xml" });
      const url = URL.createObjectURL(blob);
      const img = new Image();

      img.onload = () => {
        const canvas = canvasRef.current;
        if (!canvas) return;

        const ctx = canvas.getContext("2d");
        if (!ctx) return;
        canvas.width = img.width * 2; // scale up
        canvas.height = img.height * 2;
        ctx.scale(2, 2); // double the drawing size
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        ctx.drawImage(img, 0, 0);
        URL.revokeObjectURL(url);
      };

      img.src = url;
      mol.delete();
    } catch (e) {
      console.error("Failed to render molecule:", e);
    }
  }, [rdkitLoaded, smiles]);

  return (
    <canvas
      ref={canvasRef}
      className="w-full h-auto rounded-lg shadow-md"
      style={{ backgroundColor: "#fff", maxWidth: "100%" }}
    />
  );
}



