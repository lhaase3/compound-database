"use client";
import { useEffect, useRef, useState } from "react";

export default function DrawModal({ onClose, onSmilesSubmit }) {
  const [smiles, setSmiles] = useState("");
  const jsmeInitialized = useRef(false); // <-- 🧠 key to prevent duplicate init

  useEffect(() => {
    const container = document.getElementById("jsme_container");
    if (container) container.innerHTML = ""; // Clear just in case

    const initializeJsme = () => {
      if (window.JSApplet?.JSME && !jsmeInitialized.current) {
        // Init only once
        window.jsmeAppletInstance = new window.JSApplet.JSME(
          "jsme_container",
          "700px",
          "450px",
          { options: "oldlook,star" }
        );
        jsmeInitialized.current = true;
      } else {
        setTimeout(initializeJsme, 100);
      }
    };

    const loadJsmeScript = () => {
      const scriptId = "jsme_script";
      const existingScript = document.getElementById(scriptId);

      if (!existingScript) {
        const script = document.createElement("script");
        script.src = "/jsme/jsme.nocache.js";
        script.id = scriptId;
        script.onload = initializeJsme;
        script.onerror = () => console.error("❌ Failed to load JSME script");
        document.body.appendChild(script);
      } else {
        initializeJsme();
      }
    };

    loadJsmeScript();

    return () => {
      const container = document.getElementById("jsme_container");
      if (container) container.innerHTML = "";

      if (window.jsmeAppletList?.length) {
        // 💥 DESTROY applet completely
        window.jsmeAppletList.length = 0;
      }

      window.jsmeAppletInstance = null;
      jsmeInitialized.current = false;
    };
  }, []);

  const handleFilter = () => {
    const jsme = window.jsmeAppletInstance || window.jsmeAppletList?.[0];
    if (jsme) {
      const drawnSmiles = jsme.smiles();
      setSmiles(drawnSmiles);
      onSmilesSubmit(drawnSmiles);
      onClose();
    }
  };

return (
  <div
    className="fixed inset-0 z-50 bg-black/30 backdrop-blur-sm flex items-center justify-center"
    onClick={onClose} // Clicking background closes modal
  >
    <div
      className="bg-white p-6 rounded-2xl shadow-xl relative w-[750px] h-[600px] flex flex-col"
      onClick={(e) => e.stopPropagation()} // Prevent closing on inner click
    >
      <button
        onClick={onClose}
        className="absolute top-2 right-2 text-black text-xl font-bold"
      >
        &times;
      </button>
      <h2 className="text-2xl font-semibold mb-4 text-black">Draw Molecule</h2>
      <div id="jsme_container" className="flex-1 border border-gray-300 rounded-md" />
      <button
        onClick={handleFilter}
        className="mt-4 self-end px-4 py-2 cursor-pointer bg-blue-600 text-white rounded-md hover:bg-blue-700 transition"
      >
        Filter
      </button>
    </div>
  </div>
);

}

/*
  Copyright © 2025 Polaris Electro Optics
  This code is the property of Polaris Electro Optics and may not be reused,
  modified, or distributed without explicit permission.
*/

