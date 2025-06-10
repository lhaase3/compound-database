"use client";

import { useRef, useState, useEffect } from "react";

export default function AddStructureModal({ onClose, onSubmit }) {
  const jsmeInitialized = useRef(false);
  const [formData, setFormData] = useState({});

  const fields = [
    "id",
    "MW",
    "Lambda Max (DCM/Ac CN)",
    "Lambda Max (neat film)",
    "phase map",
    "r33",
    "CAMB3LYP SVPD CHCl3 (Cosmo)",
    "B3LYP SVPD CHCl3 dipole",
    "B3LYP SVPD CHCl3 beta",
    "beta/MW",
    "J/g DSC melt (total)",
    "kJ/mol DSC melt (total)",
    "Refractive index (ne/no)",
    "Notes",
    "lab?",
    "first PEO#",
    "registered PEO#",
    "Lab book #",
    "Max loading (%)"
  ];

    useEffect(() => {
    const container = document.getElementById("jsme_container");
    if (container) container.innerHTML = "";

    const init = () => {
        if (window.JSApplet?.JSME && !jsmeInitialized.current) {
        window.jsmeAppletInstance = new window.JSApplet.JSME(
            "jsme_container",
            "700px",
            "450px",
            { options: "oldlook,star" }
        );
        jsmeInitialized.current = true;
        } else {
        setTimeout(init, 100);
        }
    };

    const loadScript = () => {
        if (!document.getElementById("jsme_script")) {
        const script = document.createElement("script");
        script.src = "/jsme/jsme.nocache.js";
        script.id = "jsme_script";
        script.onload = init;
        document.body.appendChild(script);
        } else {
        init();
        }
    };

    loadScript();

    return () => {
        // 🧼 CLEANUP JSME instance
        const container = document.getElementById("jsme_container");
        if (container) container.innerHTML = "";
        window.jsmeAppletInstance = null;
        if (window.jsmeAppletList?.length) {
        window.jsmeAppletList.length = 0;
        }
        jsmeInitialized.current = false;
    };
    }, []);


  const handleChange = (e) => {
    setFormData({ ...formData, [e.target.name]: e.target.value });
  };

  const handleSave = async () => {
    const smiles = window.jsmeAppletInstance?.smiles?.() || "";
    const data = { ...formData, smiles };

    try {
      const res = await fetch("http://localhost:5000/add-compound", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify(data),
      });

      if (!res.ok) {
        const errData = await res.json();
        throw new Error(errData.error || "Failed to add compound");
      }

      console.log("✅ Compound added successfully");
      onClose();
    } catch (err) {
      console.error("❌ Error:", err.message);
      alert("Failed to save compound: " + err.message);
    }
  };


  return (
    <div
      className="fixed inset-0 z-50 bg-black/40 backdrop-blur-sm flex items-center justify-center"
      onClick={onClose}
    >
      <div
        className="bg-white p-6 rounded-xl w-[800px] max-h-[90vh] overflow-y-auto shadow-lg"
        onClick={(e) => e.stopPropagation()}
      >
        <h2 className="text-2xl font-bold mb-4 text-black">Add New Structure</h2>
        <div id="jsme_container" className="mb-6 border rounded-md" />
        <div className="grid grid-cols-1 sm:grid-cols-2 gap-4 mb-6">
          {fields.map((field) => (
            <div key={field} className="flex flex-col">
              <label className="text-sm text-gray-600 font-semibold">{field}</label>
              <input
                type="text"
                name={field}
                value={formData[field] || ""}
                onChange={handleChange}
                className="border px-2 py-1 rounded text-black"
              />
            </div>
          ))}
        </div>
        <button
          onClick={handleSave}
          className="bg-green-600 hover:bg-green-700 text-white px-4 py-2 cursor-pointer rounded"
        >
          ✅ Save Compound
        </button>
      </div>
    </div>
  );
}
