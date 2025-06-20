import React, { useEffect, useState } from "react";
import Link from "next/link";
import CompoundModal from "./CompoundModal";


export default function FormulationList() {
  const [formulations, setFormulations] = useState<any[]>([]);
  const [selectedFormulation, setSelectedFormulation] = useState<any | null>(null);

  useEffect(() => {
    fetch("http://localhost:5000/formulations")
      .then((res) => res.json())
      .then(setFormulations)
      .catch((err) => console.error("Failed to fetch formulations", err));
  }, []);

  return (
    <div className="p-6 text-black">
      <h1 className="text-3xl font-bold mb-4">Formulations</h1>
      {!selectedFormulation ? (
        <ul className="space-y-2">
          {formulations.map((form) => (
            <li
              key={form.id}
              className="border p-4 rounded cursor-pointer hover:bg-gray-100"
              onClick={() => setSelectedFormulation(form)}
            >
              <span className="font-semibold">{form.name || "Unnamed Formulation"}</span>
            </li>
          ))}
        </ul>
      ) : (
        <div className="border p-4 rounded">
          <button
            className="text-blue-600 underline mb-4"
            onClick={() => setSelectedFormulation(null)}
          >
            ← Back to List
          </button>

          <h2 className="text-2xl font-semibold mb-2">{selectedFormulation.name || "Unnamed"}</h2>
          <p className="text-sm text-gray-600 mb-2">
            Created: {new Date(selectedFormulation.createdAt?.seconds * 1000).toLocaleString()}
          </p>
          {selectedFormulation.imageUrls?.length > 0 && (
            <div className="mb-4">
              <h3 className="font-semibold">Images:</h3>
              <div className="flex flex-wrap gap-4 mt-2">
                {selectedFormulation.imageUrls.map((url: string, idx: number) => (
                  <img
                    key={idx}
                    src={url}
                    alt={`Formulation Image ${idx + 1}`}
                    className="w-40 h-40 object-contain border rounded shadow"
                  />
                ))}
              </div>
            </div>
          )}


          <h3 className="font-semibold mt-4">Components:</h3>
            <ul className="mb-4 list-disc pl-5">
              {selectedFormulation.components.map((comp: any, idx: number) => (
                <li key={idx}>
                  {comp.compoundId} ({comp.lotId || "original"}) – {comp.molPercent}% → {comp.mass} g
                </li>
              ))}
            </ul>



          <div className="mb-2">
            <strong>Total Moles:</strong> {selectedFormulation.totalMoles}
          </div>

          <div className="mb-2">
            <strong>Phase Map:</strong>
            <div className="whitespace-pre-line bg-gray-100 p-2 rounded mt-1">
              {selectedFormulation.phaseMap || "N/A"}
            </div>
          </div>

          <div className="mb-2">
            <strong>Notes:</strong>
            <div className="whitespace-pre-line bg-gray-100 p-2 rounded mt-1">
              {selectedFormulation.notes || "N/A"}
            </div>
          </div>
        </div>
      )}
    </div>
  );
  
}
