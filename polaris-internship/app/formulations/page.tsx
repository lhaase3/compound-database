"use client"

import React, { useEffect, useState } from "react";
import Link from "next/link";
import CompoundModal from "@/components/CompoundModal";


export default function FormulationList() {
  const [formulations, setFormulations] = useState<any[]>([]);
  const [selectedFormulation, setSelectedFormulation] = useState<any | null>(null);
  const [editMode, setEditMode] = useState(false);
  const [editData, setEditData] = useState({ name: "", phaseMap: "", notes: "" });
  const [selectedCompound, setSelectedCompound] = useState<any | null>(null);
  const [compoundSource, setCompoundSource] = useState<"main" | "lot">("main");
  const [compoundLotId, setCompoundLotId] = useState<string | null>(null);
  const [searchTerm, setSearchTerm] = useState("");




  useEffect(() => {
    fetch("http://localhost:5000/formulations")
      .then((res) => res.json())
      .then(setFormulations)
      .catch((err) => console.error("Failed to fetch formulations", err));
  }, []);

  return (
    <div className="p-6 bg-white min-h-screen text-black">
    <Link href="/">
        <button className="mb-4 px-4 py-2 bg-gray-300 hover:bg-gray-400 text-black rounded">
            ⬅ Back to Home
        </button>
    </Link>
      <h1 className="text-3xl font-bold mb-4">Formulations</h1>
      <div className="mb-4 w-full max-w-md">
        <input
          type="text"
          placeholder="🔍 Search by compound ID (e.g. PEO-0100)"
          value={searchTerm}
          onChange={(e) => setSearchTerm(e.target.value)}
          className="w-full border border-gray-300 px-3 py-2 rounded text-black"
        />
      </div>
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
        {formulations
          .filter((form) =>
            searchTerm.trim() === "" ||
            form.components?.some((comp: any) =>
              comp.compoundId?.toLowerCase().includes(searchTerm.toLowerCase())
            )
          )
          .map((form) => (
            <div
              key={form.id}
              className="bg-white border border-gray-300 rounded-lg p-4 shadow hover:shadow-md cursor-pointer"
              onClick={() => setSelectedFormulation(form)}
            >
              <h2 className="text-xl font-semibold mb-2">{form.name || "Unnamed Formulation"}</h2>
              <p className="text-sm text-gray-800 mb-1">Components:</p>
              <ul className="list-disc pl-5 text-sm text-gray-700">
                {form.components?.map((comp: any, idx: number) => (
                  <li key={idx}>{comp.compoundId} ({comp.molPercent}%)</li>
                )) || <li>No components</li>}
              </ul>
            </div>
        ))}
      </div>

      {selectedFormulation && (
        <div
          className="fixed inset-0 z-50 flex items-center justify-center bg-opacity-50 backdrop-blur-sm"
          onClick={() => setSelectedFormulation(null)}
        >
          <div
            className="bg-white p-6 rounded-lg shadow-xl max-w-2xl w-full max-h-[90vh] overflow-y-auto"
            onClick={(e) => e.stopPropagation()}
          >
            <div className="flex justify-between items-start mb-4">
            <h2 className="text-2xl font-bold text-black">
                {selectedFormulation.name || "Unnamed Formulation"}
            </h2>
            <div className="flex gap-2">
                <button
                onClick={() => {
                    setEditData({
                    name: selectedFormulation.name || "",
                    phaseMap: selectedFormulation.phaseMap || "",
                    notes: selectedFormulation.notes || "",
                    });
                    setEditMode(true);
                }}
                className="text-blue-600 border border-blue-600 px-3 py-1 text-sm rounded hover:bg-blue-100"
                >
                Edit
                </button>
                <button
                onClick={async () => {
                    try {
                    await fetch(`http://localhost:5000/delete-formulation/${selectedFormulation.id}`, {
                        method: "DELETE",
                    });
                    setFormulations((prev) =>
                        prev.filter((f) => f.id !== selectedFormulation.id)
                    );
                    setSelectedFormulation(null);
                    } catch (err) {
                    console.error("Failed to delete formulation", err);
                    alert("Failed to delete formulation.");
                    }
                }}
                className="text-red-600 border border-red-600 px-3 py-1 text-sm rounded hover:bg-red-100"
                >
                Delete
                </button>
                <button
                onClick={() => setSelectedFormulation(null)}
                className="text-black text-xl font-bold px-2"
                >
                &times;
                </button>
            </div>
            </div>


            <div className="mb-4">
                <p className="text-sm text-gray-600">
                Created: {selectedFormulation.createdAt
                    ? new Date(
                        selectedFormulation.createdAt.seconds
                        ? selectedFormulation.createdAt.seconds * 1000
                        : selectedFormulation.createdAt
                    ).toLocaleDateString()
                    : "Unknown"}
                </p>
            </div>

            <div className="mb-4">
              <h3 className="text-lg font-semibold">Components</h3>
              <ul className="list-disc pl-5">
                {selectedFormulation.components?.map((comp: any, idx: number) => (
                  <li key={idx}>
                    <button
                      className="text-blue-600 underline hover:text-blue-800"
                      onClick={async () => {
                        try {
                          const endpoint = comp.lotId
                            ? `http://localhost:5000/lot/${comp.lotId}`
                            : `http://localhost:5000/compounds/${comp.compoundId}`;

                          const res = await fetch(endpoint);
                          const data = await res.json();
                          let compound;

                          if (comp.lotId) {
                            // Assume each lot only contains one compound (or just grab the matching one)
                            compound = data.find((c: any) =>
                              c.id === comp.compoundId || c.name?.includes(comp.compoundId)
                            ) || data[0]; // fallback to first in lot
                          } else {
                            compound = data;
                          }

                          if (!compound) {
                            alert("Compound not found.");
                            return;
                          }


                          setSelectedCompound(compound);
                          setCompoundSource(comp.lotId ? "lot" : "main");
                          setCompoundLotId(comp.lotId || null);
                        } catch (err) {
                          console.error("Error loading compound:", err);
                        }
                      }}
                    >
                      {comp.compoundId} ({comp.lotId || "original"}) – {comp.molPercent}% → {comp.mass} g
                    </button>
                  </li>
                )) || <li>No components</li>}
              </ul>

            </div>

            <div className="mb-4">
              <h3 className="text-lg font-semibold">Phase Map</h3>
              <div className="bg-gray-100 p-3 rounded whitespace-pre-wrap">
                {selectedFormulation.phaseMap || "N/A"}
              </div>
            </div>

            <div className="mb-4">
              <h3 className="text-lg font-semibold">Analytical Notes</h3>
              <div className="bg-gray-100 p-3 rounded whitespace-pre-wrap">
                {selectedFormulation.notes || "N/A"}
              </div>
            </div>
            {selectedFormulation.imageUrls?.length > 0 && (
              <div className="mb-4">
                <h3 className="text-lg font-semibold">Images</h3>
                <div className="flex flex-wrap gap-4 mt-2">
                  {selectedFormulation.imageUrls.map((url: string, idx: number) => (
                    <img
                      key={idx}
                      src={url}
                      alt={`Image ${idx + 1}`}
                      className="w-160 h-100 object-contain border rounded shadow"
                    />
                  ))}
                </div>
              </div>
            )}


          </div>
        </div>
      )}
      {editMode && (
        <div
            className="fixed inset-0 backdrop-blur-sm bg-opacity-50 z-50 flex items-center justify-center"
            onClick={() => setEditMode(false)}
        >
            <div
            className="bg-white p-6 rounded-lg shadow-lg w-full max-w-xl"
            onClick={(e) => e.stopPropagation()}
            >
            <h2 className="text-xl font-bold mb-4 text-black">Edit Formulation</h2>

            <label className="block text-sm font-semibold text-gray-700 mb-1">Name</label>
            <input
                type="text"
                value={editData.name}
                onChange={(e) => setEditData((d) => ({ ...d, name: e.target.value }))}
                className="w-full mb-4 border px-2 py-1 rounded text-black"
            />

            <label className="block text-sm font-semibold text-gray-700 mb-1">Phase Map</label>
            <textarea
                value={editData.phaseMap}
                onChange={(e) => setEditData((d) => ({ ...d, phaseMap: e.target.value }))}
                className="w-full mb-4 border px-2 py-1 rounded text-black"
            />

            <label className="block text-sm font-semibold text-gray-700 mb-1">Analytical Notes</label>
            <textarea
                value={editData.notes}
                onChange={(e) => setEditData((d) => ({ ...d, notes: e.target.value }))}
                className="w-full mb-4 border px-2 py-1 rounded text-black"
            />

            <div className="flex justify-end gap-4">
                <button
                onClick={() => setEditMode(false)}
                className="px-4 py-2 bg-gray-300 rounded hover:bg-gray-400 text-black"
                >
                Cancel
                </button>
                <button
                onClick={async () => {
                    try {
                    await fetch(`http://localhost:5000/update-formulation/${selectedFormulation.id}`, {
                        method: "POST",
                        headers: { "Content-Type": "application/json" },
                        body: JSON.stringify(editData),
                    });

                    setFormulations((prev) =>
                        prev.map((f: any) =>
                            f.id === selectedFormulation.id ? { ...f, ...editData } : f
                        )
                    );
                    setSelectedFormulation((f) => f && { ...f, ...editData });
                    setEditMode(false);
                    } catch (err) {
                    console.error("Failed to update formulation", err);
                    alert("Update failed");
                    }
                }}
                className="px-4 py-2 bg-blue-600 text-white rounded hover:bg-blue-700"
                >
                Save
                </button>
            </div>
            </div>
        </div>
        )}
        {selectedCompound && (
          <CompoundModal
            compound={selectedCompound}
            onClose={() => setSelectedCompound(null)}
            onDelete={async (id: string) => {
              await fetch(`http://localhost:5000/delete-compound/${id}`, { method: "DELETE" });
              setSelectedCompound(null);
            }}
            onUpdate={async (updatedCompound) => {
              const endpoint =
                compoundSource === "lot"
                  ? "http://localhost:5000/update-lot-compound"
                  : "http://localhost:5000/update-compound";

              const payload =
                compoundSource === "lot"
                  ? { ...updatedCompound, lotId: compoundLotId }
                  : updatedCompound;

              await fetch(endpoint, {
                method: "POST",
                headers: { "Content-Type": "application/json" },
                body: JSON.stringify(payload),
              });

              setSelectedCompound(updatedCompound);
            }}
            source={compoundSource}
            lotId={compoundLotId}
            onUpdateCompoundFromLot={(compound, lotId) => {
              setSelectedCompound(compound);
              setCompoundSource(lotId ? "lot" : "main");
              setCompoundLotId(lotId);
            }}
          />
        )}

    </div>
  );
}




