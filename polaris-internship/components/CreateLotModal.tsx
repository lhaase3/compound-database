"use client";
import { useState } from "react";

type Compound = {
  id: string;
  name?: string;
  formula?: string;
  smiles: string;
  [key: string]: string | undefined;
};

type Props = {
  compounds: Compound[];
  onClose: () => void;
  onCreate: () => void;
};

export default function CreateLotModal({ compounds, onClose, onCreate }: Props) {
  const [selected, setSelected] = useState<string[]>([]);
  const [searchTerm, setSearchTerm] = useState("");
  const [lotName, setLotName] = useState(() => {
    if (compounds.length === 1) {
      return `${compounds[0].id}-`;
    }
    return "";
  });


  const toggleSelection = (id: string) => {
    setSelected((prev) =>
      prev.includes(id) ? prev.filter((x) => x !== id) : [...prev, id]
    );
  };

  const handleSubmit = async () => {
    if (selected.length !== 1) {
      alert("Please select one compound.");
      return;
    }

    try {
      const res = await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/create-lot`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ compoundIds: selected, lotName }),
      });

      const result = await res.json();
      if (!res.ok) {
        alert("Error: " + result.error);
        return;
      }

      alert(`Lot created: ${result.lotName}`);
      onCreate();
      onClose();
    } catch (error) {
      console.error("Lot creation failed", error);
    }
  };

  const filteredCompounds = compounds.filter((compound) => {
    const query = searchTerm.toLowerCase();
    return (
      (compound.id ?? "").toLowerCase().includes(query) ||
      (compound.name ?? "").toLowerCase().includes(query) ||
      (compound.formula ?? "").toLowerCase().includes(query)
    );
  });


  return (
    <div
      className="fixed inset-0 backdrop-blur-sm bg-opacity-50 z-50 flex items-center justify-center"
      onClick={onClose}
    >
      <div
        className="bg-white p-6 rounded-lg w-full max-w-lg max-h-[90vh] overflow-y-auto shadow-xl"
        onClick={(e) => e.stopPropagation()}
      >
        <h2 className="text-4xl font-extrabold text-[#002C36] uppercase tracking-wide mb-8 border-b-2 border-[#008080] pb-4">Create New Lot</h2>
        <div className="mb-3">
          <label className="block font-semibold mb-1 text-[#008080]">Lot Name</label>
          <input
            type="text"
            value={lotName}
            onChange={(e) => setLotName(e.target.value)}
            className="w-full p-2 mb-4 border border-[#008080] rounded text-black"
          />
        </div>

        <input
          type="text"
          value={searchTerm}
          onChange={(e) => setSearchTerm(e.target.value)}
          placeholder="Search compounds..."
          className="w-full p-2 mb-3 border border-[#008080] rounded text-black"
        />

        <div className="mb-4 max-h-64 overflow-y-scroll border border-gray-200 rounded p-2">
          {filteredCompounds.length > 0 ? (
            filteredCompounds.map((compound) => (
              <label key={compound.id} className="flex items-center mb-1">
                <input
                  type="checkbox"
                  checked={selected.includes(compound.id)}
                  onChange={() => toggleSelection(compound.id)}
                />
                <span className="ml-2 text-black">
                  {compound.name || compound.id}
                </span>
              </label>
            ))
          ) : (
            <p className="text-gray-500 italic">No matching compounds found.</p>
          )}
        </div>

        <div className="flex justify-end gap-4">
          <button
            onClick={onClose}
            className="px-4 py-2 rounded bg-gray-300 hover:bg-gray-400 text-black font-bold"
          >
            Cancel
          </button>
          <button
            onClick={handleSubmit}
            className="px-4 py-2 rounded bg-[#008080] hover:bg-[#006666] text-white font-bold"
          >
            Create Lot
          </button>
        </div>
      </div>
    </div>
  );
}

/*
  Copyright © 2025 Polaris Electro Optics
  This code is the property of Polaris Electro Optics and may not be reused,
  modified, or distributed without explicit permission.
*/