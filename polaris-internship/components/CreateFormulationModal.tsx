import React, { useState } from "react";

export type Compound = {
  id: string;
  name?: string;
  MW?: string; // molecular weight in g/mol
};

type Props = {
  compounds: Compound[];
  lots: Record<string, string[]>; // compoundId -> array of lotIds
  onClose: () => void;
  onCreate: (formulationData: any) => void;
};

export default function CreateFormulationModal({ compounds, lots, onClose, onCreate }: Props) {
  console.log("🧪 COMPOUND IDS:", compounds.map(c => c.id));
  console.log("🧪 LOTS OBJECT KEYS:", Object.keys(lots));

  const [formulationName, setFormulationName] = useState("");
  const [selectedComponents, setSelectedComponents] = useState<
    { compoundId: string; lotId: string; molPercent: number }[]
  >([]);

  const [totalMoles, setTotalMoles] = useState<number>(1);
  const [phaseMap, setPhaseMap] = useState("");
  const [notes, setNotes] = useState("");
  const [images, setImages] = useState<FileList | null>(null);


  const addComponent = () => {
    if (!compounds.length) return;
    setSelectedComponents((prev) => [
      ...prev,
      { compoundId: compounds[0].id, lotId: lots[compounds[0].id.toLowerCase()]?.[0] || "", molPercent: 50 }
    ]);
  };

  const updateComponent = (index: number, field: string, value: string | number) => {
    setSelectedComponents((prev) => {
      const updated = [...prev];
      updated[index] = { ...updated[index], [field]: value };
      return updated;
    });
  };

  const calculateMasses = () => {
    return selectedComponents.map((comp) => {
      const compound = compounds.find((c) => c.id === comp.compoundId);
      const MW = parseFloat(compound?.MW || "0");
      const moles = (totalMoles * comp.molPercent) / 100;
      return {
        ...comp,
        MW,
        mass: parseFloat((moles * MW).toFixed(3)),
      };
    });
  };

  const handleSubmit = async () => {
    const formulation = {
      name: formulationName,
      components: calculateMasses(),
      totalMoles,
      phaseMap,
      notes,
      createdAt: new Date().toISOString(),
      imageUrls: [] as string[]
    };

    try {
      if (images && images.length > 0) {
        for (const file of Array.from(images)) {
          const formData = new FormData();
          formData.append("file", file);
          formData.append("note", `Formulation Image: ${file.name}`);

          const res = await fetch("http://localhost:5000/upload-image-to-firebase", {
            method: "POST",
            body: formData,
          });

          const result = await res.json();
          if (res.ok && result.fileUrl) {
            formulation.imageUrls.push(result.fileUrl);
          } else {
            console.error("Upload failed:", result.error);
          }
        }
      }

      await onCreate(formulation);
      onClose();
    } catch (err) {
      console.error("Formulation creation failed:", err);
      alert("Failed to save formulation.");
    }
  };


  return (
    <div
      className="fixed inset-0 bg-black/40 z-50 flex items-center justify-center"
      onClick={onClose}
    >
      <div
        className="bg-white rounded-lg p-6 w-full max-w-2xl max-h-[90vh] overflow-y-auto"
        onClick={(e) => e.stopPropagation()}
      >
        <h2 className="text-2xl font-bold text-black mb-4">Create Formulation</h2>
        <div className="mb-4">
            <label className="text-sm font-semibold text-black">Formulation Name</label>
            <input
            type="text"
            className="border rounded w-full px-2 py-1 text-black"
            value={formulationName}
            onChange={(e) => setFormulationName(e.target.value)}
            />
        </div>
        

        <div className="mb-4">
          <label className="text-sm font-semibold text-black">Total Moles</label>
          <input
            type="number"
            step="0.01"
            className="border rounded w-full px-2 py-1 text-black"
            value={totalMoles}
            onChange={(e) => setTotalMoles(parseFloat(e.target.value))}
          />
        </div>
        
        {selectedComponents.map((comp, idx) => (
          <div key={comp.compoundId + '-' + comp.lotId + '-' + idx} className="mb-4 border p-3 rounded">
            <div className="flex gap-4 mb-2 text-black">
              <select
                className="border px-2 py-1"
                value={comp.compoundId}
                onChange={(e) => updateComponent(idx, "compoundId", e.target.value)}
              >
                {compounds.map((c) => (
                  <option key={c.id} value={c.id}>{c.name || c.id}</option>
                ))}
              </select>
              <select
                className="border px-2 py-1"
                value={comp.lotId}
                onChange={(e) => updateComponent(idx, "lotId", e.target.value)}
              >
                <option value="">(Original)</option>
                {(lots[comp.compoundId.toLowerCase()] || []).map((lot) => (
                  <option key={lot} value={lot}>{lot}</option>  // ✅ key added
                ))}
              </select>

              <input
                type="number"
                min={0}
                max={100}
                step={0.1}
                className="w-24 border px-2 py-1"
                value={comp.molPercent}
                onChange={(e) => updateComponent(idx, "molPercent", parseFloat(e.target.value))}
              />
              <span className="self-center">mol%</span>
            </div>

            <div className="text-sm text-gray-700">
              Calculated Mass: <b>{calculateMasses()[idx]?.mass ?? "-"} g</b>
            </div>
          </div>
        ))}

        <button
          className="mb-4 bg-green-600 text-white px-4 py-2 rounded hover:bg-green-700"
          onClick={addComponent}
        >
          ➕ Add Component
        </button>

        <div className="mb-4">
          <label className="text-sm font-semibold text-black">Phase Map</label>
          <textarea
            rows={3}
            className="w-full border rounded px-2 py-1 text-black"
            value={phaseMap}
            onChange={(e) => setPhaseMap(e.target.value)}
          />
        </div>

        <div className="mb-4">
          <label className="text-sm font-semibold text-black">Analytical Notes</label>
          <textarea
            rows={3}
            className="w-full border rounded px-2 py-1 text-black"
            value={notes}
            onChange={(e) => setNotes(e.target.value)}
          />
        </div>

        <div className="mb-6">
          <label className="text-sm font-semibold text-black">Upload Images</label>
            <input
              type="file"
              multiple
              onChange={(e) => setImages(e.target.files)}
              className="block w-full text-sm text-gray-700 border border-gray-300 rounded px-2 py-1 mt-2"
            />
        </div>

        <div className="flex justify-end gap-4">
          <button className="px-4 py-2 rounded bg-gray-300 hover:bg-gray-400 cursor-pointer" onClick={onClose}>Cancel</button>
          <button
            onClick={handleSubmit}
            className="px-4 py-2 rounded bg-blue-600 text-white hover:bg-blue-700 cursor-pointer"
          >
            Save Formulation
          </button>
        </div>
      </div>
    </div>
  );
}