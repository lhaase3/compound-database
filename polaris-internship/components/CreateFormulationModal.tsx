import React, { useState } from "react";
import AttachmentModal from "./AttachmentModal";

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
  const [formData, setFormData] = useState<{ [key: string]: string }>({});
  const [showAddTextField, setShowAddTextField] = useState(false);
  const [newTextFieldName, setNewTextFieldName] = useState("");
  const [newTextFieldValue, setNewTextFieldValue] = useState("");
  const [attachments, setAttachments] = useState<{ [key: string]: { note: string; imageUrl: string } }>({
    data: { note: "", imageUrl: "" }
  });
  const [showAddAttachmentField, setShowAddAttachmentField] = useState(false);
  const [newAttachmentFieldName, setNewAttachmentFieldName] = useState("");
  const [selectedAttachment, setSelectedAttachment] = useState<string | null>(null);
  




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
    ...formData,
    name: formulationName,
    components: calculateMasses(),
    totalMoles,
    phaseMap,
    notes,
    attachments,  // Include attachments
    createdAt: new Date().toISOString(),
    imageUrls: [] as string[],
  };

  try {
    await onCreate(formulation);  // ✅ Send the formulation to parent
    onClose();                    // ✅ Close the modal
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
                <option key="original" value="">(Original)</option>
                {(lots[comp.compoundId.toLowerCase()] || []).map((lot, lotIdx) => (
                  <option key={lot + '-' + lotIdx} value={lot}>{lot}</option>
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

        <div className="grid grid-cols-1 sm:grid-cols-2 gap-4 text-black mb-6">
          {Object.entries(formData).map(([key, value]) => (
            <div key={key} className="flex flex-col">
              <span className="text-sm font-semibold text-black">{key}</span>
              <input
                type="text"
                value={value}
                onChange={(e) => setFormData(prev => ({ ...prev, [key]: e.target.value }))}
                className="border rounded px-2 py-1 text-sm text-black"
              />
            </div>
          ))}
        </div>


        <div className="mb-6 flex gap-4 flex-wrap">
          {Object.keys(attachments).map((key) => (
            <button
              key={key}
              className="text-blue-600 underline hover:text-blue-800"
              type="button"
              onClick={() => setSelectedAttachment(key)}
            >
              {key.replace(/_/g, " ")}
            </button>
          ))}
        </div>

        <div className="mb-6 flex gap-4">
          <button
            className="px-4 py-2 bg-blue-500 text-white rounded hover:bg-blue-700"
            onClick={() => setShowAddTextField(true)}
          >
            ➕ Add Text Field
          </button>
          <button
            className="px-4 py-2 bg-blue-500 text-white rounded hover:bg-blue-700"
            onClick={() => setShowAddAttachmentField(true)}
          >
            ➕ Add Attachment Field
          </button>
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

        {showAddTextField && (
          <div className="fixed inset-0 bg-black/50 z-50 flex items-center justify-center p-6" onClick={() => setShowAddTextField(false)}>
            <div className="bg-white rounded-lg shadow-xl w-full max-w-md p-8" onClick={e => e.stopPropagation()}>
              <h2 className="text-xl font-bold mb-4">Add Text Field</h2>
              <input
                type="text"
                className="w-full border border-gray-300 rounded p-2 mb-4 text-black"
                placeholder="Field name"
                value={newTextFieldName}
                onChange={e => setNewTextFieldName(e.target.value)}
              />
              <textarea
                className="w-full border border-gray-300 rounded p-2 mb-4 text-black"
                placeholder="Field value"
                value={newTextFieldValue}
                onChange={e => setNewTextFieldValue(e.target.value)}
              />
              <div className="flex justify-end gap-2">
                <button
                  className="px-4 py-2 bg-gray-300 rounded hover:bg-gray-400 text-black"
                  onClick={() => setShowAddTextField(false)}
                >
                  Cancel
                </button>
                <button
                  className="px-4 py-2 bg-blue-600 text-white rounded hover:bg-blue-700"
                  onClick={() => {
                    if (!newTextFieldName.trim()) return;
                    setFormData((prev) => ({ ...prev, [newTextFieldName]: newTextFieldValue }));
                    setShowAddTextField(false);
                    setNewTextFieldName("");
                    setNewTextFieldValue("");
                  }}
                >
                  Save
                </button>
              </div>
            </div>
          </div>
        )}
        {showAddAttachmentField && (
          <div className="fixed inset-0 bg-black/50 z-50 flex items-center justify-center p-6" onClick={() => setShowAddAttachmentField(false)}>
            <div className="bg-white rounded-lg shadow-xl w-full max-w-md p-8" onClick={e => e.stopPropagation()}>
              <h2 className="text-xl font-bold mb-4">Add Attachment Field</h2>
              <input
                type="text"
                className="w-full border border-gray-300 rounded p-2 mb-4 text-black"
                placeholder="Attachment field name"
                value={newAttachmentFieldName}
                onChange={e => setNewAttachmentFieldName(e.target.value)}
              />
              <div className="flex justify-end gap-2">
                <button
                  className="px-4 py-2 bg-gray-300 rounded hover:bg-gray-400 text-black"
                  onClick={() => setShowAddAttachmentField(false)}
                >
                  Cancel
                </button>
                <button
                  className="px-4 py-2 bg-blue-600 text-white rounded hover:bg-blue-700"
                  onClick={() => {
                    if (!newAttachmentFieldName.trim()) return;
                    setAttachments((prev) => ({
                      ...prev,
                      [newAttachmentFieldName]: { note: "", imageUrl: "" },
                    }));
                    setShowAddAttachmentField(false);
                    setNewAttachmentFieldName("");
                  }}
                >
                  Save
                </button>
              </div>
            </div>
          </div>
        )}
        {selectedAttachment && (
          <AttachmentModal
            attachmentKey={selectedAttachment}
            data={attachments[selectedAttachment]}
            onClose={() => setSelectedAttachment(null)}
            onSave={(note, fileUrl) => {
              setAttachments((prev) => ({
                ...prev,
                [selectedAttachment]: { note, imageUrl: fileUrl }
              }));
              setSelectedAttachment(null);
            }}
          />
        )}


      </div>
    </div>
  );
}