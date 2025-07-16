import React, { useState } from "react";
import AttachmentModal from "./AttachmentModal";

export type Compound = {
  id: string;
  name?: string;
  MW?: string; // molecular weight in g/mol
  imageUrl?: string;
};

type Props = {
  compounds: Compound[];
  lots: Record<string, string[]>; // compoundId -> array of lotIds
  onClose: () => void;
  onCreate: (formulationData: any) => void;
};

export default function CreateFormulationModal({ compounds, lots, onClose, onCreate }: Props) {
  const [formulationName, setFormulationName] = useState("");
  const [formulationOperator, setFormulationOperator] = useState("");
  const [selectedComponents, setSelectedComponents] = useState<
    { compoundId: string; imageUrl?: string; lotId: string; molPercent: number; massPercent: number; actualMass?: number }[]
  >([]);

  const [totalMass, setTotalMass] = useState<number>(1);
  const [phaseMap, setPhaseMap] = useState("");
  const [notes, setNotes] = useState("");
  const [formData, setFormData] = useState<{ [key: string]: string }>({});
  const [showAddTextField, setShowAddTextField] = useState(false);
  const [newTextFieldName, setNewTextFieldName] = useState("");
  const [newTextFieldValue, setNewTextFieldValue] = useState("");
  // Change default attachments to thermal data and stability data
  const [attachments, setAttachments] = useState<{ [key: string]: { note: string; imageUrl: string } }>(
    {
      thermal_data: { note: "", imageUrl: "" },
      stability_data: { note: "", imageUrl: "" }
    }
  );
  const [showAddAttachmentField, setShowAddAttachmentField] = useState(false);
  const [newAttachmentFieldName, setNewAttachmentFieldName] = useState("");
  const [selectedAttachment, setSelectedAttachment] = useState<string | null>(null);
  const [compoundSearches, setCompoundSearches] = useState<string[]>([]); // Use a separate search state for each row

  const addComponent = () => {
    if (!compounds.length) return;
    setSelectedComponents((prev) => [
      ...prev,
      {
        compoundId: "",  // Don't default to compounds[0].id
        lotId: "",
        molPercent: 50,
        massPercent: 50,
        actualMass: undefined,
      },
    ]);
    setCompoundSearches((prev) => [...prev, ""]);
  };

  // Update updateComponent to accept undefined
  const updateComponent = (index: number, field: string, value: string | number | undefined) => {
    setSelectedComponents((prev) => {
      const updated = [...prev];
      updated[index] = { ...updated[index], [field]: value };
      return updated;
    });
  };

  const removeComponent = (index: number) => {
    setSelectedComponents(prev => prev.filter((_, i) => i !== index));
    setCompoundSearches(prev => prev.filter((_, i) => i !== index));
  };

  const calculateMasses = () => {
    const totalActualMass = selectedComponents.reduce((sum, c) => sum + (c.actualMass || 0), 0);
    const totalActualMoles = selectedComponents.reduce((sum, c) => {
      const compound = compounds.find(x => x.id === c.compoundId);
      const MW = parseFloat(compound?.MW || "0");
      return sum + (c.actualMass && MW > 0 ? c.actualMass / MW : 0);
    }, 0);

    // Calculate total target moles for all components (for outputMolPercent)
    const totalTargetMoles = selectedComponents.reduce((sum, c) => {
      const compound = compounds.find(x => x.id === c.compoundId);
      const MW = parseFloat(compound?.MW || "0");
      const targetMass = (totalMass * (c.massPercent || 0)) / 100;
      return sum + (MW > 0 ? targetMass / MW : 0);
    }, 0);

    return selectedComponents.map((comp) => {
      const compound = compounds.find((c) => c.id === comp.compoundId);
      const MW = parseFloat(compound?.MW || "0");
      const targetMass = (totalMass * comp.massPercent) / 100;
      const targetMoles = MW > 0 ? targetMass / MW : 0;
      const outputMolPercent = totalTargetMoles > 0 ? (targetMoles / totalTargetMoles) * 100 : 0;

      const actualMass = comp.actualMass || 0;
      const actualMassPercent = totalActualMass > 0 ? (actualMass / totalActualMass) * 100 : 0;
      const actualMoles = MW > 0 ? actualMass / MW : 0;
      const actualMolPercent = totalActualMoles > 0 ? (actualMoles / totalActualMoles) * 100 : 0;

      return {
        ...comp,
        compoundId: compound?.id || comp.compoundId,
        compoundName: compound?.name || '',
        MW,
        targetMass: parseFloat(targetMass.toFixed(3)),
        targetMoles: parseFloat(targetMoles.toFixed(4)),
        outputMolPercent: parseFloat(outputMolPercent.toFixed(2)),
        actualMass,
        actualMassPercent: parseFloat(actualMassPercent.toFixed(2)),
        actualMolPercent: parseFloat(actualMolPercent.toFixed(2)),
      };
    });
  };


  const handleSubmit = async () => {
    const formulation = {
      ...formData,
      name: formulationName,
      operator: formulationOperator,
      components: calculateMasses(),
      totalMass,
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
      style={{backdropFilter: 'blur(2px)'}}
    >
      <div
        className="bg-white border-2 border-[#008080] rounded-2xl p-10 w-full max-w-5xl max-h-[98vh] overflow-y-auto shadow-2xl relative"
        onClick={(e) => e.stopPropagation()}
      >
        <h2 className="text-4xl font-extrabold text-[#002C36] uppercase tracking-wide mb-8 border-b-2 border-[#008080] pb-4">Create Formulation</h2>
        <div className="mb-4">
          <label className="text-xs font-bold uppercase text-[#008080] mb-1 tracking-wide">Formulation Name</label>
          <input
            type="text"
            className="border border-[#008080] rounded px-2 py-1 w-full text-[#002C36] bg-white focus:ring-2 focus:ring-[#008080]"
            value={formulationName}
            onChange={(e) => setFormulationName(e.target.value)}
          />
        </div>
        {/* Operator Dropdown */}
        <div className="mb-4">
          <label className="text-xs font-bold uppercase text-[#008080] mb-1 tracking-wide">Operator</label>
          <input
            type="text"
            className="border border-[#008080] rounded px-2 py-1 w-full text-[#002C36] bg-white focus:ring-2 focus:ring-[#008080]"
            value={formulationOperator}
            onChange={(e) => setFormulationOperator(e.target.value)}
          />
        </div>
        {/* Total Desired Mass (mg) */}
        <div className="mb-4">
          <label className="text-xs font-bold uppercase text-[#008080] mb-1 tracking-wide">Total Desired Mass (mg)</label>
          <input
            type="number"
            step="0.01"
            className="border border-[#008080] rounded px-2 py-1 w-full text-[#002C36] bg-white focus:ring-2 focus:ring-[#008080]"
            value={totalMass}
            onChange={e => setTotalMass(parseFloat(e.target.value))}
          />
        </div>
        {/* Mass Calculator Table */}
        <div className="mb-8 border border-[#008080] rounded-lg p-4 bg-[#F8FAFB]">
          <h3 className="text-lg font-bold text-[#008080] mb-4">Mass Calculator</h3>
          <div className="overflow-x-auto">
            <table className="min-w-full border border-[#008080]">
              <thead>
                <tr className="text-xs text-[#008080] uppercase text-center">
                  <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Compound</th>
                  <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Desired Mass %</th>
                  <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Output Mol %</th>
                  <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Target Mass (mg)</th>
                  <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Actual Mass (mg)</th>
                  <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Actual Mass %</th>
                  <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Actual Mol %</th>
                  <th className="px-2 py-2 border-b border-[#008080]"></th>
                </tr>
              </thead>
              <tbody>
                {selectedComponents.map((comp, idx) => {
                  const compoundSearch = compoundSearches[idx] || "";
                  const filteredCompounds = compounds.filter(c =>
                    c.name?.toLowerCase().includes(compoundSearch.toLowerCase()) ||
                    c.id.toLowerCase().includes(compoundSearch.toLowerCase())
                  );
                  const compound = compounds.find(c => c.id === comp.compoundId);
                  const MW = parseFloat(compound?.MW || "0");
                  const targetMass = (totalMass * (comp.massPercent || 0)) / 100;
                  const targetMoles = MW > 0 ? targetMass / MW : 0;
                  // Calculate total target moles for all components
                  const totalTargetMoles = selectedComponents.reduce((sum, c) => {
                    const cCompound = compounds.find(x => x.id === c.compoundId);
                    const cMW = parseFloat(cCompound?.MW || "0");
                    const cTargetMass = (totalMass * (c.massPercent || 0)) / 100;
                    return sum + (cMW > 0 ? cTargetMass / cMW : 0);
                  }, 0);
                  const outputMolPercent = totalTargetMoles > 0 ? (targetMoles / totalTargetMoles) * 100 : 0;

                  // Actuals
                  const actualMass = comp.actualMass || 0;
                  const totalActualMass = selectedComponents.reduce((sum, c) => sum + (c.actualMass || 0), 0);
                  const actualMassPercent = totalActualMass > 0 ? (actualMass / totalActualMass) * 100 : 0;
                  const actualMoles = MW > 0 ? actualMass / MW : 0;
                  const totalActualMoles = selectedComponents.reduce((sum, c) => {
                    const cCompound = compounds.find(x => x.id === c.compoundId);
                    const cMW = parseFloat(cCompound?.MW || "0");
                    return sum + (c.actualMass && cMW > 0 ? c.actualMass / cMW : 0);
                  }, 0);
                  const actualMolPercent = totalActualMoles > 0 ? (actualMoles / totalActualMoles) * 100 : 0;

                  return (
                    <tr key={comp.compoundId + '-' + comp.lotId + '-' + idx} className="align-middle text-center border-b border-[#008080]">
                      {/* Compound selector */}
                      <td className="px-2 py-2 border-r border-[#008080]">
                        <div className="flex flex-col gap-1">
                          <select
                            className="border border-[#008080] px-1 py-0.5 rounded text-[#002C36] bg-white text-xs"
                            value={comp.compoundId}
                            onChange={(e) => {
                              const selectedId = e.target.value;
                              const selectedCompound = compounds.find(c => c.id === selectedId);
                              updateComponent(idx, 'compoundId', selectedId);
                              updateComponent(idx, 'imageUrl', selectedCompound?.imageUrl || '');
                            }}
                          >
                            <option value="" disabled>Select a compound...</option>
                            {compounds
                              .filter(c =>
                                c.name?.toLowerCase().includes(compoundSearches[idx]?.toLowerCase() || '') ||
                                c.id.toLowerCase().includes(compoundSearches[idx]?.toLowerCase() || '')
                              )
                              .map(c => (
                                <option key={c.id} value={c.id}>
                                  {c.name || c.id}
                                </option>
                              ))}
                          </select>
                          <input
                            type="text"
                            className="border border-[#008080] px-1 py-0.5 rounded text-[#002C36] bg-white text-xs"
                            placeholder="Search..."
                            value={compoundSearches[idx] || ""}
                            onChange={(e) => {
                              const search = e.target.value;
                              setCompoundSearches(prev => {
                                const copy = [...prev];
                                copy[idx] = search;
                                return copy;
                              });
                            }}
                          />
                          <select
                            className="border border-[#008080] px-1 py-0.5 rounded text-[#002C36] bg-white text-xs"
                            value={comp.lotId}
                            onChange={e => updateComponent(idx, 'lotId', e.target.value)}
                          >
                            <option key="original" value="">(Original)</option>
                            {(lots[comp.compoundId.toLowerCase()] || []).map((lot, lotIdx) => (
                              <option key={lot + '-' + lotIdx} value={lot}>{lot}</option>
                            ))}
                          </select>
                        </div>
                      </td>
                      {/* Desired Mass % */}
                      <td className="px-2 py-2 border-r border-[#008080]">
                        <input
                          type="number"
                          min={0}
                          max={100}
                          step={0.01}
                          className="w-20 border border-[#008080] px-2 py-1 rounded text-[#002C36] bg-white text-xs"
                          value={comp.massPercent === undefined ? '' : comp.massPercent}
                          onChange={e => {
                            const val = e.target.value;
                            updateComponent(idx, 'massPercent', val === '' ? undefined : parseFloat(val));
                          }}
                          placeholder="Mass %"
                        />
                      </td>
                      {/* Output Mol % */}
                      <td className="px-2 py-2 border-r border-[#008080] text-black">
                        <span>{outputMolPercent ? outputMolPercent.toFixed(2) : "-"}%</span>
                      </td>
                      {/* Target Mass */}
                      <td className="px-2 py-2 border-r border-[#008080] text-black">
                        <span>{targetMass ? targetMass.toFixed(2) : "-"} mg</span>
                      </td>
                      {/* Actual Mass (input) */}
                      <td className="px-2 py-2 border-r border-[#008080]">
                        <input
                          type="number"
                          min={0}
                          step={0.01}
                          className="w-20 border border-[#008080] px-2 py-1 rounded text-[#002C36] bg-white text-xs"
                          value={comp.actualMass || ''}
                          onChange={e => updateComponent(idx, 'actualMass', e.target.value === '' ? undefined : parseFloat(e.target.value))}
                          placeholder="Actual Mass"
                        />
                      </td>
                      {/* Actual Mass % */}
                      <td className="px-2 py-2 border-r border-[#008080] text-black">
                        <span>{actualMassPercent ? actualMassPercent.toFixed(2) : "-"}%</span>
                      </td>
                      {/* Actual Mol % */}
                      <td className="px-2 py-2 border-r border-[#008080] text-black">
                        <span>{actualMolPercent ? actualMolPercent.toFixed(2) : "-"}%</span>
                      </td>
                      {/* Remove button */}
                      <td className="px-2 py-2">
                        <button
                          className="text-[#008080] hover:text-red-600 font-bold text-lg px-2 py-0.5 rounded"
                          title="Remove component"
                          onClick={() => removeComponent(idx)}
                          type="button"
                        >
                          &times;
                        </button>
                      </td>
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>
          <button
            className="mt-4 bg-[#008080] hover:bg-[#006666] text-white font-bold px-4 py-2 rounded-lg uppercase tracking-wide shadow transition-all"
            onClick={addComponent}
          >
            + Add Component
          </button>
        </div>
        <div className="mb-4">
          <label className="text-xs font-bold uppercase text-[#008080] mb-1 tracking-wide">Phase Map</label>
          <textarea
            rows={3}
            className="w-full border border-[#008080] rounded px-2 py-1 text-[#002C36] bg-white focus:ring-2 focus:ring-[#008080]"
            value={phaseMap}
            onChange={(e) => setPhaseMap(e.target.value)}
          />
        </div>
        <div className="mb-4">
          <label className="text-xs font-bold uppercase text-[#008080] mb-1 tracking-wide">Analytical Notes</label>
          <textarea
            rows={3}
            className="w-full border border-[#008080] rounded px-2 py-1 text-[#002C36] bg-white focus:ring-2 focus:ring-[#008080]"
            value={notes}
            onChange={(e) => setNotes(e.target.value)}
          />
        </div>
        <div className="grid grid-cols-1 sm:grid-cols-2 gap-4 text-[#002C36] mb-6">
          {Object.entries(formData).map(([key, value]) => (
            <div key={key} className="flex flex-col">
              <span className="text-xs font-bold uppercase text-[#008080] mb-1 tracking-wide">{key}</span>
              <input
                type="text"
                value={value}
                onChange={(e) => setFormData(prev => ({ ...prev, [key]: e.target.value }))}
                className="border border-[#008080] rounded px-2 py-1 text-sm text-[#002C36] bg-white focus:ring-2 focus:ring-[#008080]"
              />
            </div>
          ))}
        </div>
        <div className="mb-6 flex gap-2 flex-wrap">
          {Object.entries(attachments).map(([key, { note }]) => (
            <button
              key={key}
              className="border border-[#008080] text-[#008080] font-bold px-2 py-0.5 rounded-md uppercase tracking-wide text-xs hover:bg-[#008080] hover:text-white transition-all mb-1"
              type="button"
              onClick={() => setSelectedAttachment(key)}
            >
              {key.replace(/_/g, " ")}
            </button>
          ))}
        </div>
        <div className="mb-6 flex gap-2 flex-wrap">
          <button
            className="px-2 py-0.5 bg-[#008080] text-white rounded-md hover:bg-[#006666] font-bold uppercase tracking-wide text-xs"
            onClick={() => setShowAddTextField(true)}
          >
            + Add Text Field
          </button>
          <button
            className="px-2 py-0.5 bg-[#008080] text-white rounded-md hover:bg-[#006666] font-bold uppercase tracking-wide text-xs"
            onClick={() => setShowAddAttachmentField(true)}
          >
            + Add Attachment Field
          </button>
        </div>
        <div className="flex justify-end gap-2 mt-6">
          {/* <button className="px-4 py-2 rounded-md bg-gray-200 hover:bg-gray-300 text-[#002C36] font-bold uppercase tracking-wide text-xs" onClick={onClose}>Cancel</button> */}
          <button
            onClick={handleSubmit}
            className="px-5 py-3 rounded-md bg-[#008080] text-white font-bold uppercase tracking-wide text-md hover:bg-[#006666]"
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

/*
  Copyright © 2025 Polaris Electro Optics
  This code is the property of Polaris Electro Optics and may not be reused,
  modified, or distributed without explicit permission.
*/