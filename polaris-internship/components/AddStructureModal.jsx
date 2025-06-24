"use client";

import { useRef, useState, useEffect } from "react";
import AttachmentModal from "./AttachmentModal";

export default function AddStructureModal({ onClose, onSubmit }) {
  const jsmeInitialized = useRef(false);
  const [formData, setFormData] = useState({});
  const [showAddTextField, setShowAddTextField] = useState(false);
  const [newTextFieldName, setNewTextFieldName] = useState("");
  const [newTextFieldValue, setNewTextFieldValue] = useState("");

  const fields = [
    "smiles",
    "id",
    "MW",
    "Lambda Max (DCM/AcCN)",
    "Lambda Max (neat film)",
    "phase map",
    "r33",
    "dipole CAMB3LYP SVPD CHCl3 (Cosmo)",
    "beta CAMB3LYP SVPD CHCl3 (Cosmo)",
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

  const handleAddTextField = () => {
    setShowAddTextField(true);
  };

  const handleSaveTextField = () => {
    if (newTextFieldName.trim()) {
      setFormData((prev) => ({ ...prev, [newTextFieldName]: newTextFieldValue }));
    }
    setShowAddTextField(false);
    setNewTextFieldName("");
    setNewTextFieldValue("");
  };


  const [attachments, setAttachments] = useState({
    uv_vis: { note: "", imageUrl: "" },
    dsc: { note: "", imageUrl: "" },
    lcms: { note: "", imageUrl: "" },
    thermal_stability: { note: "", imageUrl: "" },
    pda_detector_spectrum: { note: "", imageUrl: "" },
  });
  const [showAddAttachmentField, setShowAddAttachmentField] = useState(false);
  const [newAttachmentFieldName, setNewAttachmentFieldName] = useState("");
  const [selectedAttachment, setSelectedAttachment] = useState(null);

  const handleAddAttachmentField = () => {
    setShowAddAttachmentField(true);
  };

  const handleSave = async () => {
    const drawnSMILES = window.jsmeAppletInstance?.smiles?.() || "";
    const smiles = formData.smiles?.trim() || drawnSMILES;
    // Include all fields in fields array and all customFields, even if empty
    const data = {
      ...formData,
      smiles,
      attachments,
    };


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
        className="bg-white border-2 border-[#008080] rounded-2xl p-8 w-[800px] max-h-[92vh] overflow-y-auto shadow-2xl"
        onClick={(e) => e.stopPropagation()}
      >
        <h2 className="text-3xl font-extrabold text-[#002C36] uppercase tracking-wide mb-6 border-b-2 border-[#008080] pb-3">Add New Structure</h2>
        <div id="jsme_container" className="mb-6 border-2 border-[#008080] rounded-lg" />
        <div className="grid grid-cols-1 sm:grid-cols-2 gap-4 text-[#002C36] mb-6">
          {[...fields, ...Object.keys(formData).filter((key) => !fields.includes(key) && key !== "attachments")].map((field) => (
            <div key={field} className="flex flex-col">
              <span className="text-xs font-bold uppercase text-[#008080] mb-1 tracking-wide">{field}</span>
              <input
                type="text"
                name={field}
                value={formData[field] || ""}
                onChange={handleChange}
                className="border border-[#008080] rounded px-2 py-1 text-sm text-[#002C36] bg-white focus:ring-2 focus:ring-[#008080]"
              />
            </div>
          ))}
        </div>
        {/* Attachments Section */}
        <div className="mt-6 flex gap-2 flex-wrap justify-center">
          {Object.keys(attachments).map((key) => (
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
        {/* Add Field Buttons */}
        <div className="mt-8 flex gap-2 justify-center">
          <button
            className="px-2 py-0.5 bg-[#008080] text-white rounded-md hover:bg-[#006666] font-bold uppercase tracking-wide text-xs"
            onClick={handleAddAttachmentField}
            type="button"
          >
            + Add Attachment Field
          </button>
          <button
            className="px-2 py-0.5 bg-[#008080] text-white rounded-md hover:bg-[#006666] font-bold uppercase tracking-wide text-xs"
            onClick={handleAddTextField}
            type="button"
          >
            + Add Text Field
          </button>
        </div>
        {/* Add Text Field Modal */}
        {showAddTextField && (
          <div className="fixed inset-0 bg-black/50 z-50 flex items-center justify-center p-6" onClick={() => setShowAddTextField(false)}>
            <div className="bg-white rounded-2xl border-2 border-[#008080] shadow-2xl w-full max-w-md p-8" onClick={e => e.stopPropagation()}>
              <h2 className="text-2xl font-bold mb-4 text-[#002C36] uppercase tracking-wide">Add Text Field</h2>
              <input
                type="text"
                className="w-full border border-[#008080] rounded p-2 mb-4 text-[#002C36] bg-white focus:ring-2 focus:ring-[#008080]"
                placeholder="Field name"
                value={newTextFieldName}
                onChange={e => setNewTextFieldName(e.target.value)}
              />
              <textarea
                className="w-full border border-[#008080] rounded p-2 mb-4 text-[#002C36] bg-white focus:ring-2 focus:ring-[#008080]"
                placeholder="Field value"
                value={newTextFieldValue}
                onChange={e => setNewTextFieldValue(e.target.value)}
              />
              <div className="flex justify-end gap-2">
                <button
                  className="px-3 py-1 bg-gray-200 rounded-md hover:bg-gray-300 text-[#002C36] font-bold uppercase tracking-wide text-xs"
                  onClick={() => setShowAddTextField(false)}
                >
                  Cancel
                </button>
                <button
                  className="px-3 py-1 bg-[#008080] text-white rounded-md hover:bg-[#006666] font-bold uppercase tracking-wide text-xs"
                  onClick={() => {
                    if (!newTextFieldName.trim()) return;
                    if (!fields.includes(newTextFieldName)) {
                      fields.push(newTextFieldName);
                    }
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
        {/* Add Attachment Field Modal */}
        {showAddAttachmentField && (
          <div className="fixed inset-0 bg-black/50 z-50 flex items-center justify-center p-6" onClick={() => setShowAddAttachmentField(false)}>
            <div className="bg-white rounded-2xl border-2 border-[#008080] shadow-2xl w-full max-w-md p-8" onClick={e => e.stopPropagation()}>
              <h2 className="text-2xl font-bold mb-4 text-[#002C36] uppercase tracking-wide">Add Attachment Field</h2>
              <input
                type="text"
                className="w-full border border-[#008080] rounded p-2 mb-4 text-[#002C36] bg-white focus:ring-2 focus:ring-[#008080]"
                placeholder="Attachment field name"
                value={newAttachmentFieldName}
                onChange={e => setNewAttachmentFieldName(e.target.value)}
              />
              <div className="flex justify-end gap-2">
                <button
                  className="px-3 py-1 bg-gray-200 rounded-md hover:bg-gray-300 text-[#002C36] font-bold uppercase tracking-wide text-xs"
                  onClick={() => setShowAddAttachmentField(false)}
                >
                  Cancel
                </button>
                <button
                  className="px-3 py-1 bg-[#008080] text-white rounded-md hover:bg-[#006666] font-bold uppercase tracking-wide text-xs"
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
        {/* Attachment Modal Section */}
        {selectedAttachment && (
          <AttachmentModal
            attachmentKey={selectedAttachment}
            data={attachments[selectedAttachment]}
            onClose={() => setSelectedAttachment(null)}
            onSave={(note, fileUrl) => {
              setAttachments(prev => ({
                ...prev,
                [selectedAttachment]: { note, imageUrl: fileUrl },
              }));
              setSelectedAttachment(null);
            }}
          />
        )}
        <button
          onClick={handleSave}
          className="bg-[#008080] hover:bg-[#006666] text-white px-4 py-2 cursor-pointer rounded-lg font-bold uppercase tracking-wide text-base mt-6 float-right"
        >
          Save Compound
        </button>
      </div>
    </div>
  );
}
