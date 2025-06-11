"use client";

import React, { useState } from "react";

export default function CompoundModal({ compound, onClose, onDelete, onUpdate }) {
  const [editMode, setEditMode] = useState(false);
  const [editedCompound, setEditedCompound] = useState({ ...compound });

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
    "Max loading (%)",
    "smiles"
  ];

  const handleChange = (field, value) => {
    setEditedCompound((prev) => ({ ...prev, [field]: value }));
  };

  const handleSave = async () => {
    try {
      await onUpdate(editedCompound); // Backend update function passed as prop
      setEditMode(false);
    } catch (err) {
      console.error("Update failed:", err);
    }
  };

  const handleDelete = async () => {
    const confirmed = window.confirm("Are you sure you want to delete this compound?");
    if (confirmed) {
      try {
        await onDelete(compound.id); // Backend delete function passed as prop
        onClose();
      } catch (err) {
        console.error("Delete failed:", err);
      }
    }
  };

  return (
    <div
      className="fixed inset-0 bg-black/40 backdrop-blur-sm z-50 flex items-center justify-center"
      onClick={onClose}
    >
      <div
        className="bg-white p-6 rounded-lg shadow-xl max-w-2xl w-full max-h-[90vh] overflow-y-auto"
        onClick={(e) => e.stopPropagation()}
      >
        <div className="flex justify-between items-center mb-4">
          <h2 className="text-2xl font-bold text-black">Compound Details</h2>
          <div className="flex gap-2">
            <button
              className="text-red-600 border border-red-600 px-2 py-1 rounded"
              onClick={handleDelete}
            >
              Delete
            </button>
            <button
              className="text-blue-600 border border-blue-600 px-2 py-1 rounded"
              onClick={() => (editMode ? handleSave() : setEditMode(true))}
            >
              {editMode ? "Save" : "Edit"}
            </button>
            <button onClick={onClose} className="text-black text-xl font-bold">
              &times;
            </button>
          </div>
        </div>

        <div className="grid grid-cols-1 sm:grid-cols-2 gap-4 text-black">
          {fields.map((field) => (
            <div key={field} className="flex flex-col">
              <span className="text-sm font-semibold text-gray-600">{field}</span>
              {editMode ? (
                <input
                  className="border rounded px-2 py-1 text-sm"
                  value={editedCompound[field] ?? ""}
                  onChange={(e) => handleChange(field, e.target.value)}
                />
              ) : (
                <span>{compound[field] ?? "N/A"}</span>
              )}
            </div>
          ))}
        </div>
      </div>
    </div>
  );
}

