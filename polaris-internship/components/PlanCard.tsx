"use client";

import { motion } from "framer-motion";
import React, { useState } from "react";

export type PlannedStructure = {
  id: string;
  name: string;
  smiles: string;
  mw?: number;
  dipole?: string;
  beta?: string;
  owner?: string;
  priority: number;
  imageUrl?: string;
  b3lypDipole?: string;
  b3lypBeta?: string;
  notes?: string;
};

interface PlanCardProps {
  plan: PlannedStructure;
  onClose: () => void;
}

export default function PlanCard({ plan, onClose }: PlanCardProps) {
  const [editMode, setEditMode] = useState(false);
  const [form, setForm] = useState({ ...plan });
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Save handler
  async function handleSave() {
    setSaving(true);
    setError(null);
    try {
      const res = await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/update-plan/${plan.id}`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(form),
      });
      const data = await res.json();
      if (!data.success) throw new Error(data.message || "Failed to update plan");
      setEditMode(false);
      // Optionally, you can refresh the page or call a callback to update parent state
      window.location.reload(); // quick solution for now
    } catch (err: any) {
      setError(err.message || "Failed to update plan");
    } finally {
      setSaving(false);
    }
  }

  return (
    <div className="fixed inset-0 z-50 bg-black/40 flex items-center justify-center" onClick={onClose}>
      <motion.div
        className="relative bg-white rounded-2xl shadow-2xl p-8 w-[520px] max-w-full flex flex-col items-center border-2 border-[#00E6D2]"
        initial={{ scale: 0.9, opacity: 0 }}
        animate={{ scale: 1, opacity: 1 }}
        exit={{ scale: 0.9, opacity: 0 }}
        onClick={e => e.stopPropagation()}
      >
        <button
          className="absolute top-5 right-5 text-3xl text-[#00E6D2] hover:text-[#008080] focus:outline-none"
          onClick={onClose}
          title="Close"
        >
          ×
        </button>
        <button
          className="absolute top-5 left-5 px-4 py-2 bg-[#00E6D2] text-[#002C36] rounded font-bold shadow hover:bg-[#00bfae] transition-all"
          onClick={() => setEditMode(e => !e)}
        >
          {editMode ? "Cancel" : "Edit"}
        </button>
        <button
          className="absolute top-5 left-32 px-4 py-2 bg-red-500 text-white rounded font-bold shadow hover:bg-red-600 transition-all"
          onClick={async () => {
            if (window.confirm("Are you sure you want to delete this plan?")) {
              try {
                const res = await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/delete-plan/${plan.id}`, {
                  method: "DELETE",
                });
                const data = await res.json();
                if (!data.success) throw new Error(data.message || "Failed to delete plan");
                window.location.reload();
              } catch (err: any) {
                alert(err.message || "Failed to delete plan");
              }
            }
          }}
        >
          Delete
        </button>
        <div className="w-full flex justify-center mb-6">
          {/* Only show image once, above both view and edit modes */}
          {plan.imageUrl ? (
            <div className="flex items-center justify-center w-[260px] h-[180px] bg-white rounded border border-[#00E6D2] mt-8">
              <img
                src={plan.imageUrl}
                alt={plan.name}
                className="w-full h-full object-contain"
                style={{ display: "block", margin: "auto" }}
                onError={e => { e.currentTarget.style.display = 'none'; }}
              />
            </div>
          ) : (
            <div className="w-[260px] h-[180px] flex items-center justify-center bg-gray-100 rounded border border-[#00E6D2] text-gray-400 mt-8">No image</div>
          )}
        </div>
        {editMode ? (
          <form className="w-full max-w-lg mx-auto text-[#002C36] text-base space-y-3" onSubmit={e => { e.preventDefault(); handleSave(); }}>
            <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
              <div className="flex flex-col">
                <label className="font-bold mb-1">Name</label>
                <input className="border border-[#00E6D2] rounded px-3 py-2 focus:ring-2 focus:ring-[#00E6D2]" value={form.name} onChange={e => setForm(f => ({ ...f, name: e.target.value }))} required />
              </div>
              <div className="flex flex-col">
                <label className="font-bold mb-1">Owner</label>
                <input className="border border-[#00E6D2] rounded px-3 py-2 focus:ring-2 focus:ring-[#00E6D2]" value={form.owner || ""} onChange={e => setForm(f => ({ ...f, owner: e.target.value }))} />
              </div>
              <div className="flex flex-col">
                <label className="font-bold mb-1">MW</label>
                <input type="number" className="border border-[#00E6D2] rounded px-3 py-2 focus:ring-2 focus:ring-[#00E6D2]" value={form.mw || ""} onChange={e => setForm(f => ({ ...f, mw: Number(e.target.value) }))} />
              </div>
              <div className="flex flex-col">
                <label className="font-bold mb-1">Dipole (CAMB3LYP)</label>
                <input className="border border-[#00E6D2] rounded px-3 py-2 focus:ring-2 focus:ring-[#00E6D2]" value={form.dipole || ""} onChange={e => setForm(f => ({ ...f, dipole: e.target.value }))} />
              </div>
              <div className="flex flex-col">
                <label className="font-bold mb-1">Beta (CAMB3LYP)</label>
                <input className="border border-[#00E6D2] rounded px-3 py-2 focus:ring-2 focus:ring-[#00E6D2]" value={form.beta || ""} onChange={e => setForm(f => ({ ...f, beta: e.target.value }))} />
              </div>
              <div className="flex flex-col md:col-span-2">
                <label className="font-bold mb-1">Notes</label>
                <textarea className="border border-[#00E6D2] rounded px-3 py-2 focus:ring-2 focus:ring-[#00E6D2] min-h-[60px]" value={form.notes || ""} onChange={e => setForm(f => ({ ...f, notes: e.target.value }))} />
              </div>
              <div className="flex flex-col md:col-span-2">
                <label className="font-bold mb-1">SMILES</label>
                <input className="border border-[#00E6D2] rounded px-3 py-2 focus:ring-2 focus:ring-[#00E6D2]" value={form.smiles} onChange={e => setForm(f => ({ ...f, smiles: e.target.value }))} />
              </div>
            </div>
            <div className="flex gap-4 mt-6 justify-center">
              <button type="submit" className="bg-[#008080] text-white px-8 py-2 rounded-lg font-bold shadow hover:bg-[#00E6D2] hover:text-[#008080] transition-all text-lg" disabled={saving}>{saving ? "Saving..." : "Save"}</button>
              <button type="button" className="bg-gray-200 text-[#002C36] px-8 py-2 rounded-lg font-bold shadow hover:bg-gray-300 transition-all text-lg" onClick={() => setEditMode(false)} disabled={saving}>Cancel</button>
            </div>
            {error && <div className="text-red-600 font-bold mt-2 text-center">{error}</div>}
          </form>
        ) : (
          <>
            <h2 className="text-3xl font-bold text-[#00E6D2] mb-4 uppercase tracking-wide text-center">{plan.name}</h2>
            <div className="w-full text-[#002C36] text-lg space-y-3">
              <div><span className="font-bold">Owner:</span> {plan.owner || "-"}</div>
              {plan.mw !== undefined && <div><span className="font-bold">MW:</span> {plan.mw}</div>}
              {plan.dipole && <div><span className="font-bold">Dipole (CAMB3LYP):</span> {plan.dipole}</div>}
              {plan.beta && <div><span className="font-bold">Beta (CAMB3LYP):</span> {plan.beta}</div>}
              {plan.b3lypDipole && <div><span className="font-bold">B3LYP Dipole:</span> {plan.b3lypDipole}</div>}
              {plan.b3lypBeta && <div><span className="font-bold">B3LYP Beta:</span> {plan.b3lypBeta}</div>}
              {plan.notes && <div><span className="font-bold">Notes:</span> <span className="whitespace-pre-line">{plan.notes}</span></div>}
              <div><span className="font-bold">SMILES:</span> <span className="break-all">{plan.smiles}</span></div>
            </div>
          </>
        )}
      </motion.div>
    </div>
  );
}

/*
  Copyright © 2025 Polaris Electro Optics
  This code is the property of Polaris Electro Optics and may not be reused,
  modified, or distributed without explicit permission.
*/