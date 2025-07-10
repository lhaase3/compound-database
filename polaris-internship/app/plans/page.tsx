"use client";

import React, { useEffect, useState, useRef } from "react";
import Link from "next/link";
import DrawModal from "@/components/DrawModal";
import PlanCard from "@/components/PlanCard";
import { DragDropContext, Droppable, Draggable, DropResult } from "@hello-pangea/dnd";

// Type for a planned structure (not in main DB)
type PlannedStructure = {
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
  completed?: boolean; // <-- add this
  status?: 'Not Started' | 'In Progress' | 'Completed'; // <-- add status
};

export default function PlansPage() {
  const [plans, setPlans] = useState<PlannedStructure[]>([]);
  const [showDrawModal, setShowDrawModal] = useState(false);
  const [draft, setDraft] = useState<Partial<PlannedStructure>>({ priority: 1 });
  const [showCreatePlan, setShowCreatePlan] = useState(false);
  const [selectedPlan, setSelectedPlan] = useState<PlannedStructure | null>(null);
  const [showFullStructureModal, setShowFullStructureModal] = useState<{ imageUrl: string; name: string } | null>(null);
  const heroRef = useRef<HTMLDivElement>(null);
  const [showStickyLogo, setShowStickyLogo] = useState(false);

  useEffect(() => {
    const handleScroll = () => {
      if (!heroRef.current) return;
      const heroBottom = heroRef.current.getBoundingClientRect().bottom;
      setShowStickyLogo(heroBottom <= 0);
    };
    window.addEventListener("scroll", handleScroll, { passive: true });
    return () => window.removeEventListener("scroll", handleScroll);
  }, []);

  useEffect(() => {
    // Fetch plans from backend on mount
    fetch("http://localhost:5000/plans")
      .then((res) => res.json())
      .then((data) => {
        if (Array.isArray(data)) setPlans(data.sort((a, b) => a.priority - b.priority));
      })
      .catch((err) => console.error("Failed to fetch plans", err));
  }, []);

  // Add helper to upload structure image to backend
  async function uploadStructureImage(smiles: string, planId: string): Promise<string | undefined> {
    try {
      const res = await fetch("http://localhost:5000/upload-image-to-firebase", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ smiles, id: planId, type: "plan" }),
      });
      const data = await res.json();
      return data.fileUrl;
    } catch (err) {
      console.error("Failed to upload structure image", err);
      return undefined;
    }
  }

  // Add a new plan
  const handleAddPlan = async () => {
    if (!draft.name || !draft.smiles) return;
    try {
      // Find if the chosen priority is already used
      const chosenPriority = draft.priority || 1;
      const plansToShift = activePlans.filter(p => p.priority >= chosenPriority);
      // Increment priorities for all plans at or above the chosen priority
      for (const plan of plansToShift.sort((a, b) => b.priority - a.priority)) {
        await fetch(`http://localhost:5000/update-plan/${plan.id}`, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ priority: plan.priority + 1 }),
        });
      }
      // Generate a unique plan id for the image
      const planId = Date.now().toString();
      let imageUrl = draft.imageUrl;
      if (!imageUrl && draft.smiles) {
        imageUrl = await uploadStructureImage(draft.smiles, planId);
      }
      const planData = { ...draft, id: planId, imageUrl, priority: chosenPriority, status: draft.status || 'Not Started', completed: draft.status === 'Completed' };
      const res = await fetch("http://localhost:5000/create-plan", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(planData),
      });
      const result = await res.json();
      if (result.success) {
        // Fetch updated plans list
        const plansRes = await fetch("http://localhost:5000/plans");
        const plansData = await plansRes.json();
        setPlans(Array.isArray(plansData) ? plansData.sort((a, b) => a.priority - b.priority) : []);
        setDraft({ priority: 1 });
      }
    } catch (err) {
      console.error("Failed to create plan", err);
    }
  };

  // Add completed state to plans
  const markPlanCompleted = async (plan: PlannedStructure) => {
    try {
      // 1. Mark the plan as completed and set status
      await fetch(`http://localhost:5000/update-plan/${plan.id}`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ id: plan.id, completed: true, status: 'Completed' }),
      });
      // 2. Re-fetch plans to get the latest state
      const plansRes = await fetch("http://localhost:5000/plans");
      let plansData = await plansRes.json();
      // 3. Find the completed plan's priority
      const completedPriority = plan.priority;
      // 4. For all active plans with priority > completedPriority, decrement their priority
      const activePlansToUpdate = plansData.filter(
        (p: PlannedStructure) => !p.completed && p.priority > completedPriority
      );
      for (const p of activePlansToUpdate) {
        await fetch(`http://localhost:5000/update-plan/${p.id}`, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ priority: p.priority - 1 }),
        });
      }
      // 5. Re-fetch plans again to update UI
      const plansRes2 = await fetch("http://localhost:5000/plans");
      const plansData2 = await plansRes2.json();
      setPlans(Array.isArray(plansData2) ? plansData2.sort((a, b) => a.priority - b.priority) : []);
    } catch (err) {
      console.error("Failed to mark plan completed", err);
    }
  };

  // Split plans into active and completed
  const activePlans = plans.filter(p => p.completed === false || p.completed === undefined);
  const completedPlans = plans.filter(p => p.completed === true);

  // Helper to find the lowest available priority
  function getNextAvailablePriority() {
    const priorities = activePlans.map(p => p.priority).sort((a, b) => a - b);
    let next = 1;
    for (let i = 0; i < priorities.length; i++) {
      if (priorities[i] > next) break;
      if (priorities[i] === next) next++;
    }
    return next;
  }

  // Handle drag end for reordering priorities
  const onDragEnd = async (result: DropResult) => {
    if (!result.destination) return;
    const sourceIdx = result.source.index;
    const destIdx = result.destination.index;
    if (sourceIdx === destIdx) return;
    // Reorder the activePlans array
    const reordered = Array.from(activePlans);
    const [removed] = reordered.splice(sourceIdx, 1);
    reordered.splice(destIdx, 0, removed);
    // Update priorities in the reordered array
    for (let i = 0; i < reordered.length; i++) {
      if (reordered[i].priority !== i + 1) {
        await fetch(`http://localhost:5000/update-plan/${reordered[i].id}`, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ priority: i + 1 }),
        });
      }
    }
    // Refetch plans to update UI
    const plansRes = await fetch("http://localhost:5000/plans");
    const plansData = await plansRes.json();
    setPlans(Array.isArray(plansData) ? plansData.sort((a, b) => a.priority - b.priority) : []);
  };

  const fastScrollToTop = () => {
    const c = document.documentElement.scrollTop || document.body.scrollTop;
    if (c > 0) {
      window.scrollBy(0, -Math.max(120, Math.floor(c / 4)));
      setTimeout(fastScrollToTop, 4);
    }
  };

  // When opening the Create Plan modal, set draft priority to next available
  const handleOpenCreatePlan = () => {
    setDraft(d => ({ ...d, priority: getNextAvailablePriority() }));
    setShowCreatePlan(true);
  };

  // Status badge component
  function StatusBadge({ status }: { status: string }) {
    let color = '';
    if (status === 'Completed') color = 'bg-[#00E6D2] text-[#008080] border-[#008080]';
    else if (status === 'In Progress') color = 'bg-yellow-200 text-yellow-800 border-yellow-400';
    else color = 'bg-gray-200 text-gray-700 border-gray-400';
    return (
      <span className={`inline-block px-3 py-1 rounded-full border font-semibold text-xs ${color}`}>{status}</span>
    );
  }

  // Helper to update plan status
  const updatePlanStatus = async (plan: PlannedStructure, status: 'Not Started' | 'In Progress' | 'Completed') => {
    try {
      await fetch(`http://localhost:5000/update-plan/${plan.id}`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ status, completed: status === 'Completed' }),
      });
      // If marking completed, also update priorities
      if (status === 'Completed') {
        await markPlanCompleted(plan);
        return;
      }
      // Otherwise, just refetch
      const plansRes = await fetch("http://localhost:5000/plans");
      const plansData = await plansRes.json();
      setPlans(Array.isArray(plansData) ? plansData.sort((a, b) => a.priority - b.priority) : []);
    } catch (err) {
      console.error("Failed to update plan status", err);
    }
  }

  return (
    <div className="min-h-screen bg-[#002C36] flex flex-col items-center p-0">
      {/* Sticky Logo Taskbar */}
      <div
        className={`fixed top-0 left-0 w-full z-50 flex justify-start pointer-events-none transition-opacity duration-300 ${
          showStickyLogo ? "opacity-100" : "opacity-0"
        }`}
        style={{ height: "64px" }}
        aria-hidden={!showStickyLogo}
      >
        <button
          onClick={fastScrollToTop}
          className="pointer-events-auto bg-[#002C36]/80 rounded-full shadow-lg p-2 ml-8"
          style={{ marginTop: "8px" }}
          aria-label="Back to top"
        >
          <img
            src="/polaris-logo-only.png"
            alt="Polaris Electro-Optics Logo"
            className="w-16 h-21 drop-shadow-lg"
          />
        </button>
      </div>
      {/* Hero Section */}
      <div ref={heroRef} className="w-full bg-gradient-to-r from-[#00343F] to-[#002C36] py-12 mb-10 shadow flex flex-col items-center relative overflow-hidden">
        {/* Logo in top-left corner */}
        <img src="/polaris-logo-only.png" alt="Polaris Electro-Optics Logo" className={`w-16 h-21 absolute top-6 left-8 z-20 drop-shadow-lg transition-opacity duration-300 ${showStickyLogo ? "opacity-0" : "opacity-100"}`} />
        <div className="absolute inset-0 opacity-30 pointer-events-none select-none" style={{background: 'url(/circuit-bg.svg) center/cover no-repeat'}} />
        <h1 className="text-5xl font-extrabold mb-3 text-[#00E6D2] tracking-tight drop-shadow uppercase z-10">Plans</h1>
        <p className="text-xl text-white mb-6 max-w-2xl text-center z-10 font-semibold flex items-center justify-center gap-3">
          <img src="/white-logo.png" alt="Polaris Logo" className="w-8 h-10 inline-block" />
          Polaris Electro-Optics
        </p>
        <div className="mb-2 z-10">
          <Link href="/">
            <button className="bg-[#00E6D2] hover:bg-[#00bfae] text-[#002C36] px-6 py-2 rounded-lg shadow font-bold text-lg uppercase tracking-wide flex items-center gap-2 transition-all">
              <span role="img" aria-label="home">🏠</span> Back to Home
            </button>
          </Link>
        </div>
      </div>
      <h1 className="text-4xl font-bold text-[#00E6D2] mt-10 mb-6">Planning Workspace</h1>
      <div className="w-full max-w-3xl flex justify-center mb-6">
        <button
          className="bg-[#00E6D2] hover:bg-[#00bfae] text-[#002C36] px-6 py-2 rounded-lg shadow font-bold uppercase tracking-wide flex items-center gap-2 transition-all"
          onClick={handleOpenCreatePlan}
        >
          <span role="img" aria-label="add">➕</span> Create Plan
        </button>
      </div>
      {showCreatePlan && (
        <div className="fixed inset-0 bg-black/40 z-50 flex items-center justify-center" onClick={() => setShowCreatePlan(false)}>
          <div className="bg-white rounded-lg shadow-lg p-6 w-full max-w-xl relative max-h-[90vh] overflow-y-auto" onClick={e => e.stopPropagation()}>
            <h2 className="text-xl font-bold mb-4 text-[#008080]">Create New Plan</h2>
            <div className="mb-3">
              <label className="block font-semibold mb-1">Name</label>
              <input type="text" className="w-full border rounded p-2" value={draft.name || ""} onChange={e => setDraft(d => ({ ...d, name: e.target.value }))} />
            </div>
            <div className="mb-3">
              <label className="block font-semibold mb-1">Priority (1 = highest)</label>
              <input type="number" min={1} className="w-full border rounded p-2" value={draft.priority || 1} onChange={e => setDraft(d => ({ ...d, priority: Number(e.target.value) }))} />
            </div>
            <div className="mb-3">
              <label className="block font-semibold mb-1">Owner</label>
              <input type="text" className="w-full border rounded p-2" value={draft.owner || ""} onChange={e => setDraft(d => ({ ...d, owner: e.target.value }))} />
            </div>
            <div className="mb-3">
              <label className="block font-semibold mb-1">Draw Structure</label>
              <button className="bg-[#00E6D2] text-[#002C36] px-4 py-2 rounded font-bold" onClick={() => setShowDrawModal(true)}>
                {draft.smiles ? "Edit Structure" : "Draw Structure"}
              </button>
              {draft.smiles && (
                <div className="mt-2 text-xs text-gray-700">SMILES: {draft.smiles}</div>
              )}
            </div>
            <div className="mb-3">
              <label className="block font-semibold mb-1">SMILES</label>
              <input
                type="text"
                className="w-full border rounded p-2"
                value={draft.smiles || ""}
                onChange={e => setDraft(d => ({ ...d, smiles: e.target.value }))}
                placeholder="Enter SMILES or use Draw Structure"
              />
            </div>
            <div className="mb-3">
              <label className="block font-semibold mb-1">MW</label>
              <input type="number" className="w-full border rounded p-2" value={draft.mw || ""} onChange={e => setDraft(d => ({ ...d, mw: Number(e.target.value) }))} />
            </div>
            <div className="mb-3">
              <label className="block font-semibold mb-1">Dipole CAMB3LYP SVPD CHCl3 (Cosmo)</label>
              <input type="text" className="w-full border rounded p-2" value={draft.dipole || ""} onChange={e => setDraft(d => ({ ...d, dipole: e.target.value }))} />
            </div>
            <div className="mb-3">
              <label className="block font-semibold mb-1">Beta CAMB3LYP SVPD CHCl3 (Cosmo)</label>
              <input type="text" className="w-full border rounded p-2" value={draft.beta || ""} onChange={e => setDraft(d => ({ ...d, beta: e.target.value }))} />
            </div>
            <div className="mb-3">
              <label className="block font-semibold mb-1">Notes</label>
              <textarea className="w-full border rounded p-2" value={draft.notes || ""} onChange={e => setDraft(d => ({ ...d, notes: e.target.value }))} />
            </div>
            <div className="mb-3">
              <label className="block font-semibold mb-1">Status</label>
              <select
                className="w-full border rounded p-2"
                value={draft.status || 'Not Started'}
                onChange={e => setDraft(d => ({ ...d, status: e.target.value as any }))
                }
              >
                <option value="Not Started">Not Started</option>
                <option value="In Progress">In Progress</option>
                <option value="Completed">Completed</option>
              </select>
            </div>
            <div className="flex justify-end gap-2 mt-4">
              <button className="px-4 py-2 bg-gray-300 rounded hover:bg-gray-400 text-black" onClick={() => setShowCreatePlan(false)}>
                Cancel
              </button>
              <button className="bg-[#008080] text-white px-4 py-2 rounded font-bold" onClick={() => { handleAddPlan(); setShowCreatePlan(false); }}>
                Add Plan
              </button>
            </div>
          </div>
        </div>
      )}
      <div className="w-full max-w-6xl">
        <h2 className="text-2xl font-bold text-[#008080] mb-4">Planned Structures</h2>
        <DragDropContext onDragEnd={onDragEnd}>
          <Droppable droppableId="activePlansTable">
            {(provided) => (
              <table className="w-full min-w-full border border-[#008080] bg-white rounded-lg" ref={provided.innerRef} {...provided.droppableProps}>
                <thead>
                  <tr className="text-xs text-[#008080] uppercase text-center">
                    <th className="px-1 py-2 border-b border-r border-[#008080] w-12">Priority</th>
                    <th className="px-2 py-2 border-b border-r border-[#008080]" style={{ width: '120px', minWidth: '100px', maxWidth: '160px' }}>Name</th>
                    <th className="px-2 py-2 border-b border-r border-[#008080]">Owner</th>
                    <th className="px-2 py-2 border-b border-r border-[#008080]" style={{ width: '240px', minWidth: '200px', maxWidth: '320px' }}>Structure</th>
                    <th className="px-2 py-2 border-b border-r border-[#008080]" style={{ width: '100px', minWidth: '80px', maxWidth: '140px' }}>MW</th>
                    <th className="px-2 py-2 border-b border-r border-[#008080]" style={{ width: '100px', minWidth: '80px', maxWidth: '140px' }}>Status</th>
                  </tr>
                </thead>
                <tbody>
                  {activePlans.length === 0 && (
                    <tr><td colSpan={7} className="text-center py-4 text-gray-500">No plans yet.</td></tr>
                  )}
                  {activePlans.map((plan, idx) => (
                    <Draggable key={plan.id} draggableId={plan.id} index={idx}>
                      {(provided, snapshot) => (
                        <tr
                          ref={provided.innerRef}
                          {...provided.draggableProps}
                          {...provided.dragHandleProps}
                          className={`text-center border-b border-[#008080] hover:bg-[#e6f9f7] cursor-pointer group ${snapshot.isDragging ? 'bg-[#b2f0e9]' : ''}`}
                          onClick={() => setSelectedPlan(plan)}
                        >
                          <td className="px-1 py-2 border-r border-[#008080] font-bold w-12">{plan.priority}</td>
                          <td className="px-2 py-2 border-r border-[#008080]" style={{ width: '160px', minWidth: '100px', maxWidth: '160px' }}>{plan.name}</td>
                          <td className="px-2 py-2 border-r border-[#008080]" style={{ width: '160px', minWidth: '100px', maxWidth: '160px' }}>{plan.owner}</td>
                          {/* update this to increase structure column size */}
                          <td className="border-r border-[#008080] bg-white" style={{ padding: 0, verticalAlign: 'middle', width: '320px', minWidth: '400px', maxWidth: '320px' }}>
                            {plan.imageUrl ? (
                              <div
                                className="flex items-center justify-center bg-white rounded mx-auto overflow-hidden border border-gray-200 group relative"
                                style={{ width: '280px', height: '240px', padding: '8px', boxSizing: 'border-box' }}
                                tabIndex={0}
                                aria-label={`Structure image for ${plan.name}`}
                              >
                                <img
                                  src={plan.imageUrl}
                                  alt={plan.name}
                                  className="object-contain"
                                  style={{ width: '100%', height: '100%', display: 'block', margin: 0, padding: 0 }}
                                />
                                <button
                                  className="absolute bottom-2 right-2 bg-white/90 text-xs text-[#008080] px-3 py-1 rounded shadow opacity-0 group-hover:opacity-100 transition-opacity font-bold border border-[#00E6D2] hover:bg-[#00E6D2] hover:text-white focus:opacity-100 focus:outline-none"
                                  onClick={e => { e.stopPropagation(); setShowFullStructureModal({ imageUrl: plan.imageUrl!, name: plan.name }); }}
                                  tabIndex={0}
                                  aria-label={`View enlarged structure for ${plan.name}`}
                                  type="button"
                                >
                                  View
                                </button>
                              </div>
                            ) : (
                              <span className="text-gray-400">No image</span>
                            )}
                          </td>
                          <td className="px-2 py-2 border-r border-[#008080]" style={{ width: '160px', minWidth: '150px', maxWidth: '160px' }}>{plan.mw}</td>
                          <td className="px-2 py-2 border-r border-[#008080]" style={{ width: '100px', minWidth: '180px', maxWidth: '140px' }}>
                            {/* Status dropdown for active plans, badge for completed */}
                            {plan.status === 'Completed' || plan.completed ? (
                              <StatusBadge status="Completed" />
                            ) : (
                              <select
                                className="border rounded px-2 py-1 text-xs font-semibold bg-white text-[#008080] focus:outline-none focus:ring-2 focus:ring-[#00E6D2]"
                                value={plan.status || 'Not Started'}
                                onClick={e => e.stopPropagation()}
                                onChange={e => { e.stopPropagation(); updatePlanStatus(plan, e.target.value as any); }}
                                style={{ minWidth: '110px' }}
                              >
                                <option value="Not Started">Not Started</option>
                                <option value="In Progress">In Progress</option>
                                <option value="Completed">Completed</option>
                              </select>
                            )}
                          </td>
                        </tr>
                      )}
                    </Draggable>
                  ))}
                  {provided.placeholder}
                </tbody>
              </table>
            )}
          </Droppable>
        </DragDropContext>
        {completedPlans.length > 0 && (
          <>
            <h3 className="text-xl font-bold text-[#008080] mt-8 mb-2">Completed</h3>
            <table className="min-w-full border border-[#008080] bg-white rounded-lg">
              <thead>
                <tr className="text-xs text-[#008080] uppercase text-center">
                  <th className="px-2 py-2 border-b border-r border-[#008080]" style={{ width: '170px', minWidth: '150px', maxWidth: '200px' }}>Name</th>
                  <th className="px-2 py-2 border-b border-[#008080]" style={{ width: '170px', minWidth: '150px', maxWidth: '200px' }}>Owner</th>
                  <th className="px-2 py-2 border-b border-r border-[#008080]" style={{ width: '280px', minWidth: '420px', maxWidth: '320px' }}>Structure</th>
                  <th className="px-2 py-2 border-b border-r border-[#008080]" style={{ width: '100px', minWidth: '190px', maxWidth: '140px' }}>MW</th>
                  <th className="px-2 py-2 border-b border-[#008080]">Status</th>
                </tr>
              </thead>
              <tbody>
                {completedPlans.map((plan) => (
                  <tr key={plan.id} className="text-center border-b border-[#008080] bg-[#e6f9f7] group hover:bg-[#d2f5ee] cursor-pointer" onClick={() => setSelectedPlan(plan)}>
                    <td className="px-2 py-2 border-r border-[#008080]">{plan.name}</td>
                    <td className="px-2 py-2 border-r border-[#008080] relative">
                      {plan.owner}
                    </td>
                    <td className="border-r border-[#008080] bg-white" style={{ padding: 0, verticalAlign: 'middle', width: '280px', minWidth: '240px', maxWidth: '320px' }}>
                      {plan.imageUrl ? (
                        <div
                          className="flex items-center justify-center bg-white rounded mx-auto overflow-hidden border border-gray-200 group relative"
                          style={{ width: '280px', height: '240px', padding: '8px', boxSizing: 'border-box' }}
                          tabIndex={0}
                          aria-label={`Structure image for ${plan.name}`}
                        >
                          <img
                            src={plan.imageUrl}
                            alt={plan.name}
                            className="object-contain"
                            style={{ width: '100%', height: '100%', display: 'block', margin: 0, padding: 0 }}
                          />
                          <button
                            className="absolute bottom-2 right-2 bg-white/90 text-xs text-[#008080] px-3 py-1 rounded shadow opacity-0 group-hover:opacity-100 transition-opacity font-bold border border-[#00E6D2] hover:bg-[#00E6D2] hover:text-white focus:opacity-100 focus:outline-none"
                            onClick={e => { e.stopPropagation(); setShowFullStructureModal({ imageUrl: plan.imageUrl!, name: plan.name }); }}
                            tabIndex={0}
                            aria-label={`View enlarged structure for ${plan.name}`}
                            type="button"
                          >
                            View
                          </button>
                        </div>
                      ) : (
                        <span className="text-gray-400">No image</span>
                      )}
                    </td>
                    <td className="px-2 py-2 border-r border-[#008080] relative">
                      {plan.mw}
                    </td>
                    <td className="px-2 py-2 border-r border-[#008080]">
                      <span className="inline-flex items-center gap-1 px-2 py-1 rounded bg-[#00E6D2] text-[#008080] font-bold text-xs border border-[#008080] shadow">
                        <svg width="16" height="16" viewBox="0 0 18 18" fill="none" xmlns="http://www.w3.org/2000/svg">
                          <circle cx="9" cy="9" r="8" stroke="#008080" strokeWidth="2" fill="#00E6D2" />
                          <path d="M5 9.5L8 12.5L13 7.5" stroke="#008080" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" />
                        </svg>
                        Completed
                      </span>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </>
        )}
      </div>
      {selectedPlan && (
        <PlanCard plan={selectedPlan} onClose={() => setSelectedPlan(null)} />
      )}
      {showDrawModal && (
        <DrawModal
          onClose={() => setShowDrawModal(false)}
          onSmilesSubmit={async (smiles: string) => {
            let imageUrl = undefined;
            if (smiles) {
              // Generate a temporary plan id for image upload
              const tempPlanId = draft.id || Date.now().toString();
              imageUrl = await uploadStructureImage(smiles, tempPlanId);
            }
            setDraft(d => ({ ...d, smiles, imageUrl }));
            setShowDrawModal(false);
          }}
        />
      )}
      {showFullStructureModal && (
        <div className="fixed inset-0 bg-black/70 z-50 flex justify-center items-center" onClick={() => setShowFullStructureModal(null)}>
          <div className="bg-white p-4 rounded-2xl border-2 border-[#00E6D2] shadow-2xl max-w-6xl w-full max-h-[96vh] overflow-auto relative flex flex-col items-center" onClick={e => e.stopPropagation()}>
            <img
              src={showFullStructureModal.imageUrl}
              alt={showFullStructureModal.name}
              className="object-contain mx-auto"
              style={{ width: '900px', height: '600px', maxWidth: '100%', maxHeight: '85vh', display: 'block' }}
            />
            <div className="text-center mt-4 text-lg font-bold text-[#008080]">{showFullStructureModal.name}</div>
            <button
              className="absolute top-2 right-4 text-3xl text-[#008080] font-bold bg-white/80 rounded-full px-3 py-1 hover:bg-[#00E6D2]/80 hover:text-white transition-all"
              onClick={() => setShowFullStructureModal(null)}
              aria-label="Close enlarged structure"
            >
              &times;
            </button>
          </div>
        </div>
      )}
    </div>
  );
}
