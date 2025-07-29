"use client"

import React, { useEffect, useState, useRef } from "react";
import { useRouter } from "next/navigation";
import { auth } from "@/utils/firebase";
import { onAuthStateChanged, signInWithPopup, GoogleAuthProvider, signOut } from "firebase/auth";
import Link from "next/link";
import CompoundModal from "@/components/CompoundModal";
import CreateFormulationModal from "@/components/CreateFormulationModal";
import DrawModal from "@/components/DrawModal";


function recalculateComponents(components: any[], totalMass: number, compounds: any[]) {
  const totalTargetMoles = components.reduce((sum, c) => {
    const mw = parseFloat(c.molecularWeight || c.MW) || 1;  // ✅ Correct field now
    const targetMassLocal = (parseFloat(c.massPercent) || 0) / 100 * totalMass;
    return sum + (targetMassLocal / mw);
  }, 0);


  const updatedComponents = components.map(comp => {
    const compound = compounds.find(c => c.id === comp.compoundId);  // 🔑 THIS WAS MISSING
    const massPercent = parseFloat(comp.massPercent) || 0;
    const actualMass = parseFloat(comp.actualMass) || 0;
    const molecularWeight = parseFloat(comp.molecularWeight || comp.MW) || 1;

    const targetMass = (massPercent / 100) * totalMass;
    const targetMoles = targetMass / molecularWeight;

    const totalActualMass = components.reduce((sum, c) => sum + (parseFloat(c.actualMass) || 0), 0);
    const actualMassPercent = totalActualMass ? (actualMass / totalActualMass) * 100 : 0;

    const totalActualMoles = components.reduce((sum, c) => {
      const mw = parseFloat(c.molecularWeight || c.MW) || 1;
      return sum + ((parseFloat(c.actualMass) || 0) / mw);
    }, 0);
    const actualMolPercent = totalActualMoles ? ((actualMass / molecularWeight) / totalActualMoles) * 100 : 0;

    const outputMolPercent = totalTargetMoles ? (targetMoles / totalTargetMoles) * 100 : 0;

    return {
      ...comp,
    smiles: comp.smiles || compound?.smiles || "",
    compoundId: compound?.id || comp.compoundId,
    compoundName: compound?.name || '',
    targetMass: parseFloat(targetMass.toFixed(3)),
    targetMoles: parseFloat(targetMoles.toFixed(4)),
    outputMolPercent: parseFloat(outputMolPercent.toFixed(2)),
    actualMass,
    actualMassPercent: parseFloat(actualMassPercent.toFixed(2)),
    actualMolPercent: parseFloat(actualMolPercent.toFixed(2)),
    };
  });
  return updatedComponents;
}



// --- Types for editData state ---
type EditComponent = {
  massPercent?: number | string;
  actualMass?: number | string;
  [key: string]: any;
};
type EditData = {
  name: string;
  phaseMap: string;
  notes: string;
  totalMass?: number | string;
  components?: EditComponent[];
  [key: string]: any;
};



export default function FormulationList() {
  const router = useRouter();
  // Firebase Auth
  const [user, setUser] = useState<any>(null);
  const [authChecked, setAuthChecked] = useState(false);
  const provider = new GoogleAuthProvider();
  const [formulations, setFormulations] = useState<any[]>([]);
  const [selectedFormulation, setSelectedFormulation] = useState<any | null>(null);
  const [editMode, setEditMode] = useState(false);
  const [editData, setEditData] = useState<EditData>({ name: "", phaseMap: "", notes: "" });
  const [selectedCompound, setSelectedCompound] = useState<any | null>(null);
  const [compoundSource, setCompoundSource] = useState<"main" | "lot">("main");
  const [compoundLotId, setCompoundLotId] = useState<string | null>(null);
  const [searchTerm, setSearchTerm] = useState("");
  const [showStickyLogo, setShowStickyLogo] = useState(false);
  const [showAddTextField, setShowAddTextField] = useState(false);
  const [newTextFieldName, setNewTextFieldName] = useState("");
  const [newTextFieldValue, setNewTextFieldValue] = useState("");
  const [showAddAttachmentField, setShowAddAttachmentField] = useState(false);
  const [newAttachmentFieldName, setNewAttachmentFieldName] = useState("");
  const [showFormulationModal, setShowFormulationModal] = useState(false);
  const [compounds, setCompounds] = useState<any[]>([]);
  const [lotMapping, setLotMapping] = useState<Record<string, string[]>>({});
  const [starredFormulations, setStarredFormulations] = useState<string[]>([]);
  const [showDrawModal, setShowDrawModal] = useState(false);
  const [showOnlyStarred, setShowOnlyStarred] = useState(false);
  const heroRef = useRef<HTMLDivElement>(null);
  const [selectedAttachment, setSelectedAttachment] = useState<{
  name: string;
  data: { note: string; imageUrl: string };
} | null>(null);
  const [deleteConfirm, setDeleteConfirm] = useState(false);
  const [showFullStructureModal, setShowFullStructureModal] = useState<{
    imageUrl: string;
    smiles?: string;
    name: string;
  } | null>(null);
  // Firebase Auth: Check if user is logged in
  useEffect(() => {
    if (!auth) return;
    const unsubscribe = onAuthStateChanged(auth, (firebaseUser) => {
      setUser(firebaseUser);
      setAuthChecked(true);
    });
    return () => unsubscribe();
  }, []);

  // Load starred formulations for user
  useEffect(() => {
    if (!authChecked) return;

    const stored = localStorage.getItem("starredFormulations");

    // If user is logged in, load from backend
    if (user && user.email) {
      fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/get-starred-formulations?email=${user.email}`)
        .then(res => res.json())
        .then(data => {
          if (Array.isArray(data.starred)) {
            setStarredFormulations(data.starred);
            localStorage.setItem("starredFormulations", JSON.stringify(data.starred));
          }
        })
        .catch(err => {
          console.error("Failed to fetch from backend, falling back to localStorage");
          if (stored) setStarredFormulations(JSON.parse(stored));
        });
    } else if (stored) {
      setStarredFormulations(JSON.parse(stored));
    }
  }, [authChecked, user]);



  // Sync starred formulations to backend and localStorage
  useEffect(() => {
    localStorage.setItem("starredFormulations", JSON.stringify(starredFormulations));
    if (user && user.email) {
      fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/save-starred-formulations`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ email: user.email, starred: starredFormulations }),
      }).catch((err) => console.error("Failed to save starred formulations:", err));
    }
  }, [starredFormulations, user]);


  const handleResetFilters = async () => {
    setSearchTerm("");
    setShowOnlyStarred(false);

    try {
      const res = await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/formulations`);
      const data = await res.json();
      setFormulations(Array.isArray(data) ? data : []);
    } catch (err) {
      console.error("Error resetting filters", err);
      setFormulations([]);
    }
  };




  useEffect(() => {
    fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/formulations`)
      .then((res) => res.json())
      .then(setFormulations)
      .catch((err) => console.error("Failed to fetch formulations", err));
  }, []);

  useEffect(() => {
    fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/compounds`)
      .then((res) => res.json())
      .then((data) => setCompounds(data))
      .catch((err) => console.error("Failed to fetch compounds", err));
  }, []);

  useEffect(() => {
    fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/lots`)
      .then((res) => res.json())
      .then((data) => {
        const mapping: Record<string, string[]> = {};
        data.forEach((lot: any) => {
          if (!mapping[lot.compoundId]) {
            mapping[lot.compoundId] = [];
          }
          mapping[lot.compoundId].push(lot.lotId);
        });
        setLotMapping(mapping);
      })
      .catch((err) => console.error("Failed to fetch lots", err));
  }, []);

  useEffect(() => {
    if (typeof window !== "undefined") {
      const stored = window.localStorage.getItem("starredFormulations");
      if (stored) setStarredFormulations(JSON.parse(stored));
    }
  }, []);

  useEffect(() => {
    if (typeof window !== "undefined") {
      window.localStorage.setItem("starredFormulations", JSON.stringify(starredFormulations));
    }
  }, [starredFormulations]);

  useEffect(() => {
    const handleScroll = () => {
      if (!heroRef.current) return;
      const heroBottom = heroRef.current.getBoundingClientRect().bottom;
      setShowStickyLogo(heroBottom <= 0);
    };
    window.addEventListener("scroll", handleScroll, { passive: true });
    return () => window.removeEventListener("scroll", handleScroll);
  }, []);

  // Fast scroll-to-top function
  const fastScrollToTop = () => {
    if (typeof window !== "undefined") {
      const c = document.documentElement.scrollTop || document.body.scrollTop;
      if (c > 0) {
        window.scrollBy(0, -Math.max(120, Math.floor(c / 4)));
        setTimeout(fastScrollToTop, 4);
      }
    }
  };

  const toggleStarFormulation = (id: string) => {
    if (!user || !user.email) {
      alert("You must be signed in to star formulations.");
      return;
    }
    setStarredFormulations((prev) =>
      prev.includes(id) ? prev.filter((fid) => fid !== id) : [...prev, id]
    );
  };

  const handleDrawFilterSubmit = async (smiles: string) => {
    try {
      const res = await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/search-substructure-formulations`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ query: smiles }),
      });
      if (!res.ok) {
        throw new Error(`Request failed: ${res.status}`);
      }
      const data = await res.json();
      setFormulations(Array.isArray(data) ? data : []);
    } catch {
      // error handled by alert below
      alert("Failed to filter by drawn structure. Please check your backend and network connection.");
      setFormulations([]);
    }
  };



  // When opening edit mode, initialize editData with all editable fields
  const handleEditClick = () => {
    const editable: EditData = {
      name: selectedFormulation.name || "",
      phaseMap: selectedFormulation.phaseMap || "",
      notes: selectedFormulation.notes || "",
      totalMass: selectedFormulation.totalMass || "",
      components: selectedFormulation.components?.map((c: any) => ({
        ...c // Carry over everything: molecularWeight, compoundId, compoundName, imageUrl, lotId, etc.
      })) || [],
    };
    // Add all custom fields
    Object.entries(selectedFormulation).forEach(([key, value]) => {
      if (!["id", "name", "components", "phaseMap", "notes", "attachments", "createdAt", "imageUrls", "totalmoles"].includes(key)) {
        editable[key] = value;
      }
    });
    setEditData(editable);
    setEditMode(true);
  };

  useEffect(() => {
    if (authChecked && !user) {
      router.push("/login");
    }
  }, [authChecked, user, router]);
  if (!user) {
    return null;
  }
  return (
    <div className="min-h-screen bg-[#002C36] flex flex-col items-center p-0">
      {/* Auth Bar */}
      <div className="w-full flex justify-end items-center px-8 py-2">
        {user ? (
          <div className="flex gap-4 items-center">
            <span className="text-[#00E6D2] font-bold">{user.email}</span>
            <button
              className="bg-[#00E6D2] hover:bg-[#00bfae] text-[#002C36] px-4 py-2 rounded-lg font-bold shadow uppercase tracking-wide"
              onClick={() => signOut(auth)}
            >
              Logout
            </button>
          </div>
        ) : (
          <button
            className="bg-[#00E6D2] hover:bg-[#00bfae] text-[#002C36] px-4 py-2 rounded-lg font-bold shadow uppercase tracking-wide"
            onClick={async () => {
              if (!auth) {
                alert("Auth not initialized");
                return;
              }
              try {
                await signInWithPopup(auth, provider);
              } catch {
                alert("Sign in failed");
              }
            }}
          >
            Login
          </button>
        )}
      </div>
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
          className="pointer-events-auto p-2 ml-8"
          style={{ marginTop: "8px" }}
          aria-label="Back to top"
        >
          <img
            src="/polaris-logo-only.PNG"
            alt="Polaris Electro-Optics Logo"
            className="w-16 h-21 drop-shadow-lg"
          />
        </button>
      </div>
      {/* Hero Section */}
      <div ref={heroRef} className="w-full bg-gradient-to-r from-[#00343F] to-[#002C36] py-12 mb-10 shadow flex flex-col items-center relative overflow-hidden">
        {/* Logo in top-left corner */}
        <img src="/polaris-logo-only.PNG" alt="Polaris Electro-Optics Logo" className={`w-16 h-21 absolute top-6 left-8 z-20 drop-shadow-lg transition-opacity duration-300 ${showStickyLogo ? "opacity-0" : "opacity-100"}`} />
        <div className="absolute inset-0 opacity-30 pointer-events-none select-none"/>
        <h1 className="text-5xl font-extrabold mb-3 text-[#00E6D2] tracking-tight drop-shadow uppercase z-10">Formulations</h1>
        <p className="text-xl text-white mb-6 max-w-2xl text-center z-10 font-semibold flex items-center justify-center gap-3">
          <img src="/white-logo.PNG" alt="Polaris Logo" className="w-8 h-10 inline-block" />
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

      {/* Create Formulation Button */}
      <div className="mb-6">
        <button
          onClick={() => setShowFormulationModal(true)}
          className="bg-[#00E6D2] hover:bg-[#00bfae] text-[#002C36] px-6 py-2 rounded-lg shadow font-bold uppercase tracking-wide flex items-center gap-2 transition-all"
        >
          <span role="img" aria-label="formulation">🧪</span> Create Formulation
        </button>
      </div>


      <div className="flex flex-wrap items-center gap-4 mb-8 w-full max-w-3xl justify-center">
        {!user && (
          <button
            className="bg-[#00E6D2] hover:bg-[#00bfae] text-[#002C36] px-6 py-2 rounded-lg shadow font-bold text-lg uppercase tracking-wide flex items-center gap-2 transition-all"
        onClick={async () => {
          if (!auth) {
            alert("Auth not initialized");
            return;
          }
          try {
            await signInWithPopup(auth, provider);
          } catch {
            alert("Sign in failed");
          }
        }}
          >
            <span role="img" aria-label="google">🔒</span> Sign in with Google
          </button>
        )}
        <input
          type="text"
          placeholder="🔍 Search by compound ID (e.g. PEO-0100)"
          value={searchTerm}
          onChange={(e) => setSearchTerm(e.target.value)}
          className="flex-1 min-w-[250px] border border-[#00E6D2] bg-[#00343F] px-3 py-2 rounded text-[#00E6D2] focus:outline-none focus:ring-2 focus:ring-[#00E6D2] placeholder-[#00E6D2]/60"
        />
        
        <label className="flex items-center gap-2 cursor-pointer text-[#00E6D2] font-bold">
          <input
            type="checkbox"
            checked={showOnlyStarred}
            onChange={() => setShowOnlyStarred(prev => !prev)}
            className="accent-[#00E6D2]"
          />
          Starred Only
        </label>

        <button
          onClick={() => setShowDrawModal(true)}
          className="bg-[#00343F] hover:bg-[#00545F] text-[#00E6D2] px-6 py-2 rounded-lg shadow font-bold uppercase tracking-wide flex items-center gap-2 transition-all border border-[#00E6D2]"
        >
          <span role="img" aria-label="draw">✏️</span> Draw Structure
        </button>

        <button
          onClick={handleResetFilters}
          className="bg-[#00343F] hover:bg-[#00545F] text-[#00E6D2] px-6 py-2 rounded-lg shadow font-bold uppercase tracking-wide flex items-center gap-2 transition-all border border-[#00E6D2]"
        >
          <span role="img" aria-label="reset">🔄</span> Reset Filter
        </button>
      </div>

      <div className="w-full max-w-6xl flex items-center mb-6">
        <hr className="flex-grow border-t border-[#00E6D2]/40" />
        <span className="mx-4 text-lg text-[#00E6D2] font-bold uppercase tracking-wide">Formulations</span>
        <hr className="flex-grow border-t border-[#00E6D2]/40" />
      </div>

      {/* Formulations Grid */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-8 w-full max-w-6xl px-4">
      {formulations
        .filter((form) =>
          (searchTerm.trim() === "" ||
          form.components?.some((comp: any) =>
            (comp.compoundId?.toLowerCase().includes(searchTerm.toLowerCase()) ||
            comp.compoundName?.toLowerCase().includes(searchTerm.toLowerCase()))
          )) && 
          (!showOnlyStarred || starredFormulations.includes(form.id))
        )
        .map((form) => (
            <div
              key={form.id}
              className="relative bg-white border-2 border-[#008080] rounded-2xl p-6 shadow-lg hover:shadow-2xl cursor-pointer transition-all text-[#002C36]"
              onClick={() => setSelectedFormulation(form)}
            >
              {/* Header with title and star */}
              <div className="flex justify-between items-center mb-2">
                <h2 className="text-2xl font-bold uppercase tracking-wide text-[#008080]">
                  {form.name || "Unnamed Formulation"}
                </h2>
                <button
                  className={`absolute top-3 right-3 text-3xl z-10 focus:outline-none ${starredFormulations.includes(form.id) ? 'text-[#00E6D2]' : 'text-[#00E6D2]/50 hover:text-[#00E6D2]'}`}
                  title={starredFormulations.includes(form.id) ? 'Unstar' : 'Star'}
                  onClick={(e) => { e.stopPropagation(); toggleStarFormulation(form.id); }}
                  tabIndex={0}
                  style={{ textShadow: starredFormulations.includes(form.id) ? `0 0 8px #00E6D2` : '0 0 8px #00E6D2' }}
                >
                  {starredFormulations.includes(form.id) ? '★' : '☆'}
                </button>
              </div>

              <p className="text-xs font-bold uppercase text-[#008080] mb-1 tracking-wide">Components:</p>
              <ul className="list-disc pl-5 text-sm">
                {form.components?.map((comp: any, idx: number) => (
                  <li key={idx}>
                    <span className="font-semibold text-[#002C36]">{comp.compoundName || comp.compoundId}</span> <span className="text-[#008080]">({comp.lotId || "original"})</span> – <span className="text-[#008080]">{comp.massPercent !== undefined ? comp.massPercent + "%" : "-"}</span>
                  </li>
                )) || <li>No components</li>}
              </ul>
            </div>
          ))}
      </div>

      {selectedFormulation && (
        <div
          className="fixed inset-0 z-50 flex items-center justify-center bg-black/40 backdrop-blur-sm"
          onClick={() => setSelectedFormulation(null)}
        >
          {/* Delete Confirmation Modal */}
          {deleteConfirm && (
            <div className="fixed inset-0 bg-black/50 z-60 flex items-center justify-center p-6" onClick={() => setDeleteConfirm(false)}>
              <div className="bg-white rounded-lg shadow-xl w-full max-w-md p-8" onClick={e => e.stopPropagation()}>
                <h2 className="text-xl font-bold mb-4 text-red-700">Are you sure you want to delete '{selectedFormulation.name || "this formulation"}'?</h2>
                <div className="flex justify-end gap-2">
                  <button
                    className="px-4 py-2 bg-gray-300 rounded hover:bg-gray-400 text-black"
                    onClick={() => setDeleteConfirm(false)}
                  >
                    Cancel
                  </button>
                  <button
                    className="px-4 py-2 bg-red-600 text-white rounded hover:bg-red-700"
                    onClick={async () => {
                      try {
                        await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/delete-formulation/${selectedFormulation.id}`, { method: "DELETE" });
                        setFormulations((prev) => prev.filter((f) => f.id !== selectedFormulation.id));
                        setSelectedFormulation(null);
                        setDeleteConfirm(false);
                      } catch (err) {
                        console.error("Failed to delete formulation", err);
                        alert("Failed to delete formulation.");
                      }
                    }}
                  >
                    Delete
                  </button>
                </div>
              </div>
            </div>
          )}
          <div
            className="bg-white border-2 border-[#008080] p-8 rounded-2xl shadow-2xl max-w-7xl w-full max-h-[95vh] overflow-y-auto relative"
            onClick={(e) => e.stopPropagation()}
          >
            <div className="flex justify-between items-start mb-4">
              <h2 className="text-3xl font-extrabold text-[#002C36] uppercase tracking-wide">{selectedFormulation.name || "Unnamed Formulation"}</h2>
              <div className="flex gap-2">
                <button
                  onClick={handleEditClick}
                  className="bg-[#00E6D2] hover:bg-[#00bfae] text-[#002C36] font-bold px-2 py-0.5 rounded-md uppercase tracking-wide text-xs shadow transition-all"
                >
                  Edit
                </button>
                <button
                  onClick={() => setDeleteConfirm(true)}
                  className="bg-red-100 hover:bg-red-200 text-red-600 font-bold px-2 py-0.5 rounded-md uppercase tracking-wide text-xs shadow transition-all"
                >
                  Delete
                </button>
              </div>
            </div>

            <div className="mb-4">
              <span className="text-xs font-bold uppercase text-[#008080] mb-1 tracking-wide">Created:</span> <span className="text-[#002C36]">{selectedFormulation.createdAt
                ? new Date(
                  selectedFormulation.createdAt.seconds
                    ? selectedFormulation.createdAt.seconds * 1000
                    : selectedFormulation.createdAt
                ).toLocaleDateString()
                : "Unknown"}</span>
            </div>

            <div className="mb-4">
              <span className="text-xs font-bold uppercase text-[#008080] mb-1 tracking-wide">Operator:</span> <span className="text-[#002C36]">{selectedFormulation.operator}</span>
            </div>

            {/* Mass Calculator Table in Modal */}
            <div className="mb-8 border border-[#008080] rounded-lg p-4 bg-[#F8FAFB]">
              <div className="flex items-center justify-between mb-4">
                <h3 className="text-lg font-bold text-[#008080]">Mass Calculator</h3>
                <span className="text-sm font-semibold text-[#008080]">Total Desired Mass: <span className="text-[#002C36]">{selectedFormulation.totalMass}</span></span>
              </div>
              <div className="overflow-x-auto">
                <table className="min-w-full border border-[#008080]">
                  <thead>
                    <tr className="text-xs text-[#008080] uppercase text-center">
                      <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Compound</th>
                      <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Image</th>
                      <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Lot</th>
                      <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Desired Mass %</th>
                      <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Output Mol %</th>
                      <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Target Mass (mg)</th>
                      <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Actual Mass (mg)</th>
                      <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Actual Mass %</th>
                      <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Actual Mol %</th>
                    </tr>
                  </thead>
                  <tbody>
                    {selectedFormulation.components?.map((comp: any, idx: number) => (
                      <tr key={idx} className="align-middle text-center border-b border-[#008080]">
                        <td className="px-2 py-2 border-r border-[#008080]">
                          <button
                            className="text-blue-600 underline hover:text-blue-800 text-sm cursor-pointer"
                            onClick={async () => {
                              try {
                                const endpoint = comp.lotId
                                  ? `${process.env.NEXT_PUBLIC_API_BASE_URL}/lot/${comp.lotId}`
                                  : `${process.env.NEXT_PUBLIC_API_BASE_URL}/compounds/${comp.compoundId}`;
                                const res = await fetch(endpoint);
                                const data = await res.json();
                                let compound;
                                if (comp.lotId) {
                                  compound = data.find((c: any) =>
                                    c.id === comp.compoundId || c.name?.includes(comp.compoundId)
                                  ) || data[0];
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
                            {comp.compoundName || comp.compoundId}
                          </button>
                        </td>
                        <td className="px-2 py-2 border-r border-[#008080]">
                          {comp.imageUrl && typeof comp.imageUrl === 'string' && comp.imageUrl.trim() !== '' ? (
                            <img
                              src={comp.imageUrl}
                              alt={comp.compoundName || comp.compoundId}
                              style={{ maxWidth: '300px', maxHeight: '300px', objectFit: 'contain', background: 'white', borderRadius: '0.375rem', border: '1px solid #e5e7eb', marginBottom: '0.25rem', cursor: 'pointer' }}
                              onClick={() => setShowFullStructureModal({ imageUrl: comp.imageUrl, smiles: comp.smiles, name: comp.compoundName || comp.compoundId })}
                              onError={e => { e.currentTarget.style.display = 'none'; }}
                            />
                          ) : null}
                        </td>
                        <td className="px-2 py-2 border-r border-[#008080] text-black">{comp.lotId || "original"}</td>
                        <td className="px-2 py-2 border-r border-[#008080] text-black">{comp.massPercent !== undefined ? comp.massPercent + "%" : "-"}</td>
                        <td className="px-2 py-2 border-r border-[#008080] text-black">{comp.outputMolPercent !== undefined ? comp.outputMolPercent.toFixed(2) + "%" : "-"}</td>
                        <td className="px-2 py-2 border-r border-[#008080] text-black">{comp.targetMass !== undefined ? comp.targetMass + " mg" : "-"}</td>
                        <td className="px-2 py-2 border-r border-[#008080] text-black">{comp.actualMass !== undefined ? comp.actualMass + " mg" : "-"}</td>
                        <td className="px-2 py-2 border-r border-[#008080] text-black">{comp.actualMassPercent !== undefined ? comp.actualMassPercent + "%" : "-"}</td>
                        <td className="px-2 py-2 border-r border-[#008080] text-black">{comp.actualMolPercent !== undefined ? comp.actualMolPercent + "%" : "-"}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>

            <div className="mb-4">
              <span className="text-xs font-bold uppercase text-[#008080] mb-1 tracking-wide">Phase Map</span>
              <div className="bg-gray-100 p-3 rounded whitespace-pre-wrap text-[#002C36]">
                {selectedFormulation.phaseMap || "N/A"}
              </div>
            </div>

            <div className="mb-4">
              <span className="text-xs font-bold uppercase text-[#008080] mb-1 tracking-wide">Analytical Notes</span>
              <div className="bg-gray-100 p-3 rounded whitespace-pre-wrap text-[#002C36]">
                {selectedFormulation.notes || "N/A"}
              </div>
            </div>

            {/* Custom Fields */}
            <div className="mb-4">
              <span className="text-xs font-bold uppercase text-[#008080] mb-1 tracking-wide">Other Data</span>
              <div className="mt-2 flex flex-wrap gap-4">
                {Object.entries(selectedFormulation)
                  .filter(([key]) =>
                    !["id", "name", "components", "phaseMap", "notes", "attachments", "createdAt", "imageUrls", "totalmoles", "operator", "totalMass"].includes(key)
                  )
                  .map(([key, value], idx) => (
                    <div key={idx} className="flex flex-col mr-4 mb-2 min-w-[180px]">
                      <label className="font-semibold text-[#008080] text-xs uppercase mb-1">{key.replace(/_/g, " ")}:</label>
                      <input
                        type="text"
                        className="border border-gray-300 rounded p-2 text-[#002C36] bg-white min-w-[120px]"
                        value={typeof value === "string" || typeof value === "number" || typeof value === "boolean" ? String(value) : ""}
                        readOnly
                      />
                    </div>
                  ))}
              </div>
            </div>

            {/* Attachments with Click-to-View */}
            <div className="mb-4">
              <span className="text-xs font-bold uppercase text-[#008080] mb-1 tracking-wide">Attachments</span>
              <div className="mt-2 flex flex-wrap gap-4">
                {Object.entries(selectedFormulation.attachments || {}).length === 0 && (
                  <span className="text-[#002C36]">None</span>
                )}
                {Object.entries(selectedFormulation.attachments || {}).map(([key, { note, imageUrl }]: any, idx) => (
                  <div key={idx} className="flex items-center gap-2 border border-gray-200 rounded px-2 py-1 bg-gray-50">
                    <button
                      className="text-blue-600 underline hover:text-blue-800 cursor-pointer text-sm"
                      onClick={() => setSelectedAttachment({ name: key, data: { note, imageUrl } })}
                    >
                      {key.replace(/_/g, " ")}
                    </button>
                    {imageUrl && (
                      <img src={imageUrl} alt={key} className="w-8 h-8 object-cover rounded ml-2" />
                    )}
                  </div>
                ))}
              </div>
            </div>
            <div className="flex gap-2 mt-4">
              <button
                className="px-2 py-1 bg-[#008080] text-white rounded hover:bg-[#006666] font-bold uppercase tracking-wide text-xs"
                onClick={() => setShowAddTextField(true)}
              >
                + Add Text Field
              </button>
              <button
                className="px-2 py-1 bg-[#008080] text-white rounded hover:bg-[#006666] font-bold uppercase tracking-wide text-xs"
                onClick={() => setShowAddAttachmentField(true)}
              >
                + Add Attachment Field
              </button>
            </div>
          </div>
          {showAddTextField && (
            <div className="fixed inset-0 bg-black/50 z-50 flex items-center justify-center p-6" onClick={() => setShowAddTextField(false)}>
              <div className="bg-white rounded-lg shadow-xl w-full max-w-md p-8" onClick={e => e.stopPropagation()}>
                <h2 className="text-xl font-extrabold mb-4 text-[#002C36] uppercase tracking-wide">Add Text Field</h2>
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
                    onClick={async () => {
                      if (!newTextFieldName.trim()) return;

                      const updatedFormulation = { ...selectedFormulation, [newTextFieldName]: newTextFieldValue };

                      await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/update-formulation/${selectedFormulation.id}`, {
                        method: "POST",
                        headers: { "Content-Type": "application/json" },
                        body: JSON.stringify(updatedFormulation),
                      });

                      setFormulations((prev) =>
                        prev.map((f) => (f.id === selectedFormulation.id ? updatedFormulation : f))
                      );

                      setSelectedFormulation(updatedFormulation);
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
                <h2 className="text-xl font-extrabold mb-4 text-[#002C36] uppercase tracking-wide">Add Attachment Field</h2>
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
                    onClick={async () => {
                      if (!newAttachmentFieldName.trim()) return;

                      const updatedAttachments = {
                        ...(selectedFormulation.attachments || {}),
                        [newAttachmentFieldName]: { note: "", imageUrl: "" },
                      };

                      const updatedFormulation = { ...selectedFormulation, attachments: updatedAttachments };

                      await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/update-formulation/${selectedFormulation.id}`, {
                        method: "POST",
                        headers: { "Content-Type": "application/json" },
                        body: JSON.stringify(updatedFormulation),
                      });

                      setFormulations((prev) =>
                        prev.map((f) => (f.id === selectedFormulation.id ? updatedFormulation : f))
                      );

                      setSelectedFormulation(updatedFormulation);
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
            <div
              className="fixed inset-0 bg-black/50 z-50 flex items-center justify-center p-6"
              onClick={() => setSelectedAttachment(null)}
            >
              <div
                className="bg-white rounded-lg shadow-xl w-full max-w-2xl p-8"
                onClick={(e) => e.stopPropagation()}
              >
                <h2 className="text-xl font-bold mb-4 text-[#008080] uppercase">
                  {selectedAttachment.name.replace(/_/g, " ")}
                </h2>

                {selectedAttachment.data.imageUrl ? (
                  <img
                    src={selectedAttachment.data.imageUrl}
                    alt={selectedAttachment.name}
                    className="w-full max-h-[70vh] object-contain mb-4 border rounded shadow"
                  />
                ) : (
                  <p className="text-gray-500 mb-4">No image available.</p>
                )}

                {/* Note editing */}
                <textarea
                  className="w-full border border-gray-300 rounded p-2 mb-4 text-black"
                  placeholder="Attachment note"
                  value={selectedAttachment.data.note || ""}
                  onChange={e => {
                    setSelectedAttachment((prev) => prev && ({
                      ...prev,
                      data: { ...prev.data, note: e.target.value }
                    }));
                  }}
                />

                {/* Image upload */}
                <input
                  type="file"
                  accept="image/*"
                  className="mb-4"
                  onChange={async (e) => {
                    if (!e.target.files || e.target.files.length === 0) return;
                    const file = e.target.files[0];
                    const formData = new FormData();
                    formData.append('file', file);

                    try {
                      const res = await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/upload-image-to-firebase`, {
                        method: 'POST',
                        body: formData,
                      });
                      const data = await res.json();
                      // The backend returns fileUrl, not imageUrl
                      const imageUrl = data.fileUrl || data.imageUrl;
                      if (imageUrl) {
                        // Update image URL in attachment
                        const updatedAttachments = {
                          ...(selectedFormulation.attachments || {}),
                          [selectedAttachment.name]: {
                            ...selectedAttachment.data,
                            imageUrl,
                          },
                        };
                        const updatedFormulation = { ...selectedFormulation, attachments: updatedAttachments };
                        await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/update-formulation/${selectedFormulation.id}`, {
                          method: "POST",
                          headers: { "Content-Type": "application/json" },
                          body: JSON.stringify(updatedFormulation),
                        });
                        setFormulations((prev) =>
                          prev.map((f) => (f.id === selectedFormulation.id ? updatedFormulation : f))
                        );
                        setSelectedFormulation(updatedFormulation);
                        setSelectedAttachment({
                          name: selectedAttachment.name,
                          data: { ...selectedAttachment.data, imageUrl },
                        });
                      }
                    } catch (err) {
                      console.error("Image upload failed", err);
                      alert("Image upload failed.");
                    }
                  }}
                />

                <div className="flex justify-end gap-2">
                  <button
                    className="px-4 py-2 bg-gray-300 rounded hover:bg-gray-400 text-black"
                    onClick={async () => {
                      // Save note changes
                      const updatedAttachments = {
                        ...selectedFormulation.attachments,
                        [selectedAttachment.name]: {
                          ...selectedAttachment.data,
                          // Use the latest note
                          note: selectedAttachment.data.note,
                        },
                      };
                      const updatedFormulation = { ...selectedFormulation, attachments: updatedAttachments };
                      await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/update-formulation/${selectedFormulation.id}`, {
                        method: "POST",
                        headers: { "Content-Type": "application/json" },
                        body: JSON.stringify(updatedFormulation),
                      });
                      setFormulations((prev) =>
                        prev.map((f) => (f.id === selectedFormulation.id ? updatedFormulation : f))
                      );
                      setSelectedFormulation(updatedFormulation);
                      setSelectedAttachment(null);
                    }}
                  >
                    Save
                  </button>
                  <button
                    className="px-4 py-2 bg-[#008080] text-white rounded hover:bg-[#006666] font-bold uppercase tracking-wide"
                    onClick={() => setSelectedAttachment(null)}
                  >
                    Close
                  </button>
                </div>
              </div>
            </div>
          )}

          {showFullStructureModal && (
            <div className="fixed inset-0 bg-black/70 z-50 flex justify-center items-center" onClick={() => setShowFullStructureModal(null)}>
              <div className="bg-white p-8 rounded-2xl border-2 border-[#00E6D2] shadow-2xl max-w-4xl w-full max-h-[92vh] overflow-auto relative" onClick={e => e.stopPropagation()}>
                <div className="flex flex-col items-center">
                  {showFullStructureModal.imageUrl ? (
                    <img
                      src={showFullStructureModal.imageUrl}
                      alt={showFullStructureModal.name}
                      style={{ maxWidth: '100%', maxHeight: 560, objectFit: 'contain' }}
                    />
                  ) : (
                    <div className="text-gray-500">No image available.</div>
                  )}
                  <div className="text-lg font-bold text-[#002C36] mt-4">{showFullStructureModal.name}</div>
                  <button
                    className="mt-6 px-4 py-2 bg-[#00E6D2] text-[#002C36] rounded hover:bg-[#00bfae] font-bold"
                    onClick={() => setShowFullStructureModal(null)}
                  >
                    Close
                  </button>
                </div>
              </div>
            </div>
          )}
        </div>
      )}
      {editMode && (
          <div
            className="fixed inset-0 backdrop-blur-sm bg-opacity-50 z-50 flex items-center justify-center"
            onClick={() => setEditMode(false)}
          >
            <div
              className="bg-white p-6 rounded-lg shadow-lg w-full max-w-6xl max-h-[90vh] overflow-y-auto"
              onClick={(e) => e.stopPropagation()}
            >
              <h2 className="text-xl font-bold mb-4 text-black">Edit Formulation</h2>
              <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                {/* Name */}
                <div>
                  <label className="block text-sm font-semibold text-gray-700 mb-1">Name</label>
                  <input
                    type="text"
                    value={editData.name}
                    onChange={(e) => setEditData((d) => ({ ...d, name: e.target.value }))}
                    className="w-full mb-4 border px-2 py-1 rounded text-black"
                  />
                </div>
                {/* Total Desired Mass */}
                <div>
                  <label className="block text-sm font-semibold text-gray-700 mb-1">Total Desired Mass</label>
                  <input
                    type="number"
                    value={editData.totalMass !== undefined ? editData.totalMass : selectedFormulation.totalMass || ''}
                    onChange={(e) => setEditData((d) => ({ ...d, totalMass: e.target.value }))}
                    className="w-full mb-4 border px-2 py-1 rounded text-black"
                  />
                </div>
              </div>
              {/* Mass Calculator Table (editable desired mass % and actual mass) */}
              <div className="mb-6">
                <h3 className="text-lg font-semibold mb-2 text-[#008080]">Mass Calculator</h3>
                <div className="overflow-x-auto">
                  <table className="min-w-full border border-[#008080]">
                    <thead>
                      <tr className="text-xs text-[#008080] uppercase text-center">
                        <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Compound</th>
                        <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Image</th>
                        <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Lot</th>
                        <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Desired Mass %</th>
                        <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Actual Mass (mg)</th>
                        <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Target Mass (mg)</th>
                        <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Actual Mass %</th>
                        <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Actual Mol %</th>
                        <th className="px-2 py-2 font-bold border-b border-r border-[#008080]">Output Mol %</th>
                      </tr>
                    </thead>
                    <tbody>
                      {(editData.components || selectedFormulation.components)?.map((comp: any, idx: number) => (
                        <tr key={idx} className="align-middle text-center border-b border-[#008080]">
                          <td className="px-2 py-2 border-r border-[#008080] text-black">{comp.compoundName || comp.compoundId}</td>
                          <td className="px-2 py-2 border-r border-[#008080]">
                            {comp.imageUrl && typeof comp.imageUrl === 'string' && comp.imageUrl.trim() !== '' ? (
                              <img
                                src={comp.imageUrl}
                                alt={comp.compoundName || comp.compoundId}
                                style={{ maxWidth: '90px', maxHeight: '90px', objectFit: 'contain', background: 'white', borderRadius: '0.375rem', border: '1px solid #e5e7eb' }}
                              />
                            ) : null}
                          </td>
                          <td className="px-2 py-2 border-r border-[#008080] text-black">{comp.lotId || 'original'}</td>
                          <td className="px-2 py-2 border-r border-[#008080]">
                            <input
                              type="number"
                              value={editData.components?.[idx]?.massPercent !== undefined ? editData.components[idx].massPercent : comp.massPercent || ''}
                              onChange={e => {
                                const newComponents = [...(editData.components || selectedFormulation.components)];
                                newComponents[idx].massPercent = e.target.value;

                                const updatedComponents = recalculateComponents(newComponents, parseFloat(editData.totalMass || selectedFormulation.totalMass), compounds);

                                setEditData((d) => ({ ...d, components: updatedComponents }));
                              }}
                              className="w-20 border px-1 py-0.5 rounded text-black text-center"
                            />
                          </td>
                          <td className="px-2 py-2 border-r border-[#008080]">
                            <input
                              type="number"
                              value={editData.components?.[idx]?.actualMass !== undefined ? editData.components[idx].actualMass : comp.actualMass || ''}
                              onChange={e => {
                                const newComponents = [...(editData.components || selectedFormulation.components)];
                                if (!newComponents[idx]) newComponents[idx] = { ...comp };
                                newComponents[idx].actualMass = e.target.value;

                                const updatedComponents = recalculateComponents(newComponents, parseFloat(editData.totalMass || selectedFormulation.totalMass), compounds);

                                setEditData((d) => ({ ...d, components: updatedComponents }));
                              }}
                              className="w-20 border px-1 py-0.5 rounded text-black text-center"
                            />
                          </td>
                          <td className="px-2 py-2 border-r border-[#008080] text-black">
                            {comp.targetMass !== undefined ? comp.targetMass + " mg" : "-"}
                          </td>
                          <td className="px-2 py-2 border-r border-[#008080] text-black">
                            {comp.actualMassPercent !== undefined ? comp.actualMassPercent + "%" : "-"}
                          </td>
                          <td className="px-2 py-2 border-r border-[#008080] text-black">
                            {comp.actualMolPercent !== undefined ? comp.actualMolPercent + "%" : "-"}
                          </td>
                          <td className="px-2 py-2 border-r border-[#008080] text-black">
                            {comp.outputMolPercent !== undefined ? comp.outputMolPercent + "%" : "-"}
                          </td>

                        </tr>
                      ))}
                    </tbody>

                  </table>
                </div>
              </div>
              {/* Phase Map */}
              <label className="block text-sm font-semibold text-gray-700 mb-1">Phase Map</label>
              <textarea
                value={editData.phaseMap}
                onChange={(e) => setEditData((d) => ({ ...d, phaseMap: e.target.value }))}
                className="w-full mb-4 border px-2 py-1 rounded text-black"
              />
              {/* Analytical Notes */}
              <label className="block text-sm font-semibold text-gray-700 mb-1">Analytical Notes</label>
              <textarea
                value={editData.notes}
                onChange={(e) => setEditData((d) => ({ ...d, notes: e.target.value }))}
                className="w-full mb-4 border px-2 py-1 rounded text-black"
              />
              {/* Custom Fields */}
              <div className="mb-4">
                <h3 className="text-lg font-semibold text-gray-700 mb-2">Custom Fields</h3>
                <div className="max-h-64 overflow-y-auto border border-gray-300 rounded p-4 bg-white">
                  {Object.entries(selectedFormulation)
                    .filter(([key]) =>
                      !["id", "name", "components", "phaseMap", "notes", "attachments", "createdAt", "imageUrls", "totalmoles"].includes(key)
                    )
                    .map(([key, value], idx) => (
                      <div key={idx} className="mb-2">
                        <label className="block text-xs font-semibold text-[#008080] mb-1">{key.replace(/_/g, " ")}</label>
                        <input
                          type="text"
                          className="w-full border px-2 py-1 rounded text-[#002C36] bg-white"
                          value={editData[key] !== undefined ? editData[key] : (typeof value === "string" || typeof value === "number" || typeof value === "boolean" ? String(value) : "")}
                          onChange={e => setEditData((d) => ({ ...d, [key]: e.target.value }))}
                        />
                      </div>
                    ))}
                </div>
              </div>
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
                      // Merge updated fields and components
                      const updatedFormulation = {
                        ...selectedFormulation,
                        ...editData,
                        components: (editData.components || selectedFormulation.components).map((comp: any, idx: number) => ({
                          ...selectedFormulation.components[idx],
                          ...comp,
                        })),
                      };
                      await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/update-formulation/${selectedFormulation.id}`, {
                        method: "POST",
                        headers: { "Content-Type": "application/json" },
                        body: JSON.stringify(updatedFormulation),
                      });
                      setFormulations((prev) =>
                        prev.map((f: any) =>
                          f.id === selectedFormulation.id ? updatedFormulation : f
                        )
                      );
                      setSelectedFormulation(updatedFormulation);
                      setEditMode(false);
                    } catch (err) {
                      console.error("Failed to update formulation", err);
                      alert("Update failed");
                    }
                  }}
                  className="px-4 py-2 bg-[#008080] text-white rounded hover:bg-[#00E6D2]"
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
              await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/delete-compound/${id}`, { method: "DELETE" });
              setSelectedCompound(null);
            }}
            onUpdate={async (updatedCompound) => {
              const endpoint =
                compoundSource === "lot"
                  ? `${process.env.NEXT_PUBLIC_API_BASE_URL}/update-lot-compound`
                  : `${process.env.NEXT_PUBLIC_API_BASE_URL}/update-compound`;

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

    {showFormulationModal && (
      <CreateFormulationModal
        compounds={compounds}
        lots={lotMapping}
        onClose={() => setShowFormulationModal(false)}
        onCreate={async (data) => {
          try {
            const res = await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/create-formulation`, {
              method: "POST",
              headers: { "Content-Type": "application/json" },
              body: JSON.stringify(data),
            });
            const result = await res.json();
            if (!res.ok) throw new Error(result.error);
            alert("Formulation saved!");

            // Refresh formulation list
            const updatedFormulations = await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/formulations`).then((res) => res.json());
            setFormulations(updatedFormulations);
            setShowFormulationModal(false);
          } catch (err) {
            console.error("Failed to save formulation:", err);
            alert("Failed to save formulation.");
          }
        }}
      />
    )}

    {showDrawModal && (
      <DrawModal
        onClose={() => setShowDrawModal(false)}
        onSmilesSubmit={handleDrawFilterSubmit}
      />
    )}
    </div>
  );
}


/*
  Copyright © 2025 Polaris Electro Optics
  This code is the property of Polaris Electro Optics and may not be reused,
  modified, or distributed without explicit permission.
*/



