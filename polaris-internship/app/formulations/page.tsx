"use client"

import React, { useEffect, useState, useRef } from "react";
import Link from "next/link";
import CompoundModal from "@/components/CompoundModal";



function AttachmentViewer({ name, note, imageUrl }: { name: string, note: string, imageUrl: string }) {
  const [showImage, setShowImage] = useState(false);

  return (
    <div className="mb-4">
      <p className="font-semibold">{name}</p>
      {note && <p className="text-gray-700 mb-2">Note: {note}</p>}
      {imageUrl && (
        <p
          className="text-blue-600 underline hover:text-blue-800 cursor-pointer"
          onClick={() => setShowImage((prev) => !prev)}
        >
          {showImage ? "Hide Image" : "View Image"}
        </p>
      )}
      {imageUrl && showImage && (
        <img
          src={imageUrl}
          alt={name}
          className="w-160 h-100 object-contain border rounded shadow mt-2"
        />
      )}
    </div>
  );
}



export default function FormulationList() {
  const [formulations, setFormulations] = useState<any[]>([]);
  const [selectedFormulation, setSelectedFormulation] = useState<any | null>(null);
  const [editMode, setEditMode] = useState(false);
  const [editData, setEditData] = useState({ name: "", phaseMap: "", notes: "" });
  const [selectedCompound, setSelectedCompound] = useState<any | null>(null);
  const [compoundSource, setCompoundSource] = useState<"main" | "lot">("main");
  const [compoundLotId, setCompoundLotId] = useState<string | null>(null);
  const [searchTerm, setSearchTerm] = useState("");
  const [showStickyLogo, setShowStickyLogo] = useState(false);
  const heroRef = useRef<HTMLDivElement>(null);
  const [selectedAttachment, setSelectedAttachment] = useState<{
  name: string;
  data: { note: string; imageUrl: string };
} | null>(null);




  useEffect(() => {
    fetch("http://localhost:5000/formulations")
      .then((res) => res.json())
      .then(setFormulations)
      .catch((err) => console.error("Failed to fetch formulations", err));
  }, []);

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
    const c = document.documentElement.scrollTop || document.body.scrollTop;
    if (c > 0) {
      window.scrollBy(0, -Math.max(120, Math.floor(c / 4)));
      setTimeout(fastScrollToTop, 4);
    }
  };

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
        <h1 className="text-5xl font-extrabold mb-3 text-[#00E6D2] tracking-tight drop-shadow uppercase z-10">Formulations</h1>
        <p className="text-xl text-white mb-6 max-w-2xl text-center z-10 font-semibold">Polaris Electro-Optics</p>
        <div className="mb-2 z-10">
          <Link href="/">
            <button className="bg-[#00E6D2] hover:bg-[#00bfae] text-[#002C36] px-6 py-2 rounded-lg shadow font-bold text-lg uppercase tracking-wide flex items-center gap-2 transition-all">
              <span role="img" aria-label="home">🏠</span> Back to Home
            </button>
          </Link>
        </div>
      </div>

      {/* Search Bar */}
      <div className="mb-8 w-full max-w-md">
        <input
          type="text"
          placeholder="🔍 Search by compound ID (e.g. PEO-0100)"
          value={searchTerm}
          onChange={(e) => setSearchTerm(e.target.value)}
          className="w-full border border-[#00E6D2] bg-[#00343F] px-3 py-2 rounded text-[#00E6D2] focus:outline-none focus:ring-2 focus:ring-[#00E6D2] placeholder-[#00E6D2]/60"
        />
      </div>

      {/* Formulations Grid */}
    <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-8 w-full max-w-6xl px-4">
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
            className="bg-[#e6f9f7] border border-[#00E6D2] rounded-2xl p-6 shadow-lg hover:shadow-2xl cursor-pointer transition-all text-[#002C36]"
            onClick={() => setSelectedFormulation(form)}
          >
            <h2 className="text-2xl font-bold mb-2 uppercase tracking-wide text-[#00E6D2]">{form.name || "Unnamed Formulation"}</h2>
            <p className="text-sm mb-1 font-semibold">Components:</p>
            <ul className="list-disc pl-5 text-sm">
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
                    await fetch(`http://localhost:5000/delete-formulation/${selectedFormulation.id}`, { method: "DELETE" });
                    setFormulations((prev) => prev.filter((f) => f.id !== selectedFormulation.id));
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

          {/* Custom Fields */}
          <div className="mb-4">
            <h3 className="text-lg font-semibold">Other Data</h3>
            {Object.entries(selectedFormulation)
              .filter(([key]) =>
                !["id", "name", "components", "phaseMap", "notes", "attachments", "createdAt", "imageUrls"].includes(key)
              )
              .map(([key, value], idx) => (
                <div key={idx} className="mb-2">
                  <span className="font-semibold">{key}:</span> {value || "N/A"}
                </div>
              ))}
          </div>

          {/* Attachments with Click-to-View */}
          <div className="mb-4">
            <h3 className="text-lg font-semibold">Attachments</h3>
              {Object.entries(selectedFormulation.attachments || {}).map(([key, { note, imageUrl }]: any, idx) => (
                <div key={idx} className="mb-4">
                  <p
                    className="text-blue-600 underline hover:text-blue-800 cursor-pointer"
                    onClick={() => setSelectedAttachment({ name: key, data: { note, imageUrl } })}
                  >
                    {key.replace(/_/g, " ")}
                  </p>
                </div>
              ))}

          </div>
        </div>
        {selectedAttachment && (
  <div className="fixed inset-0 bg-black/50 z-50 flex items-center justify-center p-6" onClick={() => setSelectedAttachment(null)}>
    <div className="bg-white rounded-lg shadow-xl w-full max-w-2xl p-8" onClick={(e) => e.stopPropagation()}>
      <h2 className="text-xl font-bold mb-4">{selectedAttachment.name.replace(/_/g, " ")}</h2>
      {selectedAttachment.data.imageUrl ? (
        <img
          src={selectedAttachment.data.imageUrl}
          alt={selectedAttachment.name}
          className="w-full max-h-[70vh] object-contain mb-4 border rounded shadow"
        />
      ) : (
        <p className="text-gray-500 mb-4">No image available.</p>
      )}
      {selectedAttachment.data.note && (
        <div className="bg-gray-100 p-3 rounded whitespace-pre-wrap mb-4">
          {selectedAttachment.data.note}
        </div>
      )}
      <div className="flex justify-end">
        <button
          className="px-4 py-2 bg-blue-600 text-white rounded hover:bg-blue-700"
          onClick={() => setSelectedAttachment(null)}
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
            <div className="mb-4">
              <h3 className="text-lg font-semibold">Custom Fields</h3>
              {Object.entries(selectedFormulation)
                .filter(([key]) =>
                  !["id", "name", "components", "phaseMap", "notes", "attachments", "createdAt", "imageUrls"].includes(key)
                )
                .map(([key, value], idx) => (
                  <div key={idx} className="mb-2">
                    <span className="font-semibold">{key}:</span> {value || "N/A"}
                  </div>
                ))}
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




