"use client";
import { useEffect, useState } from "react";
import CompoundCard from "../components/CompoundCard";
import CompoundModal from "../components/CompoundModal";
import FilterBar from "../components/FilterBar";
import DrawModal from "../components/DrawModal";
import AddStructureModal from "../components/AddStructureModal";

const transitionOptions = [
  "CR - FN", "CR - I", "CR - N",  "CR - SMA",
  "FN - CR", "FN - FNG", "FN - I", "FN - N", "FN - NX",
  "FNG - FN", "I - CR", "I - FN", "I - N",
  "N - CR", "N - FN", "N - I",  "N - SMA",
  "SMA - CR", "SMA - N", "SMA - SMX",
  "FN - GLASS", "I - GLASS", "N - GLASS", "CR - GLASS"
];

type Compound = {
  id: string;
  name: string;
  formula: string;
  smiles: string;
  [key: string]: string | undefined;
};

type NewCompound = {
  id: string;
  name?: string;
  formula?: string;
  smiles: string;
  [key: string]: string | undefined;
};

export default function Home() {
  const [compounds, setCompounds] = useState<Compound[]>([]);
  const [selectedCompound, setSelectedCompound] = useState<Compound | null>(null);
  const [showDrawModal, setShowDrawModal] = useState(false);
  const [resetCount, setResetCount] = useState(0);
  const [showAddModal, setShowAddModal] = useState(false);
  const [selectedTransition, setSelectedTransition] = useState("");
  const [temperature, setTemperature] = useState("");

  useEffect(() => {
    fetch("http://localhost:5000/compounds")
      .then((res) => res.json())
      .then((data) => setCompounds(data));
  }, []);

  const handleSmilesSubmit = async (smiles: string) => {
    try {
      const res = await fetch("http://localhost:5000/search-substructure", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ query: smiles }),
      });

      const data = await res.json();
      setCompounds(Array.isArray(data) ? data : []);
    } catch (err) {
      console.error("Error filtering by drawn SMILES", err);
      setCompounds([]);
    }
  };

  const handleResetFilters = async () => {
    setResetCount((prev) => prev + 1);
    setSelectedTransition("");
    setTemperature("");

    try {
      const res = await fetch("http://localhost:5000/compounds");
      const data = await res.json();
      setCompounds(Array.isArray(data) ? data : []);
    } catch (err) {
      console.error("Error resetting filters", err);
      setCompounds([]);
    }
  };

  const handlePhaseMapFilter = async () => {
    if (!selectedTransition) return;

    try {
      const res = await fetch("http://localhost:5000/filter-phase-map", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          transition: selectedTransition,
          temperature: temperature || null,
        }),
      });

      const data = await res.json();
      setCompounds(Array.isArray(data) ? data : []);
    } catch (err) {
      console.error("Phase map filter failed", err);
      setCompounds([]);
    }
  };

  return (
    <div className="min-h-screen bg-white flex flex-col items-center p-5">
      <h1 className="text-3xl font-bold mb-6 text-black">Compound Database</h1>

      <div className="flex flex-wrap gap-4 mb-4 items-center">
        <button
          onClick={() => setShowAddModal(true)}
          className="bg-green-600 hover:bg-green-700 text-white px-4 py-2 rounded transition-colors duration-200 cursor-pointer flex items-center gap-2"
        >
          ➕ Add Structure
        </button>
        <button
          onClick={() => setShowDrawModal(true)}
          className="bg-blue-600 hover:bg-blue-700 text-white px-4 py-2 rounded transition-colors duration-200 cursor-pointer flex items-center gap-2"
        >
          ✏️ Draw Structure
        </button>

        <button
          onClick={handleResetFilters}
          className="bg-gray-300 hover:bg-gray-400 text-black px-4 py-2 rounded transition-colors duration-200 cursor-pointer"
        >
          🔄 Reset Filter
        </button>

        <select
          value={selectedTransition}
          onChange={(e) => setSelectedTransition(e.target.value)}
          className="p-2 border border-gray-300 rounded-md text-black"
        >
          <option value="">Select Transition</option>
          {transitionOptions.map((option) => (
            <option key={option} value={option}>{option}</option>
          ))}
        </select>

        {selectedTransition && (
          <input
            type="text"
            placeholder="Temperature (e.g. 87, 50+, 100-)"
            value={temperature}
            onChange={(e) => setTemperature(e.target.value)}
            className="p-2 border border-gray-300 rounded-md w-48 text-black"
          />
        )}

        <button
          onClick={handlePhaseMapFilter}
          className="bg-indigo-600 hover:bg-indigo-700 text-white px-4 py-2 rounded"
        >
          🔍 Filter
        </button>
      </div>

      <FilterBar onFilterResults={setCompounds} resetSignal={resetCount} />

      {showDrawModal && (
        <DrawModal
          onClose={() => setShowDrawModal(false)}
          onSmilesSubmit={handleSmilesSubmit}
        />
      )}

      <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-3 lg:grid-cols-4 gap-10 mx-auto max-w-screen-xl px-6">
        {Array.isArray(compounds) &&
          compounds.map((compound) => (
            <CompoundCard
              key={compound.id}
              compound={compound}
              onMoreInfo={setSelectedCompound}
            />
          ))}
      </div>

      {selectedCompound && (
        <CompoundModal
          compound={selectedCompound}
          onClose={() => setSelectedCompound(null)}
          onDelete={async (id: string) => {
            try {
              await fetch(`http://localhost:5000/delete-compound/${id}`, {
                method: "DELETE",
              });
              const updated = await fetch("http://localhost:5000/compounds").then((res) => res.json());
              setCompounds(updated);
            } catch (err) {
              console.error("Failed to delete compound:", err);
            }
          }}
          onUpdate={async (updatedCompound: Compound) => {
            try {
              await fetch(`http://localhost:5000/update-compound`, {
                method: "POST",
                headers: { "Content-Type": "application/json" },
                body: JSON.stringify(updatedCompound),
              });
              const updated = await fetch("http://localhost:5000/compounds").then((res) => res.json());
              setCompounds(updated);
              setSelectedCompound(updatedCompound);
            } catch (err) {
              console.error("Failed to update compound:", err);
            }
          }}
        />
      )}

      {showAddModal && (
        <AddStructureModal
          onClose={() => setShowAddModal(false)}
          onSubmit={async (newCompound: NewCompound) => {
            try {
              await fetch("http://localhost:5000/add-compound", {
                method: "POST",
                headers: { "Content-Type": "application/json" },
                body: JSON.stringify(newCompound),
              });

              const updated = await fetch("http://localhost:5000/compounds").then((res) => res.json());
              setCompounds(updated);
            } catch (err) {
              console.error("Failed to add compound:", err);
            }
          }}
        />
      )}
    </div>
  );
}





