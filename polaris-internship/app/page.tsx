"use client";
import { useEffect, useState, useRef } from "react";
import CompoundCard from "../components/CompoundCard";
import CompoundModal from "../components/CompoundModal";
import FilterBar from "../components/FilterBar";
import DrawModal from "../components/DrawModal";
import AddStructureModal from "../components/AddStructureModal";
import CreateLotModal from "@/components/CreateLotModal";
import { Compound } from "@/types/compound";
import CreateFormulationModal, { Compound as CompoundType } from "../components/CreateFormulationModal";
import Link from "next/link";
import FormulationList from "@/components/FormulationList";



const transitionOptions = [
  "CR - FN", "CR - I", "CR - N",  "CR - SMA",
  "FN - CR", "FN - FNG", "FN - I", "FN - N", "FN - NX",
  "FNG - FN", "I - CR", "I - FN", "I - N",
  "N - CR", "N - FN", "N - I",  "N - SMA",
  "SMA - CR", "SMA - N", "SMA - SMX",
  "FN - GLASS", "I - GLASS", "N - GLASS", "CR - GLASS"
];


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
  const [showLotModal, setShowLotModal] = useState(false);
  const [showCreateLot, setShowCreateLot] = useState(false);
  const [lotList, setLotList] = useState<string[]>([]);
  const [selectedSource, setSelectedSource] = useState<"main" | "lot" | null>(null);
  const [currentLotId, setCurrentLotId] = useState<string | null>(null);
  const [showFormulationModal, setShowFormulationModal] = useState(false);
  const [showFilterDropdown, setShowFilterDropdown] = useState(false);
  const [searchName, setSearchName] = useState("");
  const filteredCompounds = compounds.filter((compound) =>
    searchName.trim() === "" ||
    compound.id?.toLowerCase().includes(searchName.toLowerCase())
  );
  const filterRef = useRef<HTMLDivElement>(null);
  const lotRef = useRef<HTMLDivElement>(null);
  const [lotMapping, setLotMapping] = useState<Record<string, string[]>>({});







  useEffect(() => {
    fetch("http://localhost:5000/compounds")
      .then((res) => res.json())
      .then((data) => setCompounds(data));
  }, []);

  useEffect(() => {
    fetch("http://localhost:5000/lots")
      .then((res) => res.json())
      .then(setLotList)
      .catch((err) => console.error("Failed to fetch lots:", err));
  }, []);

    useEffect(() => {
    function handleClickOutside(event: MouseEvent) {
      if (
        filterRef.current &&
        !filterRef.current.contains(event.target as Node)
      ) {
        setShowFilterDropdown(false);
      }
      if (
        lotRef.current &&
        !lotRef.current.contains(event.target as Node)
      ) {
        setShowLotModal(false);
      }
    }

    document.addEventListener("mousedown", handleClickOutside);
    return () => {
      document.removeEventListener("mousedown", handleClickOutside);
    };
  }, []);

  useEffect(() => {
    const fetchLotMapping = async () => {
      const compoundIds = compounds.map(c => c.id);
      try {
        const res = await fetch("http://localhost:5000/batch-lots-for-compounds", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ compound_ids: compoundIds }),
        });

        const data = await res.json();
        setLotMapping(data);
      } catch (err) {
        console.error("Failed to fetch lots in batch:", err);
      }
    };
    if (compounds.length > 0) {
      fetchLotMapping();
    }
  }, [compounds]);



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
    setSearchName("");

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

      
      <div className="mb-6">
        <Link href="/formulations">
          <button className="bg-gray-500 hover:bg-gray-600 text-white px-4 py-2 rounded cursor-pointer">
            📋 View Formulations
          </button>
        </Link>
      </div>


      <div className="flex flex-wrap gap-4 mb-4 items-center">
        <button
        onClick={() => setShowFormulationModal(true)}
        className="bg-orange-600 hover:bg-orange-700 text-white px-4 py-2 rounded transition-colors cursor-pointer"
      >
        🧪 Create Formulation
      </button>
        <div className="relative">
          <button
            onClick={() => setShowLotModal((prev) => !prev)}
            className="bg-purple-600 hover:bg-purple-700 text-white px-4 py-2 rounded transition-colors duration-200 cursor-pointer flex items-center gap-2"
          >
            📦 Lots
          </button>
          {showLotModal && (
            <div
              ref={lotRef}
              className="absolute z-10 bg-white shadow-md mt-2 rounded w-48 max-h-80 overflow-y-auto"
            >
              <div
                className="p-2 text-blue-600 hover:bg-gray-100 cursor-pointer"
                onClick={() => {
                  setShowLotModal(false);
                  setShowCreateLot(true);
                }}
              >
                ➕ Create New Lot
              </div>
              {lotList.length > 0 &&
                lotList.map((lot) => (
                  <div
                    key={lot}
                    className="p-2 hover:bg-gray-100 cursor-pointer text-black"
                    onClick={async () => {
                      setShowLotModal(false);
                      try {
                        const res = await fetch(`http://localhost:5000/lot/${lot}`);
                        const data = await res.json();
                        setCompounds(data);
                        setSelectedSource("lot");
                        setCurrentLotId(lot);
                      } catch (err) {
                        console.error("Failed to load lot compounds:", err);
                      }
                    }}
                  >
                    {lot}
                  </div>
                ))}
            </div>
          )}
        </div>

        <button
          onClick={() => setShowAddModal(true)}
          className="bg-green-600 hover:bg-green-700 text-white px-4 py-2 rounded transition-colors duration-200 cursor-pointer flex items-center gap-2"
        >
          ➕ Add Structure
        </button>
        <div className="relative">
          <button
            onClick={() => setShowFilterDropdown((prev) => !prev)}
            className="bg-blue-700 hover:bg-blue-800 text-white px-4 py-2 rounded cursor-pointer"
          >
            🔽 Filter
          </button>

          {showFilterDropdown && (
            <div ref={filterRef} className="absolute bg-white shadow-md border mt-2 p-4 rounded z-50 w-72 text-black">
              <div className="mb-3">
                <label className="block text-sm font-semibold mb-1">Search by Name</label>
                <input
                  type="text"
                  placeholder="e.g. PEO-0100"
                  value={searchName}
                  onChange={(e) => setSearchName(e.target.value)}
                  className="w-full border px-2 py-1 rounded text-black"
                />
              </div>
              <button
                onClick={() => {
                  setShowDrawModal(true);
                  setShowFilterDropdown(false);
                }}
                className="w-full text-left hover:bg-gray-100 px-2 py-1 rounded cursor-pointer"
              >
                ✏️ Draw Structure
              </button>

              <div className="mt-2">
                <label className="block text-sm font-semibold mb-1">Select Transition</label>
                <select
                  value={selectedTransition}
                  onChange={(e) => setSelectedTransition(e.target.value)}
                  className="w-full p-1 border border-gray-300 rounded text-black"
                >
                  <option value="">-- Select --</option>
                  {transitionOptions.map((option) => (
                    <option key={option} value={option}>{option}</option>
                  ))}
                </select>
              </div>

              {selectedTransition && (
                <input
                  type="text"
                  placeholder="Temperature (e.g. 87, 50+, 100-)"
                  value={temperature}
                  onChange={(e) => setTemperature(e.target.value)}
                  className="w-full mt-2 p-1 border border-gray-300 rounded text-black"
                />
              )}

              <button
                onClick={() => {
                  handlePhaseMapFilter();
                  setShowFilterDropdown(false);
                }}
                className="mt-3 w-full bg-indigo-600 hover:bg-indigo-700 text-white px-3 py-2 rounded cursor-pointer"
              >
                Apply Filter
              </button>
            </div>
          )}
        </div>

        <button
          onClick={handleResetFilters}
          className="bg-gray-300 hover:bg-gray-400 text-black px-4 py-2 rounded transition-colors duration-200 cursor-pointer"
        >
          🔄 Reset Filter
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
          {filteredCompounds.map((compound) => (
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
          source={selectedSource || "main"}
          lotId={currentLotId}
          onClose={() => {
            setSelectedCompound(null);
            setSelectedSource(null);
            setCurrentLotId(null);
          }}
          onDelete={async (id: string) => {
            try {
              const path = selectedSource === "lot" && currentLotId
                ? `http://localhost:5000/delete-lot-compound/${currentLotId}/${id}`
                : `http://localhost:5000/delete-compound/${id}`;

              await fetch(path, { method: "DELETE" });

              const updatedLots = await fetch("http://localhost:5000/lots").then((res) => res.json());
              setLotList(updatedLots);

              const updated = await fetch(
                selectedSource === "lot" && currentLotId
                  ? `http://localhost:5000/lot/${currentLotId}`
                  : "http://localhost:5000/compounds"
              ).then((res) => res.json());

              setCompounds(updated);
            } catch (err) {
              console.error("Failed to delete compound:", err);
            }
          }}
          onUpdate={async (updatedCompound: Compound) => {
            try {
              const path = selectedSource === "lot" && currentLotId
                ? `http://localhost:5000/update-lot-compound`
                : `http://localhost:5000/update-compound`;

              await fetch(path, {
                method: "POST",
                headers: { "Content-Type": "application/json" },
                body: JSON.stringify(
                  selectedSource === "lot"
                    ? { ...updatedCompound, lotId: currentLotId }
                    : updatedCompound
                ),
              });

              const updated = await fetch(
                selectedSource === "lot" && currentLotId
                  ? `http://localhost:5000/lot/${currentLotId}`
                  : "http://localhost:5000/compounds"
              ).then((res) => res.json());

              setCompounds(updated);
              setSelectedCompound(updatedCompound);
            } catch (err) {
              console.error("Failed to update compound:", err);
            }
          }}
          onUpdateCompoundFromLot={(compound, lotId) => {
            setSelectedCompound(compound);
            setCurrentLotId(lotId);
            setSelectedSource("lot");
          }}
        />
      )}

      {showCreateLot && (
        <CreateLotModal
          compounds={compounds}
          onClose={() => setShowCreateLot(false)}
          onCreate={() => {
            fetch("http://localhost:5000/lots")
              .then((res) => res.json())
              .then(setLotList);
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

      {showFormulationModal && (
        <CreateFormulationModal
          compounds={compounds}
          lots={lotMapping}  // ✅ now uses accurate compound-to-lots map
          onClose={() => setShowFormulationModal(false)}
          onCreate={async (data) => {
            try {
              const res = await fetch("http://localhost:5000/create-formulation", {
                method: "POST",
                headers: { "Content-Type": "application/json" },
                body: JSON.stringify(data),
              });

              const result = await res.json();
              if (!res.ok) throw new Error(result.error);
              alert("Formulation saved!");
            } catch (err) {
              console.error("Failed to save formulation:", err);
              alert("Failed to save formulation.");
            }
          }}
        />
      )}
    </div>
  );
}




















