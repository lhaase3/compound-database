"use client";
import { useEffect, useState, useRef } from "react";
import CompoundCard from "../components/CompoundCard";
import CompoundModal from "../components/CompoundModal";
import DrawModal from "../components/DrawModal";
import AddStructureModal from "../components/AddStructureModal";
import CreateLotModal from "@/components/CreateLotModal";
import { Compound } from "@/types/compound";
import CreateFormulationModal from "../components/CreateFormulationModal";
import Link from "next/link";



const transitionOptions = [
  "CR - FN", "CR - I", "CR - N",  "CR - SMA",
  "FN - CR", "FN - FNG", "FN - I", "FN - N", "FN - NX",
  "I - CR", "I - FN", "I - N",
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
  const [starred, setStarred] = useState<string[]>([]);
  const [showOnlyStarred, setShowOnlyStarred] = useState(false);
  const [selectedTag, setSelectedTag] = useState<string | null>(null);
  const [selectedForComparison, setSelectedForComparison] = useState<string[]>([]);
  const [showCompareModal, setShowCompareModal] = useState(false);
  const [compareAttachment, setCompareAttachment] = useState<{compoundId: string, key: string, data: any} | null>(null);
  const [mwRange, setMwRange] = useState<[number, number]>([0, 4000]);
  const filteredCompounds = compounds.filter((compound) => {
    // Name filter
    const matchesName =
      searchName.trim() === "" || compound.id?.toLowerCase().includes(searchName.toLowerCase());
    // Starred filter
    const matchesStar = !showOnlyStarred || starred.includes(compound.id);
    // Tag filter
    const matchesTag = !selectedTag || (compound.tags || []).includes(selectedTag);
    // MW filter
    let matchesMW = true;
    if (compound.MW === undefined || compound.MW === null || compound.MW === "") {
      // If MW is empty, only show if lower bound is 0 or 1
      matchesMW = mwRange[0] <= 1;
    } else {
      const mwNum = typeof compound.MW === "number" ? compound.MW : parseFloat(compound.MW);
      matchesMW = !isNaN(mwNum) && mwNum >= mwRange[0] && mwNum <= mwRange[1];
    }
    return matchesName && matchesStar && matchesTag && matchesMW;
});
  const filterRef = useRef<HTMLDivElement>(null);
  const lotRef = useRef<HTMLDivElement>(null);
  const [lotMapping, setLotMapping] = useState<Record<string, string[]>>({});
  const [showStickyLogo, setShowStickyLogo] = useState(false);
  const heroRef = useRef<HTMLDivElement>(null);


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

  useEffect(() => {
    const handleScroll = () => {
      if (!heroRef.current) return;
      const heroBottom = heroRef.current.getBoundingClientRect().bottom;
      setShowStickyLogo(heroBottom <= 0);
    };
    window.addEventListener("scroll", handleScroll, { passive: true });
    return () => window.removeEventListener("scroll", handleScroll);
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
    setSearchName("");
    setShowOnlyStarred(false);
    setSelectedTag(null);
    setMwRange([0, 4000]);

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

  // Load starred from localStorage on mount
  useEffect(() => {
    const stored = localStorage.getItem("starredCompounds");
    if (stored) setStarred(JSON.parse(stored));
  }, []);
  // Save starred to localStorage when changed
  useEffect(() => {
    localStorage.setItem("starredCompounds", JSON.stringify(starred));
  }, [starred]);

  // Toggle star for a compound
  const toggleStar = (id: string) => {
    setStarred((prev) =>
      prev.includes(id) ? prev.filter((sid) => sid !== id) : [...prev, id]
    );
  };

  // Add/remove compound from comparison
  const toggleCompare = (id: string) => {
    setSelectedForComparison((prev) =>
      prev.includes(id)
        ? prev.filter((cid) => cid !== id)
        : prev.length < 4
          ? [...prev, id]
          : prev // max 4
    );
  };
  const clearComparison = () => setSelectedForComparison([]);

  // Fast scroll-to-top function
  const fastScrollToTop = () => {
    const c = document.documentElement.scrollTop || document.body.scrollTop;
    if (c > 0) {
      window.scrollBy(0, -Math.max(120, Math.floor(c / 4)));
      setTimeout(fastScrollToTop, 4); // Lower timeout for even faster scroll
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
      <div
        ref={heroRef}
        className="w-full bg-gradient-to-r from-[#00343F] to-[#002C36] py-12 mb-10 shadow flex flex-col items-center relative overflow-hidden"
      >
        {/* Logo in top-left corner (hide when sticky logo is visible) */}
        <img
          src="/polaris-logo-only.png"
          alt="Polaris Electro-Optics Logo"
          className={`w-16 h-21 absolute top-6 left-8 z-20 drop-shadow-lg transition-opacity duration-300 ${
            showStickyLogo ? "opacity-0" : "opacity-100"
          }`}
        />
        {/* Optional: circuit/tech background effect */}
        <div className="absolute inset-0 opacity-30 pointer-events-none select-none" style={{background: ' center/cover no-repeat'}} />
        <h1 className="text-5xl font-extrabold mb-3 text-[#00E6D2] tracking-tight drop-shadow uppercase z-10 flex items-center gap-4">
          Compound Database
        </h1>
        <p className="text-xl text-white mb-6 max-w-2xl text-center z-10 font-semibold flex items-center justify-center gap-3">
          <img src="/white-logo.png" alt="Polaris Logo" className="w-8 h-10 inline-block" />
          Polaris Electro-Optics
        </p>
        <div className="mb-2 z-10">
          <Link href="/formulations">
            <button className="bg-[#00E6D2] hover:bg-[#00bfae] text-[#002C36] px-6 py-2 rounded-lg shadow font-bold text-lg uppercase tracking-wide flex items-center gap-2 transition-all">
              <span role="img" aria-label="formulations">📋</span> View Formulations
            </button>
          </Link>
        </div>
      </div>

      {/* Action Buttons Panel */}
      <div className="flex flex-wrap gap-5 mb-10 items-center justify-center w-full max-w-3xl px-2">
        <button
          onClick={() => setShowFormulationModal(true)
          }
          className="bg-[#00E6D2] hover:bg-[#00bfae] text-[#002C36] px-6 py-2 rounded-lg shadow font-bold uppercase tracking-wide flex items-center gap-2 transition-all"
        >
          <span role="img" aria-label="formulation">🧪</span> Create Formulation
        </button>
        <div className="relative">
          <button
            onClick={() => setShowLotModal((prev) => !prev)}
            className="bg-[#00343F] hover:bg-[#00545F] text-[#00E6D2] px-6 py-2 rounded-lg shadow font-bold uppercase tracking-wide flex items-center gap-2 transition-all border border-[#00E6D2]"
          >
            <span role="img" aria-label="lots">📦</span> Lots
          </button>
          {showLotModal && (
            <div
              ref={lotRef}
              className="absolute z-20 bg-[#002C36] shadow-lg mt-2 rounded-lg w-56 max-h-80 overflow-y-auto border border-[#00E6D2]"
            >
              <div
                className="p-3 text-[#00E6D2] hover:bg-[#00343F] cursor-pointer font-bold border-b border-[#00E6D2]"
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
                    className="p-3 hover:bg-[#00343F] cursor-pointer text-white border-b border-[#00E6D2] last:border-b-0"
                    onClick={async () => {
                      setShowLotModal(false);
                      try {
                        // Fetch all compounds for this lot
                        const res = await fetch(`http://localhost:5000/lot/${lot}`);
                        const lotCompounds = await res.json();
                        if (lotCompounds.length === 1) {
                          // If only one compound in lot, open modal directly for that lot compound
                          setSelectedCompound(lotCompounds[0]);
                          setSelectedSource("lot");
                          setCurrentLotId(lot);
                        } else if (lotCompounds.length > 1) {
                          // If multiple, prompt user to select which one
                          const options = lotCompounds.map((c, i) => `${i + 1}: ${c.id}`).join("\n");
                          const idx = window.prompt(`Select which lot compound to view (enter number):\n${options}`, "1");
                          const i = Number(idx) - 1;
                          if (!isNaN(i) && i >= 0 && i < lotCompounds.length) {
                            setSelectedCompound(lotCompounds[i]);
                            setSelectedSource("lot");
                            setCurrentLotId(lot);
                          }
                        }
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
          className="bg-[#00E6D2] hover:bg-[#00bfae] text-[#002C36] px-6 py-2 rounded-lg shadow font-bold uppercase tracking-wide flex items-center gap-2 transition-all"
        >
          <span role="img" aria-label="add">➕</span> Add Structure
        </button>
        <div className="relative">
          <button
            onClick={() => setShowFilterDropdown((prev) => !prev)}
            className="bg-[#00343F] hover:bg-[#00545F] text-[#00E6D2] px-6 py-2 rounded-lg shadow font-bold uppercase tracking-wide flex items-center gap-2 transition-all border border-[#00E6D2]"
          >
            <span role="img" aria-label="filter">🔎</span> Filter
          </button>
          {showFilterDropdown && (
            <div
              ref={filterRef}
              className="absolute z-20 bg-[#002C36] shadow-lg mt-2 rounded-lg w-80 max-h-96 overflow-y-auto border border-[#00E6D2] p-4"
            >
              {/* Search by Name */}
              <div className="mb-4">
                <label className="block text-sm font-bold mb-1 uppercase tracking-wide text-[#00E6D2]">Search by Name</label>
                <input
                  type="text"
                  placeholder="e.g. PEO-0100"
                  value={searchName}
                  onChange={(e) => setSearchName(e.target.value)}
                  className="w-full border border-[#00E6D2] bg-[#00343F] px-3 py-2 rounded text-[#00E6D2] focus:outline-none focus:ring-2 focus:ring-[#00E6D2] placeholder-[#00E6D2]/60"
                />
              </div>
              {/* MW Range Slider - dual sliders, top for min, bottom for max */}
              <div className="mb-4">
                <label className="block text-[#00E6D2] font-semibold mb-1">MW Range</label>
                <div className="flex flex-col gap-2">
                  {/* Top slider for min MW */}
                  <div className="flex items-center gap-2">
                    <span className="text-white text-xs w-10 text-right">{mwRange[0]}</span>
                    <input
                      type="range"
                      min={0}
                      max={mwRange[1]}
                      step={1}
                      value={mwRange[0]}
                      onChange={e => {
                        let val = Number(e.target.value);
                        if (val > mwRange[1]) val = mwRange[1];
                        setMwRange([val, mwRange[1]]);
                      }}
                      className="w-full accent-[#00E6D2] h-2 rounded-lg appearance-none bg-[#00343F]"
                      style={{ accentColor: '#00E6D2' }}
                    />
                    <span className="text-white text-xs w-10 text-right">{mwRange[1]}</span>
                  </div>
                  {/* Bottom slider for max MW */}
                  <div className="flex items-center gap-2">
                    <span className="text-white text-xs w-10 text-right">{mwRange[0]}</span>
                    <input
                      type="range"
                      min={mwRange[0]}
                      max={4000}
                      step={1}
                      value={mwRange[1]}
                      onChange={e => {
                        let val = Number(e.target.value);
                        if (val < mwRange[0]) val = mwRange[0];
                        setMwRange([mwRange[0], val]);
                      }}
                      className="w-full accent-[#00E6D2] h-2 rounded-lg appearance-none bg-[#00343F]"
                      style={{ accentColor: '#00E6D2' }}
                    />
                    <span className="text-white text-xs w-10 text-right">{mwRange[1]}</span>
                  </div>
                </div>
              </div>
              <div className="flex gap-4 mb-4">
                {["testing", "crystals"].map((tag) => (
                  <button
                    key={tag}
                    onClick={() => setSelectedTag(tag === selectedTag ? null : tag)}
                    className={`px-3 py-1 rounded text-sm font-bold border ${
                      selectedTag === tag
                        ? "bg-[#00E6D2] text-[#002C36] border-[#00E6D2]"
                        : "bg-[#00343F] text-[#00E6D2] border-[#00E6D2]"
                    } hover:bg-[#00545F] transition-all`}
                  >
                    {tag}
                  </button>
                ))}
              </div>
              <button
                onClick={() => {
                  setShowDrawModal(true);
                  setShowFilterDropdown(false);
                }}
                className="w-full text-left hover:bg-[#00343F] px-3 py-2 rounded cursor-pointer font-bold mb-2 text-[#00E6D2]"
              >
                ✏️ Draw Structure
              </button>
              <div className="mt-2">
                <label className="block text-sm font-bold mb-1 uppercase tracking-wide text-[#00E6D2]">Select Transition</label>
                <select
                  value={selectedTransition}
                  onChange={(e) => setSelectedTransition(e.target.value)}
                  className="w-full p-2 border border-[#00E6D2] rounded bg-[#00343F] text-[#00E6D2] focus:outline-none focus:ring-2 focus:ring-[#00E6D2]"
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
                  className="w-full mt-2 p-2 border border-[#00E6D2] rounded bg-[#00343F] text-[#00E6D2] focus:outline-none focus:ring-2 focus:ring-[#00E6D2] placeholder-[#00E6D2]/60"
                />
              )}
              <div className="mb-4 flex items-center gap-2 mt-4">
                <input
                  id="starred-filter"
                  type="checkbox"
                  checked={showOnlyStarred}
                  onChange={() => setShowOnlyStarred((v) => !v)}
                  className="accent-[#00E6D2] w-4 h-4"
                />
                <label htmlFor="starred-filter" className="text-sm font-bold cursor-pointer select-none flex items-center gap-1 text-[#00E6D2]">
                  <span className="text-[#00E6D2] text-lg">★</span> Show only starred compounds
                </label>
              </div>
              <button
                onClick={() => {
                  handlePhaseMapFilter();
                  setShowFilterDropdown(false);
                }}
                className="mt-4 w-full bg-[#00E6D2] hover:bg-[#00bfae] text-[#002C36] px-4 py-2 rounded-lg font-bold shadow uppercase tracking-wide"
              >
                Apply Filter
              </button>
            </div>
          )}
        </div>
        <button
          onClick={handleResetFilters}
          className="bg-[#00343F] hover:bg-[#00545F] text-[#00E6D2] px-6 py-2 rounded-lg shadow font-bold uppercase tracking-wide flex items-center gap-2 transition-all border border-[#00E6D2]"
        >
          <span role="img" aria-label="reset">🔄</span> Reset Filter
        </button>
      </div>

      {/* Section Divider */}
      <div className="w-full max-w-6xl flex items-center mb-6">
        <hr className="flex-grow border-t border-[#00E6D2]/40" />
        <span className="mx-4 text-lg text-[#00E6D2] font-bold uppercase tracking-wide">Compounds</span>
        <hr className="flex-grow border-t border-[#00E6D2]/40" />
      </div>

      {/* Compound Grid */}
      <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-3 lg:grid-cols-4 gap-10 mx-auto max-w-screen-xl px-6">
        {filteredCompounds.map((compound, index) => (
          <CompoundCard
            key={compound.id || `compound-${index}`}
            compound={compound}
            onMoreInfo={setSelectedCompound}
            isStarred={starred.includes(compound.id)}
            onToggleStar={() => toggleStar(compound.id)}
            cardColor="#00343F"
            accentColor="#00E6D2"
            compareChecked={selectedForComparison.includes(compound.id)}
            onToggleCompare={() => toggleCompare(compound.id)}
          />
        ))}
      </div>

      {/* Floating Compare Button */}
      {selectedForComparison.length >= 2 && (
        <button
          className="fixed bottom-8 right-8 z-50 bg-[#00E6D2] hover:bg-[#00bfae] text-[#002C36] px-8 py-4 rounded-full shadow-2xl font-bold uppercase tracking-wide text-lg flex items-center gap-3 border-4 border-[#008080] animate-bounce"
          onClick={() => setShowCompareModal(true)}
        >
          <span role="img" aria-label="compare">📊</span> Compare ({selectedForComparison.length})
        </button>
      )}

      {/* Compare Modal */}
      {showCompareModal && (
        <div className="fixed inset-0 bg-black/50 z-50 flex items-center justify-center p-6" onClick={() => setShowCompareModal(false)}>
          <div className="bg-white rounded-2xl shadow-2xl max-w-5xl w-full p-8 relative overflow-x-auto max-h-[90vh] flex flex-col" onClick={e => e.stopPropagation()}>
            <h2 className="text-2xl font-extrabold text-[#008080] uppercase mb-6 tracking-wide">Compare Compounds</h2>
            <button
              className="absolute top-4 right-4 text-2xl text-[#008080] font-bold hover:text-[#00bfae]"
              onClick={() => setShowCompareModal(false)}
            >
              &times;
            </button>
            <button
              className="absolute top-4 left-4 bg-gray-200 hover:bg-gray-300 text-[#002C36] font-bold px-4 py-2 rounded-md uppercase tracking-wide text-xs"
              onClick={clearComparison}
            >
              Clear
            </button>
            <div className="overflow-x-auto flex-1 max-h-[70vh]">
              <table className="min-w-full border border-[#008080] rounded-lg">
                <tbody>
                  {/* Display only the desired fields in order */}
                  {[
                    "id", "MW", "Lambda Max (DCM/AcCN)", "Lambda Max (neat film)",
                    "phase map", "r33", "dipole CAMB3LYP SVPD CHCl3 (Cosmo)",
                    "beta CAMB3LYP SVPD CHCl3 (Cosmo)", "dipole B3LYP SVPD CHCl3",
                    "beta B3LYP SVPD CHCl3", "beta/MW", "J/g DSC melt (total)",
                    "kJ/mol DSC melt (total)", "Refractive index (ne/no)", "Notes",
                    "lab?", "first PEO#", "registered PEO#", "Lab book #", "Max loading (%)"
                  ].filter(field => {
                    // Only show if at least one compound has a non-empty value for this field
                    return selectedForComparison.some(id => {
                      const cmp = compounds.find(c => c.id === id);
                      const val = cmp?.[field];
                      return val !== undefined && val !== null && val !== "" && val !== "N/A";
                    });
                  }).map((field) => (
                    <tr key={field}>
                      <td className="p-3 font-bold uppercase text-xs text-[#008080] bg-[#f8fafb]">{field}</td>
                      {selectedForComparison.map((id) => {
                        const cmp = compounds.find(c => c.id === id);
                        const val = cmp?.[field];
                        return <td key={id} className="p-3 text-center text-[#002C36]">{val !== undefined && val !== null && val !== "" && val !== "N/A" ? val : "N/A"}</td>;
                      })}
                    </tr>
                  ))}
                  {/* Attachments row (rendered only once, after all fields) */}
                  <tr>
                    <td className="p-3 font-bold uppercase text-xs text-[#008080] bg-[#f8fafb]">Attachments</td>
                    {selectedForComparison.map((id) => {
                      const cmp = compounds.find(c => c.id === id);
                      const atts = cmp?.attachments || {};
                      const keysWithData = Object.keys(atts).filter(k => atts[k] && (atts[k].note || atts[k].imageUrl));
                      return (
                        <td key={id} className="p-3 text-center text-[#002C36] flex flex-col gap-2 items-center justify-center">
                          {keysWithData.length === 0 && <span>N/A</span>}
                          {keysWithData.map((key) => (
                            <button
                              key={key}
                              className="px-2 py-1 bg-[#008080] text-white rounded hover:bg-[#006666] font-bold uppercase tracking-wide text-xs mb-1"
                              onClick={() => setCompareAttachment({ compoundId: id, key, data: atts[key] })}
                              type="button"
                            >
                              {key.replace(/_/g, " ")}
                            </button>
                          ))}
                        </td>
                      );
                    })}
                  </tr>
                </tbody>
              </table>
            </div>
            {/* Attachment Modal for Compare */}
            {compareAttachment && (
              <div className="fixed inset-0 bg-black/60 z-50 flex items-center justify-center p-6" onClick={() => setCompareAttachment(null)}>
                <div className="bg-white rounded-lg shadow-xl w-full max-w-7xl max-h-[95vh] p-10 relative flex flex-col items-center justify-center" onClick={e => e.stopPropagation()}>
                  <h2 className="text-2xl font-bold mb-6 text-[#008080] uppercase text-center w-full">{compareAttachment.key.replace(/_/g, " ")}</h2>
                  {compareAttachment.data.imageUrl ? (
                    <div className="flex justify-center items-center w-full" style={{ minHeight: '200px' }}>
                      <img
                        src={compareAttachment.data.imageUrl}
                        alt={compareAttachment.key}
                        className="w-auto h-auto max-h-[90vh] max-w-[1600px] object-contain mb-6 border rounded shadow bg-white mx-auto"
                        style={{ display: 'block' }}
                      />
                    </div>
                  ) : null}
                  {compareAttachment.data.note && (
                    <div className="bg-gray-100 p-4 rounded whitespace-pre-wrap mb-6 text-[#002C36] max-w-4xl w-full text-center mx-auto">
                      {compareAttachment.data.note}
                    </div>
                  )}
                  {!(compareAttachment.data.imageUrl || compareAttachment.data.note) && (
                    <div className="text-gray-500 mb-6">No data available.</div>
                  )}
                  <div className="flex justify-end w-full mt-2">
                    <button
                      className="px-8 py-4 bg-[#008080] text-white rounded hover:bg-[#006666] font-bold uppercase tracking-wide text-xl"
                      onClick={() => setCompareAttachment(null)}
                    >
                      Close
                    </button>
                  </div>
                </div>
              </div>
            )}
          </div>
        </div>
      )}
      {/* Modals */}
      {showDrawModal && (
        <DrawModal
          onClose={() => setShowDrawModal(false)}
          onSmilesSubmit={handleSmilesSubmit}
        />
      )}
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
              let freshCompound: Compound;
              if (selectedSource === "lot" && currentLotId) {
                const lotCompounds = await fetch(`http://localhost:5000/lot/${currentLotId}`).then(res => res.json());
                freshCompound = lotCompounds.find((c: Compound) => c.id === updatedCompound.id) || updatedCompound;
                setCompounds(lotCompounds);
                console.log("[DEBUG] Fetched lot compound after update:", freshCompound);
              } else {
                freshCompound = await fetch(`http://localhost:5000/compounds/${updatedCompound.id}`).then(res => res.json());
                const updatedList = await fetch("http://localhost:5000/compounds").then(res => res.json());
                setCompounds(updatedList);
                console.log("[DEBUG] Fetched main compound after update:", freshCompound);
              }
              setSelectedCompound(freshCompound);
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
              // Wait for backend to finish generating the imageUrl
              await new Promise(res => setTimeout(res, 1500));
              let updated = await fetch("http://localhost:5000/compounds").then((res) => res.json());
              // If imageUrl is still missing, poll a few more times
              let tries = 0;
              while (
                tries < 4 &&
                updated.some(c => c.id === newCompound.id && (!c.imageUrl || c.imageUrl === ""))
              ) {
                await new Promise(res => setTimeout(res, 1000));
                updated = await fetch("http://localhost:5000/compounds").then((res) => res.json());
                tries++;
              }
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
          lots={lotMapping}
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




















