"use client";
import { useRouter } from "next/navigation";
import { useEffect, useState, useRef } from "react";
import { auth } from "@/utils/firebase";
import { onAuthStateChanged, signOut, User as FirebaseUser } from "firebase/auth";
import CompoundCard from "../components/CompoundCard";
import CompoundModal from "../components/CompoundModal";
import DrawModal from "../components/DrawModal";
import AddStructureModal from "../components/AddStructureModal";
import CreateLotModal from "@/components/CreateLotModal";
import { Compound } from "@/types/compound";
import CreateFormulationModal from "../components/CreateFormulationModal";
import Link from "next/link";
import { usePathname } from "next/navigation";
import CompareModal from "../components/CompareModal";



const transitionOptions = [
  "CR - FN", "CR - I", "CR - N",  "CR - SMA",
  "FN - CR", "FN - FNG", "FN - I", "FN - N", "FN - NX",
  "I - CR", "I - FN", "I - N",
  "N - CR", "N - FN", "N - I",  "N - SMA",
  "SMA - CR", "SMA - N", "SMA - SMX",
  "FN - GLASS", "I - GLASS", "N - GLASS", "CR - GLASS"
];

// Fields to show in table view (similar to CompareModal)
const tableFieldsToShow = [
  "Lambda Max (DCM/AcCN)", "Lambda Max (neat film)",
  "phase map", "r33", "dipole CAMB3LYP SVPD CHCl3 (Cosmo)",
  "beta CAMB3LYP SVPD CHCl3 (Cosmo)", "dipole B3LYP SVPD CHCl3",
  "beta B3LYP SVPD CHCl3", "beta/MW", "J/g DSC melt (total)",
  "kJ/mol DSC melt (total)", "Refractive index (ne/no)", "Notes",
  "lab?", "first PEO#", "registered PEO#", "Lab book #", "Max loading (%)"
];

const IMAGE_COL_WIDTH = 420;
const ID_COL_WIDTH = 180;
const MW_COL_WIDTH = 120;


type NewCompound = {
  id: string;
  name?: string;
  formula?: string;
  smiles: string;
  [key: string]: string | undefined;
};

export default function Home() {
  const router = useRouter();
  const [user, setUser] = useState<FirebaseUser | null>(null);
  // Listen for auth state changes
  useEffect(() => {
    const unsubscribe = onAuthStateChanged(auth, (firebaseUser) => {
      if (firebaseUser && firebaseUser.email && firebaseUser.email.endsWith("@polariseo.com")) {
        setUser(firebaseUser);
      } else {
        setUser(null);
        router.push("/login"); // Redirect to login page if not authenticated
      }
    });
    return () => unsubscribe();
  }, [router]);
  const [compounds, setCompounds] = useState<Compound[]>([]);
  const [selectedCompound, setSelectedCompound] = useState<Compound | null>(null);
  const [showDrawModal, setShowDrawModal] = useState(false);
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
  const [viewMode, setViewMode] = useState<'table' | 'cards'>('table');
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


  // On mount, fetch starred from localStorage and compounds from backend
  useEffect(() => {
    const stored = localStorage.getItem("starredCompounds");
    const starredList = stored ? JSON.parse(stored) : [];
    setStarred(starredList);
    fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/compounds`)
      .then((res) => res.json())
      .then((data) => {
        setCompounds(
          Array.isArray(data)
            ? data.map((compound: Compound) => ({
                ...compound,
                isStarred: starredList.includes(compound.id),
              }))
            : []
        );
      })
      .catch((err) => console.error("Failed to fetch compounds:", err));
  }, []);


  useEffect(() => {
    fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/lots`)
      .then((res) => res.json())
      .then(setLotList)
      .catch((err) => console.error("Failed to fetch lots:", err));
  }, []);

  useEffect(() => {
    const targetId = localStorage.getItem("compoundToOpen");
    if (targetId && compounds.length > 0) {  // ✅ Wait until compounds are loaded
      const targetCompound = compounds.find((c) => c.id === targetId);
      if (targetCompound) {
        setSelectedCompound(targetCompound);
        localStorage.removeItem("compoundToOpen");
      }
    }
  }, [compounds]);







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
        const res = await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/batch-lots-for-compounds`, {
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
      const res = await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/search-substructure`, {
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

    useEffect(() => {
    if (user && user.email) {
      fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/get-starred?email=${user.email}`)
        .then(res => res.json())
        .then(data => {
          if (Array.isArray(data.starred)) {
            setStarred(data.starred);
            localStorage.setItem("starredCompounds", JSON.stringify(data.starred));
            setCompounds(prev =>
              prev.map(compound => ({
                ...compound,
                isStarred: data.starred.includes(compound.id),
              }))
            );
          }
        })
        .catch(err => console.error("Failed to fetch user-starred compounds:", err));
    }
  }, [user]);


  const pathname = usePathname();

  useEffect(() => {
    if (pathname === "/") {
      const stored = localStorage.getItem("starredCompounds");
      const starredList = stored ? JSON.parse(stored) : [];
      setStarred(starredList);
      setCompounds((prev) =>
        prev.map((compound) => ({
          ...compound,
          isStarred: starredList.includes(compound.id),
        }))
      );
    }
  }, [pathname]);


  const handleResetFilters = async () => {
    setSelectedTransition("");
    setTemperature("");
    setSearchName("");
    setShowOnlyStarred(false);
    setSelectedTag(null);
    setMwRange([0, 4000]);

    try {
      const res = await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/compounds`);
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
      const res = await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/filter-phase-map`, {
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

  // Whenever starred or compounds change, enrich compounds and sync localStorage and backend
  useEffect(() => {
    localStorage.setItem("starredCompounds", JSON.stringify(starred));
    setCompounds((prev) =>
      prev.map((compound) => ({
        ...compound,
        isStarred: starred.includes(compound.id),
      }))
    );
    // Sync starred compounds to backend if user is logged in
    if (user && user.email) {
      fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/save-starred`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ email: user.email, starred }),
      }).catch((err) => console.error("Failed to save starred compounds:", err));
    }
  }, [starred, user]);


  // Toggle star for a compound
  // const toggleStar = (id: string) => {
  //   setStarred((prev) =>
  //     prev.includes(id) ? prev.filter((sid) => sid !== id) : [...prev, id]
  //   );
  // };

  const toggleStar = (id: string) => {
    setStarred((prev) =>
      prev.includes(id)
        ? prev.filter((sid) => sid !== id)
        : [...prev, id]
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

  useEffect(() => {
  const handleVisibilityChange = () => {
    if (document.visibilityState === "visible") {
      const stored = localStorage.getItem("starredCompounds");
      if (stored) {
        const starredList = JSON.parse(stored);
        setStarred(starredList);
        setCompounds((prev) =>
          prev.map((compound) => ({
            ...compound,
            isStarred: starredList.includes(compound.id),
          }))
        );
      }
    }
  };

  document.addEventListener("visibilitychange", handleVisibilityChange);
  return () => {
    document.removeEventListener("visibilitychange", handleVisibilityChange);
  };
}, []);

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
          onClick={() => router.push("/login")}
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
      <div
        ref={heroRef}
        className="w-full bg-gradient-to-r from-[#00343F] to-[#002C36] py-12 mb-10 shadow flex flex-col items-center relative overflow-hidden"
      >
        {/* Logo in top-left corner (hide when sticky logo is visible) */}
        <img
          src="/polaris-logo-only.PNG"
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
          <img src="/white-logo.PNG" alt="Polaris Logo" className="w-8 h-10 inline-block" />
          Polaris Electro-Optics
        </p>
        <div className="mb-2 z-10 flex flex-row gap-4 justify-center">
          <Link href="/formulations">
            <button className="bg-[#00E6D2] hover:bg-[#00bfae] text-[#002C36] px-6 py-2 rounded-lg shadow font-bold text-lg uppercase tracking-wide flex items-center gap-2 transition-all">
              <span role="img" aria-label="formulations"></span> View Formulations
            </button>
          </Link>
          <Link href="/plans">
            <button className="bg-[#00E6D2] hover:bg-[#00bfae] text-[#002C36] px-6 py-2 rounded-lg shadow font-bold text-lg uppercase tracking-wide flex items-center gap-2 transition-all">
              <span role="img" aria-label="plans"></span> View Plans
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
                lotList.map((lot: string) => (
                  <div
                    key={lot}
                    className="p-3 hover:bg-[#00343F] cursor-pointer text-white border-b border-[#00E6D2] last:border-b-0"
                    onClick={async () => {
                      setShowLotModal(false);
                      try {
                        // Fetch all compounds for this lot
                        const res = await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/lot/${lot}`);
                        const lotCompounds = await res.json();
                        if (lotCompounds.length === 1) {
                          // If only one compound in lot, open modal directly for that lot compound
                          setSelectedCompound(lotCompounds[0]);
                          setSelectedSource("lot");
                          setCurrentLotId(lot);
                        } else if (lotCompounds.length > 1) {
                          // If multiple, prompt user to select which one
                          const options = lotCompounds.map((c: any, i: number) => `${i + 1}: ${c.id}`).join("\n");
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
              className="absolute z-[100] bg-[#002C36] shadow-lg mt-2 rounded-lg w-80 max-h-96 overflow-y-auto border border-[#00E6D2] p-4"
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

        {/* View Toggle */}
        <div className="ml-4 flex bg-[#00343F] rounded-lg border border-[#00E6D2] overflow-hidden">
          <button
            onClick={() => setViewMode('table')}
            className={`px-4 py-2 text-sm font-bold transition-all ${
              viewMode === 'table'
                ? 'bg-[#00E6D2] text-[#002C36]'
                : 'bg-transparent text-[#00E6D2] hover:bg-[#00545F]'
            }`}
          >
            📊 Table
          </button>
          <button
            onClick={() => setViewMode('cards')}
            className={`px-4 py-2 text-sm font-bold transition-all ${
              viewMode === 'cards'
                ? 'bg-[#00E6D2] text-[#002C36]'
                : 'bg-transparent text-[#00E6D2] hover:bg-[#00545F]'
            }`}
          >
            📋 Cards
          </button>
        </div>
      </div>

      {/* Compound Display - Cards or Table */}
      {viewMode === 'cards' ? (
        /* Compound Grid */
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
              similarity={undefined}
            />
          ))}
        </div>
      ) : (
        /* Compound Table */
        <div className="w-full max-w-[95vw] mx-auto px-4">
          <div className="overflow-x-auto overflow-y-auto bg-white rounded-lg border border-[#008080] shadow-2xl" style={{ maxHeight: '90vh' }}>
            <table className="min-w-full border-collapse">
              <thead>
                <tr>
                  <th
                    className="p-3 text-xs font-bold uppercase text-[#008080] border-b border-r border-[#008080] sticky bg-white"
                    style={{
                      left: 0,
                      top: 0,
                      width: IMAGE_COL_WIDTH,
                      minWidth: IMAGE_COL_WIDTH,
                      boxShadow: '2px 0 0 #008080, 1px 0 0 #008080 inset',
                      zIndex: 40,
                    }}
                  >
                    Image
                  </th>
                  <th
                    className="p-3 text-xs font-bold uppercase text-[#008080] border-b border-r border-[#008080] sticky bg-white"
                    style={{
                      left: IMAGE_COL_WIDTH,
                      top: 0,
                      width: ID_COL_WIDTH,
                      minWidth: ID_COL_WIDTH,
                      boxShadow: '2px 0 0 #008080, 1px 0 0 #008080 inset',
                      zIndex: 40,
                    }}
                  >
                    ID
                  </th>
                  <th
                    className="p-3 text-xs font-bold uppercase text-[#008080] border-b border-r border-[#008080] sticky bg-white relative"
                    style={{
                      left: IMAGE_COL_WIDTH + ID_COL_WIDTH,
                      top: 0,
                      width: MW_COL_WIDTH,
                      minWidth: MW_COL_WIDTH,
                      boxShadow: '3px 0 0 #008080, 1px 0 0 #008080 inset',
                      zIndex: 40,
                      borderRight: '2px solid #008080 !important',
                    }}
                  >
                    MW
                    <div 
                      style={{
                        position: 'absolute',
                        right: 0,
                        top: 0,
                        bottom: 0,
                        width: '2px',
                        backgroundColor: '#008080',
                        zIndex: 1,
                      }}
                    />
                  </th>
                  {tableFieldsToShow.map(field => (
                    <th
                      key={field}
                      className={
                        `p-3 text-xs font-bold uppercase text-[#008080] border-b border-r-2 border-[#008080] bg-white whitespace-nowrap sticky` +
                        (field === 'r33' || field === 'R33' ? ' min-w-[60px] max-w-[80px] w-[70px] text-center' : '')
                      }
                      style={{
                        top: 0,
                        zIndex: 30,
                        ...(field === 'r33' || field === 'R33' ? { minWidth: 60, maxWidth: 80, width: 70, textAlign: 'center' } : {})
                      }}
                    >
                      {field}
                    </th>
                  ))}
                  <th className="p-3 text-xs font-bold uppercase text-[#008080] border-b border-r border-[#008080] bg-white whitespace-nowrap sticky" style={{ top: 0, zIndex: 30 }}>
                    Attachments
                  </th>
                  <th className="p-3 text-xs font-bold uppercase text-[#008080] border-b border-[#008080] bg-white whitespace-nowrap sticky" style={{ top: 0, zIndex: 30 }}>
                    Actions
                  </th>
                </tr>
              </thead>
              <tbody>
                {filteredCompounds.map((compound, index) => (
                  <tr key={compound.id || `compound-${index}`} className="border-b border-[#008080] hover:bg-[#f8fffe] cursor-pointer group" onClick={() => setSelectedCompound(compound)}>
                    <td
                      className="p-3 text-center sticky z-20 bg-white group-hover:bg-[#f8fffe] border-r border-[#008080]"
                      style={{
                        left: 0,
                        width: IMAGE_COL_WIDTH,
                        minWidth: IMAGE_COL_WIDTH,
                        maxWidth: IMAGE_COL_WIDTH,
                        boxShadow: '2px 0 0 #008080, 1px 0 0 #008080 inset',
                        zIndex: 20,
                        overflow: 'visible',
                        padding: 0,
                      }}
                    >
                      <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'center', height: '100%', padding: 24 }}>
                        {compound.imageUrl ? (
                          <img
                            src={compound.imageUrl}
                            alt={compound.id}
                            style={{ maxWidth: '370px', maxHeight: '370px', objectFit: 'contain', display: 'block', margin: '0 auto', background: '#fff', borderRadius: '0.75rem', boxShadow: '0 2px 12px 0 rgba(0,0,0,0.07)' }}
                          />
                        ) : (
                          <span className="text-gray-400">No Image</span>
                        )}
                      </div>
                    </td>
                    <td
                      className="p-3 font-bold text-[#002C36] text-center sticky z-20 bg-white group-hover:bg-[#f8fffe] border-r border-[#008080]"
                      style={{
                        left: IMAGE_COL_WIDTH,
                        width: ID_COL_WIDTH,
                        minWidth: ID_COL_WIDTH,
                        maxWidth: ID_COL_WIDTH,
                        boxShadow: '2px 0 0 #008080, 1px 0 0 #008080 inset',
                        zIndex: 20,
                        overflow: 'visible',
                      }}
                    >
                      {compound.id}
                    </td>
                    <td
                      className="p-3 font-semibold text-[#002C36] text-center sticky z-20 bg-white group-hover:bg-[#f8fffe] border-r border-[#008080] relative"
                      style={{
                        left: IMAGE_COL_WIDTH + ID_COL_WIDTH,
                        width: MW_COL_WIDTH,
                        minWidth: MW_COL_WIDTH,
                        maxWidth: MW_COL_WIDTH,
                        boxShadow: '3px 0 0 #008080, 1px 0 0 #008080 inset',
                        zIndex: 20,
                        overflow: 'visible',
                        borderRight: '2px solid #008080 !important',
                      }}
                    >
                      {compound.MW ?? "N/A"}
                      <div 
                        style={{
                          position: 'absolute',
                          right: 0,
                          top: 0,
                          bottom: 0,
                          width: '2px',
                          backgroundColor: '#008080',
                          zIndex: 1,
                        }}
                      />
                    </td>
                    {tableFieldsToShow.map(field => (
                      <td
                        key={field}
                        className={
                          `p-3 text-[#002C36] text-center border-r-2 border-[#008080]` +
                          ((field === 'phase map' || field === 'Notes') ? ' whitespace-pre-line break-words max-w-xs' : ' whitespace-nowrap') +
                          (field === 'r33' || field === 'R33' ? ' min-w-[60px] max-w-[80px] w-[70px] text-center' : '')
                        }
                        style={
                          field === 'r33' || field === 'R33'
                            ? { minWidth: 60, maxWidth: 80, width: 70, textAlign: 'center' }
                            : (field === 'phase map' || field === 'Notes')
                              ? { whiteSpace: 'pre-line', wordBreak: 'break-word', minWidth: 300, maxWidth: 600 }
                              : undefined
                        }
                      >
                        {(!compound[field] || compound[field] === 'N/A') ? (
                          <span className="text-gray-300">N/A</span>
                        ) : (
                          compound[field]
                        )}
                      </td>
                    ))}
                    <td className="p-3 text-center border-r border-[#008080]" onClick={e => e.stopPropagation()}>
                      <div className="flex flex-col gap-2 items-center">
                        {(() => {
                          const atts = compound.attachments || {};
                          const keysWithData = Object.keys(atts).filter(k => {
                            const att = atts[k];
                            if (Array.isArray(att)) {
                              return att.some(entry => entry && (entry.note || entry.imageUrl));
                            } else {
                              return att && (att.note || att.imageUrl);
                            }
                          });

                          if (keysWithData.length === 0) {
                            return <span>N/A</span>;
                          }

                          return keysWithData.map((key) => {
                            const att = atts[key];
                            if (Array.isArray(att)) {
                              return att.map((entry, idx) =>
                                entry && (entry.note || entry.imageUrl) ? (
                                  <button
                                    key={key + "-" + idx}
                                    className="px-2 py-1 bg-[#008080] text-white rounded hover:bg-[#006666] font-bold uppercase tracking-wide text-xs mb-1"
                                    onClick={() => setCompareAttachment({ compoundId: compound.id, key, data: entry })}
                                    type="button"
                                  >
                                    {key.replace(/_/g, " ")} {entry.name ? `(${entry.name})` : `#${idx + 1}`}
                                  </button>
                                ) : null
                              );
                            } else {
                              return (
                                <button
                                  key={key}
                                  className="px-2 py-1 bg-[#008080] text-white rounded hover:bg-[#006666] font-bold uppercase tracking-wide text-xs mb-1"
                                  onClick={() => setCompareAttachment({ compoundId: compound.id, key, data: att })}
                                  type="button"
                                >
                                  {key.replace(/_/g, " ")}
                                </button>
                              );
                            }
                          });
                        })()}
                      </div>
                    </td>
                    <td className="p-3 text-center border-[#008080]">
                      <div className="flex justify-center items-center" onClick={e => e.stopPropagation()}>
                        <button
                          onClick={() => toggleStar(compound.id)}
                          className={`text-4xl focus:outline-none ${
                            starred.includes(compound.id) ? 'text-[#00E6D2]' : 'text-gray-400 hover:text-[#00E6D2]'
                          }`}
                          title={starred.includes(compound.id) ? 'Unstar' : 'Star'}
                          style={{
                            textShadow: starred.includes(compound.id)
                              ? '0 0 16px #00E6D2, 0 0 8px #00E6D2'
                              : '0 0 8px #00E6D2',
                            filter: starred.includes(compound.id) ? 'drop-shadow(0 0 8px #00E6D2)' : 'none',
                            transition: 'text-shadow 0.2s, filter 0.2s',
                          }}
                        >
                          {starred.includes(compound.id) ? '★' : '☆'}
                        </button>
                      </div>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}

      {/* Floating Compare Button */}
      {selectedForComparison.length >= 2 && (
        <button
          className="fixed bottom-8 right-8 z-50 bg-[#00E6D2] hover:bg-[#00bfae] text-[#002C36] px-8 py-4 rounded-full shadow-2xl font-bold uppercase tracking-wide text-lg flex items-center gap-3 border-4 border-[#008080] animate-bounce"
          onClick={() => setShowCompareModal(true)}
        >
          <span role="img" aria-label="compare">📊</span> Compare ({selectedForComparison.length})
        </button>
      )}
      <CompareModal
        open={showCompareModal}
        onClose={() => setShowCompareModal(false)}
        compounds={compounds}
        selectedForComparison={selectedForComparison}
        compareAttachment={compareAttachment}
        setCompareAttachment={setCompareAttachment}
        clearComparison={clearComparison}
      />
      
      {/* Attachment Modal */}
      {compareAttachment && (
        <div className="fixed inset-0 bg-black/60 z-50 flex items-center justify-center p-6" onClick={() => setCompareAttachment(null)}>
          <div className="bg-white rounded-lg shadow-xl w-full max-w-7xl max-h-[95vh] p-10 relative flex flex-col items-center justify-center" onClick={e => e.stopPropagation()}>
            <h2 className="text-2xl font-bold mb-6 text-[#008080] uppercase text-center w-full">{compareAttachment.key.replace(/_/g, " ")}</h2>
            {compareAttachment.data.imageUrl && (
              <div className="flex justify-center items-center w-full" style={{ minHeight: '200px' }}>
                <img
                  src={compareAttachment.data.imageUrl}
                  alt={compareAttachment.key}
                  className="w-auto h-auto max-h-[90vh] max-w-[1600px] object-contain mb-6 border rounded shadow bg-white mx-auto"
                />
              </div>
            )}
            {compareAttachment.data.note && (
              <div className="bg-gray-100 p-4 rounded whitespace-pre-wrap mb-6 text-[#002C36] max-w-4xl w-full text-center mx-auto">
                {compareAttachment.data.note}
              </div>
            )}
            {!compareAttachment.data.imageUrl && !compareAttachment.data.note && (
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
      
      {/* Remove inline login modal and block overlay. Redirect handled above. */}
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
            localStorage.removeItem("compoundToOpen");
          }}
          onDelete={async (id: string) => {
            try {
              const path = selectedSource === "lot" && currentLotId
                ? `${process.env.NEXT_PUBLIC_API_BASE_URL}/delete-lot-compound/${currentLotId}/${id}`
                : `${process.env.NEXT_PUBLIC_API_BASE_URL}/delete-compound/${id}`;
              await fetch(path, { method: "DELETE" });
              const updatedLots = await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/lots`).then((res) => res.json());
              setLotList(updatedLots);
              const updated = await fetch(
                selectedSource === "lot" && currentLotId
                  ? `${process.env.NEXT_PUBLIC_API_BASE_URL}/lot/${currentLotId}`
                  : `${process.env.NEXT_PUBLIC_API_BASE_URL}/compounds`
              ).then((res) => res.json());
              setCompounds(updated);
            } catch (err) {
              console.error("Failed to delete compound:", err);
            }
          }}
          onUpdate={async (updatedCompound: Compound) => {
            try {
              const path = selectedSource === "lot" && currentLotId
                ? `${process.env.NEXT_PUBLIC_API_BASE_URL}/update-lot-compound`
                : `${process.env.NEXT_PUBLIC_API_BASE_URL}/update-compound`;
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
                const lotCompounds = await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/lot/${currentLotId}`).then(res => res.json());
                freshCompound = lotCompounds.find((c: Compound) => c.id === updatedCompound.id) || updatedCompound;
                setCompounds(lotCompounds);
                console.log("[DEBUG] Fetched lot compound after update:", freshCompound);
              } else {
                freshCompound = await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/compounds/${updatedCompound.id}`).then(res => res.json());
                const updatedList = await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/compounds`).then(res => res.json());
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
            fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/lots`)
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
              await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/add-compound`, {
                method: "POST",
                headers: { "Content-Type": "application/json" },
                body: JSON.stringify(newCompound),
              });
              // Wait for backend to finish generating the imageUrl
              await new Promise(res => setTimeout(res, 1500));
              let updated = await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/compounds`).then((res) => res.json());
              // If imageUrl is still missing, poll a few more times
              let tries = 0;
              while (
                tries < 4 &&
                updated.some((c: any) => c.id === newCompound.id && (!c.imageUrl || c.imageUrl === ""))
              ) {
                await new Promise(res => setTimeout(res, 1000));
                updated = await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/compounds`).then((res) => res.json());
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
              const res = await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/create-formulation`, {
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




/*
  Copyright © 2025 Polaris Electro Optics
  This code is the property of Polaris Electro Optics and may not be reused,
  modified, or distributed without explicit permission.
*/