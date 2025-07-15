"use client";
import { useEffect, useState } from "react";
import { auth } from "@/utils/firebase";
import { onAuthStateChanged } from "firebase/auth";
import CompoundCard from "@/components/CompoundCard";
import CompoundModal from "@/components/CompoundModal";
import { Compound } from "@/types/compound";

export default function SimilarCompoundsPage() {
  const [user, setUser] = useState<any>(null);
  // Listen for auth state changes (same as main page)
  useEffect(() => {
    const unsubscribe = onAuthStateChanged(auth, (firebaseUser) => {
      if (firebaseUser && firebaseUser.email && firebaseUser.email.endsWith("@polariseo.com")) {
        setUser(firebaseUser);
      } else {
        setUser(null);
      }
    });
    return () => unsubscribe();
  }, []);
  const [compounds, setCompounds] = useState<Compound[]>([]);
  const [originalCompound, setOriginalCompound] = useState<Compound | null>(null);
  const [selectedCompound, setSelectedCompound] = useState<Compound | null>(null);

  const [starred, setStarred] = useState<string[]>([]);
  const [selectedForComparison, setSelectedForComparison] = useState<string[]>([]);
  const [showCompareModal, setShowCompareModal] = useState(false);
  const [compareAttachment, setCompareAttachment] = useState<{ compoundId: string; key: string; data: any } | null>(null);

  useEffect(() => {
    const stored = localStorage.getItem("similarCompounds");
    const original = localStorage.getItem("originalCompound");
    if (stored) {
      setCompounds(JSON.parse(stored));
    }
    if (original) {
      setOriginalCompound(JSON.parse(original));
    }
    // Load starred compounds from backend if logged in, else from localStorage
    if (user && user.email) {
      fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/get-starred?email=${user.email}`)
        .then(res => res.json())
        .then(data => {
          if (Array.isArray(data.starred)) {
            setStarred(data.starred);
            localStorage.setItem("starredCompounds", JSON.stringify(data.starred));
          } else {
            setStarred([]);
            localStorage.setItem("starredCompounds", JSON.stringify([]));
          }
        })
        .catch(err => {
          console.error("Failed to fetch user-starred compounds:", err);
          // fallback to localStorage
          const starredStored = localStorage.getItem("starredCompounds");
          setStarred(starredStored ? JSON.parse(starredStored) : []);
        });
    } else {
      const starredStored = localStorage.getItem("starredCompounds");
      setStarred(starredStored ? JSON.parse(starredStored) : []);
    }
  }, [user]);

  // Star functionality matches main home page
  const toggleStar = (id: string) => {
    setStarred((prev) => {
      const updated = prev.includes(id)
        ? prev.filter((sid) => sid !== id)
        : [...prev, id];
      localStorage.setItem("starredCompounds", JSON.stringify(updated));
      return updated;
    });
  };
  // Sync starred state to localStorage and backend if logged in
  useEffect(() => {
    localStorage.setItem("starredCompounds", JSON.stringify(starred));
    if (user && user.email) {
      fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/save-starred`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ email: user.email, starred }),
      }).catch((err) => console.error("Failed to save starred compounds:", err));
    }
  }, [starred, user]);


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

    const onDelete = async (id: string) => {
    await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/delete-compound/${id}`, { method: "DELETE" });
    setCompounds((prev) => prev.filter((c) => c.id !== id));
    setSelectedCompound(null);
    };

    const onUpdate = async (updatedCompound: Compound) => {
    await fetch(`${process.env.NEXT_PUBLIC_API_BASE_URL}/update-compound`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(updatedCompound),
    });

    setCompounds((prev) =>
        prev.map((c) => (c.id === updatedCompound.id ? updatedCompound : c))
    );

    setSelectedCompound(updatedCompound);
    };


  return (
    <div className="min-h-screen bg-[#002C36] flex flex-col items-center p-6">
      <h1 className="text-4xl font-extrabold text-[#00E6D2] mb-6">Similar Compounds</h1>

      {originalCompound && (
        <div className="mb-10">
          <button
            className="mt-6 px-4 py-2 bg-[#00E6D2] text-[#002C36] rounded hover:bg-[#00bfae] font-bold uppercase tracking-wide"
            onClick={() => {
              localStorage.setItem("compoundToOpen", originalCompound.id);
              localStorage.removeItem("similarCompounds");
              localStorage.removeItem("originalCompound");
              window.location.href = "/";
            }}
          >
            Back to {originalCompound.id}
          </button>
        </div>
      )}

      {compounds.length === 0 && <p className="text-white">No similar compounds found.</p>}

      <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-3 lg:grid-cols-4 gap-10">
        {compounds.map((compound, index) => (
          <CompoundCard
            key={compound.id || `compound-${index}`}
            compound={compound}
            onMoreInfo={setSelectedCompound}
            cardColor="#00343F"
            accentColor="#00E6D2"
            isStarred={starred.includes(compound.id)}
            onToggleStar={() => toggleStar(compound.id)}
            compareChecked={selectedForComparison.includes(compound.id)}
            onToggleCompare={() => toggleCompare(compound.id)}
            similarity={compound.similarity}
          />
        ))}
      </div>

      {selectedCompound && (
        <CompoundModal
        compound={selectedCompound}
        onClose={() => setSelectedCompound(null)}
        onDelete={onDelete}
        onUpdate={onUpdate}
        source="main"
        lotId={null}
        />
      )}

      {selectedForComparison.length >= 2 && (
        <button
          className="fixed bottom-8 right-8 z-50 bg-[#00E6D2] hover:bg-[#00bfae] text-[#002C36] px-8 py-4 rounded-full shadow-2xl font-bold uppercase tracking-wide text-lg flex items-center gap-3 border-4 border-[#008080] animate-bounce"
          onClick={() => setShowCompareModal(true)}
        >
          <span role="img" aria-label="compare">📊</span> Compare ({selectedForComparison.length})
        </button>
      )}


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
              className="absolute top-4 right-10 bg-gray-200 hover:bg-gray-300 text-[#002C36] font-bold px-4 py-2 rounded-md uppercase tracking-wide text-xs"
              onClick={() => {
                clearComparison();
                setShowCompareModal(false);
              }}
            >
              Clear
            </button>
            <div className="overflow-x-auto flex-1 max-h-[70vh]">
              <table className="min-w-full border border-[#008080] rounded-lg divide-y divide-[#008080]">
                <tbody>
                  {/* Display only the desired fields in order */}
                  {["id", "MW", "Lambda Max (DCM/AcCN)", "Lambda Max (neat film)",
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
                    <tr key={field} className="border-b border-[#008080] last:border-b-0">
                      <td className="p-3 font-bold uppercase text-xs text-[#008080] bg-[#f8fafb] border-r border-[#008080]">{field}</td>
                      {selectedForComparison.map((id, idx) => {
                        const cmp = compounds.find(c => c.id === id);
                        const val = cmp?.[field];
                        return <td key={id} className={`p-3 text-center text-[#002C36] border-r border-[#008080] ${idx === selectedForComparison.length - 1 ? 'last:border-r-0' : ''}`}>{val !== undefined && val !== null && val !== "" && val !== "N/A" ? val : "N/A"}</td>;
                      })}
                    </tr>
                  ))}
                  {/* Attachments row (rendered only once, after all fields) */}
                  <tr>
                    <td className="p-3 font-bold uppercase text-xs text-[#008080] bg-[#f8fafb] border-r border-[#008080]">Attachments</td>
                    {selectedForComparison.map((id) => {
                      const cmp = compounds.find(c => c.id === id);
                      const atts = cmp?.attachments || {};
                      
                      const keysWithData = Object.keys(atts).filter(k => {
                        const att = atts[k];
                        if (Array.isArray(att)) {
                          return att.some(entry => entry && (entry.note || entry.imageUrl));
                        } else {
                          return att && (att.note || att.imageUrl);
                        }
                      });

                      return (
                        <td key={id} className="p-3 text-center text-[#002C36] border-r border-[#008080]">
                          <div className="flex flex-col gap-2 items-center justify-center min-w-[150px]">
                            {keysWithData.length === 0 && <span>N/A</span>}
                            {keysWithData.map((key) => {
                              const att = atts[key];
                              if (Array.isArray(att)) {
                                return att.map((entry, idx) => (
                                  (entry && (entry.note || entry.imageUrl)) ? (
                                    <button
                                      key={key + '-' + idx}
                                      className="px-2 py-1 bg-[#008080] text-white rounded hover:bg-[#006666] font-bold uppercase tracking-wide text-xs mb-1"
                                      onClick={() => setCompareAttachment({ compoundId: id, key, data: entry })}
                                      type="button"
                                    >
                                      {key.replace(/_/g, " ")} {entry.name ? `(${entry.name})` : `#${idx + 1}`}
                                    </button>
                                  ) : null
                                ));
                              } else {
                                return (
                                  <button
                                    key={key}
                                    className="px-2 py-1 bg-[#008080] text-white rounded hover:bg-[#006666] font-bold uppercase tracking-wide text-xs mb-1"
                                    onClick={() => setCompareAttachment({ compoundId: id, key, data: att })}
                                    type="button"
                                  >
                                    {key.replace(/_/g, " ")}
                                  </button>
                                );
                              }
                            })}
                          </div>
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

    </div>
  );
}
