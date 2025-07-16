"use client";
import { useEffect, useState } from "react";
import { auth } from "@/utils/firebase";
import { onAuthStateChanged } from "firebase/auth";
import CompoundCard from "@/components/CompoundCard";
import CompoundModal from "@/components/CompoundModal";
import { Compound } from "@/types/compound";
import CompareModal from "@/components/CompareModal";


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
      <CompareModal
        open={showCompareModal}
        onClose={() => setShowCompareModal(false)}
        compounds={compounds}
        selectedForComparison={selectedForComparison}
        compareAttachment={compareAttachment}
        setCompareAttachment={setCompareAttachment}
        clearComparison={clearComparison}
      />
    </div>
  );
}