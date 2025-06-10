"use client";
import { useRouter } from "next/navigation";
import { useEffect, useState } from "react";

export default function CompoundDetailPage({ params }) {
  const { id } = params;
  const [compound, setCompound] = useState(null);
  const router = useRouter();

  useEffect(() => {
    fetch(`http://localhost:5000/compounds`)
      .then((res) => res.json())
      .then((data) => {
        const found = data.find((c) => c.id === id);
        setCompound(found);
      });
  }, [id]);

  if (!compound) return <div className="p-10">Loading...</div>;

  return (
    <div className="min-h-screen flex flex-col items-center justify-center bg-white p-10">
      <button
        onClick={() => router.push("/")}
        className="absolute top-8 left-8 text-lg text-blue-600 hover:underline"
      >
        ← Back to Home
      </button>
      <div className="bg-gradient-to-b from-gray-100 to-gray-300 rounded-lg shadow-xl p-10 w-full max-w-2xl border">
        <h2 className="text-3xl font-bold mb-4">{compound.name}</h2>
        <p className="mb-2"><strong>Formula:</strong> {compound.formula}</p>
        {compound.description && (
          <p className="mb-2"><strong>Description:</strong> {compound.description}</p>
        )}
        {compound.molecularWeight && (
          <p className="mb-2"><strong>Molecular Weight:</strong> {compound.molecularWeight}</p>
        )}
        {/* Add more fields as needed */}
      </div>
    </div>
  );
}