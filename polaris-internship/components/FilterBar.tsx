"use client";
import { useEffect, useState } from "react";

type FilterBarProps = {
  onFilterResults: (results: any[]) => void;
  resetSignal: number; // 🔄 New prop to trigger resets
};

export default function FilterBar({ onFilterResults, resetSignal }: FilterBarProps) {
  const [inputValue, setInputValue] = useState("");

  useEffect(() => {
    // 🔁 When resetSignal changes, clear input and show all results
    setInputValue("");
    onFilterResults([]);
  }, [resetSignal]);

  const handleInputChange = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const value = e.target.value;
    setInputValue(value);

    if (value.trim() === "") {
      onFilterResults([]);
      return;
    }

    try {
      const res = await fetch("http://localhost:5000/search-substructure", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ query: value.trim() }),
      });

      const data = await res.json();

      if (Array.isArray(data)) {
        onFilterResults(data);
      } else {
        console.warn("Unexpected response from backend:", data);
        onFilterResults([]);
      }
    } catch (err) {
      console.error("Search failed", err);
      onFilterResults([]);
    }
  };
}




