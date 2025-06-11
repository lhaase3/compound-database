"use client";

import { motion } from "framer-motion";
import SMILESRenderer from "./SMILESRenderer";

export default function CompoundCard({ compound, onMoreInfo }) {
  return (
    <motion.div
      onClick={() => onMoreInfo(compound)}
      className="relative bg-gradient-to-b from-gray-100 to-gray-300 text-black rounded-xl shadow-md p-8 flex flex-col items-center gap-6 w-70 h-auto justify-between border border-gray-400 cursor-pointer hover:border-blue-500 transition-colors"
    >
      <SMILESRenderer smiles={compound.smiles || ""} />

      <div>
        <h2 className="text-xl font-bold text-black">{compound.name ?? compound.id}</h2>
      </div>
    </motion.div>

  );
}
