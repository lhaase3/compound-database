"use client";

import { motion } from "framer-motion";
import SMILESRenderer from "./SMILESRenderer";

export default function CompoundCard({ compound, onMoreInfo, isStarred, onToggleStar, cardColor = "#00343F", accentColor = "#00E6D2" }) {
  return (
    <motion.div
      className="relative rounded-2xl shadow-lg p-6 flex flex-col items-center gap-4 border transition-all duration-200 cursor-pointer w-72 group"
      style={{ background: cardColor, borderColor: accentColor }}
      onClick={() => onMoreInfo(compound)}
      whileHover={{ scale: 1.045 }}
    >
      {/* Hover border/glow */}
      <div className="absolute inset-0 rounded-2xl pointer-events-none transition-all duration-200 group-hover:shadow-[0_0_0_4px] group-hover:shadow-[var(--accent)]" style={{ '--accent': accentColor }} />
      {/* Star icon button */}
      <button
        className={`absolute top-3 right-3 text-2xl z-10 focus:outline-none ${isStarred ? 'text-[#00E6D2]' : 'text-white hover:text-[#00E6D2]'}`}
        title={isStarred ? 'Unstar' : 'Star'}
        onClick={e => { e.stopPropagation(); onToggleStar(); }}
        tabIndex={0}
        style={{ textShadow: isStarred ? `0 0 8px ${accentColor}` : 'none' }}
      >
        {isStarred ? '★' : '☆'}
      </button>
      <div className="w-full flex justify-center">
        <SMILESRenderer smiles={compound.smiles || ""} />
      </div>
      <div className="text-center">
        <h2 className="text-lg font-bold" style={{ color: accentColor, textTransform: 'uppercase', letterSpacing: '0.05em' }}>
          {compound.name ?? compound.id}
        </h2>
        {/* <p className="text-sm" style={{ color: '#B2F7F4' }}>{compound.id}</p> */}
      </div>
    </motion.div>
  );
}

