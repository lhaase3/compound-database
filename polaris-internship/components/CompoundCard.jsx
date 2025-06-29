"use client";

import { motion } from "framer-motion";
import SMILESRenderer from "./SMILESRenderer";

export default function CompoundCard({ compound, onMoreInfo, isStarred, onToggleStar, cardColor = "#00343F", accentColor = "#00E6D2", compareChecked = false, onToggleCompare, similarity }) {
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
      
      {/* Compare checkbox */}
      {onToggleCompare && (
        <label
          className="absolute top-3 left-3 z-10 flex items-center gap-1 bg-white/80 px-2 py-1 rounded shadow border border-[#008080] cursor-pointer text-xs font-bold uppercase tracking-wide text-[#008080] hover:bg-[#e6f9f7]"
          style={{ userSelect: 'none' }}
          onClick={e => e.stopPropagation()}
        >
          <input
            type="checkbox"
            checked={compareChecked}
            onChange={onToggleCompare}
            className="accent-[#008080] w-4 h-4 mr-1"
            onClick={e => e.stopPropagation()}
          />
          Compare
        </label>
      )}
      {similarity !== undefined && (
        <div className="absolute top-10 left-2 bg-[#00E6D2] text-[#002C36] font-bold text-xs px-2 py-1 rounded z-20">
          {`${(similarity * 100).toFixed(1)}% Similar`}
        </div>
      )}
      <div className="w-full flex justify-center">
        {compound.imageUrl && typeof compound.imageUrl === 'string' && compound.imageUrl.trim() !== '' ? (
          <div className="flex items-center justify-center w-[320px] h-[190px] bg-white rounded" style={{ borderRadius: "0.75rem" }}>
            <img
              src={compound.imageUrl}
              alt={compound.name || compound.id}
              className="w-full h-full object-contain"
              style={{ display: "block", margin: "auto" }}
              onError={e => { e.currentTarget.style.display = 'none'; }}
            />
          </div>
        ) : (
          <SMILESRenderer smiles={compound.smiles || ""} />
        )}
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

