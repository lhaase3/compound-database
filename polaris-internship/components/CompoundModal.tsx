import React, { useState, useEffect, useRef } from "react";
import { Compound, AttachmentData } from "../types/compound";
import AttachmentModal from "./AttachmentModal";
import CreateLotModal from "./CreateLotModal";

type Props = {
  compound: Compound;
  onClose: () => void;
  onDelete: (id: string) => void;
  onUpdate: (updatedCompound: Compound) => Promise<void>;
  source: "main" | "lot";
  lotId: string | null;
  onUpdateCompoundFromLot?: (compound: Compound, lotId: string | null) => void;
};

// --- Types for multi-attachment support ---
type MultiAttachmentEntry = { name: string; note: string; imageUrl: string };

type CompoundWithAttachments = Omit<Compound, 'attachments'> & {
  attachments: { [key: string]: AttachmentData };
};

// --- Helper to check if an attachment supports multiple entries ---
const multiEntryAttachments = ["lcms", "uv_vis", "dsc", "thermal_stability", "pda_detector_spectrum"];
function isMultiAttachment(key: string) {
  return multiEntryAttachments.includes(key);
}

// --- Add to AttachmentModal prop types (at the top of the file, or import if external) ---
// interface AttachmentModalProps extends ... {
//   isMulti?: boolean;
//   isNewEntry?: boolean;
//   renderHeaderExtra?: React.ReactNode;
// }

// Helper to merge updates with the full compound
function getFullCompoundUpdate(base: Compound, updates: Partial<CompoundWithAttachments>): Compound {
  // Always preserve id, smiles, and all required fields
  return {
    ...base,
    ...updates,
    id: base.id,
    smiles: base.smiles,
    attachments: {
      ...base.attachments,
      ...updates.attachments,
    },
  };
}

// --- NewAttachmentEntryModal: Dedicated child component for new multi-attachment entries ---
function NewAttachmentEntryModal({
  keyName,
  onAdd,
  onClose
}: {
  keyName: string;
  onAdd: (entry: any) => void;
  onClose: () => void;
}) {
  const [name, setName] = useState("");
  const [note, setNote] = useState("");
  const [imageUrl, setImageUrl] = useState("");
  const [file, setFile] = useState<File | null>(null);
  const [uploading, setUploading] = useState(false);
  const [error, setError] = useState("");

  const handleFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const selected = e.target.files?.[0];
    if (selected) {
      setFile(selected);
      setImageUrl(URL.createObjectURL(selected)); // preview
    }
  };



  const handleAdd = async () => {
    setError("");
    let finalUrl = imageUrl;
    if (file) {
      setUploading(true);
      try {
        const formData = new FormData();
        formData.append("file", file);
        formData.append("note", note);
        const res = await fetch("http://localhost:5000/upload-image-to-firebase", {
          method: "POST",
          body: formData,
        });
        if (!res.ok) {
          const errData = await res.json();
          throw new Error(errData.error || "Upload failed");
        }
        const result = await res.json();
        finalUrl = result.fileUrl;
      } catch (err: any) {
        setError(err.message || "Upload failed");
        setUploading(false);
        return;
      }
      setUploading(false);
    }
    onAdd({ name, note, imageUrl: finalUrl });
    onClose();
  };

  return (
    <div className="fixed inset-0 bg-black/50 z-50 flex items-center justify-center p-6" onClick={onClose}>
      <div className="bg-white rounded-2xl border-2 border-[#00E6D2] shadow-2xl w-full max-w-md p-8" onClick={e => e.stopPropagation()}>
        <h2 className="text-2xl font-bold mb-4 text-[#002C36] uppercase tracking-wide">Add New {keyName.replace(/_/g, " ")}</h2>
        <input
          type="text"
          className="w-full border border-[#00E6D2] rounded p-2 mb-4 text-[#002C36] bg-white focus:ring-2 focus:ring-[#00E6D2]"
          placeholder="Name (optional)"
          value={name}
          onChange={e => setName(e.target.value)}
        />
        <textarea
          className="w-full border border-[#00E6D2] rounded p-2 mb-4 text-[#002C36] bg-white focus:ring-2 focus:ring-[#00E6D2]"
          placeholder="Note"
          value={note}
          onChange={e => setNote(e.target.value)}
        />
        <input
          type="file"
          accept="image/*"
          className="mb-4"
          onChange={handleFileChange}
        />
        {imageUrl && (
          <img src={imageUrl} alt="preview" className="w-full max-h-48 object-contain mb-4 rounded border" />
        )}
        {error && <div className="text-red-600 mb-2">{error}</div>}
        <div className="flex gap-2 justify-end">
          <button className="px-4 py-2 bg-gray-200 text-gray-700 rounded hover:bg-gray-300" onClick={onClose} disabled={uploading}>Cancel</button>
          <button
            className="px-4 py-2 bg-[#00E6D2] text-[#002C36] rounded hover:bg-[#00bfae] font-bold"
            onClick={handleAdd}
            disabled={uploading}
          >
            {uploading ? "Uploading..." : "Add"}
          </button>
        </div>
      </div>
    </div>
  );
}

export default function CompoundModal({
  compound,
  onClose,
  onDelete,
  onUpdate,
  source,
  lotId,
  onUpdateCompoundFromLot,
}: Props) {
  // Removed duplicate declaration of showBWFullStructure
  const [editMode, setEditMode] = useState(false);
  const [editedCompound, setEditedCompound] = useState<CompoundWithAttachments>(() => ({
    ...compound,
    attachments: {
      uv_vis: Array.isArray(compound.attachments?.uv_vis)
        ? { note: '', imageUrl: '' }
        : { note: compound.attachments?.uv_vis?.note || '', imageUrl: compound.attachments?.uv_vis?.imageUrl || '' },
      dsc: Array.isArray(compound.attachments?.dsc)
        ? { note: '', imageUrl: '' }
        : { note: compound.attachments?.dsc?.note || '', imageUrl: compound.attachments?.dsc?.imageUrl || '' },
      lcms: Array.isArray(compound.attachments?.lcms)
        ? { note: '', imageUrl: '' }
        : { note: compound.attachments?.lcms?.note || '', imageUrl: compound.attachments?.lcms?.imageUrl || '' },
      thermal_stability: Array.isArray(compound.attachments?.thermal_stability)
        ? { note: '', imageUrl: '' }
        : { note: compound.attachments?.thermal_stability?.note || '', imageUrl: compound.attachments?.thermal_stability?.imageUrl || '' },
      pda_detector_spectrum: Array.isArray(compound.attachments?.pda_detector_spectrum)
        ? { note: '', imageUrl: '' }
        : { note: compound.attachments?.pda_detector_spectrum?.note || '', imageUrl: compound.attachments?.pda_detector_spectrum?.imageUrl || '' },
      ...Object.fromEntries(
        Object.entries(compound.attachments || {})
          .filter(([k]) => !['uv_vis', 'dsc', 'lcms', 'thermal_stability', 'pda_detector_spectrum'].includes(k))
          .map(([k, att]) => [
            k,
            Array.isArray(att)
              ? []
              : { note: (att as any)?.note || '', imageUrl: (att as any)?.imageUrl || '' }
          ])
      )
    },
  }));

  // Track user-added custom fields
  const [lotsForCompound, setLotsForCompound] = useState<string[]>([]);
  const [showLotDropdown, setShowLotDropdown] = useState(false);
  const [attachmentDropdown, setAttachmentDropdown] = useState<{key: string, anchor: HTMLElement | null} | null>(null);
  const [selectedAttachment, setSelectedAttachment] = useState<string | null>(null);
  const [selectedAttachmentIdx, setSelectedAttachmentIdx] = useState<number | null>(null);
  const [showNewEntryModal, setShowNewEntryModal] = useState<{key: string}|null>(null);
  const isLotVersion = source === "lot" && compound.original_id;
  const [showCreateLotModal, setShowCreateLotModal] = useState(false);
  const [showAddTextField, setShowAddTextField] = useState(false);
  const [newTextFieldName, setNewTextFieldName] = useState("");
  const [newTextFieldValue, setNewTextFieldValue] = useState("");
  const [showFullStructureModal, setShowFullStructureModal] = useState(false);
  const [showTagDropdown, setShowTagDropdown] = useState(false);
  const [currentLotId, setCurrentLotId] = useState<string | null>(lotId);
  const lotDropdownRef = useRef<HTMLDivElement>(null);
  const tagDropdownRef = useRef<HTMLDivElement>(null);
  const attachmentDropdownRef = useRef<HTMLDivElement | null>(null); // <-- NEW


  // List of unwanted fields to exclude from custom display
  const unwanted = [
    "attachments", "lots", "original_id", "created_at", "updated_at", "_id", "__v", "isStarred", "parsed_phase_transitions", "imageUrl", "tags", "lotId", "Lambda Max (DCM/Ac CN)", "views", "similarity", "bwImageUrl"
  ];

  // Your defined fields (always shown, in order)
  const fields = [
    "id", "MW", "Lambda Max (DCM/AcCN)", "Lambda Max (neat film)",
    "phase map", "r33", "dipole CAMB3LYP SVPD CHCl3 (Cosmo)",
    "beta CAMB3LYP SVPD CHCl3 (Cosmo)","dipole B3LYP SVPD CHCl3",
    "beta B3LYP SVPD CHCl3","beta/MW",
    "J/g DSC melt (total)", "kJ/mol DSC melt (total)", 
    "Refractive index (ne/no)", "Notes", "lab?", "first PEO#",
    "registered PEO#", "Lab book #", "Max loading (%)", "smiles"
  ];

  useEffect(() => {
    if (compound?.id) {
      fetch(`http://localhost:5000/lots-for-compound/${compound.id}`)
        .then((res) => res.json())
        .then(setLotsForCompound)
        .catch((err) => console.error("Failed to fetch lots:", err));
    }
  }, [compound.id]);

  const [showBWFullStructure, setShowBWFullStructure] = useState(false);
  const handleChange = (field: string, value: string) => {
    setEditedCompound((prev) => ({ ...prev, [field]: value }));
  };

  // --- Fix handleSave (edit/save compound) ---
  const handleSave = async () => {
    try {
      if (source === "lot" && onUpdateCompoundFromLot) {
        // Save to lot compound endpoint
        const response = await fetch("http://localhost:5000/update-lot-compound", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ ...editedCompound, lotId: currentLotId }),
        });
        if (!response.ok) {
          throw new Error("Failed to update lot compound");
        }
        // Fetch the updated lot compound from Firestore and update editedCompound
        const updated = await fetch(`http://localhost:5000/lot/${currentLotId}`)
          .then(res => res.json());
        // Find the correct compound by id
        const updatedLotCompound = updated.find((c: any) => c.id === editedCompound.id);
        if (updatedLotCompound) {
          setEditedCompound((prev) => ({ ...prev, ...updatedLotCompound }));
          onUpdateCompoundFromLot?.(updatedLotCompound, currentLotId); // ✅ Force re-sync parent
        }
        setEditMode(false);
      } else {
        await onUpdate(getFullCompoundUpdate(compound, editedCompound));
        setEditMode(false);
      }
    } catch (err) {
      console.error("Update failed:", err);
    }
  };

  const handleDelete = async () => {
    const confirmed = window.confirm("Are you sure you want to delete this compound?");
    if (confirmed) {
      try {
        await onDelete(compound.id);
        onClose();
      } catch (err) {
        console.error("Delete failed:", err);
      }
    }
  };

  // --- Normalize attachments on load to support arrays for multi ---
  useEffect(() => {
    const builtIn = multiEntryAttachments;
    const allKeys = compound.attachments
      ? Array.from(new Set([...Object.keys(compound.attachments), ...builtIn]))
      : builtIn;
    // Fix type for normalizedAttachments
    const normalizedAttachments: Record<string, AttachmentData | MultiAttachmentEntry[]> = {};
    allKeys.forEach((key) => {
      const att = compound.attachments?.[key];
      if (isMultiAttachment(key)) {
        if (Array.isArray(att)) {
          normalizedAttachments[key] = att.map((entry) => ({
            name: entry.name || '',
            note: entry.note || '',
            imageUrl: entry.imageUrl || '',
          }));
        } else if (att && (att.note || att.imageUrl)) {
          normalizedAttachments[key] = [{
            name: key.charAt(0).toUpperCase() + key.slice(1),
            note: att.note || '',
            imageUrl: att.imageUrl || '',
          }];
        } else {
          normalizedAttachments[key] = [];
              <div style={{ position: 'absolute', top: 24, right: 32, zIndex: 10 }}>
                <button
                  onClick={() => setShowBWFullStructure((prev) => !prev)}
                  className="px-4 py-2 bg-[#002C36] text-[#00E6D2] rounded-lg font-bold uppercase tracking-wide hover:bg-[#00bfae]"
                  style={{ minWidth: 120 }}
                >
                  {showBWFullStructure ? 'Show Colored' : 'Show B&W'}
                </button>
              </div>
        }
      } else {
        if (Array.isArray(att)) {
          normalizedAttachments[key] = { note: '', imageUrl: '' };
        } else {
          normalizedAttachments[key] = {
            note: att?.note || '',
            imageUrl: att?.imageUrl || '',
          };
        }
      }
    });
    setEditedCompound((prev) => ({ ...compound, attachments: normalizedAttachments }));
  }, [compound]);


  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      if (lotDropdownRef.current && !lotDropdownRef.current.contains(event.target as Node)) {
        setShowLotDropdown(false);
      }
    };

    if (showLotDropdown) {
      document.addEventListener("mousedown", handleClickOutside);
    } else {
      document.removeEventListener("mousedown", handleClickOutside);
    }

    return () => {
      document.removeEventListener("mousedown", handleClickOutside);
    };
  }, [showLotDropdown]);


  useEffect(() => {
    const handleClickOutsideTag = (event: MouseEvent) => {
      if (tagDropdownRef.current && !tagDropdownRef.current.contains(event.target as Node)) {
        setShowTagDropdown(false);
      }
    };

    if (showTagDropdown) {
      document.addEventListener("mousedown", handleClickOutsideTag);
    } else {
      document.removeEventListener("mousedown", handleClickOutsideTag);
    }

    return () => {
      document.removeEventListener("mousedown", handleClickOutsideTag);
    };
  }, [showTagDropdown]);


  useEffect(() => {
    const handleClickOutsideAttachment = (event: MouseEvent) => {
      if (
        attachmentDropdownRef.current &&
        !attachmentDropdownRef.current.contains(event.target as Node)
      ) {
        setAttachmentDropdown(null);
      }
    };
    if (attachmentDropdown) {
      document.addEventListener("mousedown", handleClickOutsideAttachment);
    } else {
      document.removeEventListener("mousedown", handleClickOutsideAttachment);
    }
    return () => {
      document.removeEventListener("mousedown", handleClickOutsideAttachment);
    };
  }, [attachmentDropdown]);



  // Helper: get all attachment keys (built-in + custom)
  const getAllAttachmentKeys = () => Object.keys(editedCompound.attachments || {});

  // Add Attachment Field
  // --- Fix handleAddAttachmentField ---
  const handleAddAttachmentField = async () => {
    const name = prompt("Enter a name for the new attachment field:");
    if (
      name &&
      !Object.keys(editedCompound.attachments || {}).includes(name)
    ) {
      const updatedCompound = {
        ...editedCompound,
        attachments: {
          ...editedCompound.attachments,
          [name]: { note: "", imageUrl: "" },
        },
      };
      setEditedCompound(updatedCompound);
      await onUpdate(getFullCompoundUpdate(compound, updatedCompound));
    }
  };


  const FullStructurePreview = ({ smiles }: { smiles: string }) => {
    if (showBWFullStructure && compound.bwImageUrl && typeof compound.bwImageUrl === 'string' && compound.bwImageUrl.trim() !== '') {
      return (
        <div className="flex justify-center items-center w-full h-[600px]">
          <img
            src={compound.bwImageUrl}
            alt={compound.name || compound.id}
            style={{ maxWidth: "100%", maxHeight: 560, objectFit: "contain", background: 'white' }}
            onError={e => { e.currentTarget.style.display = 'none'; }}
          />
        </div>
      );
    }
    if (compound.imageUrl && typeof compound.imageUrl === 'string' && compound.imageUrl.trim() !== '') {
      return (
        <div className="flex justify-center items-center w-full h-[600px]">
          <img
            src={compound.imageUrl}
            alt={compound.name || compound.id}
            style={{ maxWidth: "100%", maxHeight: 560, objectFit: "contain" }}
            onError={e => { e.currentTarget.style.display = 'none'; }}
          />
        </div>
      );
    }
    const containerRef = useRef<HTMLDivElement>(null);
    useEffect(() => {
      if (!window.RDKit || !smiles) return;
      let mol;
      try {
        mol = window.RDKit.get_mol(smiles);
      } catch (err) {
        console.error("Invalid SMILES for full preview:", err);
        return;
      }
      if (!mol) return;
      // Render a large, high-res SVG and scale it up visually
      const svg = mol.get_svg_with_highlight({}, {}, 2000, 1200);
      if (containerRef.current) {
        containerRef.current.innerHTML = svg;
        containerRef.current.style.width = "1000px";
        containerRef.current.style.height = "600px";
        containerRef.current.style.display = "flex";
        containerRef.current.style.justifyContent = "center";
        containerRef.current.style.alignItems = "center";
        const svgElem = containerRef.current.querySelector('svg');
        if (svgElem) {
          svgElem.setAttribute('width', '1000px');
          svgElem.setAttribute('height', '600px');
        }
      }
      mol.delete();
    }, [smiles, showBWFullStructure]);
    return (
      <div
        className="flex justify-center items-center overflow-auto"
        ref={containerRef}
        style={{
          width: "100%",
          height: "600px",
        }}
      />
    );
  };

  // Print handler
  const printRef = useRef<HTMLDivElement>(null);
  const handlePrint = () => {
  const printContents = document.getElementById('printable-area')?.innerHTML;

  if (!printContents) return;

  const printWindow = window.open('', 'compound-print', 'width=800,height=900');
  if (!printWindow) return;

  printWindow.document.write(`
    <html>
      <head>
        <title>Print Compound</title>
        <style>
          @media print {
            @page {
              margin: 0;
            }
            html, body {
              margin: 0;
              padding: 0;
              font-family: Arial, sans-serif;
              background: white;
              color: black;
            }
            img {
              filter: grayscale(100%) !important;
            }
          }
        </style>
      </head>
      <body onload="window.print(); window.close();">
        ${printContents}
      </body>
    </html>
  `);

  printWindow.document.close();
  printWindow.focus();
};

  return (
    <div
      className="fixed inset-0 bg-black/30 z-50 flex items-center justify-center"
      onClick={onClose}
      style={{backdropFilter: 'blur(2px)'}}
    >
      <div
        className="bg-white border-2 border-[#00E6D2] p-8 rounded-2xl shadow-2xl max-w-3xl w-full max-h-[92vh] overflow-y-auto relative"
        onClick={(e) => e.stopPropagation()}
      >
        {/* Header */}
        <div className="flex justify-between items-center mb-6 border-b-2 border-[#00E6D2] pb-3">
          <h2 className="text-3xl font-extrabold text-[#002C36] uppercase tracking-wide flex items-center gap-2">
            Compound Details
            {currentLotId && (
              <span className="text-purple-600 text-base ml-2 font-semibold">(from Lot: {currentLotId})</span>
            )}
          </h2>
          <button onClick={onClose} className="text-[#002C36] text-3xl font-bold hover:text-[#00E6D2] transition-colors ml-4">&times;</button>
        </div>

        {/* Tag and Lot Buttons */}
        <div className="flex flex-wrap gap-3 mb-6 items-center justify-between">
          <div className="flex gap-2">
            <div className="relative">
              <button
                className="border border-[#00E6D2] text-[#00E6D2] font-bold px-2 py-0.5 rounded-md uppercase tracking-wide text-xs hover:bg-[#00E6D2] hover:text-[#002C36] transition-all"
                onClick={() => setShowTagDropdown((prev) => !prev)}
              >
                🏷️ Tags
              </button>
              {showTagDropdown && (
                <div
                  ref={tagDropdownRef}
                  className="absolute left-0 bg-[#00343F] shadow-lg border border-[#00E6D2] mt-2 z-50 rounded w-44 text-[#00E6D2]"
                >
                  {["testing", "crystals"].map((tag) => (
                    <div
                      key={tag}
                      className="px-4 py-2 text-[#00E6D2] hover:bg-[#00545F] cursor-pointer font-semibold"
                      onClick={async () => {
                        const updatedTags = editedCompound.tags?.includes(tag)
                          ? editedCompound.tags.filter((t: string) => t !== tag)
                          : [...(editedCompound.tags || []), tag];
                        const updatedCompound = { ...editedCompound, tags: updatedTags };
                        setEditedCompound(updatedCompound);
                        setShowTagDropdown(false);
                        await onUpdate(getFullCompoundUpdate(compound, updatedCompound));
                      }}
                    >
                      <span className={editedCompound.tags?.includes(tag) ? "font-bold text-[#00E6D2]" : ""}>
                        {editedCompound.tags?.includes(tag) ? "✔️ " : ""}
                        {tag}
                      </span>
                    </div>
                  ))}
                </div>
              )}
            </div>
            <div className="relative">
              <button
                className="border border-purple-600 text-purple-600 font-bold px-2 py-0.5 rounded-md uppercase tracking-wide text-xs hover:bg-purple-100 transition-all"
                onClick={() => setShowLotDropdown((prev) => !prev)}
              >
                📦 Lots ({lotsForCompound.length})
              </button>
              {showLotDropdown && (
                <div ref={lotDropdownRef} className="absolute left-0 bg-[#00343F] shadow-lg border border-purple-400 mt-2 z-50 rounded w-56 text-[#00E6D2]">
                  <div
                    onClick={() => {
                      setShowLotDropdown(false);
                      setShowCreateLotModal(true);
                    }}
                    className="px-4 py-2 text-blue-300 hover:bg-[#00545F] cursor-pointer border-b border-purple-200 font-semibold"
                  >
                    ➕ Create Lot
                  </div>
                  {lotsForCompound.map((lotId) => (
                    <div
                      key={lotId}
                      className="px-4 py-2 text-[#00E6D2] hover:bg-[#00545F] cursor-pointer font-semibold"
                      onClick={async () => {
                        try {
                          const res = await fetch(`http://localhost:5000/lot/${lotId}`);
                          const lotCompounds = await res.json();
                          if (Array.isArray(lotCompounds) && lotCompounds.length > 0) {
                            let lotCompound;
                            if (lotCompounds.length === 1) {
                              lotCompound = lotCompounds[0];
                            } else {
                              const options = lotCompounds.map((c, i) => `${i + 1}: ${c.id}`).join("\n");
                              const idx = window.prompt(`Select which lot compound to view (enter number):\n${options}`, "1");
                              const i = Number(idx) - 1;
                              if (!isNaN(i) && i >= 0 && i < lotCompounds.length) {
                                lotCompound = lotCompounds[i];
                              } else {
                                return;
                              }
                            }
                            onUpdateCompoundFromLot?.(lotCompound, lotId);
                          }
                          setShowLotDropdown(false);
                        } catch (err) {
                          console.error("Failed to load lot compound:", err);
                        }
                      }}
                    >
                      {lotId}
                    </div>
                  ))}
                </div>
              )}
            </div>
            <button
              className="bg-[#00E6D2] hover:bg-[#00bfae] text-[#002C36] font-bold px-2 py-0.5 rounded-md uppercase tracking-wide text-xs shadow transition-all"
              onClick={async () => {
                try {
                  const res = await fetch("http://localhost:5000/similar-compounds", {
                    method: "POST",
                    headers: { "Content-Type": "application/json" },
                    body: JSON.stringify({ smiles: compound.smiles, id: compound.id, threshold: 0.7 }),
                  });
                  const data = await res.json();
                  localStorage.setItem("similarCompounds", JSON.stringify(data)); // ✅ Save similar compounds
                  localStorage.setItem("originalCompound", JSON.stringify(compound)); // ✅ Save the original compound
                  window.location.href = "/similar-compounds"; // ✅ Redirect
                } catch (err) {
                  console.error("Failed to fetch similar compounds:", err);
                }
              }}
            >
              Similar Compounds
            </button>
            {isLotVersion && (
              <button
                className="border border-blue-600 text-blue-600 font-bold px-2 py-0.5 rounded-md uppercase tracking-wide text-xs hover:bg-blue-100 transition-all"
                onClick={async () => {
                  try {
                    const res = await fetch(`http://localhost:5000/compounds/${compound.original_id}`);
                    const original = await res.json();
                    setEditedCompound({
                      ...original,
                      attachments: {
                        uv_vis: {
                          note: original.attachments?.uv_vis?.note || "",
                          imageUrl: original.attachments?.uv_vis?.imageUrl || "",
                        },
                        dsc: {
                          note: original.attachments?.dsc?.note || "",
                          imageUrl: original.attachments?.dsc?.imageUrl || "",
                        },
                        lcms: {
                          note: original.attachments?.lcms?.note || "",
                          imageUrl: original.attachments?.lcms?.imageUrl || "",
                        },
                      },
                    });
                    onUpdateCompoundFromLot?.(original, null);
                  } catch (err) {
                    console.error("Failed to load original compound", err);
                  }
                }}
              >
                ⬅️ Back to Original
              </button>
            )}
          </div>
          <div className="flex gap-2">
            <button
              className="bg-[#00E6D2] hover:bg-[#00bfae] text-[#002C36] font-bold px-2 py-0.5 rounded-md uppercase tracking-wide text-xs shadow transition-all"
              onClick={handlePrint}
            >
              Print
            </button>
            <button
              className="bg-[#00E6D2] hover:bg-[#00bfae] text-[#002C36] font-bold px-2 py-0.5 rounded-md uppercase tracking-wide text-xs shadow transition-all"
              onClick={() => (editMode ? handleSave() : setEditMode(true))}
            >
              {editMode ? "Save" : "Edit"}
            </button>
            <button className="bg-red-100 hover:bg-red-200 text-red-600 font-bold px-2 py-0.5 rounded-md uppercase tracking-wide text-xs shadow transition-all" onClick={handleDelete}>
              Delete
            </button>
          </div>
        </div>

        <div
          id="printable-area"
          style={{
            display: 'none',
            textAlign: 'center',
            fontSize: '14px',
            color: '#000',
            padding: 0,
            margin: 0,
            lineHeight: 1.2,       // tighter spacing
          }}
        >
          {compound.bwImageUrl && typeof compound.bwImageUrl === 'string' && compound.bwImageUrl.trim() !== '' && (
            <img
              src={compound.bwImageUrl}
              alt={compound.name || compound.id}
              style={{
                width: '60%',
                maxWidth: '500px',
                margin: '0 auto',
                height: 'auto',
                display: 'block',
                paddingTop: '15px',    // tighter gap
                paddingBottom: '0px',    // tighter gap
              }}
            />
            )}
            <div
              style={{
                fontSize: '2.4rem',
                fontWeight: '600',
                letterSpacing: '0.05em',
                marginTop: '0px',         // no extra gap
                textAlign: 'center',
              }}
            >
              {compound.id}
            </div>
        </div>
        

        {/* Section Divider */}
        <div className="w-full flex items-center mb-6">
          <hr className="flex-grow border-t border-[#00E6D2]/40" />
          <span className="mx-4 text-lg text-[#00E6D2] font-bold uppercase tracking-wide">Compound Info</span>
          <hr className="flex-grow border-t border-[#00E6D2]/40" />
        </div>

        {/* Fields Grid */}
        <div className="grid grid-cols-1 sm:grid-cols-2 gap-6 text-[#002C36] mb-8">
          {fields.map((field) => (
            <div key={field} className="flex flex-col mb-2">
              <span className="text-xs font-bold uppercase text-[#008080] mb-1 tracking-wide">{field}</span>
              {editMode ? (
                <input
                  className="border border-[#00E6D2] rounded px-2 py-1 text-sm bg-white text-[#002C36] focus:ring-2 focus:ring-[#00E6D2]"
                  value={editedCompound[field] ?? ""}
                  onChange={(e) => handleChange(field, e.target.value)}
                  style={field === "smiles" ? { wordBreak: 'break-all', whiteSpace: 'pre-wrap', maxWidth: '100%' } : {}}
                />
              ) : (
              <span
                className={editedCompound[field] ? "text-[#002C36] ..." : "text-gray-400"}
                style={field === "smiles" ? { wordBreak: 'break-all', whiteSpace: 'pre-wrap', maxWidth: '100%', display: 'block' } : {}}
              >
                {editedCompound[field] ? editedCompound[field] : "N/A"}
              </span>
              )}
            </div>
          ))}
          {Object.keys(editedCompound)
            .filter((key) => !fields.includes(key) && !unwanted.includes(key))
            .map((key) => (
              <div key={key} className="flex flex-col mb-2">
                <span className="text-xs font-bold uppercase text-[#008080] mb-1 tracking-wide">{key}</span>
                {editMode ? (
                  <input
                    className="border border-[#00E6D2] rounded px-2 py-1 text-sm bg-white text-[#002C36] focus:ring-2 focus:ring-[#00E6D2]"
                    value={editedCompound[key] ?? ""}
                    onChange={(e) => handleChange(key, e.target.value)}
                  />
                ) : (
                  <span className={editedCompound[key] ? "text-[#002C36]" : "text-gray-400"}>
                    {editedCompound[key] ? editedCompound[key] : "N/A"}
                  </span>
                )}
              </div>
            ))}
        </div>

        {/* Attachments Section */}
        <div className="w-full flex flex-col items-center mb-8">
          <div className="flex gap-2 flex-wrap justify-center mb-2">
            {Object.entries(editedCompound.attachments || {}).map(([key, value]) => (
              isMultiAttachment(key) && Array.isArray(value) ? (
                <div key={key} className="relative">
                  <button
                    className="border border-[#00E6D2] text-[#00E6D2] font-bold px-2 py-0.5 rounded-md uppercase tracking-wide text-xs hover:bg-[#00E6D2] hover:text-[#002C36] transition-all mb-1"
                    onClick={e => setAttachmentDropdown({ key, anchor: e.currentTarget })}
                  >
                    {key.replace(/_/g, " ")}
                  </button>
                  {attachmentDropdown && attachmentDropdown.key === key && (
                    <div
                      ref={attachmentDropdownRef}
                      className="absolute z-50 bg-[#00343F] border border-[#00E6D2] rounded shadow-lg mt-2 left-0 min-w-[200px] max-h-72 overflow-y-auto text-[#00E6D2]"
                    >
                      <div
                        className="px-4 py-2 hover:bg-[#00545F] cursor-pointer font-semibold border-b border-[#00E6D2]"
                        onClick={() => {
                          setAttachmentDropdown(null);
                          setShowNewEntryModal({key});
                        }}
                      >
                        + New {key.replace(/_/g, " ")}
                      </div>
                      {value.map((entry, idx) => (
                        <div
                          key={idx}
                          className="px-4 py-2 hover:bg-[#00545F] cursor-pointer font-semibold"
                          onClick={() => {
                            setAttachmentDropdown(null);
                            setSelectedAttachment(key);
                            setSelectedAttachmentIdx(idx);
                          }}
                        >
                          {entry.name || `${key} ${idx + 1}`}
                        </div>
                      ))}
                    </div>
                  )}
                </div>
              ) : (
                <button
                  key={key}
                  className="border border-[#00E6D2] text-[#00E6D2] font-bold px-2 py-0.5 rounded-md uppercase tracking-wide text-xs hover:bg-[#00E6D2] hover:text-[#002C36] transition-all mb-1"
                  onClick={() => {
                    if (isMultiAttachment(key)) {
                      // If only one entry, open directly
                      setSelectedAttachment(key);
                      setSelectedAttachmentIdx(0);
                    } else {
                      setSelectedAttachment(key);
                      setSelectedAttachmentIdx(null);
                    }
                  }}
                >
                  {key.replace(/_/g, " ")}
                </button>
              )
            ))}
          </div>
          <div className="flex gap-4 mt-2">
            <button
              className="px-4 py-2 bg-blue-500 text-white rounded hover:bg-blue-700 cursor-pointer font-bold uppercase tracking-wide"
              onClick={handleAddAttachmentField}
            >
              Add Attachment Field
            </button>
            <button
              className="px-4 py-2 bg-blue-500 text-white rounded hover:bg-blue-700 cursor-pointer font-bold uppercase tracking-wide"
              onClick={() => setShowAddTextField(true)}
            >
              Add Text Field
            </button>
          </div>
        </div>

        {/* Structure Preview */}
        <div className="mt-6 text-center">
          <div
            className="bg-white p-2 rounded-lg border-2 border-[#00E6D2] shadow hover:shadow-md transition-all cursor-pointer mx-auto flex items-center justify-center"
            style={{ width: 200, height: 150, overflow: "hidden" }}
            onClick={() => setShowFullStructureModal(true)}
          >
            {compound.imageUrl ? (
              <img
                src={compound.imageUrl}
                alt={compound.name || compound.id}
                style={{ width: "100%", height: "100%", objectFit: "contain" }}
              />
            ) : (
              typeof window !== "undefined" && window.RDKit ? (() => {
                try {
                  const mol = window.RDKit.get_mol(compound.smiles || "");
                  return mol ? (
                    <div
                      dangerouslySetInnerHTML={{ __html: mol.get_svg() }}
                      style={{ display: "flex", alignItems: "center", justifyContent: "center", height: "100%", transform: "scale(0.6)", transformOrigin: "center" }}
                    />
                  ) : null;
                } catch {
                  return null;
                }
              })() : (
                <span className="text-sm text-gray-500">Loading...</span>
              )
            )}
          </div>
          <div className="text-xs text-gray-500 mt-1">Click to view full structure</div>
        </div>

        {/* Add Text Field Modal */}
        {showAddTextField && (
          <div className="fixed inset-0 bg-black/50 z-50 flex items-center justify-center p-6" onClick={() => setShowAddTextField(false)}>
            <div className="bg-white rounded-2xl border-2 border-[#00E6D2] shadow-2xl w-full max-w-md p-8" onClick={e => e.stopPropagation()}>
              <h2 className="text-2xl font-bold mb-4 text-[#002C36] uppercase tracking-wide">Add Text Field</h2>
              <input
                type="text"
                className="w-full border border-[#00E6D2] rounded p-2 mb-4 text-[#002C36] bg-white focus:ring-2 focus:ring-[#00E6D2]"
                placeholder="Field name"
                value={newTextFieldName}
                onChange={e => setNewTextFieldName(e.target.value)}
              />
              <textarea
                className="w-full border border-[#00E6D2] rounded p-2 mb-4 text-[#002C36] bg-white focus:ring-2 focus:ring-[#00E6D2]"
                placeholder="Field value"
                value={newTextFieldValue}
                onChange={e => setNewTextFieldValue(e.target.value)}
              />
              <div className="flex justify-end gap-2">
                <button
                  className="px-4 py-2 bg-gray-200 rounded hover:bg-gray-300 text-[#002C36] font-bold uppercase tracking-wide"
                  onClick={() => setShowAddTextField(false)}
                >
                  Cancel
                </button>
                <button
                  className="px-4 py-2 bg-[#00E6D2] text-[#002C36] rounded hover:bg-[#00bfae] font-bold uppercase tracking-wide"
                  onClick={async () => {
                    if (!newTextFieldName.trim()) return;
                    const updatedCompound = {
                      ...editedCompound,
                      [newTextFieldName]: newTextFieldValue,
                    };
                    setEditedCompound(updatedCompound);
                    setShowAddTextField(false);
                    setNewTextFieldName("");
                    setNewTextFieldValue("");
                    await onUpdate(getFullCompoundUpdate(compound, updatedCompound));
                  }}
                >
                  Save
                </button>
              </div>
            </div>
          </div>
        )}

        {/* New Attachment Entry Modal (always rendered, toggled by state) */}
        {showNewEntryModal && (
          <NewAttachmentEntryModal
            keyName={showNewEntryModal.key}
            onAdd={(entry) => {
              setEditedCompound(prev => {
                const prevArr = Array.isArray(prev.attachments[showNewEntryModal.key]) ? prev.attachments[showNewEntryModal.key] : [];
                const updated = {
                  ...prev,
                  attachments: {
                    ...prev.attachments,
                    [showNewEntryModal.key]: [...(prevArr as MultiAttachmentEntry[]), entry]
                  }
                };
                // Always update with full compound
                onUpdate(getFullCompoundUpdate(compound, updated));
                return updated;
              });
            }}
            onClose={() => setShowNewEntryModal(null)}
          />
        )}

        {/* Attachment Modal */}
        {selectedAttachment && (() => {
          const key = selectedAttachment;
          const idx = selectedAttachmentIdx;
          const isMulti = isMultiAttachment(key);
          let data: MultiAttachmentEntry | { note: string; imageUrl: string } | undefined;
          if (isMulti && typeof idx === "number" && Array.isArray(editedCompound.attachments[key])) {
            data = (editedCompound.attachments[key] as MultiAttachmentEntry[])[idx];
          } else {
            data = editedCompound.attachments[key] as { note: string; imageUrl: string };
          }
          return (
            <AttachmentModal
              attachmentKey={key}
              data={data}
              isMulti={isMulti}
              onClose={() => { setSelectedAttachment(null); setSelectedAttachmentIdx(null); }}
              onSave={async (note: string, fileUrl: string, name?: string) => {
                let updatedCompound = { ...editedCompound };
                if (isMulti) {
                  let arr = Array.isArray(editedCompound.attachments[key]) ? [...(editedCompound.attachments[key] as MultiAttachmentEntry[])] : [];
                  if (typeof idx === "number") {
                    arr[idx] = { name: name || arr[idx]?.name || `Entry ${idx+1}`, note, imageUrl: fileUrl };
                  }
                  updatedCompound.attachments = { ...updatedCompound.attachments, [key]: arr };
                } else {
                  updatedCompound.attachments = { ...updatedCompound.attachments, [key]: { note, imageUrl: fileUrl } };
                }
                setEditedCompound(updatedCompound);
                await onUpdate(getFullCompoundUpdate(compound, updatedCompound));
                setSelectedAttachment(null);
                setSelectedAttachmentIdx(null);
              }}
            />
          );
        })()}

        {/* Create Lot Modal */}
        {showCreateLotModal && (
          <CreateLotModal
            compounds={[compound]}
            onClose={() => setShowCreateLotModal(false)}
            onCreate={async () => {
              const res = await fetch(`http://localhost:5000/lots-for-compound/${compound.id}`);
              const updatedLots = await res.json();
              setLotsForCompound(updatedLots);
            }}
          />
        )}

        {/* Full Structure Modal */}
        {showFullStructureModal && (
          <div className="fixed inset-0 bg-black/70 z-50 flex justify-center items-center">
            <div className="bg-white p-8 rounded-2xl border-2 border-[#00E6D2] shadow-2xl max-w-4xl w-full max-h-[92vh] overflow-auto relative">
              <h2 className="text-2xl font-bold mb-4 text-[#002C36] uppercase tracking-wide">Full Structure</h2>
              <div style={{ position: 'absolute', top: 24, right: 32, zIndex: 10 }}>
                <button
                  onClick={() => setShowBWFullStructure((prev) => !prev)}
                  className="px-4 py-2 bg-[#002C36] text-[#00E6D2] rounded-lg font-bold uppercase tracking-wide hover:bg-[#00bfae]"
                  style={{ minWidth: 120 }}
                >
                  {showBWFullStructure ? 'Show Colored' : 'Show B&W'}
                </button>
              </div>
              <FullStructurePreview smiles={compound.smiles} />
              <button
                onClick={() => setShowFullStructureModal(false)}
                className="mt-6 px-6 py-2 bg-[#00E6D2] text-[#002C36] rounded-lg font-bold uppercase tracking-wide hover:bg-[#00bfae]"
              >
                Close
              </button>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}

/*
  Copyright © 2025 Polaris Electro Optics
  This code is the property of Polaris Electro Optics and may not be reused,
  modified, or distributed without explicit permission.
*/