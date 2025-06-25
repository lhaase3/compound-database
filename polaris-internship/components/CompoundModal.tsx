import React, { useState, useEffect, useRef } from "react";
import { Compound } from "@/types/compound";
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

type AttachmentData = {
  note: string;
  imageUrl: string;
};

type CompoundWithAttachments = Compound & {
  attachments: { [key: string]: AttachmentData };
};

export default function CompoundModal({
  compound,
  onClose,
  onDelete,
  onUpdate,
  source,
  lotId,
  onUpdateCompoundFromLot,
}: Props) {
  const [editMode, setEditMode] = useState(false);
  const [editedCompound, setEditedCompound] = useState<CompoundWithAttachments>(() => ({
    ...compound,
    attachments: {
      uv_vis: {
        note: compound.attachments?.uv_vis?.note || "",
        imageUrl: compound.attachments?.uv_vis?.imageUrl || "",
      },
      dsc: {
        note: compound.attachments?.dsc?.note || "",
        imageUrl: compound.attachments?.dsc?.imageUrl || "",
      },
      lcms: {
        note: compound.attachments?.lcms?.note || "",
        imageUrl: compound.attachments?.lcms?.imageUrl || "",
      },
      thermal_stability: {
        note: compound.attachments?.thermal_stability?.note || "",
        imageUrl: compound.attachments?.thermal_stability?.imageUrl || "",
      },
      pda_detector_spectrum: {
        note: compound.attachments?.pda_detector_spectrum?.note || "",
        imageUrl: compound.attachments?.pda_detector_spectrum?.imageUrl || "",
      }
    },
  }));

  // Track user-added custom fields
  const [lotsForCompound, setLotsForCompound] = useState<string[]>([]);
  const [showLotDropdown, setShowLotDropdown] = useState(false);
  const [selectedAttachment, setSelectedAttachment] = useState<
    "uv_vis" | "dsc" | "lcms" | "thermal_stability" | "pda_detector_spectrum" | string | null
  >(null);
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


  // List of unwanted fields to exclude from custom display
  const unwanted = [
    "attachments", "lots", "original_id", "created_at", "updated_at", "_id", "__v", "parsed_phase_transitions", "imageUrl", "tags", "lotId", "Lambda Max (DCM/Ac CN)"
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

  useEffect(() => {
  setCurrentLotId(lotId);
}, [lotId]);


  const handleChange = (field: string, value: string) => {
    setEditedCompound((prev) => ({ ...prev, [field]: value }));
  };

  const handleSave = async () => {
    try {
      if (source === "lot" && currentLotId) {
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
        await onUpdate(editedCompound);
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

  useEffect(() => {
    const builtIn = [
      "uv_vis",
      "dsc",
      "lcms",
      "thermal_stability",
      "pda_detector_spectrum",
    ];
    const allKeys = compound.attachments
      ? Array.from(new Set([...Object.keys(compound.attachments), ...builtIn]))
      : builtIn;

    const normalizedAttachments: { [key: string]: AttachmentData } = {};
    allKeys.forEach((key) => {
      const att = compound.attachments?.[key] || { note: "", imageUrl: "" };
      normalizedAttachments[key] = {
        note: att.note || "",
        imageUrl: att.imageUrl || "",
      };
    });

    setEditedCompound((prev) => ({
      ...compound,
      attachments: normalizedAttachments,
    }));
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




  // Helper: get all attachment keys (built-in + custom)
  const getAllAttachmentKeys = () => Object.keys(editedCompound.attachments || {});

  // Add Attachment Field
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
      await onUpdate(updatedCompound); // This is critical to immediately save the change
    }
  };


  const FullStructurePreview = ({ smiles }: { smiles: string }) => {
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
    }, [smiles]);
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
                  className="absolute left-0 bg-white shadow-lg border border-[#00E6D2] mt-2 z-50 rounded w-44"
                >
                  {["testing", "crystals"].map((tag) => (
                    <div
                      key={tag}
                      className="px-4 py-2 text-[#002C36] hover:bg-[#00E6D2]/20 cursor-pointer font-semibold"
                      onClick={async () => {
                        const updatedTags = editedCompound.tags?.includes(tag)
                          ? editedCompound.tags.filter((t: string) => t !== tag)
                          : [...(editedCompound.tags || []), tag];
                        const updatedCompound = { ...editedCompound, tags: updatedTags };
                        setEditedCompound(updatedCompound);
                        setShowTagDropdown(false);
                        await onUpdate(updatedCompound);
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
                <div ref={lotDropdownRef} className="absolute left-0 bg-white shadow-lg border border-purple-400 mt-2 z-50 rounded w-56">
                  <div
                    onClick={() => {
                      setShowLotDropdown(false);
                      setShowCreateLotModal(true);
                    }}
                    className="px-4 py-2 text-blue-600 hover:bg-blue-50 cursor-pointer border-b border-purple-200 font-semibold"
                  >
                    ➕ Create Lot
                  </div>
                  {lotsForCompound.map((lotId) => (
                    <div
                      key={lotId}
                      className="px-4 py-2 text-[#002C36] hover:bg-purple-50 cursor-pointer font-semibold"
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
              onClick={() => (editMode ? handleSave() : setEditMode(true))}
            >
              {editMode ? "Save" : "Edit"}
            </button>
            <button className="bg-red-100 hover:bg-red-200 text-red-600 font-bold px-2 py-0.5 rounded-md uppercase tracking-wide text-xs shadow transition-all" onClick={handleDelete}>
              Delete
            </button>
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
            {getAllAttachmentKeys().map((key) => (
              <button
                key={key}
                className="border border-[#00E6D2] text-[#00E6D2] font-bold px-2 py-0.5 rounded-md uppercase tracking-wide text-xs hover:bg-[#00E6D2] hover:text-[#002C36] transition-all mb-1"
                onClick={() => setSelectedAttachment(key)}
              >
                {key.replace(/_/g, " ")}
              </button>
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
                    await onUpdate(updatedCompound); // This is what you are missing in your new version
                  }}
                >
                  Save
                </button>
              </div>
            </div>
          </div>
        )}

        {/* Attachment Modal */}
        {selectedAttachment && (
          <AttachmentModal
            attachmentKey={selectedAttachment}
            data={editedCompound.attachments?.[selectedAttachment] || undefined}
            onClose={() => setSelectedAttachment(null)}
            onSave={async (note, fileUrl) => {
              const updatedCompound = {
                ...editedCompound,
                attachments: {
                  ...(editedCompound.attachments || {}),
                  [selectedAttachment]: { note, imageUrl: fileUrl },
                },
              };
              setEditedCompound(updatedCompound);
              await onUpdate(updatedCompound);
            }}
          />
        )}

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