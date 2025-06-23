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


  // List of unwanted fields to exclude from custom display
  const unwanted = [
    "attachments", "lots", "original_id", "created_at", "updated_at", "_id", "__v", "parsed_phase_transitions", "imageUrl"
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

  const handleChange = (field: string, value: string) => {
    setEditedCompound((prev) => ({ ...prev, [field]: value }));
  };

  const handleSave = async () => {
    try {
      await onUpdate(editedCompound);
      setEditMode(false);
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

  // Normalize all attachments, including custom fields, and always include built-in keys
  useEffect(() => {
    const builtIn = [
      "uv_vis",
      "dsc",
      "lcms",
      "thermal_stability",
      "pda_detector_spectrum",
    ];
    const allKeys = compound.attachments ? Array.from(new Set([...Object.keys(compound.attachments), ...builtIn])) : builtIn;
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
    } as CompoundWithAttachments));
  }, [compound]);


  // Helper: get all attachment keys (built-in + custom)
  const getAllAttachmentKeys = () => Object.keys(editedCompound.attachments || {});

  // Add Attachment Field
  const handleAddAttachmentField = () => {
    const name = prompt("Enter a name for the new attachment field:");
    if (
      name &&
      !Object.keys(editedCompound.attachments || {}).includes(name)
    ) {
      setEditedCompound((prev) => ({
        ...prev,
        attachments: {
          ...prev.attachments,
          [name]: { note: "", imageUrl: "" },
        },
      }));
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
      className="fixed inset-0 bg-black/40 backdrop-blur-sm z-50 flex items-center justify-center"
      onClick={onClose}
    >
      <div
        className="bg-white p-6 rounded-lg shadow-xl max-w-2xl w-full max-h-[90vh] overflow-y-auto"
        onClick={(e) => e.stopPropagation()}
      >
        <div className="flex justify-between items-center mb-4">
          <h2 className="text-2xl font-bold text-black">
            Compound Details{" "}
            {source === "lot" && lotId && (
              <span className="text-purple-600 text-sm ml-2">(from Lot: {lotId})</span>
            )}
          </h2>
          <div className="flex gap-2">
            <div className="relative">
              <button
                className="text-purple-600 border border-purple-600 px-2 py-1 rounded cursor-pointer"
                onClick={() => setShowLotDropdown((prev) => !prev)}
              >
                📦 Lots ({lotsForCompound.length})
              </button>

              {showLotDropdown && (
                <div className="absolute bg-white shadow-md border border-gray-200 mt-2 right-0 z-50 rounded w-65">
                  <div
                    onClick={() => {
                      setShowLotDropdown(false);
                      setShowCreateLotModal(true);
                    }}
                    className="px-3 py-2 text-blue-600 hover:bg-gray-100 cursor-pointer border-b border-gray-200"
                  >
                    ➕ Create Lot
                  </div>

                  {lotsForCompound.map((lotId) => (

                    <div
                      key={lotId}
                      className="px-3 py-2 text-black hover:bg-gray-100 cursor-pointer"
                      onClick={async () => {
                        try {
                          const res = await fetch(`http://localhost:5000/lot/${lotId}`);
                          const data = await res.json();
                          if (Array.isArray(data) && data.length > 0) {
                            const compoundFromLot = data[0];
                            setEditedCompound({
                              ...compoundFromLot,
                              attachments: {
                                uv_vis: {
                                  note: compoundFromLot.attachments?.uv_vis?.note || "",
                                  imageUrl: compoundFromLot.attachments?.uv_vis?.imageUrl || "",
                                },
                                dsc: {
                                  note: compoundFromLot.attachments?.dsc?.note || "",
                                  imageUrl: compoundFromLot.attachments?.dsc?.imageUrl || "",
                                },
                                lcms: {
                                  note: compoundFromLot.attachments?.lcms?.note || "",
                                  imageUrl: compoundFromLot.attachments?.lcms?.imageUrl || "",
                                },
                                thermal_stability: {
                                  note: compound.attachments?.thermal_stability?.note || "",
                                  imageUrl: compound.attachments?.thermal_stability?.imageUrl || "",
                                },
                                pda_detector_spectrum: {
                                  note: compound.attachments?.pda_detector_spectrum?.note || "",
                                  imageUrl: compound.attachments?.pda_detector_spectrum?.imageUrl || "",
                                },
                              },
                            });
                            setShowLotDropdown(false);
                            onUpdateCompoundFromLot?.(compoundFromLot, lotId);
                          } else {
                            alert("No data found for that lot.");
                          }
                        } catch (err) {
                          console.error("Failed to load lot:", err);
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
              className="text-blue-600 border border-blue-600 text-sm px-2 py-1 rounded hover:bg-blue-100 cursor-pointer"
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

                    // reset context if needed
                    onUpdateCompoundFromLot?.(original, null);
                  } catch (err) {
                    console.error("Failed to load original compound", err);
                  }
                }}
              >
                ⬅️ Back to Original
              </button>
            )}

            <button
              className="text-blue-600 border border-blue-600 px-2 py-1 rounded cursor-pointer"
              onClick={() => (editMode ? handleSave() : setEditMode(true))}
            >
              {editMode ? "Save" : "Edit"}
            </button>
            <button className="text-red-600 border border-red-600 px-2 py-1 rounded cursor-pointer" onClick={handleDelete}>
              Delete
            </button>
            <button onClick={onClose} className="text-black text-xl font-bold">
              &times;
            </button>
          </div>
        </div>

        <div className="grid grid-cols-1 sm:grid-cols-2 gap-4 text-black">
          {fields.map((field) => (
            <div key={field} className="flex flex-col">
              <span className="text-sm font-semibold text-gray-600">{field}</span>
              {editMode ? (
                <input
                  className="border rounded px-2 py-1 text-sm"
                  value={editedCompound[field] ?? ""}
                  onChange={(e) => handleChange(field, e.target.value)}
                />
              ) : (
                <span className={compound[field] ? "text-black" : "text-gray-400"}>
                  {compound[field] ? compound[field] : "N/A"}
                </span>
              )}
            </div>
          ))}
            {Object.keys(editedCompound)
              .filter((key) => !fields.includes(key) && !unwanted.includes(key))
              .map((key) => (
                <div key={key} className="flex flex-col">
                  <span className="text-sm font-semibold text-gray-600">{key}</span>
                  {editMode ? (
                    <input
                      className="border rounded px-2 py-1 text-sm"
                      value={editedCompound[key] ?? ""}
                      onChange={(e) => handleChange(key, e.target.value)}
                    />
                  ) : (
                    <span className={compound[key] ? "text-black" : "text-gray-400"}>
                      {compound[key] ? compound[key] : "N/A"}
                    </span>
                  )}
                </div>
            ))}
        </div>
        {/* Attachments Section */}
        <div className="mt-6 flex gap-4 flex-wrap justify-center">
          {getAllAttachmentKeys().map((key) => (
            <button
              key={key}
              className="text-blue-600 underline hover:text-blue-800 mb-1"
              onClick={() => setSelectedAttachment(key)}
            >
              {key.replace(/_/g, " ")}
            </button>
          ))}
        </div>

        <div className="mt-6 text-center">
          <div
            className="bg-white p-2 rounded-lg border shadow hover:shadow-md transition-all cursor-pointer"
            style={{
              width: 180,
              height: 140,
              display: "flex",
              alignItems: "center",
              justifyContent: "center",
              overflow: "hidden",
              margin: "0 auto",
            }}
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
                    style={{
                      display: "flex",
                      alignItems: "center",
                      justifyContent: "center",
                      height: "100%",
                      transform: "scale(0.6)",
                      transformOrigin: "center",
                    }}
                  />
                  ) : (
                    <span className="text-sm text-gray-500">Invalid SMILES</span>
                  );
                } catch {
                  return <span className="text-sm text-gray-500">Invalid SMILES</span>;
                }
              })() : (
                <span className="text-sm text-gray-500">Loading...</span>
              )
            )}
          </div>
        </div>




        {/* Add Field Buttons */}
        <div className="mt-8 flex gap-4 justify-center">
          <button
            className="px-4 py-2 bg-blue-500 text-white rounded hover:bg-blue-700 cursor-pointer"
            onClick={handleAddAttachmentField}
          >
            Add Attachment Field
          </button>
          <button
            className="px-4 py-2 bg-blue-500 text-white rounded hover:bg-blue-700 cursor-pointer"
            onClick={() => setShowAddTextField(true)}
          >
            Add Text Field
          </button>
        </div>

        {showAddTextField && (
          <div className="fixed inset-0 bg-black/50 z-50 flex items-center justify-center p-6" onClick={() => setShowAddTextField(false)}>
            <div className="bg-white rounded-lg shadow-xl w-full max-w-md p-8" onClick={e => e.stopPropagation()}>
              <h2 className="text-xl font-bold mb-4">Add Text Field</h2>
              <input
                type="text"
                className="w-full border border-gray-300 rounded p-2 mb-4 text-black"
                placeholder="Field name"
                value={newTextFieldName}
                onChange={e => setNewTextFieldName(e.target.value)}
              />
              <textarea
                className="w-full border border-gray-300 rounded p-2 mb-4 text-black"
                placeholder="Field value"
                value={newTextFieldValue}
                onChange={e => setNewTextFieldValue(e.target.value)}
              />
              <div className="flex justify-end gap-2">
                <button
                  className="px-4 py-2 bg-gray-300 rounded hover:bg-gray-400 text-black"
                  onClick={() => setShowAddTextField(false)}
                >
                  Cancel
                </button>
                <button
                  className="px-4 py-2 bg-blue-600 text-white rounded hover:bg-blue-700"
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
                    await onUpdate(updatedCompound);
                  }}
                >
                  Save
                </button>
              </div>
            </div>
          </div>
        )}
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
              // Immediately persist the attachment update to Firestore
              await onUpdate(updatedCompound);
            }}
          />
        )}
        {showCreateLotModal && (
          <CreateLotModal
            compounds={[compound]}  // Only allow creating a lot with this compound
            onClose={() => setShowCreateLotModal(false)}
            onCreate={async () => {
              // Refresh lots after new one is created
              const res = await fetch(`http://localhost:5000/lots-for-compound/${compound.id}`);
              const updatedLots = await res.json();
              setLotsForCompound(updatedLots);
            }}
          />
        )}
        {showFullStructureModal && (
          <div className="fixed inset-0 bg-opacity-80 z-50 flex justify-center items-center">
            <div className="bg-white p-6 rounded shadow-lg max-w-4xl w-full max-h-[90vh] overflow-auto">
              <h2 className="text-lg font-bold mb-2">Full Structure</h2>
              <FullStructurePreview smiles={compound.smiles} />
              <button
                onClick={() => setShowFullStructureModal(false)}
                className="mt-4 px-4 py-2 bg-blue-600 text-white rounded"
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

