import React, { useState, useEffect } from "react";
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
  attachments: {
    uv_vis: AttachmentData;
    dsc: AttachmentData;
    lcms: AttachmentData;
    thermal_stability: AttachmentData;
    pda_detector_spectrum: AttachmentData;
    [key: string]: AttachmentData;
  };
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

  const [lotsForCompound, setLotsForCompound] = useState<string[]>([]);
  const [showLotDropdown, setShowLotDropdown] = useState(false);
  const [selectedAttachment, setSelectedAttachment] = useState<"uv_vis" | "dsc" | "lcms" | "thermal_stability" | "pda_detector_spectrum" | null>(null);
  const isLotVersion = source === "lot" && compound.original_id;
  const [showCreateLotModal, setShowCreateLotModal] = useState(false);


  const fields = [
    "id", "MW", "Lambda Max (DCM/Ac CN)", "Lambda Max (neat film)",
    "phase map", "r33", "CAMB3LYP SVPD CHCl3 (Cosmo)",
    "B3LYP SVPD CHCl3 dipole", "B3LYP SVPD CHCl3 beta", "beta/MW",
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

  useEffect(() => {
    const normalizedAttachments = {
      uv_vis: compound.attachments?.["uv-vis_spectrum"] || compound.attachments?.uv_vis || { note: "", imageUrl: "" },
      dsc: compound.attachments?.dsc || { note: "", imageUrl: "" },
      lcms: compound.attachments?.lcms || { note: "", imageUrl: "" },
      thermal_stability: compound.attachments?.thermal_stability || { note: "", imageUrl: "" },
      pda_detector_spectrum: compound.attachments?.pda_detector_spectrum || { note: "", imageUrl: ""}
    };

    setEditedCompound((prev) => ({
      ...compound,
      attachments: normalizedAttachments
    }));
  }, [compound]);



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
                  {compound[field] ?? "N/A"}
                </span>
              )}
            </div>
          ))}
        </div>

        {/* Attachments Section */}
        <div className="mt-6 flex gap-4 justify-center">
          {["UV-vis spectrum", "DSC", "LCMs", "Thermal Stability"].map((type) => (
            <button
              key={type}
              className="text-blue-600 underline hover:text-blue-800"
              onClick={() =>
                setSelectedAttachment(type.toLowerCase().replace(/\s+/g, "_") as "uv_vis" | "dsc" | "lcms" | "thermal_stability" | "pda_detector_spectrum")
              }
            >
              {type}
            </button>
          ))}
        </div>

        {selectedAttachment && (
          <AttachmentModal
            attachmentKey={selectedAttachment}
            data={editedCompound.attachments?.[selectedAttachment] || undefined}
            onClose={() => setSelectedAttachment(null)}
            onSave={(note, fileUrl) => {
              const updatedCompound = {
                ...editedCompound,
                attachments: {
                  ...(editedCompound.attachments || {}),
                  [selectedAttachment]: { note, imageUrl: fileUrl },
                },
              };
              setEditedCompound(updatedCompound);
              onUpdate(updatedCompound); // persists to Firestore
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
      </div>
    </div>
  );
}
