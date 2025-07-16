import React from "react";
import { Compound } from "@/types/compound";

interface CompareAttachment {
  compoundId: string;
  key: string;
  data: any;
}

interface CompareModalProps {
  open: boolean;
  onClose: () => void;
  compounds: Compound[];
  selectedForComparison: string[];
  compareAttachment: CompareAttachment | null;
  setCompareAttachment: (att: CompareAttachment | null) => void;
  clearComparison: () => void;
}

const fieldsToShow = [
  "MW", "Lambda Max (DCM/AcCN)", "Lambda Max (neat film)",
  "phase map", "r33", "dipole CAMB3LYP SVPD CHCl3 (Cosmo)",
  "beta CAMB3LYP SVPD CHCl3 (Cosmo)", "dipole B3LYP SVPD CHCl3",
  "beta B3LYP SVPD CHCl3", "beta/MW", "J/g DSC melt (total)",
  "kJ/mol DSC melt (total)", "Refractive index (ne/no)", "Notes",
  "lab?", "first PEO#", "registered PEO#", "Lab book #", "Max loading (%)"
];

const IMAGE_COL_WIDTH = 300;
const ID_COL_WIDTH = 160;


const CompareModal: React.FC<CompareModalProps> = ({
  open,
  onClose,
  compounds,
  selectedForComparison,
  compareAttachment,
  setCompareAttachment,
  clearComparison,
}) => {
  if (!open) return null;

  const visibleFields = fieldsToShow.filter(field =>
    selectedForComparison.some(id => {
      const cmp = compounds.find(c => c.id === id);
      const val = cmp?.[field];
      return val !== undefined && val !== null && val !== "" && val !== "N/A";
    })
  );

  return (
    <div className="fixed inset-0 bg-black/50 z-50 flex items-center justify-center p-6" onClick={onClose}>
      <div className="bg-white rounded-2xl shadow-2xl max-w-7xl w-full p-8 relative flex flex-col max-h-[90vh]" onClick={e => e.stopPropagation()}>
        {/* Fixed header bar for title and buttons */}
        <div className="sticky top-0 left-0 right-0 z-30 flex items-center justify-between bg-white pb-2" style={{ minHeight: 56 }}>
          <h2 className="text-2xl font-extrabold text-[#008080] uppercase tracking-wide">Compare Compounds</h2>
          <div className="flex gap-2">
            <button
              className="bg-gray-200 hover:bg-gray-300 text-[#002C36] font-bold px-4 py-2 rounded-md uppercase tracking-wide text-xs"
              onClick={() => {
                clearComparison();
                onClose();
              }}
              style={{ position: 'sticky', right: 60, top: 0, zIndex: 40 }}
            >
              Clear
            </button>
            <button
              className="text-2xl text-[#008080] font-bold hover:text-[#00bfae] bg-transparent border-none px-2"
              onClick={onClose}
              style={{ position: 'sticky', right: 0, top: 0, zIndex: 40 }}
            >
              &times;
            </button>
          </div>
        </div>
        <div className="overflow-x-auto overflow-y-auto flex-1 w-full" style={{ maxHeight: '70vh' }}>
          <table className="min-w-full border border-[#008080] rounded-lg mb-6">
            <thead>
              <tr>
                <th
                  className="p-3 text-xs font-bold uppercase text-[#008080] border-b border-r border-[#008080] sticky z-30 bg-white"
                  style={{ left: 0, width: IMAGE_COL_WIDTH, minWidth: IMAGE_COL_WIDTH }}
                >
                  Image
                </th>
                <th
                  className="p-3 text-xs font-bold uppercase text-[#008080] border-b border-r border-[#008080] sticky z-30 bg-white"
                  style={{ left: IMAGE_COL_WIDTH, width: ID_COL_WIDTH, minWidth: ID_COL_WIDTH }}
                >
                  ID
                </th>
                {visibleFields.map(field => (
                  <th key={field} className="p-3 text-xs font-bold uppercase text-[#008080] border-b border-r border-[#008080] bg-white whitespace-nowrap">
                    {field}
                  </th>
                ))}
                <th className="p-3 text-xs font-bold uppercase text-[#008080] border-b border-[#008080] bg-white whitespace-nowrap">Attachments</th>
              </tr>
            </thead>
            <tbody>
              {selectedForComparison.map(id => {
                const cmp = compounds.find(c => c.id === id);
                if (!cmp) return null;

                const atts = cmp.attachments || {};
                const keysWithData = Object.keys(atts).filter(k => {
                  const att = atts[k];
                  if (Array.isArray(att)) {
                    return att.some(entry => entry && (entry.note || entry.imageUrl));
                  } else {
                    return att && (att.note || att.imageUrl);
                  }
                });

                return (
                  <tr key={id} className="border-b border-[#008080]">
                    <td
                      className="p-3 text-center sticky z-20 bg-white"
                      style={{ left: 0, width: IMAGE_COL_WIDTH, minWidth: IMAGE_COL_WIDTH, borderRight: '1px solid #008080' }}
                    >
                      {cmp.imageUrl ? (
                        <img src={cmp.imageUrl} alt={cmp.id} style={{ maxWidth: '300px', maxHeight: '300px', objectFit: 'contain' }} />
                      ) : (
                        <span>N/A</span>
                      )}
                    </td>
                    <td
                      className="p-3 font-bold text-[#002C36] text-center sticky z-20 bg-white"
                      style={{ left: IMAGE_COL_WIDTH, width: ID_COL_WIDTH, minWidth: ID_COL_WIDTH, borderRight: '1px solid #008080' }}
                    >
                      {cmp.id}
                    </td>
                    {visibleFields.map(field => (
                      <td
                        key={field}
                        className={
                          `p-3 text-[#002C36] text-center border-r border-[#008080]` +
                          ((field === 'phase map' || field === 'Notes') ? ' whitespace-pre-line break-words max-w-xs' : ' whitespace-nowrap')
                        }
                        style={
                          (field === 'phase map' || field === 'Notes')
                            ? { whiteSpace: 'pre-line', wordBreak: 'break-word', maxWidth: 600 }
                            : undefined
                        }
                      >
                        {cmp[field] ?? "N/A"}
                      </td>
                    ))}
                    <td className="p-3 text-center border-[#008080]">
                      <div className="flex flex-col gap-2 items-center">
                        {keysWithData.length === 0 && <span>N/A</span>}
                        {keysWithData.map((key) => {
                          const att = atts[key];
                          if (Array.isArray(att)) {
                            return att.map((entry, idx) =>
                              entry && (entry.note || entry.imageUrl) ? (
                                <button
                                  key={key + "-" + idx}
                                  className="px-2 py-1 bg-[#008080] text-white rounded hover:bg-[#006666] font-bold uppercase tracking-wide text-xs mb-1"
                                  onClick={() => setCompareAttachment({ compoundId: id, key, data: entry })}
                                  type="button"
                                >
                                  {key.replace(/_/g, " ")} {entry.name ? `(${entry.name})` : `#${idx + 1}`}
                                </button>
                              ) : null
                            );
                          } else {
                            return (
                              <button
                                key={key}
                                className="px-2 py-1 bg-[#008080] text-white rounded hover:bg-[#006666] font-bold uppercase tracking-wide text-xs mb-1"
                                onClick={() => setCompareAttachment({ compoundId: id, key, data: att })}
                                type="button"
                              >
                                {key.replace(/_/g, " ")}
                              </button>
                            );
                          }
                        })}
                      </div>
                    </td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>
        {/* Attachment Modal */}
        {compareAttachment && (
          <div className="fixed inset-0 bg-black/60 z-50 flex items-center justify-center p-6" onClick={() => setCompareAttachment(null)}>
            <div className="bg-white rounded-lg shadow-xl w-full max-w-7xl max-h-[95vh] p-10 relative flex flex-col items-center justify-center" onClick={e => e.stopPropagation()}>
              <h2 className="text-2xl font-bold mb-6 text-[#008080] uppercase text-center w-full">{compareAttachment.key.replace(/_/g, " ")}</h2>
              {compareAttachment.data.imageUrl && (
                <div className="flex justify-center items-center w-full" style={{ minHeight: '200px' }}>
                  <img
                    src={compareAttachment.data.imageUrl}
                    alt={compareAttachment.key}
                    className="w-auto h-auto max-h-[90vh] max-w-[1600px] object-contain mb-6 border rounded shadow bg-white mx-auto"
                  />
                </div>
              )}
              {compareAttachment.data.note && (
                <div className="bg-gray-100 p-4 rounded whitespace-pre-wrap mb-6 text-[#002C36] max-w-4xl w-full text-center mx-auto">
                  {compareAttachment.data.note}
                </div>
              )}
              {!compareAttachment.data.imageUrl && !compareAttachment.data.note && (
                <div className="text-gray-500 mb-6">No data available.</div>
              )}
              <div className="flex justify-end w-full mt-2">
                <button
                  className="px-8 py-4 bg-[#008080] text-white rounded hover:bg-[#006666] font-bold uppercase tracking-wide text-xl"
                  onClick={() => setCompareAttachment(null)}
                >
                  Close
                </button>
              </div>
            </div>
          </div>
        )}
      </div>
    </div>
  );
};

export default CompareModal;
