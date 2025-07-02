// import React, { useEffect, useState } from "react";


// export default function FormulationList() {
//   const [formulations, setFormulations] = useState<any[]>([]);
//   const [selectedFormulation, setSelectedFormulation] = useState<any | null>(null);
//   const [selectedAttachment, setSelectedAttachment] = useState<{
//     name: string;
//     data: { note: string; imageUrl: string };
//   } | null>(null);
//   // New state for custom fields and attachments
//   const [customFields, setCustomFields] = useState<{ [key: string]: string }>({});
//   const [attachments, setAttachments] = useState<{ [key: string]: { note: string; imageUrl: string } }>({});
//   const [showAddTextField, setShowAddTextField] = useState(false);
//   const [newTextFieldName, setNewTextFieldName] = useState("");
//   const [newTextFieldValue, setNewTextFieldValue] = useState("");
//   const [showAddAttachmentField, setShowAddAttachmentField] = useState(false);
//   const [newAttachmentFieldName, setNewAttachmentFieldName] = useState("");
//   const [newAttachmentNote, setNewAttachmentNote] = useState("");
//   const [newAttachmentUrl, setNewAttachmentUrl] = useState("");


//   useEffect(() => {
//     fetch("http://localhost:5000/formulations")
//       .then((res) => res.json())
//       .then(setFormulations)
//       .catch((err) => console.error("Failed to fetch formulations", err));
//   }, []);


//   return (
//     <div className="p-6 text-black">
//       <h1 className="text-3xl font-bold mb-4">Formulations</h1>
//       {!selectedFormulation ? (
//         <ul className="space-y-2">
//           {formulations.map((form) => (
//             <li
//               key={form.id}
//               className="border p-4 rounded cursor-pointer hover:bg-gray-100"
//               onClick={() => setSelectedFormulation(form)}
//             >
//               <span className="font-semibold">{form.name || "Unnamed Formulation"}</span>
//             </li>
//           ))}
//         </ul>
//       ) : (
//         <div className="border p-4 rounded">
//           <button
//             className="text-blue-600 underline mb-4"
//             onClick={() => setSelectedFormulation(null)}
//           >
//             ← Back to List
//           </button>

//           <h2 className="text-2xl font-semibold mb-2">{selectedFormulation.name || "Unnamed"}</h2>
//           <p className="text-sm text-gray-600 mb-2">
//             Created: {new Date(selectedFormulation.createdAt?.seconds * 1000).toLocaleString()}
//           </p>
//           {selectedFormulation.imageUrls?.length > 0 && (
//             <div className="mb-4">
//               <h3 className="font-semibold">Images:</h3>
//               <div className="flex flex-wrap gap-4 mt-2">
//                 {selectedFormulation.imageUrls.map((url: string, idx: number) => (
//                   <img
//                     key={idx}
//                     src={url}
//                     alt={`Formulation Image ${idx + 1}`}
//                     className="w-40 h-40 object-contain border rounded shadow"
//                   />
//                 ))}
//               </div>
//             </div>
//           )}


//           <h3 className="font-semibold mt-4">Components:</h3>
//             <ul className="mb-4 list-disc pl-5">
//               {selectedFormulation.components.map((comp: any, idx: number) => (
//                 <li key={idx}>
//                   {comp.compoundId} ({comp.lotId || "original"}) – {comp.molPercent}% → {comp.mass} g
//                 </li>
//               ))}
//             </ul>



//           <div className="mb-2">
//             <strong>Total Moles:</strong> {selectedFormulation.totalMoles}
//           </div>

//           <div className="mb-2">
//             <strong>Phase Map:</strong>
//             <div className="whitespace-pre-line bg-gray-100 p-2 rounded mt-1">
//               {selectedFormulation.phaseMap || "N/A"}
//             </div>
//           </div>

//           <div className="mb-2">
//             <strong>Notes:</strong>
//             <div className="whitespace-pre-line bg-gray-100 p-2 rounded mt-1">
//               {selectedFormulation.notes || "N/A"}
//             </div>
//           </div>

//           {/* Custom Text Fields */}
//           {Object.entries(customFields).length > 0 && (
//             <div className="mb-2">
//               <strong>Custom Fields:</strong>
//               <ul className="list-disc pl-5">
//                 {Object.entries(customFields).map(([key, value]) => (
//                   <li key={key}><strong>{key}:</strong> {value}</li>
//                 ))}
//               </ul>
//             </div>
//           )}

//           {/* Attachments */}
//           {Object.keys(attachments).length > 0 && (
//             <div className="mb-2">
//               <strong>Attachments:</strong>
//               <ul className="list-disc pl-5">
//                 {Object.entries(attachments).map(([key, att]) => (
//                   <li key={key}>
//                     <span className="font-semibold cursor-pointer text-blue-600 underline" onClick={() => setSelectedAttachment({ name: key, data: att })}>{key}</span>
//                     {att.note && <span className="ml-2 text-gray-600">({att.note})</span>}
//                   </li>
//                 ))}
//               </ul>
//             </div>
//           )}

//           {/* Add Text Field Button */}
//           <button
//             className="px-2 py-1 bg-[#008080] text-white rounded-md hover:bg-[#006666] font-bold uppercase tracking-wide text-xs mr-2"
//             onClick={() => setShowAddTextField(true)}
//           >
//             + Add Text Field
//           </button>
//           {/* Add Attachment Field Button */}
//           <button
//             className="px-2 py-1 bg-[#008080] text-white rounded-md hover:bg-[#006666] font-bold uppercase tracking-wide text-xs"
//             onClick={() => setShowAddAttachmentField(true)}
//           >
//             + Add Attachment Field
//           </button>

//           {/* Add Text Field Modal */}
//           {showAddTextField && (
//             <div className="fixed inset-0 bg-black/40 z-50 flex items-center justify-center" onClick={() => setShowAddTextField(false)}>
//               <div className="bg-white p-6 rounded shadow-lg" onClick={e => e.stopPropagation()}>
//                 <h3 className="font-bold mb-2">Add Text Field</h3>
//                 <input
//                   type="text"
//                   placeholder="Field Name"
//                   className="border p-2 rounded mb-2 w-full"
//                   value={newTextFieldName}
//                   onChange={e => setNewTextFieldName(e.target.value)}
//                 />
//                 <input
//                   type="text"
//                   placeholder="Field Value"
//                   className="border p-2 rounded mb-2 w-full"
//                   value={newTextFieldValue}
//                   onChange={e => setNewTextFieldValue(e.target.value)}
//                 />
//                 <div className="flex gap-2 justify-end">
//                   <button className="px-3 py-1 bg-gray-200 rounded" onClick={() => setShowAddTextField(false)}>Cancel</button>
//                   <button
//                     className="px-3 py-1 bg-[#008080] text-white rounded"
//                     onClick={() => {
//                       if (newTextFieldName.trim()) {
//                         setCustomFields(prev => ({ ...prev, [newTextFieldName]: newTextFieldValue }));
//                         setNewTextFieldName("");
//                         setNewTextFieldValue("");
//                         setShowAddTextField(false);
//                       }
//                     }}
//                   >Add</button>
//                 </div>
//               </div>
//             </div>
//           )}

//           {/* Add Attachment Field Modal */}
//           {showAddAttachmentField && (
//             <div className="fixed inset-0 bg-black/40 z-50 flex items-center justify-center" onClick={() => setShowAddAttachmentField(false)}>
//               <div className="bg-white p-6 rounded shadow-lg" onClick={e => e.stopPropagation()}>
//                 <h3 className="font-bold mb-2">Add Attachment Field</h3>
//                 <input
//                   type="text"
//                   placeholder="Attachment Name"
//                   className="border p-2 rounded mb-2 w-full"
//                   value={newAttachmentFieldName}
//                   onChange={e => setNewAttachmentFieldName(e.target.value)}
//                 />
//                 <input
//                   type="text"
//                   placeholder="Note (optional)"
//                   className="border p-2 rounded mb-2 w-full"
//                   value={newAttachmentNote}
//                   onChange={e => setNewAttachmentNote(e.target.value)}
//                 />
//                 <input
//                   type="text"
//                   placeholder="Image URL"
//                   className="border p-2 rounded mb-2 w-full"
//                   value={newAttachmentUrl}
//                   onChange={e => setNewAttachmentUrl(e.target.value)}
//                 />
//                 <div className="flex gap-2 justify-end">
//                   <button className="px-3 py-1 bg-gray-200 rounded" onClick={() => setShowAddAttachmentField(false)}>Cancel</button>
//                   <button
//                     className="px-3 py-1 bg-[#008080] text-white rounded"
//                     onClick={() => {
//                       if (newAttachmentFieldName.trim()) {
//                         setAttachments(prev => ({ ...prev, [newAttachmentFieldName]: { note: newAttachmentNote, imageUrl: newAttachmentUrl } }));
//                         setNewAttachmentFieldName("");
//                         setNewAttachmentNote("");
//                         setNewAttachmentUrl("");
//                         setShowAddAttachmentField(false);
//                       }
//                     }}
//                   >Add</button>
//                 </div>
//               </div>
//             </div>
//           )}

//           {/* Attachment Viewer Modal */}
//           {selectedAttachment && (
//             <div className="fixed inset-0 bg-black/40 z-50 flex items-center justify-center" onClick={() => setSelectedAttachment(null)}>
//               <div className="bg-white p-6 rounded shadow-lg" onClick={e => e.stopPropagation()}>
//                 <h3 className="font-bold mb-2">{selectedAttachment.name}</h3>
//                 {selectedAttachment.data.note && <p className="mb-2">Note: {selectedAttachment.data.note}</p>}
//                 {selectedAttachment.data.imageUrl && (
//                   <img src={selectedAttachment.data.imageUrl} alt={selectedAttachment.name} className="w-80 h-80 object-contain border rounded shadow" />
//                 )}
//                 <div className="flex gap-2 justify-end mt-4">
//                   <button className="px-3 py-1 bg-gray-200 rounded" onClick={() => setSelectedAttachment(null)}>Close</button>
//                 </div>
//               </div>
//             </div>
//           )}
//         </div>
//       )}
//     </div>
//   );
  
// }
