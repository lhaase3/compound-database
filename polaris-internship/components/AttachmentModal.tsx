"use client";
import React, { useState } from "react";

type Props = {
  attachmentKey: string;
  data?: { note: string; imageUrl: string } | import("@/types/compound").MultiAttachmentEntry;
  isMulti?: boolean;
  isNewEntry?: boolean;
  onClose: () => void;
  onSave: (note: string, fileUrl: string, name?: string) => void;
  renderHeaderExtra?: React.ReactNode;
};

export default function AttachmentModal({ attachmentKey, data, isMulti, isNewEntry, onClose, onSave, renderHeaderExtra }: Props) {
  console.log("AttachmentModal opened for:", attachmentKey, data); // Debug log
  const [note, setNote] = useState(data?.note || "");
  const [imageUrl, setImageUrl] = useState(data?.imageUrl || "");
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

  const handleSubmit = async () => {
    setError("");

    try {
      let finalUrl = imageUrl;

      if (file) {
        setUploading(true);
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
        console.log("✅ Uploaded to Drive:", result.fileUrl);
        finalUrl = result.fileUrl;
      }

      onSave(note, finalUrl);
      onClose();
    } catch (err: any) {
      console.error("Upload error:", err);
      setError(err.message || "Upload failed");
    } finally {
      setUploading(false);
    }
  };

  return (
    <div
      className="fixed inset-0 bg-black/50 z-50 flex items-center justify-center p-6"
      onClick={onClose}
    >
      <div
        className="bg-white rounded-lg shadow-xl w-full max-w-7xl h-[90vh] overflow-y-auto p-8"
        onClick={(e) => e.stopPropagation()}
      >
        <div className="flex justify-between items-start mb-6">
          <h2 className="text-2xl font-bold text-black capitalize">
            {attachmentKey.replace("_", " ")}
          </h2>
          {renderHeaderExtra}
        </div>

        {imageUrl ? (
          <img
            src={imageUrl}
            alt="attachment"
            className="w-full max-h-[500px] object-contain mb-6 rounded border"
            onError={(e) => {
              e.currentTarget.style.display = "none";
              console.warn("Image failed to load:", imageUrl);
            }}
          />
        ) : (
          <p className="text-gray-600 italic mb-6">No image uploaded yet.</p>
        )}

        <input
          type="file"
          accept="image/*"
          onChange={handleFileChange}
          className="mb-6 text-black"
        />

        <textarea
          className="w-full border border-gray-300 rounded p-3 text-base text-black mb-6"
          rows={6}
          placeholder="Add a note..."
          value={note}
          onChange={(e) => setNote(e.target.value)}
        />

        {error && <p className="text-red-600 text-sm mb-4">{error}</p>}

        <div className="flex justify-end gap-3">
          <button
            onClick={onClose}
            className="px-5 py-2 rounded bg-gray-300 hover:bg-gray-400 text-black"
          >
            Cancel
          </button>
          <button
            onClick={handleSubmit}
            className="px-5 py-2 rounded bg-blue-600 hover:bg-blue-700 text-white"
            disabled={uploading}
          >
            {uploading ? "Saving..." : "Save"}
          </button>
        </div>
      </div>
    </div>
  );
}