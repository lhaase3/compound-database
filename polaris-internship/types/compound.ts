export type AttachmentData = {
  note: string;
  imageUrl: string;
};

export type Compound = {
  id: string;
  name?: string;
  formula?: string;
  smiles: string;
  [key: string]: any;
  attachments?: {
    uv_vis?: AttachmentData;
    dsc?: AttachmentData;
    lcms?: AttachmentData;
    [key: string]: AttachmentData | undefined;
  };
};
