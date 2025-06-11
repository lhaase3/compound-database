
'use client';

import { useEffect, useState } from 'react';

export default function JsmeEditor() {
  const [smiles, setSmiles] = useState('');

  useEffect(() => {
    window.jsmeOnLoad = function () {
      console.log('✅ JSME fully loaded');

      const jsmeApplet = new window.JSApplet.JSME('jsme_container', '600px', '600px', {
        options: 'oldlook,star',
      });

      setTimeout(() => {
        jsmeApplet.setCallBack('AfterStructureModified', () => {
          const s = jsmeApplet.smiles();
          console.log('SMILES:', s);
          setSmiles(s);
        });
      }, 500);
    };

    const script = document.createElement('script');
    script.src = '/jsme/jsme.nocache.js';
    script.onload = () => console.log("✅ JSME script loaded");
    script.onerror = () => console.error("❌ Failed to load JSME script");
    document.body.appendChild(script);

    return () => {
      delete window.jsmeOnLoad;
    };
  }, []);

  return (
    <div className="p-6">
      <h1 className="text-2xl font-bold mb-4">Drawing Tester-V6</h1>
      <div id="jsme_container" className="mb-6" />
      <p className="text-black font-medium">💡 Current SMILES:</p>
      <pre className="bg-gray-100 text-black p-3 rounded-md border border-gray-300">
        {smiles || 'Nothing drawn yet...'}
      </pre>
    </div>
  );
}
