import { useState } from "react";
import { useRouter } from "next/navigation";
import { auth } from "@/utils/firebase";
import { GoogleAuthProvider, signInWithPopup, signOut } from "firebase/auth";

export default function LoginModal({ onClose, onLogin }) {
  const router = useRouter();
  const [error, setError] = useState("");

  const handleGoogleLogin = async () => {
    setError("");
    const provider = new GoogleAuthProvider();
    try {
      const result = await signInWithPopup(auth, provider);
      const email = result.user.email;
      if (!email.endsWith("@polariseo.com")) {
        setError("Only polariseo.com accounts are allowed.");
        await signOut(auth);
        return;
      }
      onLogin();
      onClose();
      router.push("/"); // Redirect to home page
    } catch (err) {
      setError(err.message);
    }
  };


  return (
    <div className="fixed inset-0 bg-gradient-to-br from-[#00E6D2]/30 to-[#002C36]/80 z-50 flex items-center justify-center p-6">
      <div className="bg-white rounded-2xl shadow-2xl w-full max-w-md p-10 relative flex flex-col items-center border-2 border-[#00E6D2]">
        <button
          className="absolute top-4 right-4 text-3xl text-[#008080] font-bold hover:text-[#00bfae] focus:outline-none"
          onClick={onClose}
          aria-label="Close login modal"
        >
          &times;
        </button>
        <img src="/polaris-logo-only.png" alt="Polaris Logo" className="w-16 h-20 mb-4 drop-shadow-lg" />
        <h2 className="text-3xl font-extrabold mb-4 text-[#008080] uppercase text-center w-full tracking-wide">
          Sign in to Polaris
        </h2>
        <p className="text-md text-[#002C36] mb-6 text-center font-medium">Access the Compound Database with your <span className="font-bold text-[#00E6D2]">polariseo.com</span> account.</p>
        {error && <div className="text-red-500 text-sm mb-4 w-full text-center">{error}</div>}
        <button
          onClick={handleGoogleLogin}
          className="w-full bg-gradient-to-r from-[#00E6D2] to-[#00bfae] hover:from-[#00bfae] hover:to-[#00E6D2] text-[#002C36] px-6 py-3 rounded-xl font-bold shadow-lg uppercase tracking-wide mb-2 transition-all flex items-center justify-center gap-2"
        >
          <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 48 48" width="24" height="24" className="mr-2"><path fill="#4285F4" d="M24 9.5c3.54 0 6.73 1.22 9.24 3.22l6.93-6.93C35.64 2.34 30.13 0 24 0 14.61 0 6.27 5.7 2.44 14.02l8.51 6.62C12.7 13.13 17.89 9.5 24 9.5z"/><path fill="#34A853" d="M46.1 24.5c0-1.64-.15-3.22-.43-4.75H24v9h12.44c-.54 2.91-2.18 5.38-4.64 7.04l7.19 5.59C43.73 37.13 46.1 31.27 46.1 24.5z"/><path fill="#FBBC05" d="M12.95 28.64c-1.01-2.99-1.01-6.23 0-9.22l-8.51-6.62C1.64 16.13 0 20.89 0 24.5c0 3.61 1.64 8.37 4.44 11.7l8.51-6.62z"/><path fill="#EA4335" d="M24 46c6.13 0 11.64-2.34 15.17-6.36l-7.19-5.59c-2.01 1.35-4.59 2.15-7.98 2.15-6.11 0-11.3-3.63-13.05-8.62l-8.51 6.62C6.27 42.3 14.61 48 24 48z"/></svg>
          Sign in with Google
        </button>
        <div className="mt-4 text-xs text-[#008080] text-center w-full">Only <span className="font-bold">polariseo.com</span> accounts are allowed.</div>
      </div>
    </div>
  );
}