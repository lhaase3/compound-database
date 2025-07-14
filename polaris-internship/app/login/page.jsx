"use client";
import LoginModal from "@/components/LoginModal";

export default function TestLoginPage() {
  return (
    <div className="min-h-screen flex items-center justify-center">
      <LoginModal
        onClose={() => console.log("🛑 Modal closed")}
        onLogin={() => console.log("✅ User logged in")}
      />
    </div>
  );
}
