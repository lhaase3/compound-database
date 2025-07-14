// Firebase config and initialization
import { initializeApp } from "firebase/app";
import { getAuth } from "firebase/auth";

const firebaseConfig = {
  apiKey: "AIzaSyCfgiXJHh9ecMfpmkMTm4XOnRqud1KQNHk",
  authDomain: "polaris-test-3b8f8.firebaseapp.com",
  projectId: "polaris-test-3b8f8",
  storageBucket: "polaris-test-3b8f8.firebasestorage.app",
  messagingSenderId: "1057420182928",
  appId: "1:1057420182928:web:17e46b45f681de39af5f15",
  measurementId: "G-NHFJEQFX0X"
};

const app = initializeApp(firebaseConfig);
export const auth = getAuth(app);
// Firestore is not used, so do not export db
