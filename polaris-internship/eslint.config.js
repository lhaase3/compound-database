// // polaris-internship/eslint.config.js
// import js from "@eslint/js";
// import ts from "@typescript-eslint/eslint-plugin";
// import tsParser from "@typescript-eslint/parser";
// import react from "eslint-plugin-react";
// import path from "path";

// export default [
//   js.configs.recommended,
//   {
//     files: ["**/*.ts", "**/*.tsx"],
//     languageOptions: {
//       parser: tsParser,
//       parserOptions: {
//         project: "./tsconfig.json",
//         tsconfigRootDir: path.resolve(),
//       },
//       globals: {
//         window: true,
//         document: true,
//         alert: true,
//         fetch: true,
//         console: true,
//         Image: true,
//       }
//     },
//     plugins: {
//       "@typescript-eslint": ts,
//       "react": react,
//     },
//     rules: {
//       "no-unused-vars": "warn",
//       "react/react-in-jsx-scope": "off"
//     }
//   }
// ];


// polaris-internship/eslint.config.js
import js from "@eslint/js";
import ts from "@typescript-eslint/eslint-plugin";
import tsParser from "@typescript-eslint/parser";
import react from "eslint-plugin-react";
import path from "path";

export default [
  js.configs.recommended,
  {
    files: ["**/*.ts", "**/*.tsx"],
    languageOptions: {
      parser: tsParser,
      parserOptions: {
        project: "./tsconfig.json",
        tsconfigRootDir: path.resolve(),
      },
      globals: {
        window: true,
        document: true,
        alert: true,
        fetch: true,
        console: true,
        Image: true,
        process: true,        // ✅ Fixes 'process is not defined'
        localStorage: true,   // ✅ Fixes 'localStorage is not defined'
        setTimeout: true,      // ✅ Fixes 'setTimeout is not defined'
        prompt: true,        // ✅ Fixes 'prompt is not defined'
      }
    },
    plugins: {
      "@typescript-eslint": ts,
      "react": react,
    },
    rules: {
      "no-unused-vars": "warn",
      "react/react-in-jsx-scope": "off"
    }
  }
];
