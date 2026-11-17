// Run a Python benchmark file inside Pyodide (Node), with the project mounted at /proj.
// Usage: [PYODIDE_PACKAGES=pytest] node run.mjs <bench.py> [project_root]
import { loadPyodide } from "pyodide";
import { readFileSync } from "node:fs";
import path from "node:path";

const [, , pyFile, rootArg] = process.argv;
const root = path.resolve(rootArg ?? "../..");
const t0 = performance.now();
const py = await loadPyodide();
const extra = (process.env.PYODIDE_PACKAGES ?? "").split(",").filter(Boolean);
await py.loadPackage(["numpy", "pandas", ...extra], { messageCallback: () => {} });
console.log(`[pyodide ${py.version}] startup + numpy + pandas: ${((performance.now() - t0) / 1000).toFixed(1)} s`);
globalThis.__pyodideHeapMB = () => Math.round(py._module.HEAP8.length / 1e6);
py.FS.mkdirTree("/proj");
py.FS.mount(py.FS.filesystems.NODEFS, { root }, "/proj");
await py.runPythonAsync(readFileSync(pyFile, "utf8"));
