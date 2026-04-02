# -*- coding: utf-8 -*-
"""
hmt-cli — CLI entry point for pyHMT2D AI Tools.

Bridges the stateless-per-process Bash world with the stateful tool functions.
Each call:
  1. Loads .hmt_session.json (if present) to restore serializable state
  2. Auto-detects and loads hmt_config.json for solver paths
  3. Re-opens the project (if the tool requires live session objects)
  4. Reloads results (if the tool requires result data in memory)
  5. Calls the requested tool with provided arguments
  6. Saves updated session state back to .hmt_session.json
  7. Prints JSON result to stdout

Usage
-----
    hmt-cli <tool_name> [--args '{"key": "value"}'] [--session <path>] [--pretty]

Examples
--------
    hmt-cli get_project_info --args '{"project_file": "Muncie.srhhydro"}'
    hmt-cli get_materials
    hmt-cli set_manning_n --args '{"material_ids": [2], "values": [0.035]}'
    hmt-cli get_value_at_point --args '{"x": 411100.0, "y": 1803450.0, "variable": "Water_Elev_ft"}'
"""

from __future__ import annotations

import contextlib
import json
import os
import sys


# ── statefulness classification ───────────────────────────────────────────────

# Dynamic session attributes populated by load_config() via setattr()
# (not declared in the HydraulicSession dataclass — added at runtime)
_DYNAMIC_ATTRS = [
    "_srh_pre_path",
    "_srh_path",
    "_extra_dll_path",
    "_srh_version",
    "_hecras_version",
    "_base_case_dir",
    "_work_dir",
    "_hecras_prj_file",
]

# Tools that require live session.data / session.model objects.
# get_project_info() must be called to rebuild these before the tool runs.
_SESSION_DEPENDENT = {
    "get_materials",
    "get_boundary_conditions",
    "build_param_specs",
    "set_manning_n",
    "set_inlet_flow",
    "set_exit_wse",
    "save_modified_inputs",
    "run_preprocessing",
    "run_simulation",
    "get_simulation_status",
    "read_results",
    "get_value_at_point",
    "get_result_statistics",
    "get_flood_extent",
    "get_cross_section_profile",
    "export_to_vtk",
    "export_mesh_to_vtk",
    "load_observations",
    "evaluate_parameters",
    "run_calibration",
    "run_monte_carlo",
}

# Tools that additionally require results to be loaded in session.data.
# read_results() is called to restore them; vtk_file is restored afterward.
_RESULTS_DEPENDENT = {
    "get_value_at_point",
    "get_result_statistics",
    "get_flood_extent",
    "get_cross_section_profile",
    "export_to_vtk",
}

# Directories to search for hmt_config.json (relative to cwd)
_CONFIG_SEARCH_DIRS = [".", "..", "../.."]

# Default session file name (written in cwd)
_DEFAULT_SESSION = ".hmt_session.json"


# ── session persistence ───────────────────────────────────────────────────────

def _find_config() -> str | None:
    """Search common ancestor directories for hmt_config.json."""
    for d in _CONFIG_SEARCH_DIRS:
        p = os.path.join(d, "hmt_config.json")
        if os.path.isfile(p):
            return os.path.abspath(p)
    return None


def _load_session(session_file: str) -> dict:
    """Return saved session dict, or empty dict if file is absent/corrupt."""
    if os.path.isfile(session_file):
        try:
            with open(session_file) as fh:
                return json.load(fh)
        except Exception:
            pass
    return {}


def _save_session(session_file: str, sess, config_file: str | None = None) -> None:
    """Serialize all session state (declared fields + dynamic path attrs) to JSON.

    All paths are stored as absolute paths so the session file is
    location-independent.
    """
    def _abs(p: str | None) -> str | None:
        return os.path.abspath(p) if p else None

    data: dict = {
        "_comment": "pyHMT2D session state — auto-managed by hmt-cli",
        "model_type": sess.model_type,
        "project_file": _abs(sess.project_file),
        "result_file": _abs(sess.result_file),
        "result_loaded": sess.result_loaded,
        "result_time_steps": list(sess.result_time_steps),
        "result_variables": list(sess.result_variables),
        "vtk_file": _abs(sess.vtk_file),
        "config_file": config_file,
    }
    for attr in _DYNAMIC_ATTRS:
        data[attr] = getattr(sess, attr, None)

    try:
        with open(session_file, "w") as fh:
            json.dump(data, fh, indent=2)
    except OSError as exc:
        print(f"[hmt-cli] Warning: could not save session: {exc}", file=sys.stderr)


def _restore_session(sess, saved: dict, tool_name: str) -> None:
    """Restore session state before dispatching a tool call.

    Three-tier restoration matching the statefulness classification:

    Tier 1 (all tools): restore dynamic path attrs (``_srh_pre_path`` etc.)
    Tier 2 (session-dependent): re-call get_project_info() to rebuild live
        ``session.data`` / ``session.model`` objects.
    Tier 3 (results-dependent): re-call read_results() to reload simulation
        results into ``session.data``; then restore cached ``vtk_file``.
    """

    # Tier 1 — dynamic path attrs (set by load_config via setattr)
    for attr in _DYNAMIC_ATTRS:
        val = saved.get(attr)
        if val is not None:
            setattr(sess, attr, val)

    # Always restore scalar session metadata so stateless tools don't clobber
    # previously saved values when they call _save_session with a fresh session.
    if saved.get("model_type"):
        sess.model_type = saved["model_type"]
    if saved.get("project_file"):
        sess.project_file = saved["project_file"]
    sess.result_loaded = bool(saved.get("result_loaded", False))
    if saved.get("result_file"):
        sess.result_file = saved["result_file"]
    sess.result_time_steps = list(saved.get("result_time_steps") or [])
    sess.result_variables = list(saved.get("result_variables") or [])

    if tool_name not in _SESSION_DEPENDENT:
        return  # stateless tool — nothing more to do

    project_file = saved.get("project_file")
    if not project_file or not os.path.isfile(project_file):
        return  # no restorable project

    # Tier 2 — re-open project to rebuild live data/model objects
    from pyHMT2D.AI_Tools.tools.project_tools import get_project_info
    get_project_info(project_file)

    if tool_name not in _RESULTS_DEPENDENT:
        return

    # Tier 3 — reload results into session.data
    result_file = saved.get("result_file")
    if saved.get("result_loaded") and result_file and os.path.isfile(result_file):
        from pyHMT2D.AI_Tools.tools.result_tools import read_results
        read_results(result_file, timestep=-1)
        # read_results() clears vtk_file; restore cached copy if still on disk
        vtk_file = saved.get("vtk_file")
        if vtk_file and os.path.isfile(vtk_file):
            sess.vtk_file = vtk_file


# ── stdout capture ────────────────────────────────────────────────────────────

class _TeeStream:
    """File-like object that writes to a log file and stderr simultaneously.

    Never writes to stdout, protecting the JSON output channel from solver
    print() calls (same mechanism used in calibration_tools.py).
    """

    def __init__(self, log_path: str) -> None:
        self._log_path = log_path
        self._stderr = sys.stderr

    def write(self, text: str) -> int:
        self._stderr.write(text)
        self._stderr.flush()
        try:
            with open(self._log_path, "a") as fh:
                fh.write(text)
        except OSError:
            pass
        return len(text)

    def flush(self) -> None:
        self._stderr.flush()

    def fileno(self) -> int:
        return self._stderr.fileno()


@contextlib.contextmanager
def _capture_model_output(log_path: str):
    """Redirect sys.stdout → _TeeStream for the duration of the block."""
    tee = _TeeStream(log_path)
    saved = sys.stdout
    sys.stdout = tee
    try:
        yield
    finally:
        sys.stdout = saved


# ── entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        prog="hmt-cli",
        description="pyHMT2D CLI tool runner for Claude Code Skills",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  hmt-cli get_project_info --args '{\"project_file\": \"Muncie.srhhydro\"}'\n"
            "  hmt-cli get_materials\n"
            "  hmt-cli set_manning_n --args '{\"material_ids\": [2], \"values\": [0.035]}'\n"
            "  hmt-cli get_flood_extent --args '{\"depth_threshold\": 0.01}'"
        ),
    )
    parser.add_argument(
        "tool",
        help="Name of the AI tool function to call (e.g. get_project_info)",
    )
    def _parse_json_args(s: str) -> dict:
        # Strip surrounding single quotes added by Windows cmd.exe when the
        # user writes --args '{"key": "value"}' — single quotes are not string
        # delimiters in cmd.exe so they appear literally in the value.
        s = s.strip()
        if s.startswith("'") and s.endswith("'"):
            s = s[1:-1]
        return json.loads(s)

    parser.add_argument(
        "--args",
        dest="tool_args",
        type=_parse_json_args,
        default={},
        metavar="JSON",
        help="Tool arguments as a JSON object string",
    )
    parser.add_argument(
        "--session",
        default=_DEFAULT_SESSION,
        metavar="FILE",
        help=f"Path to session JSON file (default: {_DEFAULT_SESSION})",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Pretty-print the JSON output (useful for debugging)",
    )
    cli_args = parser.parse_args()

    # ── redirect stdout → stderr during imports so nothing corrupts JSON ──────
    _real_stdout = sys.stdout
    sys.stdout = sys.stderr

    from pyHMT2D.AI_Tools.state import get_session
    sess = get_session()

    saved = _load_session(cli_args.session)

    # Auto-detect and load hmt_config.json
    config_file = saved.get("config_file") or _find_config()
    if config_file and os.path.isfile(config_file):
        from pyHMT2D.AI_Tools.tools.project_tools import load_config
        load_config(config_file)

    # Restore session (re-opens project / results as needed for this tool)
    _restore_session(sess, saved, cli_args.tool)

    # Restore stdout before tool dispatch
    sys.stdout = _real_stdout

    # ── dispatch to tool function ─────────────────────────────────────────────
    import pyHMT2D.AI_Tools.tools as _tools_mod

    tool_fn = getattr(_tools_mod, cli_args.tool, None)

    if tool_fn is None:
        # List only public names from __all__
        available = getattr(_tools_mod, "__all__", [])
        result = {
            "status": "error",
            "data": None,
            "message": (
                f"Unknown tool '{cli_args.tool}'. "
                f"Available tools: {available}"
            ),
        }
    else:
        model_log = os.path.join(os.getcwd(), "pyHMT2D.log")
        with _capture_model_output(model_log):
            try:
                result = tool_fn(**cli_args.tool_args)
            except TypeError as exc:
                result = {
                    "status": "error",
                    "data": None,
                    "message": f"Invalid arguments for '{cli_args.tool}': {exc}",
                }
            except Exception as exc:
                result = {
                    "status": "error",
                    "data": None,
                    "message": str(exc),
                }

    # ── persist updated session ───────────────────────────────────────────────
    _save_session(cli_args.session, sess, config_file)

    # ── release HEC-RAS COM handle ────────────────────────────────────────────
    # Each hmt-cli call is a separate process. For HEC-RAS, get_project_info()
    # (called by _restore_session) creates a COM instance that must be explicitly
    # closed via QuitRas() — Python's GC does not do this on exit, leaving
    # HEC-RAS windows open. Close after every call, including run_simulation.
    if sess.model_type == "HEC-RAS" and sess.model is not None:
        try:
            sess.model.exit_model()
        except Exception:
            pass  # best-effort — never block the result

    # ── emit JSON result to stdout ────────────────────────────────────────────
    indent = 2 if cli_args.pretty else None
    print(json.dumps(result, indent=indent, default=str))


if __name__ == "__main__":
    main()
