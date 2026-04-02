# -*- coding: utf-8 -*-
"""
Group 3 — Model execution tools.

run_preprocessing / run_simulation require the SRH-2D or HEC-RAS solver to be
installed on the machine. Paths to the solver executables must be provided at
run time (the session stores them after the first successful call).
"""

from __future__ import annotations

import os
import time
from typing import Optional

from pyHMT2D.AI_Tools.state import get_session


def _ok(data, message: str) -> dict:
    return {"status": "ok", "data": data, "message": message}


def _err(message: str) -> dict:
    return {"status": "error", "data": None, "message": message}


# ── tool 10 ───────────────────────────────────────────────────────────────────

def run_preprocessing(
    srh_pre_path: Optional[str] = None,
    version: str = "3.3",
    extra_dll_path: Optional[str] = None,
) -> dict:
    """Run the SRH-2D pre-processor (srh_2d_pre) for the open project.

    The pre-processor converts the .srhhydro / SIF control file and .srhgeom
    mesh into the binary format that the SRH-2D solver expects. Must be called
    before run_simulation() for a new or modified SRH-2D case.

    Only relevant for SRH-2D; HEC-RAS does not have a separate pre-processing
    step managed by pyHMT2D.

    Parameters
    ----------
    srh_pre_path : str
        Full path to srh_2d_pre.exe (e.g. "C:/SRH-2D/srh_2d_pre_3.3.exe").
        Required on the first call; remembered for subsequent calls in the
        same session.
    version : str
        SRH-2D version string, e.g. "3.3". Defaults to "3.3".
    extra_dll_path : str, optional
        Directory that contains HDF5/XMDL DLLs required by SRH-2D on Windows.

    Returns
    -------
    dict with "data" containing {"elapsed_seconds": float, "log": str}
    """
    session = get_session()
    try:
        session.require_open()
    except RuntimeError as e:
        return _err(str(e))

    if session.model_type != "SRH-2D":
        return _err(
            "run_preprocessing() is only applicable to SRH-2D projects. "
            "HEC-RAS pre-processing is handled internally by HEC-RAS."
        )

    # Resolve srh_pre_path
    pre_path = srh_pre_path or getattr(session, "_srh_pre_path", None)
    if not pre_path:
        return _err(
            "srh_pre_path is required. Provide the full path to srh_2d_pre.exe."
        )
    if not os.path.isfile(pre_path):
        return _err(f"SRH-2D pre-processor not found: {pre_path}")

    dll_path = extra_dll_path or getattr(session, "_extra_dll_path", "")
    srh_path = getattr(session, "_srh_path", pre_path)  # placeholder

    try:
        import pyHMT2D
        model = pyHMT2D.SRH_2D.SRH_2D_Model(
            version, pre_path, srh_path, dll_path, faceless=True
        )
        model.init_model()
        model.set_simulation_case(session.data)

        # Cache paths on session for run_simulation
        session._srh_pre_path = pre_path
        session._extra_dll_path = dll_path
        session._srh_version = version
        session.model = model

        # SRH-2D Pre must run from the case directory so that it receives only
        # the filename (not a long absolute path) and can resolve relative refs.
        case_dir = os.path.dirname(os.path.abspath(session.project_file))
        _orig_dir = os.getcwd()
        os.chdir(case_dir)
        t0 = time.time()
        try:
            model.run_pre_model()
        finally:
            os.chdir(_orig_dir)
        elapsed = round(time.time() - t0, 1)

        return _ok(
            {"elapsed_seconds": elapsed},
            f"SRH-2D pre-processing completed in {elapsed}s.",
        )

    except Exception as exc:
        return _err(f"Pre-processing failed: {exc}")


# ── tool 11 ───────────────────────────────────────────────────────────────────

def run_simulation(
    timeout_minutes: float = 60,
    srh_path: Optional[str] = None,
    srh_pre_path: Optional[str] = None,
    extra_dll_path: Optional[str] = None,
    srh_2d_version: str = "3.7",
    hecras_version: str = "6.6",
    faceless: bool = False,
) -> dict:
    """Run the hydraulic simulation for the currently open project.

    For SRH-2D: executes srh_2d.exe via subprocess. The pre-processor must
    have been run first (call run_preprocessing()).
    For HEC-RAS: launches HEC-RAS via COM automation (Windows only).

    Parameters
    ----------
    timeout_minutes : float
        Maximum time to wait for the simulation to finish. Default 60 min.
    srh_path : str
        (SRH-2D only) Full path to srh_2d.exe.
        Required on the first call; remembered for subsequent calls.
    srh_pre_path : str
        (SRH-2D only) Full path to srh_2d_pre.exe (needed to build the model
        object even if pre-processing already ran).
    extra_dll_path : str
        (SRH-2D only) Directory containing HDF5/XMDL DLLs.
    srh_2d_version : str
        SRH-2D version string. Default "3.7".
    hecras_version : str
        HEC-RAS version string used by the COM interface. Default "6.6".
    faceless : bool
        Run without showing the solver GUI window. Default False.

    Returns
    -------
    dict with "data" containing
        {"elapsed_seconds": float, "model_type": str, "status": str}
    """
    session = get_session()
    try:
        session.require_open()
    except RuntimeError as e:
        return _err(str(e))

    import pyHMT2D
    t0 = time.time()

    try:
        if session.model_type == "SRH-2D":
            # Resolve paths
            exe_path = srh_path or getattr(session, "_srh_path", None)
            pre_path = srh_pre_path or getattr(session, "_srh_pre_path", None)
            dll_path = extra_dll_path or getattr(session, "_extra_dll_path", "")
            ver = srh_2d_version or getattr(session, "_srh_version", "3.7")

            if not exe_path:
                return _err(
                    "srh_path is required. Provide the full path to srh_2d.exe."
                )
            if not os.path.isfile(exe_path):
                return _err(f"SRH-2D solver not found: {exe_path}")

            if pre_path is None:
                pre_path = exe_path  # use same dir as fallback

            model = session.model
            if model is None:
                model = pyHMT2D.SRH_2D.SRH_2D_Model(
                    ver, pre_path, exe_path, dll_path, faceless=faceless
                )
                model.init_model()
                model.set_simulation_case(session.data)
                session.model = model
            else:
                # Update exe path if supplied
                model._srh_path = exe_path

            session._srh_path = exe_path
            session._srh_pre_path = pre_path
            session._extra_dll_path = dll_path

            sleep_time = min(10.0, timeout_minutes * 60 / 10)
            model.run_model(sleepTime=sleep_time, bShowProgress=False)

        else:  # HEC-RAS
            model = session.model
            cached_faceless = getattr(session, "_hecras_faceless", None)
            if model is None or cached_faceless != faceless:
                # Close the existing instance before opening a new one with
                # a different faceless setting (avoids orphaned HEC-RAS windows)
                if model is not None:
                    try:
                        model.exit_model()
                    except Exception:
                        pass
                model = pyHMT2D.RAS_2D.HEC_RAS_Model(
                    version=hecras_version, faceless=faceless
                )
                model.init_model()
                prj_file = session._hecras_prj_file or session.project_file
                model.open_project(prj_file)
                session.model = model
                session._hecras_faceless = faceless

            model.run_model()

        elapsed = round(time.time() - t0, 1)
        return _ok(
            {"elapsed_seconds": elapsed, "model_type": session.model_type, "status": "completed"},
            f"{session.model_type} simulation completed in {elapsed}s.",
        )

    except Exception as exc:
        elapsed = round(time.time() - t0, 1)
        return _err(f"Simulation failed after {elapsed}s: {exc}")


# ── tool 12 ───────────────────────────────────────────────────────────────────

def exit_model() -> dict:
    """Quit the HEC-RAS GUI and release the COM object.

    Only applicable to HEC-RAS projects. For SRH-2D this is a no-op because
    SRH-2D runs as a plain subprocess that exits on its own when the simulation
    finishes.

    Call this after a simulation with ``faceless=False`` to close the HEC-RAS
    window, or any time you want to release the HEC-RAS COM handle without
    clearing the rest of the session.

    Returns
    -------
    dict with "data" containing {"model_type": str, "exited": bool}
    """
    session = get_session()

    if session.model_type == "SRH-2D":
        return _ok(
            {"model_type": "SRH-2D", "exited": False},
            "SRH-2D runs as a subprocess and exits automatically — nothing to do.",
        )

    if session.model is None:
        return _ok(
            {"model_type": session.model_type or "unknown", "exited": False},
            "No active model instance found (already exited or never started).",
        )

    try:
        session.model.exit_model()
        session.model = None
        session._hecras_faceless = None
        return _ok(
            {"model_type": "HEC-RAS", "exited": True},
            "HEC-RAS has been shut down and the COM handle released.",
        )
    except Exception as exc:
        return _err(f"Failed to exit HEC-RAS: {exc}")


# ── tool 13 ───────────────────────────────────────────────────────────────────

def close_session() -> dict:
    """Exit the solver (if HEC-RAS) and reset the session to a blank state.

    Use this to start fresh — e.g. before opening a different project or
    switching between SRH-2D and HEC-RAS cases.  For HEC-RAS this also calls
    QuitRas() to close the GUI window and release the COM object before
    clearing state.  For SRH-2D only the in-memory state is cleared (no
    process to kill).

    Returns
    -------
    dict with "data" containing {"model_type": str, "session_cleared": bool}
    """
    session = get_session()
    model_type = session.model_type or "none"

    # Gracefully shut down HEC-RAS if running
    if session.model_type == "HEC-RAS" and session.model is not None:
        try:
            session.model.exit_model()
        except Exception:
            pass  # best-effort — clear the session regardless

    from pyHMT2D.AI_Tools.state import reset_session
    reset_session()

    return _ok(
        {"model_type": model_type, "session_cleared": True},
        f"Session closed (was: {model_type}). Ready to open a new project.",
    )


def kill_all_hec_ras() -> dict:
    """Force-kill all running HEC-RAS processes on this machine.

    Uses TASKLIST to find all ``ras.exe`` processes and TASKKILL /F to
    terminate each one.  No confirmation is asked — calling this tool is
    itself the confirmation.  Use with caution: any unsaved work in open
    HEC-RAS windows will be lost.

    Returns
    -------
    dict with "data" containing {"n_killed": int, "pids": [int, ...]}
    """
    import os
    import subprocess

    try:
        proc = subprocess.Popen(
            'TASKLIST /FO "CSV"', stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        tasklist = proc.stdout.read().split(b'\n')
    except Exception as exc:
        return _err(f"Failed to enumerate processes: {exc}")

    pids = []
    for line in tasklist:
        if b'ras.exe' in line.lower():
            parts = line.split(b',')
            try:
                pids.append(int(eval(parts[1])))
            except Exception:
                pass

    if not pids:
        return _ok({"n_killed": 0, "pids": []}, "No running HEC-RAS instances found.")

    errors = []
    killed = []
    for pid in pids:
        rc = os.system(f'TASKKILL /PID {pid} /F >nul 2>&1')
        if rc == 0:
            killed.append(pid)
        else:
            errors.append(pid)

    msg = f"Killed {len(killed)} HEC-RAS instance(s) (PIDs: {killed})."
    if errors:
        msg += f" Failed to kill PIDs: {errors}."

    return _ok({"n_killed": len(killed), "pids": killed}, msg)


# ── tool 14 ───────────────────────────────────────────────────────────────────

def get_simulation_status() -> dict:
    """Return the current status of the model and session.

    Reports whether a project is open, whether results have been loaded, and
    the last known state of the solver process (if run_simulation was called).

    Returns
    -------
    dict with "data" containing
        {"project_open": bool, "model_type": str, "results_loaded": bool,
         "project_file": str, "result_file": str}
    """
    session = get_session()

    status_data = {
        "project_open": session.is_open(),
        "model_type": session.model_type,
        "project_file": session.project_file,
        "results_loaded": session.result_loaded,
        "result_file": session.result_file,
        "result_variables": session.result_variables,
        "n_result_time_steps": len(session.result_time_steps),
        "vtk_file": session.vtk_file,
    }

    if not session.is_open():
        msg = "No project open. Call get_project_info() to open a project."
    elif not session.result_loaded:
        msg = (
            f"{session.model_type} project open: {session.project_file}. "
            "Results not yet loaded. Run the model or call read_results()."
        )
    else:
        msg = (
            f"{session.model_type} project open with results loaded. "
            f"{len(session.result_variables)} variable(s), "
            f"{len(session.result_time_steps)} time step(s)."
        )

    return _ok(status_data, msg)
