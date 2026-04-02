# -*- coding: utf-8 -*-
"""
Session state for pyHMT2D AI Tools.

All tools share a single HydraulicSession singleton so the model type and
open project are detected once (at get_project_info / open time) and reused
by every subsequent tool call without requiring the caller to re-specify them.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, Any


@dataclass
class HydraulicSession:
    """Holds the currently open hydraulic model and its loaded data.

    Only one project is active at a time. Calling get_project_info() or
    open_project() on a new file replaces the current session.
    """

    # ── identity ──────────────────────────────────────────────────────────────
    model_type: Optional[str] = None          # "SRH-2D" | "HEC-RAS"
    project_file: Optional[str] = None        # path to the control/project file (.hdf for HEC-RAS)
    _hecras_prj_file: Optional[str] = None    # path to .prj (HEC-RAS only; needed by open_project())

    # ── live objects (set by project_tools.open_project) ─────────────────────
    model: Optional[Any] = None               # SRH_2D_Model | HEC_RAS_Model
    data: Optional[Any] = None                # SRH_2D_Data  | RAS_2D_Data

    # ── result state (set by result_tools.read_results) ──────────────────────
    result_file: Optional[str] = None         # path to loaded result file
    result_loaded: bool = False
    result_time_steps: list = field(default_factory=list)
    result_variables: list = field(default_factory=list)

    # ── vtk reader (cached after export / read) ───────────────────────────────
    vtk_file: Optional[str] = None            # last exported / used VTK file

    def clear(self) -> None:
        """Reset the session to an empty state."""
        self.model_type = None
        self.project_file = None
        self._hecras_prj_file = None
        self.model = None
        self.data = None
        self.result_file = None
        self.result_loaded = False
        self.result_time_steps = []
        self.result_variables = []
        self.vtk_file = None

    def is_open(self) -> bool:
        """Return True if a project has been opened."""
        return self.model_type is not None and self.data is not None

    def require_open(self) -> None:
        """Raise RuntimeError if no project is open."""
        if not self.is_open():
            raise RuntimeError(
                "No project is currently open. "
                "Call get_project_info() with a model file first."
            )

    def require_results(self) -> None:
        """Raise RuntimeError if results have not been loaded."""
        self.require_open()
        if not self.result_loaded:
            raise RuntimeError(
                "No results loaded. Call read_results() first."
            )


# Module-level singleton — all tools import and mutate this object.
_session = HydraulicSession()


def get_session() -> HydraulicSession:
    """Return the global session singleton."""
    return _session


def reset_session() -> None:
    """Clear the global session (useful in tests)."""
    _session.clear()
