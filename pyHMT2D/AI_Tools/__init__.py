# -*- coding: utf-8 -*-
"""
pyHMT2D AI Tools
================

A collection of 25 AI-agent-callable tools for automating 2D hydraulic
modelling workflows with SRH-2D and HEC-RAS.

Quick start
-----------
Python API::

    from pyHMT2D.AI_Tools.tools import get_project_info, get_materials, set_manning_n

MCP server (stdio, for Claude Desktop / Claude Code)::

    python -m pyHMT2D.AI_Tools.mcp_server

Tool groups
-----------
1. Discovery      : list_model_files, get_project_info, get_materials,
                    get_boundary_conditions, get_result_variables
2. Parameters     : set_manning_n, set_inlet_flow, set_exit_wse,
                    save_modified_inputs
3. Execution      : run_preprocessing, run_simulation, get_simulation_status
4. Results        : read_results, get_value_at_point, get_result_statistics,
                    get_flood_extent, get_cross_section_profile
5. Export         : export_to_vtk, export_mesh_to_vtk
6. Calibration    : load_observations, evaluate_parameters, run_calibration
7. Monte Carlo    : generate_mc_samples, run_monte_carlo, get_mc_statistics
"""

from pyHMT2D.AI_Tools.state import get_session, reset_session
from pyHMT2D.AI_Tools.tools import *  # noqa: F401, F403

__all__ = [
    "get_session",
    "reset_session",
]
