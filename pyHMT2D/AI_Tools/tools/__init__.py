# -*- coding: utf-8 -*-
"""
pyHMT2D AI Tools — tool modules.

Import individual functions or use the flat namespace below.
"""

from .project_tools import (
    list_model_files,
    get_project_info,
    get_materials,
    get_boundary_conditions,
    get_result_variables,
    build_param_specs,
    load_config,
)
from .parameter_tools import (
    set_manning_n,
    set_inlet_flow,
    set_exit_wse,
    save_modified_inputs,
)
from .execution_tools import (
    run_preprocessing,
    run_simulation,
    exit_model,
    close_session,
    kill_all_hec_ras,
    get_simulation_status,
)
from .result_tools import (
    read_results,
    get_value_at_point,
    get_result_statistics,
    get_flood_extent,
    get_cross_section_profile,
)
from .export_tools import (
    export_to_vtk,
    export_mesh_to_vtk,
)
from .calibration_tools import (
    check_observation_format,
    load_observations,
    evaluate_parameters,
    run_calibration,
)
from .monte_carlo_tools import (
    generate_mc_samples,
    run_monte_carlo,
    get_mc_statistics,
)
from .conversion_tools import (
    ras_to_srh,
    srh_to_vtk,
)

__all__ = [
    # Discovery
    "list_model_files",
    "get_project_info",
    "get_materials",
    "get_boundary_conditions",
    "get_result_variables",
    "build_param_specs",
    "load_config",
    # Parameters
    "set_manning_n",
    "set_inlet_flow",
    "set_exit_wse",
    "save_modified_inputs",
    # Execution
    "run_preprocessing",
    "run_simulation",
    "exit_model",
    "close_session",
    "kill_all_hec_ras",
    "get_simulation_status",
    # Results
    "read_results",
    "get_value_at_point",
    "get_result_statistics",
    "get_flood_extent",
    "get_cross_section_profile",
    # Export
    "export_to_vtk",
    "export_mesh_to_vtk",
    # Calibration
    "check_observation_format",
    "load_observations",
    "evaluate_parameters",
    "run_calibration",
    # Monte Carlo
    "generate_mc_samples",
    "run_monte_carlo",
    "get_mc_statistics",
    # Conversion
    "ras_to_srh",
    "srh_to_vtk",
]
