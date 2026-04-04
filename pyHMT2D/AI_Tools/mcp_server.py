# -*- coding: utf-8 -*-
"""
pyHMT2D MCP Server

Exposes 32 AI tools as Model Context Protocol (MCP) tools so any
MCP-compatible AI agent (Claude Desktop, Claude Code, custom agents) can
call them via stdio.

Usage
-----
Run directly:
    python -m pyHMT2D.AI_Tools.mcp_server

Or via the installed entry point (after `pip install pyHMT2D[ai]`):
    hmt-mcp-server

Configure in Claude Desktop (~/.config/claude/claude_desktop_config.json):
    {
      "mcpServers": {
        "pyHMT2D": {
          "command": "python",
          "args": ["-m", "pyHMT2D.AI_Tools.mcp_server"]
        }
      }
    }
"""

from __future__ import annotations

from mcp.server.fastmcp import FastMCP

from pyHMT2D.AI_Tools.tools import (
    # Discovery / Config
    load_config,
    list_model_files,
    get_project_info,
    get_materials,
    get_boundary_conditions,
    get_result_variables,
    build_param_specs,
    # Parameters
    set_manning_n,
    set_inlet_flow,
    set_exit_wse,
    save_modified_inputs,
    # Execution
    run_preprocessing,
    run_simulation,
    exit_model,
    close_session,
    get_simulation_status,
    # Results
    read_results,
    get_value_at_point,
    get_result_statistics,
    get_flood_extent,
    get_cross_section_profile,
    # Export
    export_to_vtk,
    export_mesh_to_vtk,
    # Calibration
    check_observation_format,
    load_observations,
    evaluate_parameters,
    run_calibration,
    # Monte Carlo
    generate_mc_samples,
    run_monte_carlo,
    get_mc_statistics,
    # Conversion
    ras_to_srh,
    srh_to_vtk,
)

mcp = FastMCP(
    name="pyHMT2D",
    instructions=(
        "Tools for controlling and automating 2D hydraulic simulations "
        "using SRH-2D (USBR) and HEC-RAS (USACE). Supports model parameter "
        "inspection and modification, simulation execution, result queries, "
        "calibration, and Monte Carlo uncertainty analysis."
    ),
)

# ── Register all tools ────────────────────────────────────────────────────────
# FastMCP reads each function's type annotations and docstring to auto-generate
# the JSON schema exposed to the AI agent.

# Group 1: Discovery / Config
mcp.tool()(load_config)
mcp.tool()(list_model_files)
mcp.tool()(get_project_info)
mcp.tool()(get_materials)
mcp.tool()(get_boundary_conditions)
mcp.tool()(get_result_variables)
mcp.tool()(build_param_specs)

# Group 2: Parameters
mcp.tool()(set_manning_n)
mcp.tool()(set_inlet_flow)
mcp.tool()(set_exit_wse)
mcp.tool()(save_modified_inputs)

# Group 3: Execution
mcp.tool()(run_preprocessing)
mcp.tool()(run_simulation)
mcp.tool()(exit_model)
mcp.tool()(close_session)
mcp.tool()(get_simulation_status)

# Group 4: Results
mcp.tool()(read_results)
mcp.tool()(get_value_at_point)
mcp.tool()(get_result_statistics)
mcp.tool()(get_flood_extent)
mcp.tool()(get_cross_section_profile)

# Group 5: Export
mcp.tool()(export_to_vtk)
mcp.tool()(export_mesh_to_vtk)

# Group 6: Calibration
mcp.tool()(check_observation_format)
mcp.tool()(load_observations)
mcp.tool()(evaluate_parameters)
mcp.tool()(run_calibration)

# Group 7: Monte Carlo
mcp.tool()(generate_mc_samples)
mcp.tool()(run_monte_carlo)
mcp.tool()(get_mc_statistics)

# Group 8: Conversion
mcp.tool()(ras_to_srh)
mcp.tool()(srh_to_vtk)


def _preload_heavy_modules() -> None:
    """Import vtk and other heavy submodules before mcp.run() takes over stdout.

    The MCP server communicates via stdio (stdout = JSON messages).  Any library
    that prints to stdout during import (e.g. vtk banner) would corrupt the
    protocol.  By importing here — with stdout redirected to stderr — we ensure
    all startup noise is safely discarded before the protocol begins.
    """
    import sys
    _stdout = sys.stdout
    sys.stdout = sys.stderr
    try:
        import pyHMT2D.Hydraulic_Models_Data.SRH_2D  # noqa: F401
        import pyHMT2D.Hydraulic_Models_Data.RAS_2D  # noqa: F401
        import pyHMT2D.Misc                          # noqa: F401
    except Exception:
        pass  # missing optional deps (e.g. pywin32 on non-Windows) are OK
    finally:
        sys.stdout = _stdout


def _cleanup_session() -> None:
    """Best-effort cleanup: close the solver and session on process exit."""
    try:
        from pyHMT2D.AI_Tools.state import get_session
        session = get_session()
        if session.model is not None:
            try:
                session.model.close_project()
            except Exception:
                pass
            try:
                session.model.exit_model()
            except Exception:
                pass
        session.clear()
    except Exception:
        pass


def main() -> None:
    import atexit
    atexit.register(_cleanup_session)
    _preload_heavy_modules()
    mcp.run()


if __name__ == "__main__":
    main()
