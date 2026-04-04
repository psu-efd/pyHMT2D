# AGENTS.md

This file provides guidance to Codex (Codex.ai/code) when working with code in this repository.

## Commands

```bash
# Install in editable mode (development)
pip install -e .

# Run all tests
pytest

# Run a specific test file
pytest tests/01_SRH_2D_Data/01_01_SRH_2D_Data_test.py

# Run tests with verbose output (already enabled via pytest.ini addopts = -s)
pytest -v
```

**CLI tools** (available after install):
```bash
# Convert HEC-RAS 2D mesh to SRH-2D format
hmt-cli ras_to_srh --args '{"ras_hdf_file": "input.hdf", "terrain_tif_file": "terrain.tif", "srh_case_name": "output"}'

# Convert SRH-2D results to VTK (last time step)
hmt-cli srh_to_vtk --args '{"srhhydro_file": "case.srhhydro", "output_file": "result.h5"}'

# Export SRH-2D mesh only to VTK
hmt-cli srh_to_vtk --args '{"srhhydro_file": "case.srhhydro", "output_file": "mesh.vtk", "mesh_only": true}'
```

## Architecture

pyHMT2D automates and post-processes 2D hydraulic simulations for two external solvers: **SRH-2D** (USBR) and **HEC-RAS 2D** (USACE). It is **Windows-only** due to solver dependencies (SRH-2D via subprocess, HEC-RAS via COM interface).

### Package Layout

```
pyHMT2D/
├── Hydraulic_Models_Data/         # Model-specific data and control
│   ├── Hydraulic_Models_Data_Base/  # Abstract base classes
│   │   ├── HydraulicModel.py        # Base for all model controllers
│   │   └── HydraulicData.py         # Base for all data readers
│   ├── SRH_2D/                      # SRH-2D implementation
│   ├── RAS_2D/                      # HEC-RAS 2D implementation
│   └── Backwater_1D/                # Simple 1D demo model
├── Calibration/                   # Model calibration (scipy-based)
├── Parametric_Study/              # Batch/Monte Carlo parametric studies
├── Misc/                          # Utilities shared across models
│   ├── Terrain.py                   # DEM/raster processing (rasterio)
│   ├── RAS_to_SRH_Converter.py      # Mesh/BC format conversion
│   ├── vtk_utilities.py             # VTK export and manipulation
│   └── tools.py                     # Result sampling and probing
└── __common__.py                  # Global config (verbose flag, etc.)
```

### Key Design Patterns

**Model + Data split**: Each solver has a `*_Model` class (controls solver execution) and a `*_Data` class (reads/writes files). They are used independently or together.

```python
# SRH-2D: runs via subprocess
model = pyHMT2D.SRH_2D.SRH_2D_Model(version, pre_path, exe_path, dll_path)
data  = pyHMT2D.SRH_2D.SRH_2D_Data("case.srhhydro")

# HEC-RAS: runs via Windows COM automation (pywin32)
model = pyHMT2D.RAS_2D.HEC_RAS_Model(version="6.6")
data  = pyHMT2D.RAS_2D.RAS_2D_Data("case.prj")
```

**Data flow**: Solver output files → `*_Data` classes (parse HDF5/binary) → NumPy arrays → optionally exported to VTK for visualization.

**Calibration**: `Calibration/` wraps scipy optimizers around the Model+Data cycle to tune Manning's n or other parameters by minimizing error against observed data.

**Parametric studies**: `Parametric_Study/` iterates over parameter combinations (single-param sweeps, multi-param grids, Monte Carlo) using the same Model+Data cycle.

**Format conversion**: `Misc/RAS_to_SRH_Converter.py` translates HEC-RAS mesh (HDF5), terrain, and boundary conditions into SRH-2D input format. `Misc/vtk_utilities.py` handles VTK output for both models.

### AI Tools (`pyHMT2D/AI_Tools/`)

27 AI-agent-callable tools exposed as an MCP server. Self-contained directory.

```
AI_Tools/
├── mcp_server.py          # FastMCP server entry point (hmt-mcp-server)
├── state.py               # HydraulicSession singleton (model type auto-detected)
├── tools/                 # One module per tool group
│   ├── project_tools.py   # Discovery: list files, open project, get materials/BCs
│   ├── parameter_tools.py # set_manning_n, set_inlet_flow, set_exit_wse, save
│   ├── execution_tools.py # run_preprocessing, run_simulation, exit_model, close_session, get_status
│   ├── result_tools.py    # read_results, get_value_at_point, flood_extent, profile
│   ├── export_tools.py    # export_to_vtk, export_mesh_to_vtk
│   ├── calibration_tools.py  # load_observations, evaluate_parameters, run_calibration
│   ├── monte_carlo_tools.py  # generate_mc_samples, run_monte_carlo, get_mc_statistics
│   └── conversion_tools.py   # ras_to_srh, srh_to_vtk (stateless format conversion)
└── schemas/
    └── generate_schemas.py   # generates openai_tool_schemas.json
```

**Key design**: model type (SRH-2D vs HEC-RAS) is auto-detected from the file
passed to `get_project_info()` and stored in the session — no tool after that
needs `model_type` specified. All tools return `{"status", "data", "message"}`.

**MCP setup** (Codex Desktop / Codex):
```bash
pip install pyHMT2D[ai]
Codex mcp add pyHMT2D -- python -m pyHMT2D.AI_Tools.mcp_server
```

### Tests

Tests are under `tests/` organized by module. `tests/AI_Tools/` contains 42
unit tests that run without any solver installation. Other tests require actual
solver installations (SRH-2D or HEC-RAS) and real input files.
