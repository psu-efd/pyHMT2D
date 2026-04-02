# pyHMT2D AI Tools

AI-agent-callable tools for automating 2D hydraulic modelling workflows with
SRH-2D (USBR) and HEC-RAS (USACE), exposed as a
[Model Context Protocol (MCP)](https://modelcontextprotocol.io/) server.

## Installation

```bash
pip install pyHMT2D[ai]      # includes the MCP dependency
```

Or install the MCP package separately:

```bash
pip install mcp
```

## Quick start — Python API

```python
from pyHMT2D.AI_Tools.tools import (
    list_model_files, get_project_info,
    get_materials, set_manning_n,
    run_simulation, read_results, get_flood_extent,
)

# 1. Discover files
list_model_files("C:/project/Muncie")

# 2. Open project — auto-detects SRH-2D or HEC-RAS
get_project_info("C:/project/Muncie/Muncie.srhhydro")

# 3. Inspect and modify
get_materials()
set_manning_n(material_ids=[2, 4], values=[0.035, 0.055])

# 4. Run and query
run_simulation(srh_path="C:/SRH-2D/srh_2d_3.3.exe", ...)
read_results("C:/project/Muncie/Muncie_XMDFC.h5")
get_flood_extent(depth_threshold=0.1)
```

All tools return a consistent dict:
```python
{"status": "ok" | "error", "data": <payload>, "message": "<human summary>"}
```

## MCP server — connect an AI agent

### Run the server

```bash
python -m pyHMT2D.AI_Tools.mcp_server
# or, after pip install pyHMT2D[ai]:
hmt-mcp-server
```

### Configure Claude Desktop

Add to `~/.config/claude/claude_desktop_config.json`
(macOS/Linux) or `%APPDATA%\Claude\claude_desktop_config.json` (Windows):

```json
{
  "mcpServers": {
    "pyHMT2D": {
      "command": "python",
      "args": ["-m", "pyHMT2D.AI_Tools.mcp_server"]
    }
  }
}
```

Then restart Claude Desktop. The AI agent will have access to all 25 tools.

### Configure Claude Code (this CLI)

```bash
claude mcp add pyHMT2D -- python -m pyHMT2D.AI_Tools.mcp_server
```

## OpenAI / Anthropic function-calling schemas

Generate portable JSON schemas for use with any LLM API:

```bash
python -m pyHMT2D.AI_Tools.schemas.generate_schemas
# writes pyHMT2D/AI_Tools/schemas/openai_tool_schemas.json
```

Use in your own agent code:

```python
import json
from openai import OpenAI

tools = json.load(open("openai_tool_schemas.json"))
client = OpenAI()
response = client.chat.completions.create(
    model="gpt-4o",
    messages=[{"role": "user", "content": "What is the flood extent?"}],
    tools=tools,
)
```

## Tool reference

### Group 1 — Discovery (no solver needed)

| Tool | Description |
|------|-------------|
| `list_model_files(directory)` | Find `.srhhydro`, `_SIF.dat`, `.prj`, `.hdf` files |
| `get_project_info(project_file)` | Open project, auto-detect model type, return metadata |
| `get_materials()` | List material zones with Manning's n values |
| `get_boundary_conditions()` | List BCs with IDs, types, and current values |
| `get_result_variables(result_file)` | List variable names and time steps in a result file |

### Group 2 — Parameters

| Tool | Description |
|------|-------------|
| `set_manning_n(material_ids, values)` | Modify Manning's n (SRH-2D and HEC-RAS) |
| `set_inlet_flow(bc_ids, values)` | Modify inlet discharge, m³/s (SRH-2D only) |
| `set_exit_wse(bc_ids, values)` | Modify exit water surface elevation (SRH-2D only) |
| `save_modified_inputs(output_path)` | Write changes back to the input files |

### Group 3 — Execution

| Tool | Description |
|------|-------------|
| `run_preprocessing(srh_pre_path, ...)` | Run SRH-2D pre-processor |
| `run_simulation(srh_path, ...)` | Run SRH-2D or HEC-RAS |
| `get_simulation_status()` | Report session state and last run status |

### Group 4 — Results

| Tool | Description |
|------|-------------|
| `read_results(result_file)` | Load XMDF / HDF results into session |
| `get_value_at_point(x, y, variable)` | Interpolate a variable at (x, y) |
| `get_result_statistics(variable)` | Domain-wide min/max/mean/std |
| `get_flood_extent(depth_threshold)` | Flooded area above threshold |
| `get_cross_section_profile(x1,y1,x2,y2, variable)` | Profile along a transect |

### Group 5 — Export

| Tool | Description |
|------|-------------|
| `export_to_vtk(timestep, output_dir)` | Export results to VTK (ParaView) |
| `export_mesh_to_vtk(output_path)` | Export mesh geometry to VTK |

### Group 6 — Calibration

| Tool | Description |
|------|-------------|
| `load_observations(csv_file)` | Load high water marks / field measurements |
| `evaluate_parameters(param_specs, observation_csv)` | Single run → RMSE |
| `run_calibration(param_specs, observation_csv, n_iterations, method)` | Full automated calibration (GP or Nelder-Mead) |

`param_specs` structure:
```python
[
  {"type": "manning_n", "material_id": 2, "material_name": "channel",
   "initial": 0.04, "min": 0.02, "max": 0.08},
]
```

### Group 7 — Monte Carlo

| Tool | Description |
|------|-------------|
| `generate_mc_samples(param_specs, n_samples)` | Sample from truncated-normal or uniform distributions |
| `run_monte_carlo(base_case_dir, param_specs, n_samples, n_processes)` | Run N-sample ensemble |
| `get_mc_statistics(results_json, variable, observation_points)` | Exceedance probabilities at points and spatially |

`param_specs` for MC (extends calibration schema with distribution info):
```python
[
  {"type": "manning_n", "material_id": 2,
   "distribution": "truncated_normal",
   "mean": 0.04, "std": 0.005, "min": 0.03, "max": 0.05},
]
```

## Model type detection

The AI agent **never needs to specify the model type**. It is auto-detected:

| File passed to `get_project_info` | Detected model |
|---|---|
| `*.srhhydro`, `*_SIF.dat` | SRH-2D |
| `*.hdf`, `*.hdf5` | HEC-RAS |
| `*.prj` (with HEC-RAS content) | HEC-RAS |

Once a project is opened, all subsequent tools route to the correct code path
via the session singleton (`pyHMT2D.AI_Tools.state.get_session()`).

## Examples

See `examples/AI_Tools/`:
```
examples/AI_Tools/
├── SRH_2D/
│   ├── 01_open_inspect/demo.py
│   ├── 02_modify_and_run/demo.py
│   ├── 03_query_results/demo.py
│   ├── 04_calibration/demo.py
│   └── 05_monte_carlo/demo.py
└── RAS_2D/
    ├── 01_open_inspect/demo.py
    ├── 02_modify_and_run/demo.py
    ├── 03_query_results/demo.py
    ├── 04_calibration/demo.py
    └── 05_monte_carlo/demo.py
```

## Running tests

```bash
pytest tests/AI_Tools/ -v
```

No solver installation is required for the unit tests — they test input
validation, session state, error handling, and sample generation only.
