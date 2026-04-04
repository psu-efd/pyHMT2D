## Demonstrate the use of `HEC_RAS_Model` and `HEC_RAS_Project` classes

This directory contains example scripts that demonstrate how to use **pyHMT2D**'s
`HEC_RAS_Model` and `HEC_RAS_Project` classes to control HEC-RAS runs and post-process results.

---

## `Muncie2D/` — single plan, single 2D area

A compact example using the Muncie 2D case (one plan, one 2D flow area). The case files are
included in the repository under `examples/assets/RAS_models/Muncie2D/` and are copied into
this directory automatically when the script runs.

```bash
cd Muncie2D
python demo_HEC_RAS_Model.py
```

The script performs two independent tasks:

1. **Run HEC-RAS** via the `HEC_RAS_Model` class (requires HEC-RAS to be installed):
   ```python
   my_hec_ras_model = pyHMT2D.RAS_2D.HEC_RAS_Model(version="6.6", faceless=False)
   my_hec_ras_model.init_model()
   my_hec_ras_model.open_project("Muncie2D.prj")
   my_hec_ras_model.run_model()
   my_hec_ras_model.close_project()
   my_hec_ras_model.exit_model()
   ```

2. **Convert results to VTK** using the project-centric API (no HEC-RAS installation needed):
   ```python
   from pyHMT2D.Hydraulic_Models_Data.RAS_2D.HEC_RAS_Model import HEC_RAS_Project
   project = HEC_RAS_Project("Muncie2D.prj")
   plan = project.get_plan("p01")
   my_ras_2d_data = plan.load_results()
   my_ras_2d_data.saveHEC_RAS2D_results_to_VTK(lastTimeStep=True)
   ```

`HEC_RAS_Project` parses the `.prj` file and auto-discovers plans, geometries, and terrain — no
need to specify HDF or terrain file paths manually. Terrain is resolved automatically from the
geometry HDF and the `.rasmap` file.

---

## `BaldEagleCrkMulti2D/` — multiple plans, multiple 2D areas

A more advanced example using the **Bald Eagle Creek Dam Break** case, which ships with
HEC-RAS 6.x as part of the *2D Unsteady Flow Hydraulics* example set. This case contains:

- **11 plans** covering different modelling configurations (1D/2D linked, pure 2D, gridded
  precipitation, dam break, etc.)
- **Multiple 2D flow areas** within a single plan — plan `p13` ("PMF with Multi 2D Areas")
  has three areas: `193`, `194`, and `LockHaven`

The demo script `BaldEagleCrkMulti2D/demo_HEC_RAS_Model.py` shows how pyHMT2D handles these
cases, including:

- Listing all plans and geometries from the project file
- Selecting and running a specific plan by ID or name
- Inspecting the 2D area names and cell counts after loading results
- Exporting each 2D area to a **separate VTK file** (default, `combined_areas=False`)
- Exporting all 2D areas **merged into one VTK file** (`combined_areas=True`)
- Batch-processing all plans that have result HDF files

```bash
cd BaldEagleCrkMulti2D
python demo_HEC_RAS_Model.py
```

### Case files not included

The Bald Eagle Creek case is approximately **338 MB** and is therefore too large to include
in the GitHub repository. The case files are excluded via `.gitignore`; only the demo script
is tracked.

To use this example, copy the case files from the HEC-RAS 6.x installation directory:

```
<HEC-RAS install>\Example Projects\2D Unsteady Flow Hydraulics\BaldEagleCrkMulti2D\
```

into this directory (`examples/HEC_RAS_Model/BaldEagleCrkMulti2D/`) before running the script.
The HEC-RAS installer places the example projects under:

```
C:\Users\<username>\Documents\HEC\HEC-RAS\6.6\Example Projects\
```

or under the HEC-RAS program directory depending on the installation type. Once the files are
in place, run the script as shown above.
