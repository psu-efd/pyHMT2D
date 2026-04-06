"""
Demo: HEC-RAS project with multiple plans and multiple 2D flow areas.

Case: Bald Eagle Creek Example Dam Break Study (BaldEagleDamBrk.prj)

This script demonstrates pyHMT2D's handling of HEC-RAS projects that contain:
  - Multiple plans (11 plans with different geometry/flow configurations)
  - Multiple 2D flow areas within a single plan (plan p13 has three areas:
    '193', '194', and 'LockHaven')

Run from this directory:
    cd examples/HEC_RAS_Model/BaldEagleCrkMulti2D
    python demo_HEC_RAS_Model.py
"""

import pyHMT2D
from pyHMT2D.Hydraulic_Models_Data.RAS_2D.HEC_RAS_Model import HEC_RAS_Project

PRJ_FILE = "BaldEagleDamBrk.prj"

# ── 1. Inspect the project (no COM / no HEC-RAS installation required) ────────

def inspect_project():
    """Open the project file and list all plans and geometries."""

    print("\n=== Project inspection ===")
    project = HEC_RAS_Project(PRJ_FILE)

    print(f"Project title : {project.title}")
    print(f"Current plan  : {project.current_plan_id}")
    print(f"\nAll plans ({len(project.plans)}):")
    for pid, plan in project.plans.items():
        marker = " <-- current" if pid == project.current_plan_id else ""
        hdf_status = "result HDF exists" if plan.hdf_exists() else "not yet run"
        print(f"  {pid}: \"{plan.plan_name}\"  geom={plan.geometry_id}  [{hdf_status}]{marker}")

    print(f"\nAll geometries ({len(project.geometries)}):")
    for gid, geom in project.geometries.items():
        terrain = geom.terrain_tiff_file or "(not found)"
        print(f"  {gid}: {geom.geom_file}  terrain={terrain}")

    return project


# ── 2. Run a selected plan via HEC-RAS COM ────────────────────────────────────

def run_plan(plan_id: str):
    """Run the specified plan using the HEC-RAS COM interface.

    Requires HEC-RAS 6.x to be installed. The plan is set as the current plan
    before running.
    """
    print(f"\n=== Running plan {plan_id} ===")

    model = pyHMT2D.RAS_2D.HEC_RAS_Model(version="6.6", faceless=False)
    model.init_model()

    project = model.open_project(PRJ_FILE)
    print(f"Opened project: {project.title}")

    # Switch to the requested plan (accepts plan ID like 'p13' or plan title)
    model.set_current_plan(plan_id)
    current = project.current_plan
    print(f"Active plan   : {current.plan_id} — \"{current.plan_name}\"")
    print(f"Geometry      : {current.geometry_id}")

    success = model.run_model()
    print(f"Simulation {'completed successfully' if success else 'failed'}.")

    model.close_project()
    model.exit_model()
    return success


# ── 3. Load results and explore multiple 2D areas ────────────────────────────

def load_and_inspect_results(plan_id: str):
    """Load plan results and report on each 2D flow area.

    Demonstrates how pyHMT2D handles plans with multiple 2D flow areas.
    Plan p13 ('PMF with Multi 2D Areas') contains three areas:
        '193', '194', and 'LockHaven'.
    """
    print(f"\n=== Loading results for plan {plan_id} ===")

    project = HEC_RAS_Project(PRJ_FILE)
    plan = project.get_plan(plan_id)

    if not plan.hdf_exists():
        print(f"  Result HDF not found for plan {plan_id}.")
        print(f"  Run the simulation first (call run_plan('{plan_id}')).")
        return None

    ras_data = plan.load_results()

    print(f"Plan          : {plan.plan_name}")
    print(f"Units         : {ras_data.units}")
    print(f"Simulation    : {ras_data.start_time}  ->  {ras_data.end_time}")
    print(f"\n2D flow areas ({len(ras_data.TwoDAreaNames)}):")
    for i, name in enumerate(ras_data.TwoDAreaNames):
        area_name = name.decode("utf-8") if isinstance(name, bytes) else name
        cell_count = ras_data.TwoDAreaCellCounts[i]
        print(f"  [{i}] \"{area_name}\"  —  {cell_count} cells")

    return ras_data


# ── 4. Export results to VTK — separate files per area (default) ──────────────

def export_vtk_separate(ras_data, plan_id: str):
    """Export each 2D area to its own VTK file (default behaviour).

    For a plan with N areas and T time steps this produces N × T files,
    named  RAS2D_<AreaName>_<timestep>.vtk.
    The last time step is exported here for brevity.
    """
    print(f"\n=== VTK export — separate file per area (plan {plan_id}) ===")

    vtk_files = ras_data.saveHEC_RAS2D_results_to_VTK(
        lastTimeStep=True,
        combined_areas=False,   # default: one VTK per 2D area
    )

    print(f"Wrote {len(vtk_files)} VTK file(s):")
    for f in vtk_files:
        print(f"  {f}")

    return vtk_files


# ── 5. Export results to VTK — all areas merged into one file ─────────────────

def export_vtk_combined(ras_data, plan_id: str):
    """Export all 2D areas into a single VTK file.

    Useful for visualising the complete domain at once in ParaView.
    """
    print(f"\n=== VTK export — all areas combined (plan {plan_id}) ===")

    vtk_files = ras_data.saveHEC_RAS2D_results_to_VTK(
        lastTimeStep=True,
        combined_areas=True,    # merge all areas into one VTK
    )

    print(f"Wrote {len(vtk_files)} combined VTK file(s):")
    for f in vtk_files:
        print(f"  {f}")

    return vtk_files


# ── 6. Iterate over all plans that have result HDF files ─────────────────────

def export_all_available_plans():
    """Load results and export VTK for every plan whose result HDF exists.

    This shows how to batch-process multiple plans in a project.
    """
    print("\n=== Batch VTK export for all available plans ===")

    project = HEC_RAS_Project(PRJ_FILE)
    exported = []

    for pid, plan in project.plans.items():
        if not plan.hdf_exists():
            print(f"  {pid}: skipped (no result HDF)")
            continue

        print(f"  {pid}: \"{plan.plan_name}\" — loading results ...")
        ras_data = plan.load_results()
        n_areas = len(ras_data.TwoDAreaNames)
        vtk_files = ras_data.saveHEC_RAS2D_results_to_VTK(lastTimeStep=True)
        print(f"       {n_areas} area(s), wrote {len(vtk_files)} VTK file(s)")
        exported.append((pid, vtk_files))

    if not exported:
        print("  No plans with result HDF found. Run at least one plan first.")
    return exported


# ── main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":

    import os
    # Ensure we run from the directory containing the project files
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # Step 1: inspect the project structure (no HEC-RAS needed)
    project = inspect_project()

    # Step 2: run plan p13 ("PMF with Multi 2D Areas") — the plan with three
    # 2D flow areas.  Comment this out if HEC-RAS is not installed or the
    # simulation has already been run.
    # run_plan("p13")

    # Step 3: load and inspect the multi-area results
    # (plan p13 must have been run first)
    DEMO_PLAN = "p13"   # "PMF with Multi 2D Areas" — 3 areas: 193, 194, LockHaven
    ras_data = load_and_inspect_results(DEMO_PLAN)

    if ras_data is not None:
        # Step 4: export each 2D area to a separate VTK file (default)
        export_vtk_separate(ras_data, DEMO_PLAN)

        # Step 5: export all areas merged into one VTK file
        export_vtk_combined(ras_data, DEMO_PLAN)

    # Step 6: batch-process all plans that have been run
    export_all_available_plans()

    print("\nAll done!")
