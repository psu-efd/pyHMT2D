"""
Example test runner for pyHMT2D.

Runs all examples, records pass/fail, and saves a structured
JSON log to examples/example_runs_results.json.

Usage:
    python run_all_examples.py                          # run everything
    python run_all_examples.py --dry-run                # list what would run
    python run_all_examples.py --solver ras              # HEC-RAS examples only
    python run_all_examples.py --category ai_tools       # AI Tools examples only
    python run_all_examples.py --name "01_open"          # substring filter

Requirements by category:
    traditional & cli:
        - pyHMT2D installed:  pip install -e .
        - Solver(s) installed (HEC-RAS 6.6 and/or SRH-2D) for examples
          that run simulations.

    ai_tools (currently uses Claude Code CLI; see note below):
        - Claude Code CLI installed and on PATH (npm install -g @anthropic-ai/claude-code)
        - Authentication: Anthropic API key (ANTHROPIC_API_KEY) or Claude Max/Pro subscription
        - pyHMT2D MCP server registered:
              claude mcp add pyHMT2D -- python -m pyHMT2D.AI_Tools.mcp_server
        - pyHMT2D installed with AI extras:  pip install -e ".[ai]"
        - Solver(s) installed for examples that run simulations

    Note: The AI Tools examples are designed around the MCP (Model Context
    Protocol) server, which is an open standard. Any AI coding assistant that
    supports MCP (e.g., Claude Code, Cursor, Windsurf, Codex, etc.) can run
    these examples interactively. The automated test runner currently invokes
    Claude Code CLI for non-interactive testing, but the examples themselves
    are not tied to any specific AI provider.
"""

import argparse
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

EXAMPLES_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = EXAMPLES_DIR.parent

DEFAULT_TIMEOUT = 600       # 10 min
LONG_TIMEOUT = 1800         # 30 min
NO_TIMEOUT = None           # unlimited – calibration/MC may run for hours
AI_TOOLS_TIMEOUT = 900      # 15 min

STDOUT_TAIL_LINES = 50
STDERR_TAIL_LINES = 50

# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class ExampleEntry:
    name: str
    category: str               # "traditional" | "ai_tools" | "cli"
    solver: str                 # "ras" | "srh" | "none"
    script: Optional[str]       # relative path to .py from EXAMPLES_DIR (None for cli / ai_tools)
    cwd: str                    # working directory relative to EXAMPLES_DIR
    args: list = field(default_factory=list)
    timeout: int = DEFAULT_TIMEOUT
    tags: list = field(default_factory=list)
    precondition_files: list = field(default_factory=list)
    cli_args: Optional[list] = None   # for cli category: [tool_name, json_args_string]


@dataclass
class TestResult:
    name: str
    category: str
    solver: str
    status: str                 # PASS | FAIL | SKIP | TIMEOUT | WARN
    exit_code: Optional[int]
    duration_seconds: float
    started_at: str
    error_message: str = ""
    stdout_tail: str = ""
    stderr_tail: str = ""
    command: str = ""

# ---------------------------------------------------------------------------
# Example registry
# ---------------------------------------------------------------------------

EXAMPLE_REGISTRY: List[ExampleEntry] = [
    # ---- Traditional: HEC-RAS ----
    ExampleEntry(
        name="HEC_RAS_Model/Muncie2D",
        category="traditional",
        solver="ras",
        script="HEC_RAS_Model/Muncie2D/demo_HEC_RAS_Model.py",
        cwd="HEC_RAS_Model/Muncie2D",
        tags=["ras", "model"],
    ),
    ExampleEntry(
        name="HEC_RAS_Model/BaldEagleCrkMulti2D",
        category="traditional",
        solver="ras",
        script="HEC_RAS_Model/BaldEagleCrkMulti2D/demo_HEC_RAS_Model.py",
        cwd="HEC_RAS_Model/BaldEagleCrkMulti2D",
        tags=["ras", "model", "multi2d"],
        precondition_files=["*.prj"],
    ),
    ExampleEntry(
        name="calibration/RAS-2D",
        category="traditional",
        solver="ras",
        script="calibration/RAS-2D/Munice2D_ManningN_calibration/demo_calibration_RAS_2D.py",
        cwd="calibration/RAS-2D/Munice2D_ManningN_calibration",
        tags=["ras", "calibration"],
        timeout=NO_TIMEOUT,
    ),
    ExampleEntry(
        name="Monte_Carlo/RAS_2D",
        category="traditional",
        solver="ras",
        script="Monte_Carlo/RAS_2D/Munice2D_ManningN_Monte_Carlo/demo_HEC_RAS_Monte_Carlo.py",
        cwd="Monte_Carlo/RAS_2D/Munice2D_ManningN_Monte_Carlo",
        args=["sampledManningN_2022_03_10-10_24_34_PM.dat"],
        tags=["ras", "monte_carlo"],
        timeout=NO_TIMEOUT,
        precondition_files=["sampledManningN_*.dat"],
    ),

    # ---- Traditional: SRH-2D ----
    ExampleEntry(
        name="SRH_2D_Model",
        category="traditional",
        solver="srh",
        script="SRH_2D_Model/demo_SRH_2D_Model.py",
        cwd="SRH_2D_Model",
        tags=["srh", "model"],
    ),
    ExampleEntry(
        name="SRH_2D_Model_SIF",
        category="traditional",
        solver="srh",
        script="SRH_2D_Model_SIF/demo_SRH_2D_Model_w_SIF.py",
        cwd="SRH_2D_Model_SIF",
        tags=["srh", "model"],
    ),
    ExampleEntry(
        name="calibration/SRH-2D (GP)",
        category="traditional",
        solver="srh",
        script="calibration/SRH-2D/Munice2D_ManningN_calibration/demo_calibration_SRH_2D_gp.py",
        cwd="calibration/SRH-2D/Munice2D_ManningN_calibration",
        tags=["srh", "calibration"],
        timeout=NO_TIMEOUT,
    ),
    ExampleEntry(
        name="Monte_Carlo/SRH_2D",
        category="traditional",
        solver="srh",
        script="Monte_Carlo/SRH_2D/Munice2D_ManningN_Monte_Carlo/demo_SRH_2D_Monte_Carlo.py",
        cwd="Monte_Carlo/SRH_2D/Munice2D_ManningN_Monte_Carlo",
        args=["sampledManningN_2022_03_10-10_24_34_PM.dat"],
        tags=["srh", "monte_carlo"],
        timeout=NO_TIMEOUT,
        precondition_files=["sampledManningN_*.dat"],
    ),

    # ---- Traditional: post-processing (no solver) ----
    ExampleEntry(
        name="compare/process_RAS_2D_Data",
        category="traditional",
        solver="none",
        script="compare_SRH_2D_RAS_2D/Muncie2D/process_RAS_2D_Data.py",
        cwd="compare_SRH_2D_RAS_2D/Muncie2D",
        tags=["postprocess"],
        precondition_files=["HEC-RAS/*.hdf"],
    ),
    ExampleEntry(
        name="compare/process_SRH_2D_Data",
        category="traditional",
        solver="none",
        script="compare_SRH_2D_RAS_2D/Muncie2D/process_SRH_2D_Data.py",
        cwd="compare_SRH_2D_RAS_2D/Muncie2D",
        tags=["postprocess"],
        precondition_files=["SRH-2D/*.srhhydro"],
    ),
    ExampleEntry(
        name="compare/compare_SRH_2D_HEC_RAS_2D",
        category="traditional",
        solver="none",
        script="compare_SRH_2D_RAS_2D/Muncie2D/compare_SRH_2D_HEC_RAS_2D.py",
        cwd="compare_SRH_2D_RAS_2D/Muncie2D",
        tags=["postprocess"],
        precondition_files=["HEC-RAS/*.vtk", "SRH-2D/*.vtk"],
    ),

    # ---- CLI examples ----
    ExampleEntry(
        name="cli/ras_to_srh",
        category="cli",
        solver="none",
        script=None,
        cwd="cli/ras_to_srh/Muncie",
        tags=["cli", "conversion"],
        precondition_files=["*.hdf"],
        cli_args=["ras_to_srh", '{"ras_hdf_file": "Muncie2D.p01.hdf", "srh_case_name": "srh_Muncie"}'],
    ),
    ExampleEntry(
        name="cli/srh_to_vtk",
        category="cli",
        solver="none",
        script=None,
        cwd="cli/srh_to_vtk",
        tags=["cli", "conversion"],
        precondition_files=["Muncie_XMDFC.h5"],
        cli_args=["srh_to_vtk", '{"srhhydro_file": "Muncie.srhhydro", "output_file": "Muncie_XMDFC.h5"}'],
    ),
    ExampleEntry(
        name="cli/srh_mesh_to_vtk",
        category="cli",
        solver="none",
        script=None,
        cwd="cli/srh_mesh_to_vtk",
        tags=["cli", "conversion"],
        cli_args=["srh_to_vtk", '{"srhhydro_file": "Muncie.srhhydro", "output_file": "case_mesh.vtk", "mesh_only": true}'],
    ),

    # ---- AI Tools: RAS_2D ----
    ExampleEntry(
        name="AI_Tools/RAS_2D/01_open_inspect",
        category="ai_tools",
        solver="ras",
        script=None,
        cwd="AI_Tools/RAS_2D/01_open_inspect",
        tags=["ai_tools", "ras"],
        timeout=AI_TOOLS_TIMEOUT,
    ),
    ExampleEntry(
        name="AI_Tools/RAS_2D/02_modify_and_run",
        category="ai_tools",
        solver="ras",
        script=None,
        cwd="AI_Tools/RAS_2D/02_modify_and_run",
        tags=["ai_tools", "ras"],
        timeout=AI_TOOLS_TIMEOUT,
    ),
    ExampleEntry(
        name="AI_Tools/RAS_2D/03_query_results",
        category="ai_tools",
        solver="ras",
        script=None,
        cwd="AI_Tools/RAS_2D/03_query_results",
        tags=["ai_tools", "ras"],
        timeout=AI_TOOLS_TIMEOUT,
    ),
    ExampleEntry(
        name="AI_Tools/RAS_2D/04_calibration",
        category="ai_tools",
        solver="ras",
        script=None,
        cwd="AI_Tools/RAS_2D/04_calibration",
        tags=["ai_tools", "ras", "calibration"],
        timeout=NO_TIMEOUT,
    ),
    ExampleEntry(
        name="AI_Tools/RAS_2D/05_monte_carlo",
        category="ai_tools",
        solver="ras",
        script=None,
        cwd="AI_Tools/RAS_2D/05_monte_carlo",
        tags=["ai_tools", "ras", "monte_carlo"],
        timeout=NO_TIMEOUT,
    ),

    # ---- AI Tools: SRH_2D ----
    ExampleEntry(
        name="AI_Tools/SRH_2D/01_open_inspect",
        category="ai_tools",
        solver="srh",
        script=None,
        cwd="AI_Tools/SRH_2D/01_open_inspect",
        tags=["ai_tools", "srh"],
        timeout=AI_TOOLS_TIMEOUT,
    ),
    ExampleEntry(
        name="AI_Tools/SRH_2D/02_modify_and_run",
        category="ai_tools",
        solver="srh",
        script=None,
        cwd="AI_Tools/SRH_2D/02_modify_and_run",
        tags=["ai_tools", "srh"],
        timeout=AI_TOOLS_TIMEOUT,
    ),
    ExampleEntry(
        name="AI_Tools/SRH_2D/03_query_results",
        category="ai_tools",
        solver="srh",
        script=None,
        cwd="AI_Tools/SRH_2D/03_query_results",
        tags=["ai_tools", "srh"],
        timeout=AI_TOOLS_TIMEOUT,
    ),
    ExampleEntry(
        name="AI_Tools/SRH_2D/04_calibration",
        category="ai_tools",
        solver="srh",
        script=None,
        cwd="AI_Tools/SRH_2D/04_calibration",
        tags=["ai_tools", "srh", "calibration"],
        timeout=NO_TIMEOUT,
    ),
    ExampleEntry(
        name="AI_Tools/SRH_2D/05_monte_carlo",
        category="ai_tools",
        solver="srh",
        script=None,
        cwd="AI_Tools/SRH_2D/05_monte_carlo",
        tags=["ai_tools", "srh", "monte_carlo"],
        timeout=NO_TIMEOUT,
    ),
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def tail(text: str, n: int) -> str:
    """Return the last *n* lines of *text*."""
    lines = text.splitlines()
    return "\n".join(lines[-n:]) if len(lines) > n else text


def extract_readme_prompt(readme_path: Path, index: int = 0) -> str:
    """Extract the *index*-th ```text``` fenced block from a README."""
    content = readme_path.read_text(encoding="utf-8")
    blocks = re.findall(r"```text\s*\n(.*?)```", content, re.DOTALL)
    if not blocks:
        raise ValueError(f"No ```text``` block found in {readme_path}")
    if index >= len(blocks):
        raise ValueError(
            f"Only {len(blocks)} ```text``` block(s) in {readme_path}, "
            f"but index {index} requested"
        )
    return blocks[index].strip()


def convert_claude_json_to_markdown(json_text: str) -> str:
    """Convert Claude CLI JSON output to a human-readable markdown conversation."""
    try:
        data = json.loads(json_text)
    except (json.JSONDecodeError, TypeError):
        return json_text  # fallback: return raw text

    if not isinstance(data, list):
        return json_text

    lines = ["# AI Tool Conversation Log\n"]

    # Extract summary from the result entry (last item)
    for item in data:
        if item.get("type") == "result":
            cost = item.get("total_cost_usd")
            duration = item.get("duration_ms")
            turns = item.get("num_turns")
            parts = []
            if turns is not None:
                parts.append(f"**Turns:** {turns}")
            if duration is not None:
                parts.append(f"**Duration:** {duration / 1000:.1f}s")
            if cost is not None:
                parts.append(f"**Cost:** ${cost:.4f}")
            if parts:
                lines.append(" | ".join(parts) + "\n")
            lines.append("---\n")
            break

    for item in data:
        msg_type = item.get("type")

        if msg_type == "system":
            continue
        if msg_type == "rate_limit_event":
            continue
        if msg_type == "result":
            continue

        message = item.get("message")
        if not message:
            continue

        role = message.get("role", "unknown")
        content_blocks = message.get("content", [])

        for block in content_blocks:
            btype = block.get("type")

            if btype == "thinking":
                # Skip thinking blocks to keep the log concise
                continue

            elif btype == "text":
                text = block.get("text", "")
                if role == "assistant":
                    lines.append(f"## Assistant\n\n{text}\n")
                else:
                    # user text blocks are usually system-injected context
                    lines.append(f"## System Context\n\n{text}\n")

            elif btype == "tool_use":
                tool_name = block.get("name", "unknown")
                tool_input = block.get("input", {})
                lines.append(f"### Tool Call: `{tool_name}`\n")
                if tool_input:
                    input_str = json.dumps(tool_input, indent=2, ensure_ascii=False)
                    # Truncate very long inputs
                    if len(input_str) > 2000:
                        input_str = input_str[:2000] + "\n... (truncated)"
                    lines.append(f"```json\n{input_str}\n```\n")

            elif btype == "tool_result":
                content_parts = block.get("content", "")
                is_error = block.get("is_error", False)
                label = "Tool Error" if is_error else "Tool Result"

                # content can be a string or a list of content blocks
                if isinstance(content_parts, list):
                    text_parts = []
                    for part in content_parts:
                        if isinstance(part, dict) and part.get("type") == "text":
                            text_parts.append(part.get("text", ""))
                        elif isinstance(part, str):
                            text_parts.append(part)
                    result_text = "\n".join(text_parts)
                else:
                    result_text = str(content_parts)

                # Truncate very long results
                if len(result_text) > 3000:
                    result_text = result_text[:3000] + "\n... (truncated)"
                lines.append(f"**{label}:**\n```\n{result_text}\n```\n")

    return "\n".join(lines)


def check_preconditions(entry: ExampleEntry) -> tuple:
    """Return (ok: bool, reason: str). Checks glob patterns in cwd."""
    import glob as globmod

    cwd = EXAMPLES_DIR / entry.cwd
    for pattern in entry.precondition_files:
        matches = globmod.glob(str(cwd / pattern))
        if not matches:
            return False, f"Missing files matching '{pattern}' in {cwd}"
    return True, ""


def check_claude_cli(claude_path: str) -> bool:
    """Return True if the Claude CLI is reachable."""
    return shutil.which(claude_path) is not None

# ---------------------------------------------------------------------------
# Runners
# ---------------------------------------------------------------------------

def run_traditional(entry: ExampleEntry, verbose: bool = False) -> TestResult:
    """Run a traditional Python demo script via subprocess."""
    script_path = EXAMPLES_DIR / entry.script
    cwd = EXAMPLES_DIR / entry.cwd
    cmd = [sys.executable, str(script_path)] + entry.args
    cmd_str = " ".join(cmd)

    started = datetime.now(timezone.utc).isoformat()
    t0 = time.monotonic()

    try:
        proc = subprocess.run(
            cmd,
            cwd=str(cwd),
            timeout=entry.timeout,
            capture_output=True,
            text=True,
        )
        elapsed = time.monotonic() - t0

        if verbose:
            if proc.stdout:
                print(proc.stdout)
            if proc.stderr:
                print(proc.stderr, file=sys.stderr)

        status = "PASS" if proc.returncode == 0 else "FAIL"
        error_msg = ""
        if proc.returncode != 0:
            error_msg = f"Exit code {proc.returncode}"

        return TestResult(
            name=entry.name,
            category=entry.category,
            solver=entry.solver,
            status=status,
            exit_code=proc.returncode,
            duration_seconds=round(elapsed, 2),
            started_at=started,
            error_message=error_msg,
            stdout_tail=tail(proc.stdout, STDOUT_TAIL_LINES),
            stderr_tail=tail(proc.stderr, STDERR_TAIL_LINES),
            command=cmd_str,
        )
    except subprocess.TimeoutExpired:
        elapsed = time.monotonic() - t0
        return TestResult(
            name=entry.name,
            category=entry.category,
            solver=entry.solver,
            status="TIMEOUT",
            exit_code=None,
            duration_seconds=round(elapsed, 2),
            started_at=started,
            error_message=f"Timed out after {entry.timeout}s",
            command=cmd_str,
        )


def run_cli(entry: ExampleEntry, verbose: bool = False) -> TestResult:
    """Run an hmt-cli command directly (no shell needed)."""
    cwd = EXAMPLES_DIR / entry.cwd
    tool_name, json_args = entry.cli_args
    cmd = ["hmt-cli", tool_name, "--args", json_args]
    cmd_str = f"hmt-cli {tool_name} --args '{json_args}'"

    started = datetime.now(timezone.utc).isoformat()
    t0 = time.monotonic()

    try:
        proc = subprocess.run(
            cmd,
            cwd=str(cwd),
            timeout=entry.timeout,
            capture_output=True,
            text=True,
        )
        elapsed = time.monotonic() - t0

        if verbose:
            if proc.stdout:
                print(proc.stdout)
            if proc.stderr:
                print(proc.stderr, file=sys.stderr)

        status = "PASS" if proc.returncode == 0 else "FAIL"
        error_msg = ""
        if proc.returncode != 0:
            error_msg = f"Exit code {proc.returncode}"

        return TestResult(
            name=entry.name,
            category=entry.category,
            solver=entry.solver,
            status=status,
            exit_code=proc.returncode,
            duration_seconds=round(elapsed, 2),
            started_at=started,
            error_message=error_msg,
            stdout_tail=tail(proc.stdout, STDOUT_TAIL_LINES),
            stderr_tail=tail(proc.stderr, STDERR_TAIL_LINES),
            command=cmd_str,
        )
    except subprocess.TimeoutExpired:
        elapsed = time.monotonic() - t0
        return TestResult(
            name=entry.name,
            category=entry.category,
            solver=entry.solver,
            status="TIMEOUT",
            exit_code=None,
            duration_seconds=round(elapsed, 2),
            started_at=started,
            error_message=f"Timed out after {entry.timeout}s",
            command=cmd_str,
        )


def run_ai_tools(
    entry: ExampleEntry,
    claude_path: str = "claude",
    verbose: bool = False,
) -> TestResult:
    """Run an AI Tools example by invoking Claude CLI with the README prompt."""
    cwd = EXAMPLES_DIR / entry.cwd
    readme_path = cwd / "README.md"

    started = datetime.now(timezone.utc).isoformat()

    # Extract the prompt
    try:
        prompt = extract_readme_prompt(readme_path, index=0)
    except (ValueError, FileNotFoundError) as exc:
        return TestResult(
            name=entry.name,
            category=entry.category,
            solver=entry.solver,
            status="SKIP",
            exit_code=None,
            duration_seconds=0.0,
            started_at=started,
            error_message=str(exc),
        )

    # Build the allowed tools list: Bash, core tools, Skill, and all MCP tools
    allowed_tools = [
        "Bash", "Read", "Edit", "Write", "Glob", "Grep",
        "Skill", "ToolSearch",
        "mcp__pyHMT2D__list_model_files",
        "mcp__pyHMT2D__get_project_info",
        "mcp__pyHMT2D__get_materials",
        "mcp__pyHMT2D__get_boundary_conditions",
        "mcp__pyHMT2D__get_result_variables",
        "mcp__pyHMT2D__set_manning_n",
        "mcp__pyHMT2D__set_inlet_flow",
        "mcp__pyHMT2D__set_exit_wse",
        "mcp__pyHMT2D__save_modified_inputs",
        "mcp__pyHMT2D__run_preprocessing",
        "mcp__pyHMT2D__run_simulation",
        "mcp__pyHMT2D__get_simulation_status",
        "mcp__pyHMT2D__exit_model",
        "mcp__pyHMT2D__close_session",
        "mcp__pyHMT2D__get_status",
        "mcp__pyHMT2D__read_results",
        "mcp__pyHMT2D__get_result_statistics",
        "mcp__pyHMT2D__get_value_at_point",
        "mcp__pyHMT2D__get_flood_extent",
        "mcp__pyHMT2D__get_cross_section_profile",
        "mcp__pyHMT2D__export_to_vtk",
        "mcp__pyHMT2D__export_mesh_to_vtk",
        "mcp__pyHMT2D__load_observations",
        "mcp__pyHMT2D__evaluate_parameters",
        "mcp__pyHMT2D__run_calibration",
        "mcp__pyHMT2D__generate_mc_samples",
        "mcp__pyHMT2D__run_monte_carlo",
        "mcp__pyHMT2D__get_mc_statistics",
        "mcp__pyHMT2D__ras_to_srh",
        "mcp__pyHMT2D__srh_to_vtk",
    ]

    cmd = [
        claude_path,
        "-p", prompt,
        "--output-format", "json",
        "--verbose",
        "--allowedTools", ",".join(allowed_tools),
    ]
    cmd_str = f'{claude_path} -p "<prompt from README>" (cwd: {cwd})'

    t0 = time.monotonic()
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(cwd),
            timeout=entry.timeout,
            capture_output=True,
            text=True,
        )
        elapsed = time.monotonic() - t0

        if verbose:
            if proc.stdout:
                print(proc.stdout)
            if proc.stderr:
                print(proc.stderr, file=sys.stderr)

        # Save full conversation to files in the example directory
        raw_output = proc.stdout or ""
        try:
            (cwd / "ai_tool_output.json").write_text(raw_output, encoding="utf-8")
            md_text = convert_claude_json_to_markdown(raw_output)
            (cwd / "ai_tool_output.md").write_text(md_text, encoding="utf-8")
        except OSError:
            pass  # don't fail the test over a log write error

        # Determine status from exit code and output content
        output_combined = (proc.stdout or "") + (proc.stderr or "")
        has_error_indicators = bool(
            re.search(r'"status"\s*:\s*"error"', output_combined)
            or re.search(r"Traceback \(most recent call last\)", output_combined)
        )

        if proc.returncode != 0:
            status = "FAIL"
            error_msg = f"Exit code {proc.returncode}"
        elif has_error_indicators:
            status = "WARN"
            error_msg = "Exit code 0 but error indicators found in output"
        else:
            status = "PASS"
            error_msg = ""

        return TestResult(
            name=entry.name,
            category=entry.category,
            solver=entry.solver,
            status=status,
            exit_code=proc.returncode,
            duration_seconds=round(elapsed, 2),
            started_at=started,
            error_message=error_msg,
            stdout_tail=tail(proc.stdout, STDOUT_TAIL_LINES),
            stderr_tail=tail(proc.stderr, STDERR_TAIL_LINES),
            command=cmd_str,
        )
    except subprocess.TimeoutExpired:
        elapsed = time.monotonic() - t0
        return TestResult(
            name=entry.name,
            category=entry.category,
            solver=entry.solver,
            status="TIMEOUT",
            exit_code=None,
            duration_seconds=round(elapsed, 2),
            started_at=started,
            error_message=f"Timed out after {entry.timeout}s",
            command=cmd_str,
        )


def run_one(
    entry: ExampleEntry,
    claude_path: str = "claude",
    verbose: bool = False,
) -> TestResult:
    """Run a single example, checking preconditions first."""
    # Check preconditions
    ok, reason = check_preconditions(entry)
    if not ok:
        return TestResult(
            name=entry.name,
            category=entry.category,
            solver=entry.solver,
            status="SKIP",
            exit_code=None,
            duration_seconds=0.0,
            started_at=datetime.now(timezone.utc).isoformat(),
            error_message=reason,
        )

    # Check Claude CLI for ai_tools
    if entry.category == "ai_tools" and not check_claude_cli(claude_path):
        return TestResult(
            name=entry.name,
            category=entry.category,
            solver=entry.solver,
            status="SKIP",
            exit_code=None,
            duration_seconds=0.0,
            started_at=datetime.now(timezone.utc).isoformat(),
            error_message=f"Claude CLI not found at '{claude_path}'",
        )

    if entry.category == "traditional":
        return run_traditional(entry, verbose=verbose)
    elif entry.category == "cli":
        return run_cli(entry, verbose=verbose)
    elif entry.category == "ai_tools":
        return run_ai_tools(entry, claude_path=claude_path, verbose=verbose)
    else:
        return TestResult(
            name=entry.name,
            category=entry.category,
            solver=entry.solver,
            status="SKIP",
            exit_code=None,
            duration_seconds=0.0,
            started_at=datetime.now(timezone.utc).isoformat(),
            error_message=f"Unknown category: {entry.category}",
        )

# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

def filter_entries(
    entries: List[ExampleEntry],
    solver: str = "all",
    category: str = "all",
    tags: Optional[List[str]] = None,
    name_pattern: Optional[str] = None,
) -> List[ExampleEntry]:
    """Filter the registry by solver, category, tags, and name substring."""
    result = entries

    if solver != "all":
        result = [e for e in result if e.solver == solver]

    if category != "all":
        result = [e for e in result if e.category == category]

    if tags:
        result = [e for e in result if any(t in e.tags for t in tags)]

    if name_pattern:
        result = [e for e in result if name_pattern.lower() in e.name.lower()]

    return result

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def print_summary_table(results: List[TestResult]) -> None:
    """Print a human-readable summary table to stdout."""
    # Header
    print()
    print("=" * 80)
    print(f"{'Example':<45} {'Status':<8} {'Time':>8}  Notes")
    print("-" * 80)

    for r in results:
        time_str = f"{r.duration_seconds:.1f}s" if r.duration_seconds > 0 else "-"
        notes = r.error_message[:30] if r.error_message else ""
        print(f"{r.name:<45} {r.status:<8} {time_str:>8}  {notes}")

    # Summary counts
    total = len(results)
    counts = {}
    for r in results:
        counts[r.status] = counts.get(r.status, 0) + 1

    print("-" * 80)
    parts = [f"Total: {total}"]
    for status in ["PASS", "FAIL", "SKIP", "TIMEOUT", "WARN"]:
        if counts.get(status, 0) > 0:
            parts.append(f"{status}: {counts[status]}")
    print("  ".join(parts))
    print("=" * 80)
    print()


def save_log(results: List[TestResult], output_path: Path, filters: dict) -> None:
    """Write structured JSON log to disk."""
    summary = {}
    for r in results:
        summary[r.status] = summary.get(r.status, 0) + 1
    summary["total"] = len(results)

    log = {
        "run_id": datetime.now(timezone.utc).isoformat(),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "filters": filters,
        "summary": summary,
        "results": [asdict(r) for r in results],
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(log, f, indent=2, ensure_ascii=False)

    print(f"Log saved to {output_path}")

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run all pyHMT2D examples and record results.",
        epilog=(
            "AI Tools requirements:\n"
            "  Claude Code CLI installed and on PATH\n"
            "  Authentication via ANTHROPIC_API_KEY or Claude subscription\n"
            "  MCP server registered: claude mcp add pyHMT2D -- python -m pyHMT2D.AI_Tools.mcp_server\n"
            "  pyHMT2D installed with AI extras: pip install -e '.[ai]'\n"
            "\n"
            "The AI Tools examples use the MCP open standard and can also be\n"
            "run interactively with any MCP-compatible AI assistant (Cursor,\n"
            "Windsurf, Codex, etc.). The automated runner uses Claude Code CLI."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--solver",
        choices=["ras", "srh", "none", "all"],
        default="all",
        help="Filter by solver requirement (default: all)",
    )
    parser.add_argument(
        "--category",
        choices=["traditional", "ai_tools", "cli", "all"],
        default="all",
        help="Filter by example category (default: all)",
    )
    parser.add_argument(
        "--tag",
        action="append",
        dest="tags",
        help="Only run examples with this tag (repeatable)",
    )
    parser.add_argument(
        "--name",
        dest="name_pattern",
        help="Only run examples whose name contains this substring",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="List examples that would run without executing them",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=None,
        help="Override default timeout (seconds) for all examples",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=EXAMPLES_DIR / "example_runs_results.json",
        help="Output JSON log path (default: examples/example_runs_results.json)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print stdout/stderr to console as tests run",
    )
    parser.add_argument(
        "--claude-path",
        default="claude",
        help="Path to the Claude CLI executable (default: claude)",
    )
    parser.add_argument(
        "-y", "--yes",
        action="store_true",
        help="Skip confirmation prompt for AI Tools permissions",
    )
    parser.add_argument(
        "--cleanup-ras",
        action="store_true",
        help="Kill lingering Ras.exe processes after each RAS example",
    )
    return parser

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = build_parser()
    args = parser.parse_args()

    # Filter the registry
    entries = filter_entries(
        EXAMPLE_REGISTRY,
        solver=args.solver,
        category=args.category,
        tags=args.tags,
        name_pattern=args.name_pattern,
    )

    if not entries:
        print("No examples match the given filters.")
        sys.exit(0)

    # Apply global timeout override
    if args.timeout is not None:
        for e in entries:
            e.timeout = args.timeout

    # Warn about AI Tools permissions
    ai_entries = [e for e in entries if e.category == "ai_tools"]
    if ai_entries and not args.dry_run and not args.yes:
        print(f"\nWARNING: {len(ai_entries)} AI Tools example(s) will be run.")
        print("This invokes the Claude CLI with --allowedTools, which auto-approves:")
        print("  - Bash commands (read/write)")
        print("  - File edits (Read, Edit, Write)")
        print("  - All pyHMT2D MCP tools (model control, simulation, calibration)")
        print("  - Skill invocations")
        print()
        answer = input("Continue? [y/N] ").strip().lower()
        if answer not in ("y", "yes"):
            print("Aborted.")
            sys.exit(0)

    # Dry run
    if args.dry_run:
        print(f"\n{'#':<4} {'Name':<45} {'Category':<12} {'Solver':<6} Tags")
        print("-" * 90)
        for i, e in enumerate(entries, 1):
            tags_str = ", ".join(e.tags)
            print(f"{i:<4} {e.name:<45} {e.category:<12} {e.solver:<6} {tags_str}")
        print(f"\nTotal: {len(entries)} example(s)")
        return

    # Run all examples
    results: List[TestResult] = []
    for i, entry in enumerate(entries, 1):
        print(f"[{i}/{len(entries)}] Running: {entry.name} ...", end=" ", flush=True)
        result = run_one(
            entry,
            claude_path=args.claude_path,
            verbose=args.verbose,
        )
        results.append(result)
        print(result.status)

        # Optional RAS cleanup
        if args.cleanup_ras and entry.solver == "ras":
            subprocess.run(
                ["taskkill", "/F", "/IM", "Ras.exe"],
                capture_output=True,
            )

    # Output
    filters_info = {
        "solver": args.solver,
        "category": args.category,
        "tags": args.tags,
        "name_pattern": args.name_pattern,
    }
    print_summary_table(results)
    save_log(results, args.output, filters_info)

    # Exit with non-zero if any failed
    has_failures = any(r.status in ("FAIL", "TIMEOUT") for r in results)
    sys.exit(1 if has_failures else 0)


if __name__ == "__main__":
    main()
