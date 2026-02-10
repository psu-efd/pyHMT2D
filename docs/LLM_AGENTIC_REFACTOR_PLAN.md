# pyHMT2D LLM Agentic Integration Refactor Plan (Experimental)

## Purpose
This plan proposes an incremental refactor to make pyHMT2D suitable for safe, repeatable LLM-driven automation of 2D hydraulic modeling workflows.

## Guiding principles
- Preserve current scientific/modeling behavior while adding integration layers around it.
- Keep the core solver/model adapters independent from any specific LLM vendor.
- Make every automated action schema-validated, auditable, and resumable.
- Prefer additive changes first, then selective refactoring where needed.

---

## Phase 0 — Discovery, architecture baseline, and safety envelope (1–2 weeks)

### Goals
- Establish clear boundaries between existing core modules and new agent-facing layers.
- Define non-goals and risk constraints (e.g., destructive file edits, unsupported solver versions).

### Work items
1. **Codebase mapping and dependency graph**
   - Map stable extension points in:
     - `pyHMT2D/Calibration/`
     - `pyHMT2D/Hydraulic_Models_Data/`
     - `pyHMT2D/Misc/`
     - `pyHMT2D/cli/`
2. **Workflow inventory**
   - Document canonical workflows to support first:
     - run model
     - calibrate model
     - convert outputs/formats
     - postprocess and summarize
3. **Safety policy draft**
   - Define file-access policy, max runtime defaults, retry/backoff policy, and sandbox assumptions.
4. **Architecture Decision Records (ADRs)**
   - ADR-001: tool API style (HTTP vs MCP vs both)
   - ADR-002: job lifecycle model
   - ADR-003: typed config system

### Exit criteria
- Approved target architecture diagram.
- Prioritized workflow list with acceptance criteria.
- Risk register with mitigation owners.

---

## Phase 1 — Stabilize core execution contracts (2–4 weeks)

### Goals
- Make core library behavior predictable for autonomous callers.
- Remove abrupt process exits and ambiguous side effects.

### Work items
1. **Exception normalization**
   - Replace `sys.exit()` in library paths with typed exceptions.
   - Introduce common error hierarchy (e.g., `PyHMT2DError`, `ConfigError`, `ModelLaunchError`, `UnsupportedPlatformError`).
2. **Structured result objects**
   - Return dataclass/Pydantic-style result objects for key operations (`init_model`, `open_project`, `run_model`, calibration steps).
3. **Logging consistency**
   - Standardize to structured JSON logging with operation IDs.
4. **Idempotency review**
   - Ensure repeated calls do not corrupt state where avoidable.

### Exit criteria
- Core operations expose stable typed exceptions.
- Logs include operation IDs and machine-readable event payloads.
- No library-level hard exits on recoverable errors.

---

## Phase 2 — Typed schemas and validation-first job specs (2–3 weeks)

### Goals
- Ensure LLM-generated inputs are validated before execution.

### Work items
1. **Define canonical schemas**
   - `RunModelRequest`
   - `CalibrationRequest`
   - `ConvertRequest`
   - `ResultSummaryRequest`
2. **Validation service layer**
   - Add preflight validation methods/CLI (`--validate-only`).
3. **Backward compatibility adapters**
   - Add translators from existing JSON examples to canonical typed schema.
4. **Config normalization**
   - Convert string booleans/None-like values to strict JSON booleans/nulls.

### Exit criteria
- All agent-facing workflows fail fast on invalid schema.
- Existing examples can be auto-normalized or produce actionable validation errors.

---

## Phase 3 — Agent-facing tool interface (MCP/HTTP) (3–5 weeks)

### Goals
- Provide a stable tool-call surface for LLM agents.

### Work items
1. **Tool façade module**
   - New package (proposed): `pyHMT2D/agent_tools/`
   - Stateless wrappers over core functions.
2. **Tool catalog (initial)**
   - `validate_job_spec`
   - `run_calibration`
   - `run_single_simulation`
   - `convert_ras_to_srh`
   - `extract_result_metrics`
   - `summarize_artifacts`
3. **Transport adapters**
   - Option A: FastAPI server (`pyHMT2D/services/http_api/`)
   - Option B: MCP server adapter (`pyHMT2D/services/mcp_server/`)
4. **Deterministic output contract**
   - Every tool returns: `status`, `error_code`, `message`, `artifacts`, `metrics`, `trace_id`.

### Exit criteria
- End-to-end agent tool calls work with deterministic JSON responses.
- Tool manifest/spec is documented for model prompting.

---

## Phase 4 — Job orchestration, persistence, and resumability (3–4 weeks)

### Goals
- Support long-running hydraulic tasks robustly.

### Work items
1. **Job runner subsystem**
   - Create queued execution model with states:
     - `pending`, `running`, `succeeded`, `failed`, `cancelled`
2. **Run manifests**
   - Persist `run_manifest.json` for each run with:
     - normalized input
     - environment metadata
     - timestamps
     - artifact inventory and hashes
3. **Cancellation/timeouts**
   - Add max runtime, cancellation hooks, graceful shutdown semantics.
4. **Retry policy**
   - Explicit retry classes (transient vs non-transient errors).

### Exit criteria
- Jobs can be resumed/inspected post-failure.
- Artifacts and metrics are reproducibly discoverable.

---

## Phase 5 — Security and governance hardening (2–3 weeks)

### Goals
- Prevent unsafe automation and improve enterprise readiness.

### Work items
1. **Path and command restrictions**
   - Enforce workspace boundary and allowlists for executable paths.
2. **Policy guardrails**
   - Validate file mutation scopes and protected outputs.
3. **Audit trail and provenance**
   - Add immutable audit events per tool invocation.
4. **Secrets/credentials handling**
   - Document secure handling for environment-bound dependencies.

### Exit criteria
- Security review checklist passes.
- High-risk actions require explicit policy approval config.

---

## Phase 6 — Test expansion and reliability gates (ongoing; first cut 2–4 weeks)

### Goals
- Raise confidence for unattended LLM operation.

### Work items
1. **Contract tests** for all tool endpoints.
2. **Schema fuzz tests** for malformed LLM outputs.
3. **Golden-path integration tests** against sample datasets.
4. **Failure-mode tests** (missing files, wrong model version, unsupported OS, timeout).
5. **CI quality gates**
   - Lint + type checks + unit tests + contract tests.

### Exit criteria
- CI blocks merges on contract/schema regressions.
- Regression suite covers top workflows and critical failure modes.

---

## Phase 7 — Documentation and rollout (1–2 weeks)

### Goals
- Make integration easy for downstream agent developers.

### Work items
1. **Agent integration guide**
   - Required tool sequence patterns.
   - Prompting best practices for deterministic specs.
2. **Reference examples**
   - End-to-end notebook and CLI/API examples.
3. **Migration notes**
   - Legacy CLI to tool API mapping.
4. **Versioning and deprecation policy**
   - Define semantic versioning for tool contracts.

### Exit criteria
- External integrator can complete first automated workflow using docs only.

---

## Recommended implementation order (first 90 days)
1. Phase 0
2. Phase 1
3. Phase 2
4. Phase 3 (minimal tool set)
5. Phase 6 (minimum reliability gate)

Phases 4, 5, and 7 should run in parallel once Phase 3 is stable.

## Proposed success metrics
- ≥95% of invalid job specs rejected at validation stage (before solver launch).
- 100% of tool calls return structured JSON with traceable `trace_id`.
- ≥90% automated workflow success on supported environments for golden datasets.
- Mean time to diagnose failed run reduced by ≥50% via run manifests/logs.

## Risks and mitigations
- **Risk:** Solver/version coupling on Windows-only tooling.
  - **Mitigation:** Isolate platform-specific adapters and clearly mark capability matrix.
- **Risk:** LLM-generated malformed configs.
  - **Mitigation:** strict schema validation + normalization + friendly error hints.
- **Risk:** Long-running job instability.
  - **Mitigation:** job state persistence, timeout policies, resumable manifests.

## Branching and release strategy
- Keep this effort in experimental branches until Phase 3 MVP is complete.
- Use feature flags for agent-facing APIs.
- Cut a preview release once Phase 3 + minimum Phase 6 gates pass.
