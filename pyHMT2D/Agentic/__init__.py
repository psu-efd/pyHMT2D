"""Agent-facing contracts and tool wrappers for pyHMT2D."""

from .errors import (
    PyHMT2DError,
    ConfigError,
    ValidationError,
    ModelLaunchError,
    UnsupportedPlatformError,
    JobExecutionError,
)
from .models import ToolResponse, ArtifactRecord
from .schemas import (
    validate_run_model_request,
    validate_calibration_request,
    validate_convert_request,
    validate_result_summary_request,
    validate_modify_srh_case_request,
    validate_monte_carlo_srh_request,
)
from .tools import (
    validate_job_spec,
    run_single_simulation,
    run_calibration,
    convert_ras_to_srh,
    summarize_artifacts,
    modify_srh_case_parameters,
    run_monte_carlo_srh,
)
from .job_runner import run_job, JobHandle

__all__ = [
    "PyHMT2DError",
    "ConfigError",
    "ValidationError",
    "ModelLaunchError",
    "UnsupportedPlatformError",
    "JobExecutionError",
    "ToolResponse",
    "ArtifactRecord",
    "validate_run_model_request",
    "validate_calibration_request",
    "validate_convert_request",
    "validate_result_summary_request",
    "validate_modify_srh_case_request",
    "validate_monte_carlo_srh_request",
    "validate_job_spec",
    "run_single_simulation",
    "run_calibration",
    "convert_ras_to_srh",
    "summarize_artifacts",
    "modify_srh_case_parameters",
    "run_monte_carlo_srh",
    "run_job",
    "JobHandle",
]
