"""Custom exception hierarchy for pyHMT2D agent-facing workflows."""


class PyHMT2DError(Exception):
    """Base error class for pyHMT2D agent-facing utilities."""


class ConfigError(PyHMT2DError):
    """Raised for invalid or missing configuration content."""


class ValidationError(ConfigError):
    """Raised when request schemas fail validation."""


class ModelLaunchError(PyHMT2DError):
    """Raised when a hydraulic model cannot be initialized or run."""


class UnsupportedPlatformError(PyHMT2DError):
    """Raised when a request needs an unsupported platform capability."""


class JobExecutionError(PyHMT2DError):
    """Raised for runtime failures in tool/job execution."""
