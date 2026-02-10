"""Dataclasses for structured request and response contracts."""

from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional
import uuid


def utc_now_iso() -> str:
    """Return current UTC time in ISO 8601 format with Z suffix."""
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def new_trace_id() -> str:
    """Generate a unique trace identifier for a tool call."""
    return str(uuid.uuid4())


@dataclass
class ArtifactRecord:
    """A generated or discovered file artifact."""

    path: str
    exists: bool
    size_bytes: Optional[int] = None
    sha256: Optional[str] = None


@dataclass
class ToolResponse:
    """Standard response shape for agent-facing tool calls."""

    status: str
    error_code: Optional[str] = None
    message: str = ""
    artifacts: List[ArtifactRecord] = field(default_factory=list)
    metrics: Dict[str, Any] = field(default_factory=dict)
    trace_id: str = field(default_factory=new_trace_id)
    started_at: str = field(default_factory=utc_now_iso)
    completed_at: Optional[str] = None

    def finalize(self) -> "ToolResponse":
        """Set completion timestamp and return self for chaining."""
        self.completed_at = utc_now_iso()
        return self

    def to_dict(self) -> Dict[str, Any]:
        """Convert response to a plain dictionary."""
        return {
            "status": self.status,
            "error_code": self.error_code,
            "message": self.message,
            "artifacts": [a.__dict__ for a in self.artifacts],
            "metrics": self.metrics,
            "trace_id": self.trace_id,
            "started_at": self.started_at,
            "completed_at": self.completed_at,
        }
