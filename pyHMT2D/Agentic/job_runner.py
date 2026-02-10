"""Simple local job runner with persisted status manifests."""

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Any
import json
import uuid

from .models import utc_now_iso


@dataclass
class JobHandle:
    """Small handle containing a job id and manifest path."""

    job_id: str
    manifest_path: Path


def run_job(job_dir: str, job_type: str, fn: Callable[[Dict[str, Any]], Any], payload: Dict[str, Any]) -> JobHandle:
    """Run a synchronous job and persist state transitions to a manifest."""
    root = Path(job_dir)
    root.mkdir(parents=True, exist_ok=True)

    job_id = str(uuid.uuid4())
    manifest_path = root / f"job_{job_id}.json"

    manifest = {
        "job_id": job_id,
        "job_type": job_type,
        "status": "running",
        "created_at": utc_now_iso(),
        "payload": payload,
    }
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    try:
        result = fn(payload)
        manifest["status"] = "succeeded" if result.status == "succeeded" else "failed"
        manifest["result"] = result.to_dict()
    except Exception as exc:  # noqa: BLE001
        manifest["status"] = "failed"
        manifest["error"] = str(exc)
    finally:
        manifest["completed_at"] = utc_now_iso()
        with open(manifest_path, "w", encoding="utf-8") as f:
            json.dump(manifest, f, indent=2)

    return JobHandle(job_id=job_id, manifest_path=manifest_path)
