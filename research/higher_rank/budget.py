"""Explicit per-rank computation policy, retaining every earlier timed run."""
import json
from pathlib import Path


def time_limit(directory, default=3600):
    path = Path(directory) / 'computation-policy.json'
    return json.loads(path.read_text())['time_limit_seconds'] if path.exists() else default


def remaining_seconds(directory, spent, default=3600):
    limit = time_limit(directory, default)
    return float('inf') if limit is None else max(0, limit - spent)
