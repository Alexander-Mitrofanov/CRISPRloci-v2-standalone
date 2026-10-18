"""Optional, structured installer progress from operator commands."""
import json
import os

PREFIX = "CRISPR_PROGRESS "


def emit(label, completed=0, total=None, unit="items"):
    if os.environ.get("CRISPR_INSTALL_PROGRESS") == "1":
        print(PREFIX + json.dumps({"label": str(label), "completed": completed,
                                   "total": total, "unit": unit}), flush=True)
