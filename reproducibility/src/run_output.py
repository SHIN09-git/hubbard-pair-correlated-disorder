from pathlib import Path
from datetime import datetime, timezone
from uuid import uuid4
import os

def evidence_path(filename: str) -> Path:
    """Use caller-supplied run folder or a fresh one. Never overwrite output."""
    folder=os.environ.get('PHY_PAIR_OUTPUT_DIR')
    root=Path(folder) if folder else Path(__file__).resolve().parents[1]/'results'/'runs'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid4().hex[:8])
    root.mkdir(parents=True,exist_ok=True)
    path=root/filename
    if path.exists():raise FileExistsError(f'Refusing to overwrite {path}')
    return path
