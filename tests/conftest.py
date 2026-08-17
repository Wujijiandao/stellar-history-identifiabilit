"""Repository-wide pytest path bootstrap.

All tests must be independently runnable; no test module may rely on another
module having mutated sys.path earlier in collection order.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'src'
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
