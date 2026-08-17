"""Re-run E5-C5 and E5-C6 using the frozen E5-C2/E5-C3 derived history tables.

This path does not require redistribution of the third-party MORS track archive.
"""
from pathlib import Path
import shutil
import subprocess
import sys
import os

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
(RESULTS / "E5C2").mkdir(parents=True, exist_ok=True)
(RESULTS / "E5C3").mkdir(parents=True, exist_ok=True)

shutil.copy2(ROOT / "data/derived/toi700_rotation_conditioned_histories.csv",
             RESULTS / "E5C2/e5c2_mors_rotation_conditioned_histories.csv")
shutil.copy2(ROOT / "data/derived/lhs1140_rotation_compatible_histories.csv",
             RESULTS / "E5C3/e5c3_rotation_compatible_histories.csv")

env = dict(os.environ)
env['PYTHONPATH'] = str(ROOT / 'src') + os.pathsep + env.get('PYTHONPATH', '')

for script in [
    ROOT / "experiments/e5c5_activity_mapping_crosscalibration.py",
    ROOT / "experiments/e5c6_multiobservable_history_information.py",
]:
    subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, env=env)

print("Reproduced E5-C5 and E5-C6 into results/E5C5 and results/E5C6.")
