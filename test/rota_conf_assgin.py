"""
rota_conf_assgin.py
-------------------
Alias script for test/rota_conf_assign.py.
"""

from pathlib import Path
import sys

TEST_DIR = Path(__file__).resolve().parent
if str(TEST_DIR) not in sys.path:
    sys.path.insert(0, str(TEST_DIR))

from rota_conf_assign import run_test

if __name__ == "__main__":
    run_test()
