#!/usr/bin/env python3
"""
run_all_tests.py
------------------
Master test runner for the rota_assign package.
Executes both end-to-end torsion calculation and rotamer/RNA conformer assignment tests.
All output datasets are generated under test/test_results/.
"""

import sys
import subprocess
from pathlib import Path

TEST_DIR = Path(__file__).resolve().parent
RESULTS_DIR = TEST_DIR / "test_results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def main():
    print("==========================================================================")
    print("           rota_assign Package Test Suite Execution                       ")
    print("==========================================================================")
    print(f"Output Directory: {RESULTS_DIR}\n")

    # 1. Run test_torsion_calculation.py
    print(">>> [Test 1/3] Running Complex Torsion, DSSP, and SASA Calculation Test...")
    cmd1 = [sys.executable, str(TEST_DIR / "test_torsion_calculation.py")]
    res1 = subprocess.run(cmd1, capture_output=False)

    if res1.returncode != 0:
        print("\n❌ Test 1 (test_torsion_calculation.py) failed!")
        sys.exit(res1.returncode)

    print("\n--------------------------------------------------------------------------")

    # 2. Run rota_conf_assign.py
    print(">>> [Test 2/3] Running Rotamer & RNA Conformer Assignment Test...")
    cmd2 = [sys.executable, str(TEST_DIR / "rota_conf_assign.py")]
    res2 = subprocess.run(cmd2, capture_output=False)

    if res2.returncode != 0:
        print("\n❌ Test 2 (rota_conf_assign.py) failed!")
        sys.exit(res2.returncode)

    print("\n--------------------------------------------------------------------------")

    # 3. Run test_regressions.py (unit tests for fixed defects)
    print(">>> [Test 3/3] Running Regression Unit Tests...")
    cmd3 = [sys.executable, str(TEST_DIR / "test_regressions.py")]
    res3 = subprocess.run(cmd3, capture_output=False)

    if res3.returncode != 0:
        print("\n❌ Test 3 (test_regressions.py) failed!")
        sys.exit(res3.returncode)

    print("\n==========================================================================")
    print("  ALL ROTA_ASSIGN UTILITY TESTS PASSED SUCCESSFULLY!                      ")
    print(f"  All outputs generated in: {RESULTS_DIR}")
    print("==========================================================================")


if __name__ == "__main__":
    main()
