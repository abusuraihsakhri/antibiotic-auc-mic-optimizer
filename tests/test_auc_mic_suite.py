#!/usr/bin/env python3
"""
Test Suite for Antibiotic AUC/MIC Optimizer.
"""

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import unittest
from test_auc_mic import (
    TestRenalAndAnthropometrics,
    TestTrapezoidalAUC,
    TestVancomycinOptimizer,
    TestAminoglycosideOptimizer,
    TestBetaLactamOptimizer,
    TestMonteCarloPTA,
    TestCSVBatchAssessor,
)

if __name__ == "__main__":
    unittest.main()
