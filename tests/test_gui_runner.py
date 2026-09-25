from __future__ import annotations

import unittest

from PySide6.QtCore import QCoreApplication

from rna_mod_pipeline.gui.runner import ProcessRunner
from rna_mod_pipeline.gui.state import RunState


class ProcessRunnerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.application = QCoreApplication.instance() or QCoreApplication([])

    def test_finalizing_cancelled_run_disarms_force_stop_timer(self) -> None:
        runner = ProcessRunner()
        runner._force_stop_timer.start(5_000)
        self.assertTrue(runner._force_stop_timer.isActive())

        runner._finalize(RunState.INTERRUPTED, -1)

        self.assertFalse(runner._force_stop_timer.isActive())


if __name__ == "__main__":
    unittest.main()
