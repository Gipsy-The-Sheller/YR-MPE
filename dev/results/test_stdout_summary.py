"""Tests for the machine-readable terminal-analysis summaries."""

from __future__ import annotations

import io
import json
import sys
import unittest
from contextlib import redirect_stdout
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from YR_MPE.plugins.components.stdout_summary import (  # noqa: E402
    SUMMARY_PREFIX, emit_summary,
)


class StdoutSummaryTests(unittest.TestCase):
    def test_line_is_prefixed_json(self):
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            line = emit_summary("codon_ztest", pairs=3, significant=1)
        output = buffer.getvalue().strip()
        self.assertTrue(output.startswith(SUMMARY_PREFIX))
        payload = json.loads(output[len(SUMMARY_PREFIX):])
        self.assertEqual(payload["schema"], "yrmpe.summary.v1")
        self.assertEqual(payload["tool"], "codon_ztest")
        self.assertEqual(payload["pairs"], 3)
        self.assertEqual(json.loads(line), payload)


if __name__ == "__main__":
    unittest.main(verbosity=2)
