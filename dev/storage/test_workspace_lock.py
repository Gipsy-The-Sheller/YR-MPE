"""Tests for the single-writer workspace lock."""

from __future__ import annotations

import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from YR_MPE.platforms.methods.workspace_lock import (  # noqa: E402
    LOCK_NAME, WorkspaceLock,
)


class WorkspaceLockTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix="yrmpe-lock-test-"))

    def tearDown(self):
        for item in self.root.glob(LOCK_NAME):
            try:
                item.unlink()
            except OSError:
                pass
        try:
            self.root.rmdir()
        except OSError:
            pass

    def test_one_writer_at_a_time(self):
        first = WorkspaceLock()
        second = WorkspaceLock()
        self.assertEqual(first.acquire(str(self.root))[0], True)
        ok, message = second.acquire(str(self.root))
        self.assertFalse(ok)
        self.assertIn("already open", message)
        first.release()
        self.assertEqual(second.acquire(str(self.root))[0], True)
        second.release()

    def test_same_lock_releases_and_reacquires(self):
        lock = WorkspaceLock()
        self.assertTrue(lock.acquire(str(self.root))[0])
        lock.release()
        self.assertFalse((self.root / LOCK_NAME).exists())
        self.assertTrue(lock.acquire(str(self.root))[0])
        lock.release()

    def test_stale_owner_is_recovered(self):
        marker = self.root / LOCK_NAME
        marker.write_text(json.dumps({
            "pid": 99999999,
            "host": "stale-host",
            "token": "stale",
        }), encoding="utf-8")
        lock = WorkspaceLock()
        self.assertTrue(lock.acquire(str(self.root))[0])
        lock.release()

    def test_release_does_not_remove_another_owners_lock(self):
        marker = self.root / LOCK_NAME
        marker.write_text(json.dumps({
            "pid": os.getpid(),
            "host": "other",
            "token": "another-instance",
        }), encoding="utf-8")
        lock = WorkspaceLock()
        lock.release()
        self.assertTrue(marker.exists())
        marker.unlink()


if __name__ == "__main__":
    unittest.main(verbosity=2)
