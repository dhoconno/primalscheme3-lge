"""Subprocess lifecycle checks for parallel observed-only discovery."""
import multiprocessing
import os
from pathlib import Path
import signal
import subprocess
import sys
import threading
import time
import unittest

import numpy as np

from primalscheme3.core.config import Config
from primalscheme3.core.parallel_discovery import pool_results


def _array(rows=3, columns=500):
    row = np.resize(np.array(list("ACGT")), columns)
    return np.stack([row.copy() for _ in range(rows)])


def _active_pids():
    return sorted(child.pid for child in multiprocessing.active_children())


def _fixture_consumer_exception():
    generator = pool_results(
        [[("f", 20)], [("r", 30)]],
        _array(),
        Config(terminal_gap_policy="observed-only", ncores=2),
        2,
    )
    try:
        next(generator)
        raise RuntimeError("synthetic consumer failure")
    except RuntimeError:
        pass
    finally:
        generator.close()
    print(f"ACTIVE_AFTER:{_active_pids()}", flush=True)


def _fixture_worker_exception():
    generator = pool_results(
        [[("invalid", 0)]],
        _array(),
        Config(terminal_gap_policy="observed-only", ncores=2),
        2,
    )
    try:
        list(generator)
    except ValueError as error:
        print(f"WORKER_ERROR:{error}", flush=True)
    finally:
        generator.close()
    print(f"ACTIVE_AFTER:{_active_pids()}", flush=True)


def _fixture_sigterm():
    # Repeated synthetic jobs keep next() blocked while the supervisor signals
    # the parent. A daemon reporter exposes the exact spawned worker PIDs.
    chunks = [[("f", 40)] * 4_000, [("r", 40)] * 4_000]
    generator = pool_results(
        chunks,
        _array(rows=24, columns=2_000),
        Config(terminal_gap_policy="observed-only", ncores=2),
        2,
    )

    def report_workers():
        deadline = time.monotonic() + 10
        while time.monotonic() < deadline:
            pids = _active_pids()
            if pids:
                print("PIDS:" + ",".join(map(str, pids)), flush=True)
                return
            time.sleep(0.01)
        print("PIDS:", flush=True)

    threading.Thread(target=report_workers, daemon=True).start()
    try:
        next(generator)
    except KeyboardInterrupt:
        generator.close()
        print(f"ACTIVE_AFTER:{_active_pids()}", flush=True)
        return
    finally:
        generator.close()
    raise AssertionError("synthetic discovery completed before SIGTERM")


class ParallelCancellationTests(unittest.TestCase):
    @staticmethod
    def _fixture_environment():
        root = str(Path(__file__).resolve().parents[2])
        environment = os.environ.copy()
        environment["PYTHONPATH"] = root + (os.pathsep + environment["PYTHONPATH"] if environment.get("PYTHONPATH") else "")
        return environment

    def _run_fixture(self, name, timeout=30):
        result = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "fixture", name],
            cwd=Path(__file__).resolve().parents[2],
            env=self._fixture_environment(),
            text=True,
            capture_output=True,
            timeout=timeout,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("ACTIVE_AFTER:[]", result.stdout)
        return result

    def test_consumer_exception_closes_and_joins_workers(self):
        self._run_fixture("consumer")

    def test_worker_exception_closes_and_joins_workers(self):
        result = self._run_fixture("worker")
        self.assertIn("WORKER_ERROR:Invalid discovery task direction", result.stdout)

    def test_parent_sigterm_terminates_and_joins_workers(self):
        process = subprocess.Popen(
            [sys.executable, str(Path(__file__).resolve()), "fixture", "sigterm"],
            cwd=Path(__file__).resolve().parents[2],
            env=self._fixture_environment(),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.assertIsNotNone(process.stdout)
        line = process.stdout.readline().strip()
        self.assertTrue(line.startswith("PIDS:"), line)
        worker_pids = [int(value) for value in line.removeprefix("PIDS:").split(",") if value]
        self.assertTrue(worker_pids, "fixture did not start worker processes")
        process.send_signal(signal.SIGTERM)
        stdout, stderr = process.communicate(timeout=30)
        self.assertEqual(process.returncode, 0, line + "\n" + stdout + stderr)
        self.assertIn("ACTIVE_AFTER:[]", stdout)
        deadline = time.monotonic() + 5
        while time.monotonic() < deadline and any(_pid_exists(pid) for pid in worker_pids):
            time.sleep(0.02)
        self.assertFalse([pid for pid in worker_pids if _pid_exists(pid)])


def _pid_exists(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "fixture":
        {"consumer": _fixture_consumer_exception,
         "worker": _fixture_worker_exception,
         "sigterm": _fixture_sigterm}[sys.argv[2]]()
    else:
        unittest.main()
