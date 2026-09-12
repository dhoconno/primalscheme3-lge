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
from types import SimpleNamespace
from unittest.mock import patch

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


def _fixture_startup_signal(signum, sigint_policy='default'):
    """Deliver through a preexisting unblocked thread before pool ownership."""
    context = multiprocessing.get_context('spawn')
    pool_ready = threading.Event()
    signal_sent = threading.Event()
    owned_pools = []
    constructor_returned = False
    handler_states = []
    if sigint_policy == 'custom':
        def prior_sigint_handler(signum, frame):
            handler_states.append(constructor_returned)
        signal.signal(signal.SIGINT, prior_sigint_handler)
    elif sigint_policy == 'ignore':
        signal.signal(signal.SIGINT, signal.SIG_IGN)
    previous = {sig: signal.getsignal(sig) for sig in (signal.SIGTERM, signal.SIGINT)}

    def interrupt_constructor():
        if not pool_ready.wait(10):
            return
        # raise_signal targets this unblocked thread, while Python dispatches
        # the installed handler on the main thread inside the constructor.
        signal.raise_signal(signum)
        signal_sent.set()

    reporter = threading.Thread(target=interrupt_constructor, daemon=True)
    reporter.start()

    def construct_pool(*args, **kwargs):
        nonlocal constructor_returned
        pool = context.Pool(*args, **kwargs)
        owned_pools.append(pool)
        pool_ready.set()
        if not signal_sent.wait(10):
            raise AssertionError('startup signal was not delivered')
        constructor_returned = True
        return pool

    generator = pool_results(
        [[('f', 40)]], _array(),
        Config(terminal_gap_policy='observed-only', ncores=2), 2,
    )
    try:
        with patch('primalscheme3.core.parallel_discovery.multiprocessing.get_context',
                   return_value=SimpleNamespace(Pool=construct_pool)):
            try:
                next(generator)
            except KeyboardInterrupt:
                print('INTERRUPTED', flush=True)
            finally:
                generator.close()
        print(f'CONSTRUCTOR_RETURNED:{constructor_returned}', flush=True)
        print(f'PRIOR_HANDLER_STATES:{handler_states}', flush=True)
        print(f'ACTIVE_AFTER:{_active_pids()}', flush=True)
        print('HANDLERS_RESTORED:' + str(all(signal.getsignal(sig) == handler
                                           for sig, handler in previous.items())), flush=True)
    finally:
        # A failing regression must not leave its test-owned pools behind.
        for pool in owned_pools:
            pool.terminate()
            pool.join()
        reporter.join(timeout=10)


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

    def test_startup_sigterm_waits_for_pool_ownership_before_cleanup(self):
        result = self._run_fixture('startup_sigterm')
        self.assertIn('INTERRUPTED', result.stdout)
        self.assertIn('CONSTRUCTOR_RETURNED:True', result.stdout)
        self.assertIn('HANDLERS_RESTORED:True', result.stdout)

    def test_startup_sigint_waits_for_pool_ownership_before_cleanup(self):
        result = self._run_fixture('startup_sigint')
        self.assertIn('INTERRUPTED', result.stdout)
        self.assertIn('CONSTRUCTOR_RETURNED:True', result.stdout)
        self.assertIn('HANDLERS_RESTORED:True', result.stdout)

    def test_custom_sigint_handler_runs_after_ownership_and_is_restored(self):
        result = self._run_fixture('startup_custom_sigint')
        self.assertIn('PRIOR_HANDLER_STATES:[True]', result.stdout)
        self.assertIn('HANDLERS_RESTORED:True', result.stdout)

    def test_ignored_sigint_remains_ignored_and_is_restored(self):
        result = self._run_fixture('startup_ignored_sigint')
        self.assertNotIn('INTERRUPTED', result.stdout)
        self.assertIn('CONSTRUCTOR_RETURNED:True', result.stdout)
        self.assertIn('HANDLERS_RESTORED:True', result.stdout)

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
         "sigterm": _fixture_sigterm,
         "startup_sigterm": lambda: _fixture_startup_signal(signal.SIGTERM),
         "startup_sigint": lambda: _fixture_startup_signal(signal.SIGINT),
         "startup_custom_sigint": lambda: _fixture_startup_signal(signal.SIGINT, 'custom'),
         "startup_ignored_sigint": lambda: _fixture_startup_signal(signal.SIGINT, 'ignore')}[sys.argv[2]]()
    else:
        unittest.main()
