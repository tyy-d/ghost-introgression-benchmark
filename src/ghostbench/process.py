"""Process-group cancellation and process-tree resource sampling."""

import os
import signal
import subprocess
import time


def execute(command, log, cwd=None, timeout=None):
    import psutil

    from .io import terminate_group, write_json

    timeout = timeout if timeout is not None else float(os.environ.get("GHOSTBENCH_TIMEOUT", "14400"))
    environment = os.environ | {key: "1" for key in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]}
    start = time.monotonic()
    peak = 0
    cpu_by_pid = {}
    sampling_available = False
    status = "failed"
    cleanup_error = None
    old_handler = None

    def interrupted(signum, frame):
        raise KeyboardInterrupt("Process terminated")

    try:
        old_handler = signal.signal(signal.SIGTERM, interrupted)
    except ValueError:
        pass
    with log.open("w") as stream:
        process = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT,
                                   cwd=cwd, env=environment, start_new_session=True)
        try:
            while process.poll() is None:
                try:
                    root = psutil.Process(process.pid)
                    children = root.children(recursive=True)
                    rss = 0
                    for child in [root, *children]:
                        try:
                            rss += child.memory_info().rss
                            times = child.cpu_times()
                            cpu_by_pid[(child.pid, child.create_time())] = times.user + times.system
                            sampling_available = True
                        except (psutil.NoSuchProcess, psutil.AccessDenied, PermissionError):
                            continue
                    peak = max(peak, rss)
                except (psutil.NoSuchProcess, psutil.AccessDenied, PermissionError):
                    pass
                if time.monotonic() - start > timeout:
                    status = "timeout"
                    raise TimeoutError(f"Job exceeded {timeout:g} seconds: {log}")
                time.sleep(0.1)
            status = "complete" if process.returncode == 0 else "failed"
        except BaseException:
            if status != "timeout":
                status = "interrupted"
            try:
                terminate_group(process)
            except OSError as error:
                # Preserve cancellation/timeout as the primary outcome. A denied
                # signal must remain visible for an explicit process audit.
                cleanup_error = f"{type(error).__name__}: {error}"
            raise
        finally:
            telemetry = {"command": [str(x) for x in command], "cwd": str(cwd) if cwd else None,
                         "seconds": time.monotonic()-start, "returncode": process.poll(),
                         "status": status, "cpu_seconds_sampled": sum(cpu_by_pid.values()) if sampling_available else None,
                         "peak_tree_rss_bytes_sampled": peak if sampling_available else None, "sampling_seconds": 0.1,
                         "sampling_available": sampling_available,
                         "cleanup_error": cleanup_error,
                         "resource_note": "Sampled process-tree lower bounds; very short-lived processes may be missed"}
            write_json(log.with_suffix(".command.json"), telemetry)
            if old_handler is not None:
                signal.signal(signal.SIGTERM, old_handler)
    return telemetry["seconds"], process.returncode, telemetry
