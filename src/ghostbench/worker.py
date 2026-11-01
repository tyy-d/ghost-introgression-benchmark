"""Private worker entry point for the bounded local process queue."""

import json
import os
import signal
import sys


def main():
    from .pipeline import run

    def cancelled(signum, frame):
        raise KeyboardInterrupt("Worker cancelled")

    signal.signal(signal.SIGTERM, cancelled)
    spec = json.loads(open(sys.argv[1]).read())
    os.environ["GHOSTBENCH_TIMEOUT"] = str(spec["timeout"])
    run(spec["config"], spec["output"], spec["model"], thresholds=spec["thresholds"])


if __name__ == "__main__":
    main()
