"""Run two disjoint, single-thread workers and finalize local study artifacts."""
from datetime import datetime, timezone
from pathlib import Path
import fcntl
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

STUDY = Path(__file__).resolve().parent
ROOT = STUDY.parent.parent


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def main():
    lock = (STUDY / "supervisor.lock").open("a")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise SystemExit("Another supervisor already owns this study.")
    environment = os.environ.copy()
    environment.update({name: "1" for name in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
        "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS")})
    environment.update(R_PROFILE_USER="/dev/null", R_ENVIRON_USER="/dev/null")
    rscript = shutil.which("Rscript")
    if not rscript:
        raise SystemExit("Rscript is unavailable.")
    state = dict(status="starting", started_at=timestamp(), supervisor_pid=os.getpid(),
                 worker_count=2, cpu_threads_per_worker=1, planned_datasets=90,
                 workers=[], phase="fitting")
    files = list(STUDY.glob("*.R")) + list(STUDY.glob("*.Rmd")) + [Path(__file__).resolve()]
    files += [STUDY / "design.json", STUDY / "manifest.csv", STUDY / "source/MPCurver_0.3.2.tar.gz"]
    hashes = {str(p.relative_to(STUDY)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    (STUDY / "execution_checksums.json").write_text(json.dumps(hashes, indent=2) + "\n")

    def save():
        state["updated_at"] = timestamp()
        state["completed_method_outcomes"] = len(list((STUDY / "results").glob("main_*.rds")))
        temporary = STUDY / "run_status.json.tmp"
        temporary.write_text(json.dumps(state, indent=2) + "\n")
        temporary.replace(STUDY / "run_status.json")

    processes = []
    logs = []
    try:
        for worker in (1, 2):
            log = (STUDY / f"worker_{worker}.log").open("a", buffering=1)
            logs.append(log)
            command = [rscript, "--vanilla", str(STUDY / "run_study.R"),
                       "--phase=main", f"--worker={worker}", "--workers=2"]
            process = subprocess.Popen(command, cwd=ROOT, env=environment,
                                       stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT)
            processes.append(process)
            state["workers"].append(dict(worker=worker, pid=process.pid, exit_code=None))
        state["status"] = "running"
        save()
        while any(p.poll() is None for p in processes):
            for entry, process in zip(state["workers"], processes):
                entry["exit_code"] = process.poll()
            save()
            time.sleep(30)
        for entry, process in zip(state["workers"], processes):
            entry["exit_code"] = process.returncode
        if any(p.returncode != 0 for p in processes):
            raise RuntimeError("A worker failed; inspect worker and dataset logs before resuming.")
        state.update(status="finalizing", phase="validation_and_report")
        save()
        with (STUDY / "finalize.log").open("a", buffering=1) as log:
            for script, arguments in [("validate_results.R", []), ("summarize.R", ["main"]),
                                      ("render_examples.R", []), ("render_report.R", [])]:
                subprocess.run([rscript, "--vanilla", str(STUDY / script), *arguments],
                               cwd=ROOT, env=environment, stdin=subprocess.DEVNULL,
                               stdout=log, stderr=subprocess.STDOUT, check=True)
        state.update(status="completed", phase="completed", finished_at=timestamp())
    except Exception as error:
        state.update(status="failed", error=str(error), finished_at=timestamp())
        save()
        raise
    finally:
        save()
        for log in logs:
            log.close()


if __name__ == "__main__":
    main()
