"""Run the replay commands of ci/replays.json.

Each entry comes from one of the per-version workflows that the single CI workflow replaced, and it
keeps that workflow's file filter. An entry runs when a changed file matches one of its paths
(GitHub path-filter globs: `*` stays inside a directory, `**` crosses directories), or always with
--all. Unit tests are not listed here: the CI tests job runs the whole suite.

    python ci/run_replays.py --kind python --all
    python ci/run_replays.py --kind sage --changed-files changed.txt
    python ci/run_replays.py --kind python --list --changed-files changed.txt
"""
import argparse
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def glob_regex(pattern):
    out, i = "", 0
    while i < len(pattern):
        c = pattern[i]
        if pattern.startswith("**", i):
            out += ".*"
            i += 2
            if i < len(pattern) and pattern[i] == "/":
                out = out[:-2]+"(?:.*/)?"
                i += 1
            continue
        if c == "*":
            out += "[^/]*"
        elif c == "?":
            out += "[^/]"
        elif c == "[":
            j = pattern.index("]", i)
            out += pattern[i:j+1]
            i = j
        else:
            out += re.escape(c)
        i += 1
    return re.compile(out+r"\Z")


def selected(entries, changed, run_all):
    if run_all:
        return entries
    regs = {p: glob_regex(p) for e in entries for p in e["paths"]}
    return [e for e in entries if not e["paths"] or any(regs[p].match(f) for p in e["paths"] for f in changed)]


def run(entry):
    failures = []
    for cmd in entry["commands"]:
        env = dict(os.environ, **{k: str(v) for k, v in cmd.get("env", {}).items()})
        cwd = ROOT/cmd.get("working_directory", ".")
        print(f"::group::{entry['source']}: {cmd['run'].splitlines()[0][:100]}", flush=True)
        start = time.time()
        try:
            proc = subprocess.run(["bash", "-eo", "pipefail", "-c", cmd["run"]], cwd=cwd, env=env,
                                  timeout=60*entry["timeout_minutes"])
            ok = proc.returncode == 0
        except subprocess.TimeoutExpired:
            ok = False
        print(f"::endgroup::{'ok' if ok else 'FAILED'} in {time.time()-start:.0f}s", flush=True)
        if not ok:
            failures.append(cmd["run"])
    return failures


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", choices=["python", "sage"], required=True)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--changed-files", help="file with one changed path per line")
    ap.add_argument("--list", action="store_true", help="print the selection without running it")
    args = ap.parse_args()
    entries = [e for e in json.load(open(ROOT/"ci"/"replays.json"))["entries"] if e["kind"] == args.kind]
    changed = []
    if args.changed_files:
        changed = [l.strip() for l in open(args.changed_files) if l.strip()]
    elif not args.all:
        ap.error("give --all or --changed-files")
    chosen = selected(entries, changed, args.all)
    print(f"{len(chosen)} of {len(entries)} {args.kind} replay entries selected")
    for e in chosen:
        print("  ", e["source"])
    if args.list:
        return 0
    failed = {}
    for e in chosen:
        f = run(e)
        if f:
            failed[e["source"]] = f
    if failed:
        print("FAILED:", json.dumps(failed, indent=1))
        return 1
    print("all selected replays passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
