#!/usr/bin/env bash
# Choose which replays to run: all of them on a schedule or manual run (or when the base commit is
# unknown), otherwise those whose files changed since BASE. Writes the runner arguments to
# $GITHUB_OUTPUT as `args`.
set -euo pipefail
if [[ "$EVENT" == "schedule" || "$EVENT" == "workflow_dispatch" || -z "${BASE:-}" || "$BASE" =~ ^0+$ ]] \
   || ! git cat-file -e "$BASE^{commit}" 2>/dev/null; then
  echo "args=--all" >> "$GITHUB_OUTPUT"
else
  git diff --name-only "$BASE" HEAD > "$RUNNER_TEMP/changed.txt"
  echo "changed files:"; cat "$RUNNER_TEMP/changed.txt"
  echo "args=--changed-files $RUNNER_TEMP/changed.txt" >> "$GITHUB_OUTPUT"
fi
