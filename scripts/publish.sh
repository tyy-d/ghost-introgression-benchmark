#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
expected=$(pwd -P)
actual=$(git rev-parse --show-toplevel)
if [ "$actual" != "$expected" ]; then
    echo "Refusing to operate outside the project's own repository" >&2
    exit 1
fi
# Staging and committing are deliberate separate steps; this script publishes
# only existing commits, never course materials or untracked output directories.
if gh repo view tyy-d/ghost-introgression-benchmark >/dev/null 2>&1; then
    if ! git remote get-url origin >/dev/null 2>&1; then
        git remote add origin https://github.com/tyy-d/ghost-introgression-benchmark.git
    fi
    git push --set-upstream origin main
else
    gh repo create tyy-d/ghost-introgression-benchmark --public --source . --remote origin --push
fi
local_head=$(git rev-parse HEAD)
remote_head=$(git ls-remote origin refs/heads/main | cut -f1)
if [ "$local_head" != "$remote_head" ]; then
    echo "Remote HEAD does not match the tested local commit" >&2
    exit 1
fi
echo "Verified remote commit: $remote_head"
