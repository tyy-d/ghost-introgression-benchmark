#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p vendor .tools .cache
export PIP_CACHE_DIR="$PWD/.cache/pip"
export CONDA_PKGS_DIRS="$PWD/.cache/conda"

# Prefer a Conda environment because the original ArchIE model is fitted in R.
conda env create --prefix "$PWD/.tools/environment" --file environment.yml
PYTHON="$PWD/.tools/environment/bin/python"
PIP="$PWD/.tools/environment/bin/pip"
export PATH="$PWD/.tools/environment/bin:$PATH"

if [ ! -d vendor/ArchIE/.git ]; then
    git clone https://github.com/sriramlab/ArchIE vendor/ArchIE
fi
if [ ! -d vendor/trace/.git ]; then
    git clone https://github.com/YulinZhang9806/trace vendor/trace
fi
if [ ! -d vendor/SINGER-source/.git ]; then
    git clone --depth 1 --filter=blob:none --sparse https://github.com/popgenmethods/SINGER vendor/SINGER-source
fi
git -C vendor/ArchIE checkout 00820bc61beca8a958f04fb17291ff66e03b9e94
git -C vendor/trace checkout 70cab6c6687900826eb0b2ded69c39554b1b4146
git -C vendor/SINGER-source fetch --depth 1 origin eb8e39b1a15be4a9a4df4fdaab61847bf73515d7
git -C vendor/SINGER-source sparse-checkout set --no-cone \
    '/README.md' '/SINGER/SINGER/*.cpp' '/SINGER/SINGER/*.hpp' \
    '/SINGER/SINGER/singer_master' '/SINGER/SINGER/convert_to_tskit'
git -C vendor/SINGER-source checkout eb8e39b1a15be4a9a4df4fdaab61847bf73515d7
"$PIP" install ./vendor/trace

mkdir -p .tools/singer
"${CXX:-c++}" -std=c++17 -O3 vendor/SINGER-source/SINGER/SINGER/*.cpp -o .tools/singer/singer
cp vendor/SINGER-source/SINGER/SINGER/singer_master .tools/singer/
cp vendor/SINGER-source/SINGER/SINGER/convert_to_tskit .tools/singer/
chmod +x .tools/singer/singer_master .tools/singer/convert_to_tskit
"$PYTHON" -m pip freeze > .tools/python-freeze.txt
echo "Environment ready. Activate it with: conda activate $PWD/.tools/environment"
echo "Run: ghostbench smoke --config configs/smoke.yaml --output results/smoke"
