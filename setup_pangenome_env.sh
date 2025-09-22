#!/bin/bash

set -euo pipefail

# CONFIG
VENV_DIR="${1:-$HOME/venvs/pangenome-env}"   # Use first argument or fallback
HTSLIB_VERSION="1.19"
HTSLIB_DIR="$HOME/htslib-$HTSLIB_VERSION"
HTSLIB_PREFIX="$HOME/htslib-install"

# STEP 1: Create virtualenv
echo "Creating virtual environment in: $VENV_DIR"
python3 -m venv "$VENV_DIR"
source "$VENV_DIR/bin/activate"

# STEP 2: Upgrade pip and install biopython
echo "Installing biopython..."
pip install --upgrade pip
pip install biopython

# STEP 3: Install htslib from source (for bgzip)
echo "Installing HTSlib $HTSLIB_VERSION..."
cd "$HOME"
curl -LO "https://github.com/samtools/htslib/releases/download/$HTSLIB_VERSION/htslib-$HTSLIB_VERSION.tar.bz2"
tar -xjf "htslib-$HTSLIB_VERSION.tar.bz2"
cd "$HTSLIB_DIR"
./configure --prefix="$HTSLIB_PREFIX"
make -j4
make install

# STEP 4: Add bgzip to path (per session)
echo "Setup complete."
echo "To activate the environment in future, run:"
echo "   source $VENV_DIR/bin/activate"
echo "   export PATH=$HTSLIB_PREFIX/bin:\$PATH"
