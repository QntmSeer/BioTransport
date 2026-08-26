#!/usr/bin/env bash
# Script to synchronize the codebase and dispatch high-performance GROMACS/Rust compute on workstation (agni@omarchy).

set -e

REMOTE_HOST="192.168.1.112"
REMOTE_USER="agni"
REMOTE_DIR="~/multiscale-bioparticle-transport"

echo "================================================================="
echo "  Dispatching Multiscale Bioparticle Transport Compute to Agni  "
echo "  Target: ${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_DIR}            "
echo "================================================================="

# 1. Synchronize repository files (excluding local build artifacts)
echo "[1/3] Synchronizing repository files..."
rsync -avz --delete \
    --exclude '.git' \
    --exclude 'target' \
    --exclude '__pycache__' \
    --exclude '*.egg-info' \
    --exclude '.pytest_cache' \
    --exclude '*.xtc' \
    --exclude '*.trr' \
    ./ ${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_DIR}/

# 2. Execute Rust build & test suite on workstation
echo "[2/3] Building and running Rust continuum transport engine..."
ssh ${REMOTE_USER}@${REMOTE_HOST} "bash -l -c '
    export PATH=\"\$HOME/.cargo/bin:\$PATH\"
    [ -f \"\$HOME/.cargo/env\" ] && source \"\$HOME/.cargo/env\"
    cd ${REMOTE_DIR}/02_continuum_transport/biotransport-rs && \
    cargo build --release && \
    cargo test
'"

# 3. Execute sample continuum parameter simulation
echo "[3/3] Running high-speed simulation on Agni..."
ssh ${REMOTE_USER}@${REMOTE_HOST} "bash -l -c '
    export PATH=\"\$HOME/.cargo/bin:\$PATH\"
    [ -f \"\$HOME/.cargo/env\" ] && source \"\$HOME/.cargo/env\"
    cd ${REMOTE_DIR} && \
    ./02_continuum_transport/biotransport-rs/target/release/biotransport-cli run \
        --params data/sample_md_params.json \
        --tmp 200000 \
        --time 3600 \
        --out data/remote_simulation_results.json
'"

echo "================================================================="
echo "  Compute successfully completed on agni!                        "
echo "================================================================="
