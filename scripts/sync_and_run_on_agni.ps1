# PowerShell Remote Execution Script for Agni (omarchy)
# Usage: .\scripts\sync_and_run_on_agni.ps1

$ErrorActionPreference = "Stop"

$RemoteHost = "192.168.1.112"
$RemoteUser = "agni"
$RemoteDir  = "~/multiscale-bioparticle-transport"

Write-Host "=================================================================" -ForegroundColor Cyan
Write-Host "  Dispatching Multiscale Bioparticle Transport Compute to Agni  " -ForegroundColor Green
Write-Host "  Target: ${RemoteUser}@${RemoteHost}:${RemoteDir}               " -ForegroundColor Cyan
Write-Host "=================================================================" -ForegroundColor Cyan

# 1. Create directory on remote
Write-Host "[1/3] Ensuring remote directory structure..." -ForegroundColor Yellow
ssh "${RemoteUser}@${RemoteHost}" "mkdir -p ${RemoteDir}"

# 2. Copy files via scp
Write-Host "[2/3] Transferring project files..." -ForegroundColor Yellow
scp -r 01_microscale_md 02_continuum_transport scripts data notebooks Makefile pyproject.toml README.md "${RemoteUser}@${RemoteHost}:${RemoteDir}/"

# 3. Execute build and test suite on Agni
Write-Host "[3/3] Executing Rust engine build & full master multiscale pipeline on Agni..." -ForegroundColor Yellow
ssh "${RemoteUser}@${RemoteHost}" "bash -l -c '
    export PATH=\"\$HOME/.cargo/bin:\$PATH\"
    [ -f \"\$HOME/.cargo/env\" ] && source \"\$HOME/.cargo/env\"
    cd ${RemoteDir}/02_continuum_transport/biotransport-rs && \
    cargo build --release && \
    cargo test && \
    cd ${RemoteDir} && \
    [ -d \".venv\" ] && source .venv/bin/activate
    python3 scripts/master_run_pipeline.py
'"

Write-Host "=================================================================" -ForegroundColor Cyan
Write-Host "  Full Multiscale Pipeline successfully executed on Agni!        " -ForegroundColor Green
Write-Host "=================================================================" -ForegroundColor Cyan
