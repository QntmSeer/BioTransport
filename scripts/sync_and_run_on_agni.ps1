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
scp -r 01_microscale_md 02_continuum_transport data notebooks Makefile pyproject.toml README.md "${RemoteUser}@${RemoteHost}:${RemoteDir}/"

# 3. Execute build and test suite on Agni
Write-Host "[3/3] Executing Rust engine build & simulation on Agni..." -ForegroundColor Yellow
ssh "${RemoteUser}@${RemoteHost}" "bash -l -c '
    cd ${RemoteDir}/02_continuum_transport/biotransport-rs && \
    cargo build --release && \
    cargo test && \
    cd ${RemoteDir} && \
    ./02_continuum_transport/biotransport-rs/target/release/biotransport-cli run \
        --params data/sample_md_params.json \
        --tmp 200000 \
        --time 3600 \
        --out data/remote_simulation_results.json
'"

Write-Host "=================================================================" -ForegroundColor Cyan
Write-Host "  Compute successfully completed on agni!                        " -ForegroundColor Green
Write-Host "=================================================================" -ForegroundColor Cyan
