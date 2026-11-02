#!/usr/bin/env bash
# Start PlantXpresso by hand (for development/testing).
# The always-on public server is the systemd service (plantxpresso.service).
# Both can't use the same port at once: either stop the service first
#   sudo systemctl stop plantxpresso
# or test on another port through your SSH tunnel:
#   PLANTX_PORT=8000 ./run.sh --reload
set -euo pipefail
cd "$(dirname "$0")"

export PLANTX_DATA_DIR="${PLANTX_DATA_DIR:-/vol/data/datasets/PanXpresso}"
export PLANTX_METADATA_TSV="${PLANTX_METADATA_TSV:-/vol/data/datasets/PlantXpresso_metadata.tsv}"
export PLANTX_METADATA_WEBSERVER_TSV="${PLANTX_METADATA_WEBSERVER_TSV:-/vol/data/datasets/PlantXpresso_metadata_webserver.tsv}"
export PLANTX_LOGO="${PLANTX_LOGO:-/vol/data/datasets/PlantXpresso_logo.png}"
PORT="${PLANTX_PORT:-3838}"      # 3838 = the port the SimpleVM web service forwards to

# --reload restarts the app automatically when you edit a .py/.html file (development only)
exec ./venv/bin/uvicorn app.main:app --host 0.0.0.0 --port "$PORT" --forwarded-allow-ips='*' ${1:-}
