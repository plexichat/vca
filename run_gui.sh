#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")"
exec python3 gui/vector_router_gui.py "$@"