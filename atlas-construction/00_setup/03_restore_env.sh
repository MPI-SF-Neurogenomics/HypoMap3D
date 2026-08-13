#!/usr/bin/env bash

# This script should only be run whilst inside the development container (for maximum reproducibility)

## Restore python environment (via uv)
uv sync --frozen --all-extras