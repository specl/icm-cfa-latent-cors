#!/usr/bin/env bash
set -euo pipefail

# run render
R -q -f srp-render.R

# open the pdf (macOS)
open -a Preview "_out/p.pdf"