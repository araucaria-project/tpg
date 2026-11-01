#!/usr/bin/env bash
# Replace ~/.config/tpg/local_config.yaml with the OCM default from this repo (local_config.yaml).
# The old file is kept as local_config.yaml.<timestamp>.bak and the differences are printed.
set -euo pipefail

src="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/local_config.yaml"
dst="$HOME/.config/tpg/local_config.yaml"

mkdir -p "$(dirname "$dst")"
if [ -f "$dst" ]; then
    if cmp -s "$src" "$dst"; then
        echo "$dst is already up to date"
        exit 0
    fi
    bak="$dst.$(date +%Y%m%d_%H%M%S).bak"
    cp -p "$dst" "$bak"
    echo "old config saved as $bak; differences (old -> new):"
    diff "$bak" "$src" || true
fi
cp "$src" "$dst"
echo "$dst reset from $src"
