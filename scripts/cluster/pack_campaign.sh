#!/usr/bin/env bash
# Create a portable archive from immutable raw artifacts without deleting them.
set -euo pipefail
[[ $# -eq 1 ]] || { echo "Usage: $0 /path/to/campaign-directory" >&2; exit 2; }
campaign="${1%/}"
[[ -d "$campaign" ]] || { echo "not a directory: $campaign" >&2; exit 2; }
parent="$(cd "$campaign/.." && pwd)"; name="$(basename "$campaign")"
archive="$parent/$name.tar.gz"
[[ ! -e "$archive" ]] || { echo "refusing to overwrite: $archive" >&2; exit 2; }
(cd "$parent" && tar -czf "$archive" "$name")
if command -v sha256sum >/dev/null; then
    (cd "$parent" && sha256sum "$(basename "$archive")") > "$archive.sha256"
else
    (cd "$parent" && shasum -a 256 "$(basename "$archive")") > "$archive.sha256"
fi
echo "$archive"
