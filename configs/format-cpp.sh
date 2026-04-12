#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="${1:-.}"
FORMAT_FILE="$ROOT_DIR/configs/.clang-format"
DRY_RUN="${DRY_RUN:-0}"

if ! command -v clang-format >/dev/null 2>&1; then
  echo "Error: clang-format is not installed or not in PATH." >&2
  exit 1
fi

if [[ ! -f "$FORMAT_FILE" ]]; then
  echo "Error: format file not found: $FORMAT_FILE" >&2
  exit 1
fi

echo "Scanning: $ROOT_DIR"
echo "Using format file: $FORMAT_FILE"

found=0

while IFS= read -r -d '' file; do
  found=1
  if [[ "$DRY_RUN" == "1" ]]; then
    echo "Would format: $file"
  else
    echo "Formatting: $file"
    clang-format -i -style="file:$FORMAT_FILE" "$file"
  fi
done < <(
  find "$ROOT_DIR" \
    \( -type d \( -name .git -o -name build -o -name out -o -name cmake-build-debug -o -name cmake-build-release \) -prune \) -o \
    \( -type f \( -name '*.cpp' -o -name '*.hpp' \) ! -name 'toml.hpp' -print0 \)
)

if [[ "$found" == "0" ]]; then
  echo "No matching .cpp or .hpp files found."
else
  echo "Done."
fi