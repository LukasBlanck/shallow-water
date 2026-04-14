#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="${1:-.}"
MODE="${2:-check}"   # check | fix
FORMAT_FILE="$ROOT_DIR/configs/.clang-format"
REQUIRED_MAJOR="${CLANG_FORMAT_REQUIRED_MAJOR:-21}"

if ! command -v clang-format >/dev/null 2>&1; then
  echo "Error: clang-format is not installed or not in PATH." >&2
  exit 1
fi

if [[ ! -f "$FORMAT_FILE" ]]; then
  echo "Error: format file not found: $FORMAT_FILE" >&2
  exit 1
fi

VERSION_OUTPUT="$(clang-format --version)"
INSTALLED_MAJOR="$(
  printf '%s\n' "$VERSION_OUTPUT" |
  sed -nE 's/.*version[[:space:]]+([0-9]+)(\.[0-9]+.*)?/\1/p' |
  head -n1
)"

if [[ -z "$INSTALLED_MAJOR" ]]; then
  echo "Error: could not detect clang-format major version." >&2
  echo "clang-format --version output: $VERSION_OUTPUT" >&2
  exit 1
fi

if [[ "$INSTALLED_MAJOR" != "$REQUIRED_MAJOR" ]]; then
  echo "Error: clang-format major version $REQUIRED_MAJOR is required, but found:" >&2
  echo "  $VERSION_OUTPUT" >&2
  exit 1
fi

if [[ "$MODE" != "check" && "$MODE" != "fix" ]]; then
  echo "Error: mode must be 'check' or 'fix', got '$MODE'" >&2
  exit 1
fi

echo "Using format file: $FORMAT_FILE"
echo "$VERSION_OUTPUT"

SEARCH_DIRS=()
for dir in apps src include configs tests; do
  if [[ -d "$ROOT_DIR/$dir" ]]; then
    SEARCH_DIRS+=("$ROOT_DIR/$dir")
  fi
done

if [[ ${#SEARCH_DIRS[@]} -eq 0 ]]; then
  echo "No source directories found under $ROOT_DIR." >&2
  exit 0
fi

found=0
failed=0

while IFS= read -r -d '' file; do
  found=1
  if [[ "$MODE" == "fix" ]]; then
    echo "Formatting: $file"
    clang-format -i -style="file:$FORMAT_FILE" "$file"
  else
    echo "Checking: $file"
    if ! clang-format --dry-run --Werror -style="file:$FORMAT_FILE" "$file"; then
      failed=1
    fi
  fi
done < <(
  find "${SEARCH_DIRS[@]}" \
    -type f \( -name '*.cpp' -o -name '*.hpp' -o -name '*.h' \) \
    -print0
)

if [[ "$found" == "0" ]]; then
  echo "No matching source files found."
  exit 0
fi

if [[ "$MODE" == "check" && "$failed" == "1" ]]; then
  echo "Formatting check failed."
  exit 1
fi

echo "Done."