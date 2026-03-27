#!/bin/bash

BASE_DIR="/alice/sim/2026/LHC26b8/"
DEST_BASE="/mnt/hdd3TB/UPCcandMuon/MC/sim/2026/LHC26b8/"
FILE_PATTERN="*AO2D*.root"
MAX_JOBS=6

export ERROR_LOG="${DEST_BASE}/failed_downloads_$(date +%Y%m%d_%H%M%S).txt"

download_file() {
  local SRC="$1"
  local BASE_DIR="$2"
  local DEST_BASE="$3"

  if [[ ! "$SRC" == /* ]]; then
    return 0
  fi

  local RELATIVE_PATH="${SRC#$BASE_DIR}"
  local DEST_FILE="${DEST_BASE}${RELATIVE_PATH}"
  local DEST_DIR=$(dirname "$DEST_FILE")

  if [ -s "$DEST_FILE" ]; then
    echo "-> Already exists. Skipping: $DEST_FILE"
    return 0
  fi

  mkdir -p "$DEST_DIR"
  
  echo "-> Copying $SRC ..."
  
  local n=0
  until [ $n -ge 3 ]; do
    if alien_cp "$SRC" "file://${DEST_FILE}"; then
      break
    fi
    n=$((n+1))
    echo "   Failed. Retrying ($n/3)... $SRC"
    sleep 10
  done
  
# 3 times retry
  if [ $n -ge 3 ]; then
    echo "   [ERROR] Failed after 3 attempts: $SRC"
    echo "$SRC" >> "$ERROR_LOG"
  fi
}

export -f download_file

echo ">>> Searching for ${FILE_PATTERN} recursively in ${BASE_DIR} ..."
echo ">>> Error log will be saved to: ${ERROR_LOG}"

mkdir -p "$DEST_BASE"
touch "$ERROR_LOG"

alien_find "$BASE_DIR" "$FILE_PATTERN" | \
  xargs -I {} -P "$MAX_JOBS" bash -c 'download_file "$@"' _ "{}" "$BASE_DIR" "$DEST_BASE"

echo ">>> All processes completed."
echo ">>> Please check ${ERROR_LOG} for any failed downloads."