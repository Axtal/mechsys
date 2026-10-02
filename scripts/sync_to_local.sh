#!/bin/bash
# sync_to_local.sh – watch for DONE files and copy results to local PC
# Runs on the LOCAL PC (uses NFS-mounted HPC storage directly — no SSH needed)

# ===== ADAPT THIS =====
SHARED_BASE="/home/msys/HPC/storage/folder_watcher/shearbox_pressure_steady/shearbox_pressure_sweep"
LOCAL_DEST="/mnt/hdd_data/HPC_data/shearbox_pressure_steady"
SLEEP_SECONDS=10
# ======================

echo "Watcher started. Using base: $SHARED_BASE"
echo "Destination: $LOCAL_DEST"

while true; do
    for simdir in "$SHARED_BASE"/shearbox-*; do
        [ -d "$simdir" ] || continue
        if [ -f "$simdir/DONE" ] && [ ! -f "$simdir/COPIED" ]; then
            safename="$(basename "$simdir")"
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying $simdir to local PC..."

            # Capture rsync stdout+stderr to a temp file for debugging
            RSYNC_LOG="/tmp/rsync_${safename}.log"
            rsync -avP "$simdir/" "$LOCAL_DEST/$safename/" > "$RSYNC_LOG" 2>&1
            RSYNC_EXIT=$?

            if [ $RSYNC_EXIT -eq 0 ]; then
                touch "$simdir/COPIED"
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] SUCCESS: $safename copied."
                echo "Deleting $simdir ..."
                rm -rf "$simdir"
                rm -f "$RSYNC_LOG"
                echo "Done."
            else
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: rsync failed for $safename (exit code: $RSYNC_EXIT)"
                echo "--- rsync output (last 30 lines) ---"
                tail -30 "$RSYNC_LOG"
                echo "--- end rsync output ---"
                echo "Full log saved to: $RSYNC_LOG"
                echo "Will retry later."
            fi
        fi
    done
    sleep "$SLEEP_SECONDS"
done