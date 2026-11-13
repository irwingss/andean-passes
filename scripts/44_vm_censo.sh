#!/bin/bash
# 44_vm_censo.sh - Startup script of the census machine created by 44_censo_nube.sh (PRD O2.5, phase R3). Not run by
# hand. It installs the environment, reads FABDEM from Earth Engine (39), runs 40 and 41 on every 1-degree block in
# parallel, copies each block's results to the bucket as soon as it ends, and deletes the machine when all are done.
# Guards: each block has 30 minutes for 40 and 10 for 41; a block that fails or runs out of time is marked and the
# rest go on. The machine itself is deleted by Compute Engine at its maximum run time whatever happens.
set -u
M() { curl -sf -H Metadata-Flavor:Google "http://metadata.google.internal/computeMetadata/v1/instance/attributes/$1"; }
export B=$(M bucket) EP=$(M ee_project) PAR=$(M par) ZONE=$(M zone_name) NAME=$(hostname)
log() { echo "$(date -u +%FT%TZ) $*" >> /var/log/censo.log; gsutil -q cp /var/log/censo.log "$B/censo/registro.log"; }
finish() { log "end: $1"; gcloud compute instances delete "$NAME" --zone "$ZONE" --quiet || shutdown -h now; exit; }
export -f log

log "start, $(nproc) vCPU, $(free -g | awk '/Mem/ {print $2}') GB, $PAR blocks at a time"
mkdir -p /censo && cd /censo
curl -LsSf https://astral.sh/uv/install.sh | sh > /dev/null 2>&1
export PATH=/root/.local/bin:$PATH
gsutil -q cp "$B/entrada.tar.gz" . && tar xzf entrada.tar.gz || finish "inputs not found in $B"
uv venv -q -p 3.13 .venv && . .venv/bin/activate && uv pip install -q -r requirements.txt || finish "pip install failed"
python -c "import whitebox; whitebox.WhiteboxTools().version()" > /dev/null || finish "WhiteboxTools failed"
log "environment ready"
timeout 40m python -P scripts/39_descarga_fabdem.py "$EP" >> /var/log/censo_39.log 2>&1 || {
  gsutil -q cp /var/log/censo_39.log "$B/censo/"; finish "FABDEM from Earth Engine failed (see censo_39.log)"; }
log "FABDEM ready: $(ls data/raw/fabdem | wc -l) tiles"

bloque() {
  north=$1
  tag=$(python -c "n=$north; print(f\"{'S' if n <= 0 else 'N'}{abs(n):05.2f}\".replace('.', 'p'))")
  t0=$(date +%s)
  gsutil -q rm "$B/censo/estado/$tag.*" 2>/dev/null       # a rerun replaces the earlier state
  export NUMBA_CACHE_DIR=/censo/numba/$tag      # one cache per block: no race between parallel compilations
  if timeout -k 1m 30m python -P scripts/40_censo_bloque.py "$north" > "log_$tag.txt" 2>&1 &&
     timeout -k 1m 10m python -P scripts/41_etiquetas_passes.py "$tag" >> "log_$tag.txt" 2>&1; then
    st=ok
  else
    st=failed; pkill -f "whitebox_tools.*_${tag}\.tif"
  fi
  o=data/interim/censo
  rm -f $o/tmp/*_"$tag".tif
  for f in $o/passes_"$tag".csv $o/passes_"$tag".gpkg $o/passes_"$tag"_etiquetas.csv $o/sumideros_"$tag".csv $o/divisorias_"$tag".tif; do
    [ -f "$f" ] && gsutil -q cp "$f" "$B/censo/"
  done
  gsutil -q cp "log_$tag.txt" "$B/censo/logs/"
  echo "$(( ($(date +%s) - t0) / 60 )) min" | gsutil -q cp - "$B/censo/estado/$tag.$st"
  log "$tag $st in $(( ($(date +%s) - t0) / 60 )) min"
}
export -f bloque
BL=$(M blocks || true)
{ [ -n "$BL" ] && echo "$BL" | tr ';' '\n' || seq 12 -1 -56; } | xargs -P "$PAR" -I{} bash -c 'bloque {}'
finish "census done: $(gsutil ls "$B/censo/estado/" | grep -c '\.ok$') blocks ok, $(gsutil ls "$B/censo/estado/" | grep -c '\.failed$') failed"
