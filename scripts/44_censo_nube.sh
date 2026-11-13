#!/usr/bin/env bash
# 44_censo_nube.sh - Runs the whole census (40 and 41 on every block, PRD R3) on one Google Cloud machine of the
# user's own account. Run from the repository root with the Google Cloud CLI logged in:
#   bash scripts/44_censo_nube.sh [gcp_project] [earth_engine_project]
# Both default to the active gcloud project. The Earth Engine project must be registered for Earth Engine.
# Settings by environment variable: ZONE (us-central1-a), MACHINE (n2-highmem-32), PAR (16 blocks at a time),
# MAXRUN (3h), BLOCKS (only these blocks, by north edge of the core, e.g. "-32 -40", to rerun the failed ones).
# Expected cost with the defaults: about 3 USD (1.2 h of n2-highmem-32 plus disk); at most MAXRUN of the
# machine, because Compute Engine deletes it at that time whatever happens.
# Results go to gs://<project>-censo-l1/censo/; fetch them with:
#   gcloud storage cp -r "gs://<project>-censo-l1/censo/*" data/interim/censo/
set -euo pipefail
P=${1:-$(gcloud config get-value project 2>/dev/null)}
EP=${2:-$P}
ZONE=${ZONE:-us-central1-a}; MACHINE=${MACHINE:-n2-highmem-32}; PAR=${PAR:-16}; MAXRUN=${MAXRUN:-3h}
B=gs://$P-censo-l1
echo "project $P, Earth Engine $EP, $MACHINE in $ZONE, $PAR blocks at a time, deleted after $MAXRUN at most"

gcloud services enable compute.googleapis.com storage.googleapis.com earthengine.googleapis.com --project "$P"
gcloud storage buckets describe "$B" --project "$P" > /dev/null 2>&1 || gcloud storage buckets create "$B" --project "$P" --location "${ZONE%-*}"
if gcloud storage objects describe "$B/censo/registro.log" > /dev/null 2>&1; then
  if [ -z "${BLOCKS:-}" ]; then echo "$B/censo already holds a run; set BLOCKS to rerun some blocks" >&2; exit 1; fi
  gcloud storage mv "$B/censo/registro.log" "$B/censo/registro_$(date +%Y%m%d%H%M%S).log"
fi
# the machine's service account reads Earth Engine and deletes the machine when it ends
SA=$(gcloud compute project-info describe --project "$P" --format="value(defaultServiceAccount)")
for r in roles/earthengine.viewer roles/serviceusage.serviceUsageConsumer; do
  gcloud projects add-iam-policy-binding "$EP" --member "serviceAccount:$SA" --role "$r" --condition=None -q > /dev/null
done

# code and the small inputs (FABDEM is read on the machine from Earth Engine)
tar czf /tmp/entrada.tar.gz config.yml requirements.txt scripts/39_descarga_fabdem.py scripts/40_censo_bloque.py \
  scripts/41_etiquetas_passes.py data/refs/pasos_control.csv data/interim/M1_gmba_andes.gpkg data/interim/M1_domains.gpkg \
  data/raw/hybas_sa_lev01-12_v1c/hybas_sa_lev12_v1c.* data/raw/GMBA_Inventory_v2.0_standard_300
gcloud storage cp -q /tmp/entrada.tar.gz "$B/entrada.tar.gz"

gcloud compute instances create censo-l1 --project "$P" --zone "$ZONE" --machine-type "$MACHINE" \
  --image-family debian-12 --image-project debian-cloud --boot-disk-size 300GB --boot-disk-type pd-balanced \
  --scopes cloud-platform --max-run-duration "$MAXRUN" --instance-termination-action DELETE \
  --metadata "bucket=$B,ee_project=$EP,par=$PAR,zone_name=$ZONE,blocks=$(echo ${BLOCKS:-} | tr ' ' ';')" --metadata-from-file startup-script=scripts/44_vm_censo.sh
echo "running. Progress: gcloud storage cat $B/censo/registro.log    blocks: gcloud storage ls $B/censo/estado/"
