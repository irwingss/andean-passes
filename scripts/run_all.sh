#!/usr/bin/env bash
# Rebuilds the census and the present connectivity from the inputs, in order, from the repository root.
# Inputs first: see "Reproducing the analysis" in README.md. The census (40 to 46) ran on Google Cloud with
# 44_censo_nube.sh; the lines below run the same scripts locally, block by block, which takes many hours.
set -euo pipefail
Rscript scripts/02_dominios_divisoria.R --force                  # study area and sea domains
EP=${EP:-ee-irwingssaldu}                                          # Earth Engine project (FABDEM)
python -P scripts/39_descarga_fabdem.py "$EP"                      # FABDEM tiles for every block window
for n in $(seq 12 -1 -56); do                                      # census: one block of 1 degree per call
  tag=$(python -c "n=$n; print(f\"{'S' if n <= 0 else 'N'}{abs(n):05.2f}\".replace('.', 'p'))")
  python -P scripts/40_censo_bloque.py "$n"
  python -P scripts/41_etiquetas_passes.py "$tag"
done
python -P scripts/45_union_aislamiento.py                          # blocks joined, isolation
python -P scripts/46_validacion_O5.py                              # validation of the census
python -P scripts/47_conectividad_insumos.py "$EP" principal 11 -23.5 all "0.1 0.05 0.2" env
for v in $(seq 1 19); do                                           # one variable per call keeps memory flat
  for s in 200 100 400; do python -P scripts/50_conectividad_teselas.py principal "$s" "$v" env; done
  for q in 0.05 0.2; do python -P scripts/50_conectividad_teselas.py principal 200 "$v" "env_q$q"; done
done
python -P scripts/49_passes_corriente.py principal 11 -23.5          # values at every pass, profiles, control
python -P scripts/51_sintesis_conectividad.py principal            # database, tables
python -P scripts/52_figuras.py principal                          # figures
python -P scripts/53_suplemento.py                                 # Supplementary Information
