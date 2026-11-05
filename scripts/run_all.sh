#!/usr/bin/env bash
# Rebuilds everything from the inputs, in order, from the repository root. --force redoes steps whose outputs are in the repository.
# Inputs first: see "Reproducing the analysis" in README.md. The flood (09) skips blocks already in data/interim/M1_glo30_divide;
# delete that folder to run it.
set -euo pipefail
Rscript scripts/02_dominios_divisoria.R --force
Rscript scripts/03_referencias_y_grilla.R
Rscript scripts/08_marcadores.R --force
python -P scripts/07_descarga_glo30.py
python -P scripts/09_divisoria_glo30.py
python -P scripts/09b_divisoria_altiplano.py
Rscript scripts/10_perfil_glo30.R
Rscript scripts/04_collados.R --force
Rscript scripts/06_validacion_huancabamba.R
Rscript scripts/11_figuras_M1.R
Rscript scripts/12_kml_M1.R
Rscript scripts/14_flancos_M2.R --force
Rscript scripts/15_metricas_clima_M2.R --force
Rscript scripts/17_janzen_M2.R --force
Rscript scripts/31_janzen_mecanismo_M2.R
Rscript scripts/20_permeabilidad_M3.R --force
Rscript scripts/29_corredor_montano_M3.R
Rscript scripts/30_proteccion_M3.R
Rscript scripts/26_figuras_manuscrito.R
Rscript scripts/27_tablas_manuscrito.R
Rscript scripts/28_APD.R
Rscript scripts/32_fig5_sintesis.R
Rscript scripts/33_suplementos.R
python -P scripts/34_build_supp.py
python -P scripts/visor_datos.py
