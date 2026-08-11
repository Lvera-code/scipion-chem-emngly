# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Enzo Sierra (enzogael57@gmail.com)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

DEFAULT_VERSION = '1.0'

EMNGLY_DIC = {
    'name': 'EMNGly',
    'version': DEFAULT_VERSION,
    'home': 'EMNGLY_HOME',
    'activation': 'EMNGLY_ACTIVATION_CMD',
    'esm_checkpoint': 'EMNGLY_ESM_CHECKPOINT',
    'svm_checkpoint': 'EMNGLY_SVM_CHECKPOINT',
}

READ_URL = 'https://github.com/Lvera-code/scipion-chem-emngly'
UPSTREAM_URL = 'https://github.com/StellaHxy/EMNgly'

# Confirmado leyendo scripts/emngly_runner.py: 'torch.device("cuda" if
# torch.cuda.is_available() else "cpu")' -- decidido en codigo, sin ningun
# flag de CLI. No se agrega parametro useGPU al protocolo.
GPU_REQUIRED = True

# Licencia de EMNGly (upstream): MIT (resuelto 2026-08-11). Los autores de correspondencia (Yaojun Wang, Shiwei Sun) anadieron el LICENSE MIT al repo tras el correo del 2026-08-10 y autorizaron la integracion en Scipion. Paper CC BY 4.0 (PMC10627407).

# MIF (Microsoft protein-sequence-models, licencia BSD-2 verificada contra
# el repo original -- la copia vendorizada dentro de EMNgly perdio su propio
# LICENSE, ver docstring del runner) viene BUNDLED con el clon del repo
# (model/MIF/weights/mif.pt) -- no requiere descarga aparte.
#
# ESM-1b (checkpoint + companero '-contact-regression.pt') y el SVM
# (N-GlyDE.pickle, ~36MB, descarga manual desde Google Drive segun el README
# real de EMNgly) SI son manuales -- mismo patron que DEEPMVP_MODEL_DIR/
# DEEPPTMPRED_ESM_CHECKPOINT.
ESM_CHECKPOINT_FILENAME = 'esm1b_t33_650M_UR50S.pt'
ESM_CONTACT_REGRESSION_FILENAME = 'esm1b_t33_650M_UR50S-contact-regression.pt'
ESM_DOWNLOAD_URL = 'https://dl.fbaipublicfiles.com/fair-esm/models/esm1b_t33_650M_UR50S.pt'
ESM_CONTACT_REGRESSION_URL = (
    'https://dl.fbaipublicfiles.com/fair-esm/regression/'
    'esm1b_t33_650M_UR50S-contact-regression.pt'
)
SVM_CHECKPOINT_FILENAME = 'N-GlyDE.pickle'

NOINSTALL_WARNING = (
    "EMNGly no esta instalado correctamente. Revisa que el repo se haya clonado (EMNGLY_HOME, "
    "trae MIF/mif.pt bundled), que EMNGLY_ESM_CHECKPOINT apunte al checkpoint ESM-1b (+ su "
    "companero '-contact-regression.pt' en el mismo directorio), y que EMNGLY_SVM_CHECKPOINT "
    "apunte a 'N-GlyDE.pickle' (descarga manual desde Google Drive, ver README.rst). Motor de "
    "consenso OPCIONAL: su ausencia degrada el consenso de N-glicosilacion a solo DeepMVP, no "
    'bloquea el resto del pipeline.'
)
