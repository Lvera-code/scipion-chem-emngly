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

UPSTREAM_URL = 'https://github.com/StellaHxy/EMNgly'

# Confirmed by reading scripts/emngly_runner.py: 'torch.device("cuda" if
# torch.cuda.is_available() else "cpu")' -- decided in code, with no CLI
# flag. The protocol's USE_GPU/GPU_LIST hidden params act on that decision
# indirectly, via CUDA_VISIBLE_DEVICES (see runEMNGly in __init__.py).
GPU_REQUIRED = True

# EMNGly license (upstream): MIT. The LICENSE file was added to the
# original repository afterward, and the Scipion integration was
# explicitly authorized. Paper CC BY 4.0 (PMC10627407).

# MIF (Microsoft protein-sequence-models, BSD-2 license verified against
# the original repo -- the vendorized copy inside EMNgly lost its own
# LICENSE, see the runner's docstring) comes BUNDLED with the repo clone
# (model/MIF/weights/mif.pt) -- no separate download needed.
#
# ESM-1b (checkpoint + '-contact-regression.pt' companion, ~7.8GB) and the
# SVM (N-GlyDE.pickle, ~36MB) are now auto-downloaded at install time into
# '<EMNGLY_HOME>/checkpoints/' (see addEMNGlyPackage in __init__.py). The
# SVM's real direct URL is a Google Drive file-id link, confirmed via a
# real HTTP Range request against the live file ('curl -sIL': real 200,
# 'content-disposition: attachment; filename="N-GlyDE.pickle"').
ESM_CHECKPOINT_FILENAME = 'esm1b_t33_650M_UR50S.pt'
ESM_CONTACT_REGRESSION_FILENAME = 'esm1b_t33_650M_UR50S-contact-regression.pt'
ESM_DOWNLOAD_URL = 'https://dl.fbaipublicfiles.com/fair-esm/models/esm1b_t33_650M_UR50S.pt'
ESM_CONTACT_REGRESSION_URL = (
    'https://dl.fbaipublicfiles.com/fair-esm/regression/'
    'esm1b_t33_650M_UR50S-contact-regression.pt'
)
SVM_CHECKPOINT_FILENAME = 'N-GlyDE.pickle'
SVM_DOWNLOAD_URL = (
    'https://drive.usercontent.google.com/download?'
    'id=1hbnEtHHXTGnQAFm-cCHMj3pWQiAYAUsw&export=download&confirm=t'
)

NOINSTALL_WARNING = (
    "EMNGly is not installed correctly. Check that the repo has been cloned (EMNGLY_HOME, "
    "ships MIF/mif.pt bundled) and that the ESM-1b checkpoint + 'N-GlyDE.pickle' were "
    "auto-downloaded into '<EMNGLY_HOME>/checkpoints/' -- re-run 'scipion3 installb emngly' or, "
    "if that keeps failing, download them manually and point EMNGLY_ESM_CHECKPOINT/"
    "EMNGLY_SVM_CHECKPOINT at them. OPTIONAL consensus engine: its absence degrades the "
    "N-glycosylation consensus to DeepMVP alone, it does not block the rest of the pipeline."
)
