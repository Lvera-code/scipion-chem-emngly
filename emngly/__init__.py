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
"""
This package contains a protocol for N-linked glycosylation consensus
corroboration using a local EMNGly installation (ESM-1b + MIF -> SVM).
"""

import os
import subprocess

from scipion.install.funcs import InstallHelper

from pwchem import Plugin as pwchemPlugin

from .constants import (
    EMNGLY_DIC, ESM_CONTACT_REGRESSION_FILENAME, NOINSTALL_WARNING, UPSTREAM_URL,
)

# 'Hou2023' (Bioinformatics 39(11):btad650, 2023, PMC10627407) deliberately
# NOT listed here yet: bibtex.py still has the citation marked TODO (exact
# title/authors not verified from this machine) -- listing the key without
# the real bibtex would break citation resolution instead of just omitting
# it. Complete both at the same time once the exact text is verified.
_references = []


class Plugin(pwchemPlugin):
    """EMNGly (StellaHxy/EMNgly, MIT -- see constants.py) is installed by
    cloning the upstream repo (it ships the MIF weights bundled,
    ``model/MIF/weights/mif.pt``) and building a dedicated conda
    environment (Python 3.10, fair-esm/torch/scikit-learn), torch
    installed from the CPU-only index + purge of stray nvidia/triton
    packages to avoid the same real SIGSEGV already documented in
    scipion-chem-stackglyembed on a machine with no GPU. TWO pieces remain
    manual: the ESM-1b checkpoint (+ its contact-regression companion) and
    the trained SVM (``N-GlyDE.pickle``, downloaded from Google Drive)."""

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(EMNGLY_DIC['home'], cls.getEnvName(EMNGLY_DIC))
        cls._defineVar(EMNGLY_DIC['activation'], cls.getEnvActivationCommand(EMNGLY_DIC))
        cls._defineVar(EMNGLY_DIC['esm_checkpoint'], '')
        cls._defineVar(EMNGLY_DIC['svm_checkpoint'], '')

    @classmethod
    def defineBinaries(cls, env):
        cls.addEMNGlyPackage(env)

    @classmethod
    def addEMNGlyPackage(cls, env, default=True):
        home = cls.getVar(EMNGLY_DIC['home'])

        installer = InstallHelper(EMNGLY_DIC['name'], packageHome=home,
                                  packageVersion=EMNGLY_DIC['version'])

        # Clone BEFORE the conda environment (same rule already documented
        # in netcleave/deepmvp/deepptmpred/stackglyembed).
        #
        # CPU-only torch + purge of nvidia-*/triton (a real bug already
        # found and fixed for exactly this class of problem in
        # scipion-chem-stackglyembed/stackglyembed/__init__.py: a plain
        # 'pip install torch' resolves the CUDA build by default on Linux
        # even with no GPU/drivers, causing a real SIGSEGV inside
        # torch._dynamo when importing transformers/fair-esm) -- applied
        # here preventively, not just reactively: a default pip install of
        # fair-esm/torch pulls in the CUDA build even on a machine with no
        # GPU/drivers.
        #
        # Versions pinned to a known-good combination (fair-esm, numpy,
        # pandas, scikit-learn).
        installer.addCommand(
            f"git clone --depth 1 {UPSTREAM_URL} {home}",
            'EMNGLY_CLONED'
        ).getCondaEnvCommand(
            EMNGLY_DIC['name'], binaryVersion=EMNGLY_DIC['version'], pythonVersion='3.10'
        ).addCommand(
            f"{cls.getEnvActivationCommand(EMNGLY_DIC)} && "
            "pip install --index-url https://download.pytorch.org/whl/cpu torch && "
            "pip install 'numpy==1.23.5' 'pandas==2.3.3' 'scikit-learn==1.1.1' 'fair-esm==2.0.0' wget && "
            "pip uninstall -y cuda-bindings cuda-pathfinder cuda-toolkit nvidia-cublas "
            "nvidia-cuda-cupti nvidia-cuda-nvrtc nvidia-cuda-runtime nvidia-cudnn-cu13 "
            "nvidia-cufft nvidia-cufile nvidia-curand nvidia-cusolver nvidia-cusparse "
            "nvidia-cusparselt-cu13 nvidia-nccl-cu13 nvidia-nvjitlink nvidia-nvshmem-cu13 "
            "nvidia-nvtx triton || true",
            'EMNGLY_INSTALLED'
        ).addPackage(env, dependencies=['conda', 'git'], default=default)

    @classmethod
    def validateInstallation(cls):
        """Check that this plugin's requirements are met. Returns a list of
        actionable error messages, empty if the installation is correct."""
        errors = []

        mifWeights = os.path.join(cls.getEMNGlyDir(), 'model', 'MIF', 'weights', 'mif.pt')
        if not os.path.isfile(mifWeights):
            errors.append(f"Could not find MIF weights under EMNGLY_HOME: '{mifWeights}'.")
        elif not cls.checkCallEnv(EMNGLY_DIC):
            errors.append("Activation of the EMNGly conda environment failed.")

        esmCheckpoint = cls.getEsmCheckpointPath()
        if not esmCheckpoint or not os.path.isfile(esmCheckpoint):
            errors.append(f"EMNGLY_ESM_CHECKPOINT ('{esmCheckpoint}') not found.")
        elif not os.path.isfile(os.path.join(os.path.dirname(esmCheckpoint), ESM_CONTACT_REGRESSION_FILENAME)):
            errors.append(
                f"'{ESM_CONTACT_REGRESSION_FILENAME}' (required companion file) not found next to "
                'EMNGLY_ESM_CHECKPOINT.'
            )

        svmCheckpoint = cls.getSvmCheckpointPath()
        if not svmCheckpoint or not os.path.isfile(svmCheckpoint):
            errors.append(f"EMNGLY_SVM_CHECKPOINT ('{svmCheckpoint}') not found.")

        if errors:
            errors.append(NOINSTALL_WARNING)
        return errors

    @classmethod
    def checkCallEnv(cls, packageDic):
        actCommand = cls.getVar(packageDic['activation'])
        try:
            if 'conda' in actCommand and 'shell.bash hook' not in actCommand:
                actCommand = f'{cls.getCondaActivationCmd()}{actCommand}'
            subprocess.check_output(f'{actCommand} && python -c "import torch, esm, sklearn"', shell=True)
            return True
        except subprocess.CalledProcessError:
            return False

    # ---------------------------------- Utils -----------------------------------

    @classmethod
    def getEMNGlyDir(cls):
        return cls.getVar(EMNGLY_DIC['home'])

    @classmethod
    def getEsmCheckpointPath(cls):
        return cls.getVar(EMNGLY_DIC['esm_checkpoint'])

    @classmethod
    def getSvmCheckpointPath(cls):
        return cls.getVar(EMNGLY_DIC['svm_checkpoint'])

    @classmethod
    def getMifWeightsPath(cls):
        return os.path.join(cls.getEMNGlyDir(), 'model', 'MIF', 'weights', 'mif.pt')

    @classmethod
    def getRunnerScriptPath(cls):
        pluginDir = os.path.dirname(os.path.abspath(__file__))
        return os.path.join(pluginDir, 'scripts', 'emngly_runner.py')

    # ---------------------------------- Protocol functions-----------------------

    @classmethod
    def runEMNGly(cls, protocol, args, cwd=None):
        activation = cls.getVar(EMNGLY_DIC['activation'])
        scriptPath = cls.getRunnerScriptPath()
        fullProgram = f'MPLBACKEND=Agg {activation} && python {scriptPath}'
        protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
