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
    EMNGLY_DIC, ESM_CHECKPOINT_FILENAME, ESM_CONTACT_REGRESSION_FILENAME,
    ESM_CONTACT_REGRESSION_URL, ESM_DOWNLOAD_URL, NOINSTALL_WARNING, SVM_CHECKPOINT_FILENAME,
    SVM_DOWNLOAD_URL, UPSTREAM_URL,
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
    environment (Python 3.11, fair-esm/torch/scikit-learn), torch
    installed from the CPU-only index + purge of stray nvidia/triton
    packages to avoid the same real SIGSEGV already documented in
    scipion-chem-stackglyembed on a machine with no GPU. The ESM-1b
    checkpoint (+ its contact-regression companion) and the trained SVM
    (``N-GlyDE.pickle``) are auto-downloaded at install time into
    ``<EMNGLY_HOME>/checkpoints/``."""

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(EMNGLY_DIC['home'], cls.getEnvName(EMNGLY_DIC))
        cls._defineVar(EMNGLY_DIC['activation'], cls.getEnvActivationCommand(EMNGLY_DIC))
        # Empty by default: 'getEsmCheckpointPath()'/'getSvmCheckpointPath()'
        # below fall back to where addEMNGlyPackage auto-downloads them
        # ('<EMNGLY_HOME>/checkpoints/') when unset.
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
        # Installed FROM the repo's own real 'environment.yml' (via
        # 'conda env update -f'), not a hand-reconstructed package list --
        # pythonVersion bumped to '3.11' to match its real 'python=3.11.0'
        # pin (previously '3.10', never actually read from the file).
        # 'pytorch'/'numpy'/'pandas'/'scikit-learn' are filtered out of the
        # file before that update, then installed separately with the
        # EXACT versions already production-validated for this plugin
        # (MCC=0.82 real go/no-go run, matches the published benchmark) --
        # NOT the file's own newer pins, which would regress a validated
        # result:
        #   * scikit-learn==1.1.1 (not the file's 1.5.1): the downloaded
        #     SVM pickle's own '_sklearn_version' is 1.1.1 -- confirmed by
        #     inspecting the pickle directly, not assumed.
        #   * numpy==1.23.5/pandas==2.3.3 (not the file's 1.26.4/2.2.3):
        #     the already-tested combination compatible with that
        #     scikit-learn==1.1.1 build.
        #   * torch (CPU-only wheel, not the file's conda-channel
        #     'pytorch=2.3.1' -- avoids pulling a CUDA build via conda
        #     with no matching NVIDIA drivers, same reasoning as
        #     scipion-chem-deepptmpred).
        # 'fair-esm'/'wget' are NOT in the file at all (verified) --
        # genuinely additional packages, not an override.
        #
        # Purge of nvidia-*/triton kept as a defensive no-op-if-absent
        # safety net (the real SIGSEGV this originally fixed, in
        # scipion-chem-stackglyembed, came from a default 'pip install
        # torch' pulling in a CUDA build -- the explicit CPU-only index
        # above should already prevent that here, but this costs nothing
        # to keep).
        installer.addCommand(
            f"git clone --depth 1 {UPSTREAM_URL} {home}",
            'EMNGLY_CLONED'
        ).getCondaEnvCommand(
            EMNGLY_DIC['name'], binaryVersion=EMNGLY_DIC['version'], pythonVersion='3.11'
        ).addCommand(
            f"grep -vE '^[[:space:]]*-[[:space:]]*(pytorch|numpy|numpy-base|pandas|scikit-learn)"
            f"([[:space:]]*=|$)' {home}/environment.yml > {home}/environment_filtered.yml && "
            f"{cls.getCondaActivationCmd()}conda env update -n {cls.getEnvName(EMNGLY_DIC)} "
            f"-f {home}/environment_filtered.yml",
            'EMNGLY_BASE_ENV_UPDATED'
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
        ).addCommand(
            # ESM-1b + SVM checkpoint auto-download (see constants.py for
            # the SVM Google Drive URL's real verification history).
            f"mkdir -p {home}/checkpoints && "
            f"curl -fsSL --retry 3 -o {home}/checkpoints/{ESM_CHECKPOINT_FILENAME} {ESM_DOWNLOAD_URL} && "
            f"curl -fsSL --retry 3 -o {home}/checkpoints/{ESM_CONTACT_REGRESSION_FILENAME} "
            f"{ESM_CONTACT_REGRESSION_URL} && "
            f"curl -fsSL --retry 3 -o {home}/checkpoints/{SVM_CHECKPOINT_FILENAME} \"{SVM_DOWNLOAD_URL}\"",
            'EMNGLY_CHECKPOINTS_DOWNLOADED'
        ).addPackage(env, dependencies=['conda', 'git', 'curl'], default=default)

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
        configured = cls.getVar(EMNGLY_DIC['esm_checkpoint'])
        if configured:
            return configured
        return os.path.join(cls.getEMNGlyDir(), 'checkpoints', ESM_CHECKPOINT_FILENAME)

    @classmethod
    def getSvmCheckpointPath(cls):
        configured = cls.getVar(EMNGLY_DIC['svm_checkpoint'])
        if configured:
            return configured
        return os.path.join(cls.getEMNGlyDir(), 'checkpoints', SVM_CHECKPOINT_FILENAME)

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
