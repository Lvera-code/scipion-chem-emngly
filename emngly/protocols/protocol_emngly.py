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
This protocol corroborates N-linked glycosylation candidate sites (from
scipion-chem-deepmvp's output) with a local EMNGly installation.
"""

import csv
import os

from pwchem.objects import SetOfSequenceROIs
from pwem.protocols import EMProtocol
from pyworkflow.object import Float
from pyworkflow.protocol import params

from .. import Plugin as emnglyPlugin
from ..utils.structure_position_mapping import extractSequenceAndMapping, writePositionMappingCsv

NGLYCO_TYPES = {'n_linked_glycosylation', 'glycosylation_n'}


class ProtEMNGlyPrediction(EMProtocol):
    """
    AI Generated:

    Corroborates N-linked glycosylation candidate sites reported by
    scipion-chem-deepmvp with a local EMNGly installation (ESM-1b sequence
    embeddings + MIF structural embeddings -> SVM). Unlike
    scipion-chem-metoken/Kinase Library (purely informative), EMNGly is a
    REAL consensus engine for this ONE PTM type -- it replaces DeepPTMPred
    here (DeepPTMPred has no real discriminative power for
    n_linked_glycosylation, AUROC~=0.5 in its own per-type calibration) --
    but this protocol itself follows the SAME "annotate, never decide"
    mechanics as every other corroboration step: it only ADDS a score,
    downstream ``ProtPTMAnnotation`` (run on this
    protocol's output instead of DeepMVP's raw output) is what actually
    decides the N-glycosylation consensus.

    Only rows with ``_type`` in {'n_linked_glycosylation',
    'glycosylation_n'} are sent to EMNGly; every other row (and every row if
    EMNGly failed/isn't installed) gets ``_scoreEmngly=None`` -- schema
    stays uniform across the whole output set either way.

    Output
    ------
    outputROIs: the same ROIs as ``inputROIs``, each with a new
    ``_scoreEmngly`` (Float, ``None`` where not applicable/not computed).
    Does not filter.
    """

    _label = 'emngly n-glycosylation corroboration'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       label='DeepMVP candidates: ',
                       help='Output of scipion-chem-deepmvp. Only N-glycosylation rows '
                            "('n_linked_glycosylation'/'glycosylation_n') get a real score; "
                            'every other row is passed through unchanged.')
        form.addParam('inputStructure', params.PointerParam, pointerClass='AtomStruct',
                       label='Input structure (single chain): ',
                       help='MUST be derived from the SAME single-chain structure that produced '
                            "the Sequence fed into scipion-chem-deepmvp (otherwise positions won't "
                            'align) -- same convention as scipion-chem-discotope/-scannet.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.emnglyStep)
        self._insertFunctionStep(self.createOutputStep)

    # ---------------------------------- Steps -----------------------------------

    def _getRois(self):
        return [roi.clone() for roi in self.inputROIs.get()]

    def emnglyStep(self):
        rois = self._getRois()
        nglycoPositions = sorted({roi.getROIIdx() for roi in rois if roi.getType() in NGLYCO_TYPES})

        # ABSOLUTE paths are mandatory (same pattern as scipion-chem-deepmvp):
        # the subprocess runs with cwd=EMNGLY_HOME.
        outCsv = os.path.abspath(self._getExtraPath('emngly_scores.csv'))
        if not nglycoPositions:
            return

        pdbPath = os.path.abspath(self.inputStructure.get().getFileName())
        sequence, mappingRows = extractSequenceAndMapping(pdbPath)
        accession = os.path.splitext(os.path.basename(pdbPath))[0]
        mappingCsv = os.path.abspath(self._getExtraPath('position_mapping.csv'))
        writePositionMappingCsv(mappingRows, mappingCsv)

        cacheDir = os.path.abspath(self._getExtraPath('emngly_cache'))
        os.makedirs(cacheDir, exist_ok=True)

        args = (
            f'--emngly-home {emnglyPlugin.getEMNGlyDir()} --accession {accession} '
            f'--sequence {sequence} --pdb-path {pdbPath} --position-mapping-csv {mappingCsv} '
            f'--positions {" ".join(str(p) for p in nglycoPositions)} '
            f'--mif-weights {emnglyPlugin.getMifWeightsPath()} '
            f'--esm-checkpoint {emnglyPlugin.getEsmCheckpointPath()} '
            f'--svm-checkpoint {emnglyPlugin.getSvmCheckpointPath()} '
            f'--cache-dir {cacheDir} --out-csv {outCsv}'
        )
        try:
            emnglyPlugin.runEMNGly(self, args, cwd=emnglyPlugin.getEMNGlyDir())
        except Exception as exc:  # noqa: BLE001 -- OPTIONAL consensus engine, degrades without failing the protocol
            self.warning(f'EMNGly failed (non-fatal, degrades to no corroboration): {exc}')

    def createOutputStep(self):
        rois = self._getRois()
        if not rois:
            return

        scores = {}
        outCsv = self._getExtraPath('emngly_scores.csv')
        if os.path.isfile(outCsv):
            with open(outCsv, newline='') as fh:
                for row in csv.DictReader(fh):
                    scores[int(row['position'])] = float(row['probability'])

        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for roi in rois:
            score = scores.get(roi.getROIIdx()) if roi.getType() in NGLYCO_TYPES else None
            roi._scoreEmngly = Float(score) if score is not None else Float(None)
            outROIs.append(roi)

        if len(outROIs) > 0:
            self._defineOutputs(outputROIs=outROIs)
            self._defineSourceRelation(self.inputROIs, outROIs)
            self._defineSourceRelation(self.inputStructure, outROIs)

    # ---------------------------------- Validation -------------------------------

    def _validate(self):
        # EMNGly is OPTIONAL (degrades, see docstring) -- its absence never
        # blocks launching, unlike DeepMVP/DeepPTMPred.
        return []

    def _summary(self):
        summary = []
        if self.isFinished():
            outROIs = getattr(self, 'outputROIs', None)
            if outROIs is not None:
                nScored = sum(1 for roi in outROIs if roi._scoreEmngly.get() is not None)
                summary.append(f'{nScored} N-glycosylation site(s) corroborated by EMNGly.')
        return summary
