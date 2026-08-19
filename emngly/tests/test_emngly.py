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

import os

from pwchem.objects import SetOfSequenceROIs
from pwchem.protocols import ProtChemPrepareReceptor, ProtDefineSeqROI
from pwem.protocols import ProtImportPdb, ProtImportSequence
from pyworkflow.object import Boolean, Float, String
from pyworkflow.tests import BaseTest, setupTestProject

from ..protocols import ProtEMNGlyPrediction

# Same real fixture already used in scipion-chem-discotope/scipion-chem-deepptmpred
# (7c4s, mmCIF label_asym_id 'C' == author chain 'A' gotcha, 283 residues).
_TEST_PDB_ID = '7c4s'
_TEST_CHAIN = 'C'


class TestEMNGlyPrediction(BaseTest):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        setupTestProject(cls)

        protImportPdb = cls.newProtocol(ProtImportPdb, inputPdbData=0, pdbId=_TEST_PDB_ID)
        cls.proj.launchProtocol(protImportPdb, wait=True)

        cls.protPrepareReceptor = cls.newProtocol(
            ProtChemPrepareReceptor, inputAtomStruct=protImportPdb.outputPdb,
            usePDBFixer=True, addRes=False, HETATM=False, rchains=True,
            chain_name='{"model": 0, "chain": "%s"}' % _TEST_CHAIN,
        )
        cls.proj.launchProtocol(cls.protPrepareReceptor, wait=True)

        # Synthetic 'DeepMVP-shaped' ROIs -- one N-glycosylation (to be
        # corroborated) and one acetylation (must pass through untouched,
        # '_scoreEmngly=None').
        protImportSeq = cls.newProtocol(
            ProtImportSequence, inputSequenceName='EMNGLY_TEST_SEQ',
            inputSequenceDescription='placeholder, not used by ProtEMNGlyPrediction',
            inputRawSequence='X' * 20,
        )
        cls.proj.launchProtocol(protImportSeq, wait=True)

        inROIs = '1) Residues: {"index": "5-5", "residues": "X", "desc": "None"}\n' \
                 '2) Residues: {"index": "10-10", "residues": "X", "desc": "None"}'
        protDefSeqROIs = cls.newProtocol(ProtDefineSeqROI, chooseInput=0, inROIs=inROIs)
        protDefSeqROIs.inputSequence.set(protImportSeq)
        protDefSeqROIs.inputSequence.setExtended('outputSequence')
        cls.proj.launchProtocol(protDefSeqROIs, wait=True)

        sqlitePath = protDefSeqROIs.outputROIs.getFileName()
        oldItems = [roi.clone() for roi in protDefSeqROIs.outputROIs]
        os.remove(sqlitePath)
        rebuilt = SetOfSequenceROIs(filename=sqlitePath)
        for i, roi in enumerate(oldItems):
            roi.setType('n_linked_glycosylation' if i == 0 else 'acetylation_k')
            roi._scoreDeepmvp = Float(0.8)
            roi._fpr = Float(0.02)
            roi._passesThreshold = Boolean(True)
            roi._residueWt = String('N' if i == 0 else 'K')
            rebuilt.append(roi)
        rebuilt.write()
        cls.protDeepmvpROIs = protDefSeqROIs

    def _runEMNGlyPrediction(self):
        protEMNGly = self.newProtocol(ProtEMNGlyPrediction)
        protEMNGly.inputROIs.set(self.protDeepmvpROIs)
        protEMNGly.inputROIs.setExtended('outputROIs')
        protEMNGly.inputStructure.set(self.protPrepareReceptor)
        protEMNGly.inputStructure.setExtended('outputStructure')
        self.launchProtocol(protEMNGly, wait=True)
        return protEMNGly

    def test(self):
        protEMNGly = self._runEMNGlyPrediction()

        outROIs = getattr(protEMNGly, 'outputROIs', None)
        self.assertIsNotNone(outROIs)
        self.assertEqual(len(outROIs), 2)

        # .clone() is mandatory (same real gotcha already found in
        # scipion-chem-ptmannotation via a real 'scipion3 test' run:
        # iterating a SetOfXXX without cloning reuses the same Python
        # object -- underlying sqlite cursor -- for every row).
        byType = {roi.getType(): roi.clone() for roi in outROIs}
        # 'acetylation_k' is never sent to EMNGly -- always None.
        self.assertIsNone(byType['acetylation_k']._scoreEmngly.get())
        # 'n_linked_glycosylation' IS attempted -- with EMNGly installed
        # (ESM-1b/SVM/MIF checkpoints present), it produces a real
        # probability, not None.
        score = byType['n_linked_glycosylation']._scoreEmngly.get()
        self.assertIsNotNone(score)
        self.assertGreaterEqual(score, 0.0)
        self.assertLessEqual(score, 1.0)
