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

        cls.protImportPdb = cls._runImportPdb()
        cls.protPrepareReceptor = cls._runPrepareReceptorChainC(cls.protImportPdb)
        cls.protDeepmvpROIs = cls._buildSyntheticDeepmvpRois()
        # Run once here (real conda subprocess: ESM-1b + MIF + SVM), not per
        # test method -- the two test_ methods below only assert on its
        # already-computed output.
        cls.protEMNGly = cls._runEMNGlyPrediction(cls.protDeepmvpROIs, cls.protPrepareReceptor)

    @classmethod
    def _runImportPdb(cls):
        protImportPdb = cls.newProtocol(ProtImportPdb, inputPdbData=0, pdbId=_TEST_PDB_ID)
        cls.proj.launchProtocol(protImportPdb, wait=True)
        return protImportPdb

    @classmethod
    def _runPrepareReceptorChainC(cls, protImportPdb):
        protPrepareReceptor = cls.newProtocol(
            ProtChemPrepareReceptor, inputAtomStruct=protImportPdb.outputPdb,
            usePDBFixer=True, addRes=False, HETATM=False, rchains=True,
            chain_name='{"model": 0, "chain": "%s"}' % _TEST_CHAIN,
        )
        cls.proj.launchProtocol(protPrepareReceptor, wait=True)
        return protPrepareReceptor

    @classmethod
    def _buildSyntheticDeepmvpRois(cls):
        """Builds two synthetic 'DeepMVP-shaped' ROIs: one N-glycosylation
        (must be corroborated by EMNGly) and one acetylation (must pass
        through untouched, '_scoreEmngly=None').

        ``ProtDefineSeqROI`` has no way to set DeepMVP-specific attributes
        (``_type``, ``_scoreDeepmvp``, ``_fpr``, ``_passesThreshold``,
        ``_residueWt``) on its own output -- a real ``SetOfSequenceROIs``
        cannot be mutated in place once written (its rows are backed by an
        append-only sqlite file), so the only way to inject them is to
        materialize the items, delete that sqlite file, and rebuild a new
        set from scratch with the extra attributes added per item. This is
        the same rebuild trick this project already uses elsewhere for the
        same reason (see scipion-chem-ptmannotation).
        """
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
        return protDefSeqROIs

    @classmethod
    def _runEMNGlyPrediction(cls, protDeepmvpROIs, protPrepareReceptor):
        protEMNGly = cls.newProtocol(ProtEMNGlyPrediction)
        protEMNGly.inputROIs.set(protDeepmvpROIs)
        protEMNGly.inputROIs.setExtended('outputROIs')
        protEMNGly.inputStructure.set(protPrepareReceptor)
        protEMNGly.inputStructure.setExtended('outputStructure')
        cls.proj.launchProtocol(protEMNGly, wait=True)
        return protEMNGly

    def _getRoisByType(self):
        outROIs = getattr(self.protEMNGly, 'outputROIs', None)
        self.assertIsNotNone(outROIs)
        self.assertEqual(len(outROIs), 2)
        # .clone() is mandatory (same real gotcha already found in
        # scipion-chem-ptmannotation via a real 'scipion3 test' run:
        # iterating a SetOfXXX without cloning reuses the same Python
        # object -- underlying sqlite cursor -- for every row).
        return {roi.getType(): roi.clone() for roi in outROIs}

    def test_untouchedForOtherType(self):
        """A PTM type EMNGly never scores ('acetylation_k') must pass
        through the protocol with '_scoreEmngly' left unset (None)."""
        byType = self._getRoisByType()
        self.assertIsNone(byType['acetylation_k']._scoreEmngly.get())

    def test_scoreAssignedForGlycosylation(self):
        """'n_linked_glycosylation' IS sent to EMNGly -- with EMNGly
        installed (ESM-1b/SVM/MIF checkpoints present), it must produce a
        real probability, not None."""
        byType = self._getRoisByType()
        score = byType['n_linked_glycosylation']._scoreEmngly.get()
        self.assertIsNotNone(score)
        self.assertGreaterEqual(score, 0.0)
        self.assertLessEqual(score, 1.0)
