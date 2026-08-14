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
Extraction of the ATMSEQ sequence + position mapping table
(fasta_position <-> pdb_seqid) from a pwchem ``AtomStruct``, already
assumed to be a SINGLE chain. Each plugin in this project keeps its own
minimal copy of this logic (same policy as
scipion-chem-deepptmpred/deepptmpred/utils/structure_sequence.py) rather
than a shared dependency.

The mapping table is essential for EMNGly: the vendorized runner
(``scripts/emngly_runner.py``) needs to translate the 1-based position of
``sequence`` (ATMSEQ) to the REAL PDB residue number (``pdb_seqid``, which
can have gaps/not start at 1) before indexing ``structure_emb`` -- see that
runner's docstring.
"""

import csv

import gemmi

POSITION_MAPPING_COLUMNS = ['fasta_position', 'pdb_seqid']


class StructurePositionMappingError(Exception):
    pass


def _resolveResidueLetter(resname):
    info = gemmi.find_tabulated_residue(resname)
    code = info.one_letter_code.strip().upper() if info is not None else ''
    return code if len(code) == 1 and code.isalpha() else 'X'


def extractSequenceAndMapping(pdbPath):
    """Returns ``(sequence, positionMappingRows)`` for the first polymer chain in ``pdbPath``.

    ``positionMappingRows`` is a list of dicts with
    ``POSITION_MAPPING_COLUMNS`` (one row per residue, same order as
    ``sequence``) -- exactly what
    ``emngly_runner.py::_fasta_position_to_pdb_seqid`` expects to read from
    a CSV.

    Raises:
        StructurePositionMappingError: if model 1 has no chain with a
            valid amino acid polymer.
    """
    structure = gemmi.read_structure(str(pdbPath))
    structure.setup_entities()

    if len(structure) == 0:
        raise StructurePositionMappingError(f"'{pdbPath}' does not contain any parseable model (MODEL).")

    model = structure[0]
    chain = None
    for candidate in model:
        if candidate.get_polymer().length() > 0:
            chain = candidate
            break
    if chain is None:
        raise StructurePositionMappingError(
            f"Model 1 of '{pdbPath}' has no chain with at least one valid amino acid "
            "residue in its polymer."
        )

    residues = list(chain.get_polymer())
    letters = []
    mappingRows = []
    for fastaPosition, residue in enumerate(residues, start=1):
        letter = _resolveResidueLetter(residue.name)
        letters.append(letter)
        mappingRows.append({'fasta_position': fastaPosition, 'pdb_seqid': residue.seqid.num})

    return ''.join(letters), mappingRows


def writePositionMappingCsv(mappingRows, csvPath):
    with open(csvPath, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=POSITION_MAPPING_COLUMNS)
        writer.writeheader()
        writer.writerows(mappingRows)
