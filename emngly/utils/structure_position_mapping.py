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
Extraccion de secuencia ATMSEQ + tabla de mapeo de posiciones
(fasta_position <-> pdb_seqid) desde un ``AtomStruct`` de pwchem, ya asumido
de UNA sola cadena. Vendorizado (misma politica que
scipion-chem-deepptmpred/deepptmpred/utils/structure_sequence.py) a partir
de ``PTM-Prediction/src/utils/structure_parser.py::parse_structure``, ya
validado end-to-end en el pipeline standalone.

La tabla de mapeo es indispensable para EMNGly: el runner vendorizado
(``scripts/emngly_runner.py``) necesita traducir la posicion 1-based de
``sequence`` (ATMSEQ) al numero de residuo REAL del PDB
(``pdb_seqid``, puede tener huecos/no arrancar en 1) antes de indexar
``structure_emb`` -- ver docstring de ese runner.
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
    """Devuelve ``(sequence, positionMappingRows)`` de la primera cadena polimero de ``pdbPath``.

    ``positionMappingRows`` es una lista de dicts con
    ``POSITION_MAPPING_COLUMNS`` (una fila por residuo, mismo orden que
    ``sequence``) -- exactamente lo que
    ``emngly_runner.py::_fasta_position_to_pdb_seqid`` espera leer de un CSV.

    Raises:
        StructurePositionMappingError: si el modelo 1 no tiene ninguna
            cadena con un polimero de aminoacidos valido.
    """
    structure = gemmi.read_structure(str(pdbPath))
    structure.setup_entities()

    if len(structure) == 0:
        raise StructurePositionMappingError(f"'{pdbPath}' no contiene ningun modelo (MODEL) parseable.")

    model = structure[0]
    chain = None
    for candidate in model:
        if candidate.get_polymer().length() > 0:
            chain = candidate
            break
    if chain is None:
        raise StructurePositionMappingError(
            f"El modelo 1 de '{pdbPath}' no tiene ninguna cadena con al menos un residuo de "
            "aminoacido valido en su polimero."
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
