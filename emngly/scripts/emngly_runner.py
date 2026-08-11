#!/usr/bin/env python
"""Standalone runner for EMNGly (real consensus engine for 'n_linked_glycosylation', PDB path).

VENDORIZED byte-for-byte from
``PTM-Prediction/src/engines/_emngly_runner.py`` -- same policy as
``scipion-chem-deepptmpred``: the patches it contains (scoped
weights_only=False, structure_emb alignment via position_mapping) are
never rewritten from memory, they are synced from the standalone project.

NEVER imported from the ``src`` package -- it requires fair-esm/torch/
scikit-learn, dependencies of the dedicated venv ``Settings.EMNGLY_PYTHON_BIN``
(never shared with DEEPMVP_PYTHON_BIN/DEEPPTMPRED_PYTHON_BIN/others). It is
invoked EXCLUSIVELY via subprocess from ``src/engines/emngly_engine.py``,
same pattern as ``_deepptmpred_runner.py``.

## Why it exists (role in the pipeline)

Replaces CoNglyPred (Decision 2's original candidate, with no weights
published anywhere verifiable -- see STATUS.md). EMNGly
(``github.com/StellaHxy/EMNgly``, Hou et al., Bioinformatics 39(11):btad650,
2023) DOES have real weights, verified at the byte level
(see ``src/config/settings.py``, EMNGLY_* block): ESM-1b (sequence,
``site_emb``+``local_emb``, 1280+1280) + Microsoft's MIF (real structure
over the PDB's N/CA/C backbone, ``structure_emb``, 256) -> SVM (2816
features). Preserves the design property "the second engine of this type
uses real 3D structure", already decided when MTPrompt-PTM was ruled out
as a DeepPTMPred replacement for being sequence-only.

This runner directly imports the ``MIF`` package vendorized inside the
EMNgly clone (``EMNgly/model/MIF/``, Microsoft's ``protein-sequence-models``,
permissive BSD-2 license verified by reading the original at
``github.com/microsoft/protein-sequence-models`` -- EMNgly's copy lost the
LICENSE file when vendorizing it, which is why this project documents the
real license here instead of assuming it from the copy) for computing
``structure_emb``, and REIMPLEMENTS ``model/get_esm_embedding.py::ESMEmbeddingExtractor``
for ``site_emb``/``local_emb`` -- the EXACT same chunking logic (needed to
reproduce bit-for-bit what the SVM learned), but loading the ESM-1b
checkpoint from a LOCAL ``.pt`` path instead of
``torch.hub.load("facebookresearch/esm:main", ...)`` (the original script
hits the network on every run, violating this project's local-only
policy, same pattern already resolved in ``_deepptmpred_runner.py``).

## ``site_emb``/``local_emb`` indexing convention (verified by reading
   ``model/get_esm_embedding.py`` line by line, not assumed)

``ESMEmbeddingExtractor.extract()`` NEVER drops the BOS/start token from
the raw ESM representation for the first chunk (``delta=0`` when
``i==0``, the slice ``[0:j-i+1]`` includes index 0) -- the resulting
array ends up with the BOS token at index 0 and the 1-based residue ``k``
at index ``k`` (NOT at ``k-1``). This makes ``get_site_features``
(``emb = extract([seq])[pos]``, ``pos`` 1-based) CORRECT exactly as
written -- there is no real off-by-one despite the initial appearance (a
careful read of the original code confirms it is not a bug).
``get_local_features`` (``emb = extract([local_seq])[0]``) is a
deliberately different design: it takes the BOS token's embedding from the
local window (not a residue), as a pooled summary of the context -- same
pattern as using the [CLS] token in BERT as a representation of the whole
segment. This runner reproduces BOTH conventions exactly, without
"fixing" either: the SVM was trained on exactly these features, changing
them would create a new train/inference mismatch (the same class of bug
that already cost 3 real iterations in ``_deepptmpred_runner.py``, see
its docstring).

## ``structure_emb`` indexing convention (a different, real risk, see below)

``model/MIF/sequence_models/pdb_utils.py::parse_PDB`` builds the output
array indexed by PDB RESIDUE NUMBER minus 1 (``resn = int(resid) - 1``,
with ``min_resn`` forced to ``min(min_resn, 0)`` and gaps filled with
'X'/NaN) -- NEVER by the sequential order of observed residues. EMNgly's
own ``predict.py`` indexes ``structure_emb[pos-1]``, implicitly assuming
that ``pos`` (the position in its training CSV) matches the real PDB
numbering -- true for the AlphaFold2 PDBs in its datasets (continuous 1..N
numbering with no gaps), but FALSE in general for crystallographic
structures with gaps or numbering that does not start at 1. This project
uses the mapping table ALREADY built in Phase 1.5
(``src/utils/structure_parser.py::StructureRecord.position_mapping``,
``fasta_position``/``pdb_seqid`` columns) to translate the 1-based
position of ``sequence`` (ATMSEQ, the one the rest of the pipeline uses)
to the REAL PDB residue number before indexing ``structure_emb`` --
correctly generalizes the original script's implicit assumption instead
of blindly assuming it (a proactive fix for the risk identified during
the CoNglyPred-replacement investigation, before a real case broke it).

``parse_PDB`` with ``chain=None`` (default, same as
``model/get_mif_embedding.py::run``) reads ALL chains in the file
unfiltered -- only safe if ``--pdb-path`` is a SINGLE-chain PDB. This
runner requires ``record.chain_pdb_path`` (Phase 1.5, never the raw
``record.pdb_path``), same criterion already established for
MeToken/DeepPTMPred.
"""

import argparse
import hashlib
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

OUTPUT_COLUMNS = ["position", "probability"]

# --- Reimplementation of model/get_esm_embedding.py::ESMEmbeddingExtractor ---
# Constants identical to the original script (see the module docstring).
_ESM_EMBED_LAYER = 33
_ESM_EMBED_DIM = 1280
_LOCAL_WINDOW_LEFT = 21
_LOCAL_WINDOW_RIGHT = 20


class _ESMEmbeddingExtractor:
    """Local-only port of ``ESMEmbeddingExtractor`` (same chunking logic, LOCAL checkpoint)."""

    def __init__(self, checkpoint_path: Path, device):
        import torch
        from esm import pretrained

        # ``pretrained.load_model_and_alphabet`` takes the local branch
        # (``load_model_and_alphabet_local``, only ``torch.load()`` on
        # disk) as soon as the name ends in '.pt' -- same mechanism already
        # verified in ``_deepptmpred_runner.py`` for ESM-2. Requires the
        # companion file ``<checkpoint>-contact-regression.pt`` in the same
        # directory (fair-esm's internal heuristic, does not exclude
        # esm1b_t33_650M_UR50S), see README.md - Installation section.
        #
        # Verified by running the real load (not just reading the code):
        # ``esm/pretrained.py::load_model_and_alphabet_local``
        # calls ``torch.load(path, map_location="cpu")`` WITHOUT
        # ``weights_only=False`` -- since PyTorch 2.6 the default for
        # ``weights_only`` changed to True, and fair-esm's checkpoint
        # (an ``argparse.Namespace`` in its serialized state) does not pass
        # the strict allowlist. fair-esm (a pip package, not vendored in
        # this repo) was never updated for the new default. Fix scoped to
        # this class's choke point (the runner's only place that loads the
        # checkpoint), not a global ``torch.load`` patch for the whole
        # process -- same philosophy as the monkeypatch in
        # ``_deepptmpred_runner.py::_load_predict_module``. Safe because the
        # checkpoint comes from Meta's official URL
        # (dl.fbaipublicfiles.com), not an untrusted source.
        _original_torch_load = torch.load

        def _torch_load_weights_only_false(*args, **kwargs):
            kwargs.setdefault("weights_only", False)
            return _original_torch_load(*args, **kwargs)

        torch.load = _torch_load_weights_only_false
        try:
            self.model, alphabet = pretrained.load_model_and_alphabet(str(checkpoint_path))
        finally:
            torch.load = _original_torch_load
        self.model.eval()
        self.model = self.model.to(device)
        self.device = device
        self.batch_converter = alphabet.get_batch_converter()
        self.max_input_len = 1022
        self.max_step_len = 511

    def extract(self, seqs):
        """Identical to the original script's ``ESMEmbeddingExtractor.extract`` (see the module docstring)."""
        import torch

        max_seq_len = len(seqs[0])
        representations = []
        with torch.no_grad():
            for i in range(0, max_seq_len, self.max_step_len):
                j = min(i + self.max_input_len, max_seq_len)
                delta = 0 if i == 0 else self.max_input_len - self.max_step_len
                if i > 0 and j < i + self.max_input_len:
                    delta += i + self.max_input_len - max_seq_len
                    i = max_seq_len - self.max_input_len

                batch_seqs = [("", seq[i:j]) for seq in seqs]
                _, _, batch_tokens = self.batch_converter(batch_seqs)
                batch_tokens = batch_tokens.to(self.device)
                results = self.model(batch_tokens, repr_layers=[_ESM_EMBED_LAYER])
                representations.append(results["representations"][_ESM_EMBED_LAYER][..., delta : j - i + 1, :])

                if j >= max_seq_len:
                    break

        return torch.cat(representations, dim=1)[0]


def _sequence_cache_key(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()[:12]


def _load_or_compute_esm_full(sequence: str, extractor: "_ESMEmbeddingExtractor", cache_dir: Path, accession: str):
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{accession}_{_sequence_cache_key(sequence)}_esm1b_full.npy"
    if cache_path.is_file():
        return np.load(cache_path)
    full = extractor.extract([sequence]).cpu().numpy()
    np.save(cache_path, full)
    return full


def _load_or_compute_structure_emb(pdb_path: Path, mif_weights: Path, cache_dir: Path, accession: str) -> np.ndarray:
    """Full ``structure_emb`` (256-dim per residue), indexed by PDB RESIDUE NUMBER - 1.

    Cached by accession + hash of the PDB's content (not just the
    accession): a different PDB under the same accession never silently
    reuses the old embedding -- same reasoning as
    ``_deepptmpred_runner.py``'s cache.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    pdb_bytes = pdb_path.read_bytes()
    pdb_hash = hashlib.sha256(pdb_bytes).hexdigest()[:12]
    cache_path = cache_dir / f"{accession}_{pdb_hash}_mif_structure.npy"
    if cache_path.is_file():
        return np.load(cache_path)

    import torch
    from MIF.sequence_models.pdb_utils import parse_PDB, process_coords
    from MIF.sequence_models.pretrained import load_model_and_alphabet

    model, collater = load_model_and_alphabet(str(mif_weights))
    model.eval()

    coords, wt, _valid_resn = parse_PDB(str(pdb_path))
    coords_dict = {"N": coords[:, 0], "CA": coords[:, 1], "C": coords[:, 2]}
    dist, omega, theta, phi = process_coords(coords_dict)

    batch = [[
        wt,
        torch.tensor(dist, dtype=torch.float),
        torch.tensor(omega, dtype=torch.float),
        torch.tensor(theta, dtype=torch.float),
        torch.tensor(phi, dtype=torch.float),
    ]]
    src, nodes, edges, connections, edge_mask = collater(batch)
    with torch.no_grad():
        outputs = model(src, nodes, edges, connections, edge_mask)

    structure_emb = outputs[0][: len(wt)].detach().cpu().numpy()
    np.save(cache_path, structure_emb)
    return structure_emb


def _fasta_position_to_pdb_seqid(position_mapping_csv: Path) -> dict:
    df = pd.read_csv(position_mapping_csv)
    return dict(zip(df["fasta_position"].astype(int), df["pdb_seqid"].astype(int)))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Standalone EMNGly runner (n_linked_glycosylation consensus, PDB path)."
    )
    parser.add_argument("--emngly-home", required=True, help="Path to the root of the EMNgly clone")
    parser.add_argument("--accession", required=True)
    parser.add_argument("--sequence", required=True, help="Full ATMSEQ sequence (Phase 1.5)")
    parser.add_argument("--pdb-path", required=True, help="Single-chain PDB (record.chain_pdb_path)")
    parser.add_argument(
        "--position-mapping-csv", required=True,
        help="Phase 1.5 CSV (fasta_position <-> pdb_seqid) to correctly align structure_emb.",
    )
    parser.add_argument("--positions", required=True, type=int, nargs="+", help="1-based positions (sequon Asn).")
    parser.add_argument("--mif-weights", required=True)
    parser.add_argument("--esm-checkpoint", required=True)
    parser.add_argument("--svm-checkpoint", required=True)
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--out-csv", required=True)
    args = parser.parse_args()

    # model/MIF/__init__.py makes 'MIF' a valid top-level package by
    # adding 'EMNgly/model' to sys.path (used by the 'from
    # MIF.sequence_models...' import below). Also -- verified by running
    # the real import, not just reading the code -- 'MIF/sequence_models/*.py'
    # (pdb_utils.py, pretrained.py, etc.) do BARE imports of
    # 'sequence_models.xxx' (neither relative nor 'MIF.sequence_models.xxx'),
    # same as the original script 'model/get_mif_embedding.py' which,
    # besides 'sys.path.append("./MIF")', also depends on its own
    # directory ('EMNgly/model') already being in sys.path (automatic
    # interpreter behavior when running a script, not present when
    # importing this runner via subprocess) -- that is why
    # 'EMNgly/model/MIF' is ALSO needed in sys.path, not just
    # 'EMNgly/model'.
    sys.path.insert(0, str(Path(args.emngly_home) / "model"))
    sys.path.insert(0, str(Path(args.emngly_home) / "model" / "MIF"))

    import torch

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    cache_dir = Path(args.cache_dir)

    seq_length = len(args.sequence)
    valid_positions = [p for p in args.positions if 1 <= p <= seq_length]
    skipped = sorted(set(args.positions) - set(valid_positions))
    if skipped:
        print(
            f"[emngly_runner] {len(skipped)} position(s) out of range (sequence of "
            f"{seq_length} residues), skipped: {skipped}",
            file=sys.stderr,
        )
    if not valid_positions:
        pd.DataFrame(columns=OUTPUT_COLUMNS).to_csv(args.out_csv, index=False)
        return 0

    fasta_to_pdb_seqid = _fasta_position_to_pdb_seqid(Path(args.position_mapping_csv))
    structure_emb = _load_or_compute_structure_emb(
        Path(args.pdb_path), Path(args.mif_weights), cache_dir, args.accession
    )

    extractor = _ESMEmbeddingExtractor(Path(args.esm_checkpoint), device)
    site_full = _load_or_compute_esm_full(args.sequence, extractor, cache_dir, args.accession)

    with open(args.svm_checkpoint, "rb") as f:
        svm = pickle.load(f)

    rows = []
    skipped_alignment = []
    for pos in valid_positions:
        pdb_seqid = fasta_to_pdb_seqid.get(pos)
        structure_idx = pdb_seqid - 1 if pdb_seqid is not None else None
        if (
            pdb_seqid is None
            or structure_idx < 0
            or structure_idx >= structure_emb.shape[0]
            or pos >= site_full.shape[0]
        ):
            skipped_alignment.append(pos)
            continue

        site_vec = site_full[pos]
        left = max(0, pos - _LOCAL_WINDOW_LEFT)
        right = min(pos + _LOCAL_WINDOW_RIGHT, seq_length)
        local_window = args.sequence[left:right]
        local_vec = extractor.extract([local_window]).cpu().numpy()[0]
        struct_vec = structure_emb[structure_idx]

        feature_row = np.concatenate([site_vec, local_vec, struct_vec]).reshape(1, -1)
        probability = float(svm.predict_proba(feature_row)[0, 1])
        rows.append({"position": pos, "probability": probability})

    if skipped_alignment:
        print(
            f"[emngly_runner] {len(skipped_alignment)} position(s) with no valid PDB mapping or out "
            f"of range for structure_emb/site_emb, skipped: {skipped_alignment}",
            file=sys.stderr,
        )

    pd.DataFrame(rows, columns=OUTPUT_COLUMNS).to_csv(args.out_csv, index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
