=========
CHANGES
=========

0.2.0
=====
- Real protocol (``ProtEMNGlyPrediction``): corroborates N-glycosylation
  candidates from scipion-chem-deepmvp (``_scoreEmngly``, None on
  non-applicable rows). Runner (ESM-1b+MIF->SVM) vendorized byte-for-byte
  from the upstream inference script. Automatic installation of the repo
  (MIF bundled) + conda environment (CPU-only torch + nvidia/triton purge,
  same real fix as scipion-chem-stackglyembed); ESM-1b checkpoint and SVM
  manual. Real test on 7c4s (same fixture as
  scipion-chem-discotope/-deepptmpred).

0.1.0
=====
- Initial scaffolding: Scipion plugin structure generated following the
  same one-plugin-per-tool pattern used across this project's other
  plugins. No installation or protocol logic yet.
