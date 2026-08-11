=========
CHANGES
=========

0.2.0
=====
- Real protocol (``ProtEMNGlyPrediction``): corroborates N-glycosylation
  candidates from scipion-chem-deepmvp (``_scoreEmngly``, None on
  non-applicable rows). Runner (ESM-1b+MIF->SVM) vendorized byte-for-byte
  from the standalone project. Automatic installation of the repo (MIF
  bundled) + conda environment (CPU-only torch + nvidia/triton purge, same
  real fix as scipion-chem-stackglyembed); ESM-1b checkpoint and SVM
  manual. Real test on 7c4s (same fixture as
  scipion-chem-discotope/-deepptmpred).

0.1.0
=====
- Initial scaffolding: Scipion plugin structure generated following the
  same pattern as the BCell-Epitope-Prediction plugins (one plugin per
  tool). No installation or protocol logic yet -- pending end-to-end
  validation of the pipeline in Colab, see the ``PTM-Prediction`` project's
  STATUS.md.
