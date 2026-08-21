=========
CHANGES
=========

0.4.1
=====
- Reverted the 0.3.0 Python bump (3.10 -> 3.11) found broken via an
  actual end-to-end fresh install on a Colab GPU session (Tesla T4,
  2026-08-21): ``scikit-learn==1.1.1`` (the pinned version required for
  the SVM pickle's ABI, see 0.3.0) has NO prebuilt wheel for Python 3.11
  (confirmed via ``pip download --python-version 311``, real "No
  matching distribution" error) and fails to build from source. Back to
  Python 3.10 (the file's own ``python=3.11.0`` pin is now also filtered
  out, not just installed-over). Verified after the fix: full env
  creation + real GPU torch install + all 3 pinned packages install
  cleanly, ``torch.cuda.is_available()`` True, SVM pickle unpickles
  correctly (``sklearn.svm._classes.SVC``, ``n_features_in_=2816``,
  matching ESM-1b+MIF+structural dims exactly).

0.4.0
=====
- GPU support: ``USE_GPU``/``GPU_LIST`` hidden params added to
  ``ProtEMNGlyPrediction``, wired to ``CUDA_VISIBLE_DEVICES`` in
  ``runEMNGly`` (the runner decides GPU/CPU itself via
  ``torch.cuda.is_available()``, no native CLI flag). Torch install +
  nvidia/triton purge are now GPU-conditional: default (CUDA-capable)
  wheel and no purge when a GPU is detected; without one (this dev
  machine's case, the only branch verified here) stays exactly the
  already-verified CPU-only-wheel + purge behavior. The
  ``CUDA_VISIBLE_DEVICES`` lever itself was verified for real against
  torch on a Colab GPU session (Tesla T4): ``torch.cuda.is_available()``
  flips False/True exactly as expected.

0.3.0
=====
- Installed from the repo's own real ``environment.yml`` (Python bumped to
  the file's real 3.11.0; GPU/pickle-incompatible conda entries filtered
  out, kept overriding numpy/pandas/scikit-learn to the already
  production-validated versions). ESM-1b checkpoint and the SVM
  (``N-GlyDE.pickle``) now auto-downloaded at install time into
  ``<EMNGLY_HOME>/checkpoints/``. Removed unused ``READ_URL`` constant.
  Test file split into per-behavior methods instead of one ``setUpClass``
  blob.

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
