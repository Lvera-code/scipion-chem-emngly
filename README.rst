================================
EMNGly Scipion plugin
================================

Scipion framework plugin wrapping EMNGly (PMC10627407 (CC BY 4.0)) --
second consensus engine for structural N-glycosylation
(n_linked_glycosylation, PDB path), alongside DeepMVP.

``ProtEMNGlyPrediction`` corroborates N-glycosylation candidates from
``scipion-chem-deepmvp`` with a vendorized runner (identical to the one
already validated end-to-end in the standalone pipeline).

Original repo: https://github.com/StellaHxy/EMNgly

Citation: PMC10627407 (CC BY 4.0)

**EMNGly license (upstream)**: MIT. The original repository did not
include a LICENSE at integration time; the MIT file was added afterward
and the Scipion integration was explicitly authorized. Paper CC BY 4.0
(PMC10627407).

===================
Install this plugin
===================

**Developer's version**

.. code-block::

            git clone https://github.com/Lvera-code/scipion-chem-emngly.git
            cd scipion-chem-emngly
            scipion3 installp -p . --devel
            scipion3 installb EMNGly

The repo (with the MIF weights bundled) and the conda environment are
installed automatically. TWO pieces remain **manual**:

- ESM-1b checkpoint (``esm1b_t33_650M_UR50S.pt``) + its required companion
  (``esm1b_t33_650M_UR50S-contact-regression.pt``) in the SAME directory --
  point ``EMNGLY_ESM_CHECKPOINT`` (in ``scipion.conf``) at the former.
- Trained SVM (``N-GlyDE.pickle``, downloaded from Google Drive, see
  EMNgly's real README) -- point ``EMNGLY_SVM_CHECKPOINT`` at it.

.. code-block::

            scipion3 tests emngly.tests
