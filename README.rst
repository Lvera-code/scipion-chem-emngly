================================
EMNGly Scipion plugin
================================

Scipion framework plugin wrapping EMNGly (PMC10627407 (CC BY 4.0)) --
second consensus engine for structural N-glycosylation
(n_linked_glycosylation, PDB path), alongside DeepMVP.

``ProtEMNGlyPrediction`` corroborates N-glycosylation candidates from
``scipion-chem-deepmvp`` with a vendorized runner, a maintained
byte-for-byte copy of the upstream inference script.

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

The repo (with the MIF weights bundled), the conda environment (installed
from the repo's own ``environment.yml``), the ESM-1b checkpoint
(``esm1b_t33_650M_UR50S.pt`` + its ``-contact-regression.pt`` companion) and
the trained SVM (``N-GlyDE.pickle``) are all installed automatically, into
``<EMNGLY_HOME>/checkpoints/``. Only set ``EMNGLY_ESM_CHECKPOINT``/
``EMNGLY_SVM_CHECKPOINT`` in ``scipion.conf`` if you want to point at
different files instead.

.. code-block::

            scipion3 tests emngly.tests
