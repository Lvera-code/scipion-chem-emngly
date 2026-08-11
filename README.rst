================================
EMNGly Scipion plugin
================================

Scipion framework plugin wrapping EMNGly (PMC10627407 (CC BY 4.0)) --
segundo motor de consenso para N-glicosilacion estructural (n_linked_glycosylation, Camino PDB), junto con DeepMVP.

**Estado: protocolo real implementado, pendiente de instalacion+test real**
(ver ``PTM-Prediction/STATUS.md``, entrada 2026-08-11). ``ProtEMNGlyPrediction``
corrobora candidatos de N-glicosilacion de ``scipion-chem-deepmvp`` con un
runner vendorizado (identico al ya validado end-to-end en el pipeline
standalone).

Repo original: https://github.com/StellaHxy/EMNgly

Cita: PMC10627407 (CC BY 4.0)

**Licencia de EMNGly (upstream)**: MIT (resuelto 2026-08-11). Los autores de correspondencia (Yaojun Wang, Shiwei Sun) anadieron el LICENSE MIT al repo tras el correo del 2026-08-10 y autorizaron explicitamente la integracion en Scipion. Paper CC BY 4.0 (PMC10627407).

===================
Install this plugin
===================

**Developer's version**

.. code-block::

            git clone https://github.com/Lvera-code/scipion-chem-emngly.git
            cd scipion-chem-emngly
            scipion3 installp -p . --devel
            scipion3 installb EMNGly

El repo (con los pesos MIF bundled) y el entorno conda se instalan
automaticamente. DOS piezas quedan **manuales**:

- Checkpoint ESM-1b (``esm1b_t33_650M_UR50S.pt``) + su companero obligatorio
  (``esm1b_t33_650M_UR50S-contact-regression.pt``) en el MISMO directorio --
  apunta ``EMNGLY_ESM_CHECKPOINT`` (en ``scipion.conf``) al primero.
- SVM entrenado (``N-GlyDE.pickle``, descarga desde Google Drive, ver el
  README real de EMNgly) -- apunta ``EMNGLY_SVM_CHECKPOINT`` a el.

.. code-block::

            scipion3 tests emngly.tests
