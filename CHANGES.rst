=========
CHANGES
=========

0.2.0
=====
- Protocolo real (``ProtEMNGlyPrediction``): corrobora candidatos de
  N-glicosilacion de scipion-chem-deepmvp (``_scoreEmngly``, None en filas
  no aplicables). Runner (ESM-1b+MIF->SVM) vendorizado byte-a-byte desde el
  proyecto standalone. Instalacion automatica del repo (MIF bundled)+entorno
  conda (torch CPU-only + purga nvidia/triton, mismo fix real que
  scipion-chem-stackglyembed); checkpoint ESM-1b y SVM manuales. Test real
  sobre 7c4s (mismo fixture que scipion-chem-discotope/-deepptmpred).

0.1.0
=====
- Scaffolding inicial: estructura de plugin de Scipion generada siguiendo el
  mismo patron que los plugins de BCell-Epitope-Prediction (un plugin por
  herramienta). Sin logica de instalacion ni de protocolo todavia -- pendiente
  de la validacion end-to-end del pipeline en Colab, ver STATUS.md del
  proyecto ``PTM-Prediction``.
