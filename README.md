# PESOL Package for Power Systems Analysis 
# Install the corresponding folders in the latest version of the Matpower for running PESOL Toolbox. 

The PESOL toolbox contains following blocks,

1) Load Flow Solver (LFS) with optimal step size strategy (Rectangular coordinates formulation).
2) Continuation Power Flow (CPF) solver for tracing power flow solution curves (PV-Curves).
3) Direct-Methods (DMs) for calculating loadability and transfer capability limits based on Transversality Enforced Newton Raphson (TENR) Algorithm [1] .
  * Without technical constraints.
  * With voltage magnitude constraints on load buses [2]. 
5) An adaptive Homotopy Continuation (HC) method for tracing power flow solution boundaries (1-manifold curves) [3].

Kindly, readme.text file in each of the folder to run these packages. 
The input datacase files are standard case files from MATPOWER.

For any inquiry please drop me an email:  mazhar.ali@skolkovotech.ru


You may cite following papers for using the PESOL Toolbox for power systems analysis.

* [1]Ali, Mazhar, Anatoly Dymarsky, and Konstantin Turitsyn. "Transversality enforced Newton–Raphson algorithm for fast calculation of maximum loadability." IET Generation, Transmission & Distribution 12.8 (2018): 1729-1737.

* [2] Ali, Mazhar, Elena Gryazina, and Konstantin S. Turitsyn. "Fast calculation of the transfer capability margins." 2019 IEEE Milan PowerTech. IEEE, 2019.

* [3] Ali, Mazhar, et al. "Voltage Feasibility Boundaries for Power System Security Assessment." arXiv preprint arXiv:2103.00168 (2021).
