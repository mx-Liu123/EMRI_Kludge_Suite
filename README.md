# Augmented Analytic Kludge Waveform With Quadrupole Moment Correction

This is a mod of the AAK waveform in EMRI Kludge Suite (version 0.5.0) [1,2,3] https://github.com/alvincjk/EMRI_kludge_Suite. It computes quadrapole included augmented analytic kludge (QAAK) [4] waveforms for extreme-mass-ratio inspirals (EMRIs), where the system could possess arbitrary quadrapole moment deviation from Kerr value. The orbital frequency deviations caused by arbitrary quadrupoles are evaluated by the method in [5].
Note that the current version has removed the usage of AAK TDIs, AAK phases, Python interface of EMRI Kludge Suite.

&mdash; Miaoxin Liu, July 2021

## QAAK Mod files
Main:
- `src/exec/AAK_Waveform.cc`  Calculate the orbital frequency deviations from Kerr value, which correct the AAK mapping.
- `src/inclecc/GKR.cc`  Quadrupole corrections in NK.
- `src/suite/KSParMap.cc`  Quadrupole corrections in 3pn AK; The quadrupole won't be mapped to an unphysical value like spin in AAK mapping.
- `src/suite/AAK.cc`  Quadrupole is included when calling the 3pn AK functions in KSParMap.cc.
- `src/suite/GKTrajFast.cc`  Quadrupole is included when calling the GKR.
- `src/suite/KSTools.cc`  The quadrupole deviation from Kerr value q_q is added as an additional input parameter.  

These header files are also modified:
- `include/AAK.h` 
- `include/GKR.h` 
- `include/GKTrajFast.h` 
- `include/KSParMap.h` 
- `include/KSTools.h` 

The following files have nothing to do with QAAK waveform, modified to avoid errors:
- `src/exec/AAK_Phase.cc` 
- `src/exec/AAK_TDI.cc` 
- `src/inclecc/GKTrajcc.cc` 
- `src/suite/AAKPhase.cc` 
- `src/suite/AAKpy.cc` 
- `src/suite/AAKTDI.cc` 

## Bugs of QAAK

- Orbital frequency deviations may be NaNs when p is too small, but it is safe for the position where they are needed to correct the mappings.

## Authors

**Miaoxin Liu**  
School of Physics and Astronomy, Sun Yat-sen University  
`liumx37@mail2.sysu.edu.cn`

**Jian-dong Zhang**  
TianQin Research Center for Gravitational Physics  
`zhangjd9@mail.sysu.edu.cn`

## Acknowledgments

We thank Alvin J. K. Chua and Yiming Hu for helpful discussion. This code makes use of the Black Hole Perturbation Toolkit.

## References

[1] A. J. K. Chua & J. R. Gair. Improved analytic extreme-mass-ratio inspiral model for scoping out eLISA data analysis. *Class. Quantum Grav.* 32:232002, 2015.

[2] A. J. K. Chua, C. J. Moore & J. R. Gair. Augmented kludge waveforms for detecting extreme-mass-ratio inspirals. *Physical Review D* 96:044005, 2017.

[3] S. Babak, H. Fang, J. R. Gair, K. Glampedakis & S. A. Hughes. "Kludge" gravitational waveforms for a test-body orbiting a Kerr black hole. *Physical Review D* 75:024005, 2007.

[4] M. Liu & J. Zhang. Augmented analytic kludge waveform with quadrupole moment correction. e-Print: 2008.11396.

[5] S. J. Vigeland & S. A. Hughes. Spacetime and orbits of bumpy black holes. *Physical Review D* 81:024030, 2010.
