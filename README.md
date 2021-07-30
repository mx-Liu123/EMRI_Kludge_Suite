# Augmented Analytic Kludge Waveform With Quadrupole Moment Correction

This is a C/C++ code that computes quadrapole included augmented analytic kludge (QAAK) [1] waveforms for extreme-mass-ratio inspirals (EMRIs), where the system could possess arbitrary quadrapole moment deviation from Kerr value. The code is modified from EMRI Kludge Suite (version 0.5.0) [2,3,4] https://github.com/alvincjk/EMRI_kludge_Suite. And the orbital frequency deviations caused by arbitrary quadrupole is evaluated by the method in [5].

&mdash; Miaoxin Liu, July 2021

## Installation

The Installation is same to EMRI Kludge Suite: The GSL and FFTW libraries are required for compilation. Clean up any previous installation first with `make clean`. Running `make` will create the following executables in the folder `./bin`:

- `AAK_Waveform`
- `AK_Waveform`
- `NK_Waveform`
- `AAK_TDI`
- `AK_TDI`
- `AAK_Phase`

Note that the current version has removed the usage of AAK TDIs, AAK phases, Python interface of EMRI Kludge Suite.

## Usage

The file `./examples/SetPar_Waveform` is a template for a formatted settings/parameters file that contains the output file path, various toggles, and the EMRI parameters. More details are provided in the template file itself. In QAAK Mod, the quadrupole deviation from Kerr value q_q is included as an additional parameter.

For example, running `bin/AAK_Waveform examples/SetPar_Waveform` will generate a QAAK waveform with default settings and parameters. Three files will be created in `./bin`:

- `example_wave.dat` contains waveform data (t, h_I, h_II)
- `example_traj.dat` contains inspiral trajectory data (t, p/M, e, iota, E, L_z, Carter Constant: K), where quadrupole correction do not affect the output (p,e,iota)->(E,Lz,K) since K is not well defined in a spacetime with arbitrary quadrupole moment. 
- `example_info.txt` contains additional information such as signal-to-noise ratio and waveform timing

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

## List of (important) known bugs

Known bugs of AAK:
- The approximate LISA response functions h_I/h_II for the AAK/AK and the NK do not match up, likely due to differing conventions when implementing the Doppler shift
- The NK may produce NaNs at times that coincide with specific fractions of the LISA orbital period (when using integer values of dt), or near plunge; this needs to be fixed, but in the meantime a band-aid solution for dealing with isolated NaNs is included
- All waveforms may have difficulties with zero or epsilon values for certain parameters such as spin, eccentricity and inclination

Bugs of QAAK:
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

[1] M. Liu & J. Zhang. Augmented analytic kludge waveform with quadrupole moment correction. e-Print: 2008.11396.

[2] A. J. K. Chua & J. R. Gair. Improved analytic extreme-mass-ratio inspiral model for scoping out eLISA data analysis. *Class. Quantum Grav.* 32:232002, 2015.

[3] A. J. K. Chua, C. J. Moore & J. R. Gair. Augmented kludge waveforms for detecting extreme-mass-ratio inspirals. *Physical Review D* 96:044005, 2017.

[4] S. Babak, H. Fang, J. R. Gair, K. Glampedakis & S. A. Hughes. "Kludge" gravitational waveforms for a test-body orbiting a Kerr black hole. *Physical Review D* 75:024005, 2007.

[5] S. J. Vigeland & S. A. Hughes. Spacetime and orbits of bumpy black holes. *Physical Review D* 81:024030, 2010.
