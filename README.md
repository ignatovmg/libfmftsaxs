# libfmftsaxs
##### Provides tools for manifold-FFT SAXS calculation.

---

## Requirements
* libmol2
* gmp
* fftw3
* Check
* Doxygen
* CMake (version > 3.0)

## Installation
After cloning the repository, do
```bash
$ cd libfmftsaxs
$ mkdir build && cd build
$ cmake -DBUILD_TESTS=True -DBUILD_DOCS=True -DMPI=False -DCMAKE_BUILD_TYPE=Release  ..
$ make
$ make test && make doc
$ make install
```
If you are planning to use MPI, you can optionally add `-DMPI=True`. To add verbosity `-DVERBOSE=True`.

---
## Capabilities
This is a library for calculation of SAXS on a single protein molecule as well as on a dimer. It utilizes expansion of scattering amplitudes into spherical harmonics to efficiently compute scattering intensity as described  [here](https://www.researchgate.net/publication/244634749_CRYSOL-_a_Program_to_Evaluate_X-ray_Solution_Scattering_of_Biological_Macromolecules_from_Atomic_Coordinates). The novelty and the main feature of the library, is that it provides functions for SAXS calculation on the heterodimers using FFT (Fast Fourier Transform),  representing the dimer's scattering intensity as a correlation function, which allows to score millions of dimer conformations very fast.

---
## HOWTO
There are several tools and examples included in the directory `examples`. Each of the scripts there must be launched from this directory. They do not require you to install the library, just to build. 
Most scripts generally start from specifing mapping file and parameter file, where the former is needed to create the form-factor table and the latter to add parameters of the atoms. They reside in `prms` directory.

1. Computing SAXS profile

   `run_single_saxs.sh` provides you with an example of SAXS profile calculation. It uses files stored in `4g9s`. You can specify **c1**, **c2** and **lmax** - expansion depth (see documentation for the details). After the script is run, the profile is written into `4g9s_saxs_profile` and you can plot it using 
   ```bash
   $ python plot_curves.py 4g9s_saxs_profile
   ```
   
2. Scoring SAXS profiles
   `run_silly_scoring.sh` gives you an example of how to score dimer conformations using `score_ft_silly` binary. The conformations are written in the form of ft-file and respective file with rotation matrices. Each conformation is written, transformed into Euler coordinates, then SAXS profile is created for this conformation and minimized through L-BFGS-B algorithm with respect to the parameters **c1**, **c2**, provided the reference (experimental) profile. You can specify the output, where the scores will be written. They are written the form of a table with the following columns: 
| FT index | SAXS-score | c1 | c2 |.

3. Converting FT-file to Euler angles.
   `run_ft2euler.sh` converts ft-file from Cartesian to Euler coordinates and writes it into `euler_list` in the form:
| $z$ | $\beta_{rec}$ | $\gamma_{rec}$ | $\alpha_{lig}$ | $\beta_{lig}$ | $\gamma_{lig}$ |. This means, that receptor is rotated by $(0.0, \beta_{rec}, \gamma_{rec})$ and ligand is rotated by $(\alpha_{lig}, \beta_{lig}, \gamma_{lig})$ and translated along z-axis by $z$.

4. Ultra-fast FFT-SAXS scoring.
   `run_correlate.sh` demonstrates the main feature of this library. It uses the protein complex `1a2k`. First SAXS curve is calculated for the native conformation of the dimer for a chosen pair of parameters **c1**, **c2** (1.0, 1.0 by default). Then we pretend, that this curve is experimental and feed it to the executable `correlate` together with the 3 ft-files, which in total contain 210000 conformations. The ft-files are concatenated and fed to `correlate` as a single list. It uses FFT to score the conformations in super-fast fashion and provides a happy user with the output of the form | Serial number | FT index | SAXS-score | c1 | c2 |, which is then sorted by the serial number (line index in the concatenated list) and split into three lists with SAXS-scores, each of which corresponds to the one of the ft-files. 
By default the utility `correlate` is launched with `mpirun`, but you can use the serial version, if you did not specify MPI flag during the compilation (which will make it proportinally slower).

You can `plot_curves.py` to build as many curves as you want of a single plot just like this: 
```bash 
$ python plot_curves.py curve1 curve2 ... 
```
It also computes $\chi$-score between every pair of them.

---
## License ##
See License.TXT for details.

---
## Authors ##
Mikhail Ignatov, Andrey Kazennov and Dima Kozakov.

---
## References ##
1. Svergun, D., C. Barberato, and M. H. J. Koch. "CRYSOL - a Program to Evaluate X-ray Solution Scattering of Biological Macromolecules from Atomic Coordinates." Journal of Applied Crystallography 28.6 (1995): 768-73. Web.
2. Konarev, P. V., M. V. Petoukhov, and D. I. Svergun. "MASSHA - a Graphics System for Rigid-body Modelling of Macromolecular Complexes against Solution Scattering Data." Journal of Applied Crystallography 34.4 (2001): 527-32. Web.
3. Xia, Bing, Artem Mamonov, Seppe Leysen, Karen N. Allen, Sergei V. Strelkov, Ioannis Ch. Paschalidis, Sandor Vajda, and Dima Kozakov. "Accounting for Observed Small Angle X-ray Scattering Profile in the Protein-protein Docking Server Cluspro." Journal of Computational Chemistry 36.20 (2015): 1568-572. Web.


