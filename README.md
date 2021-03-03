# scspin
A repository for the scspin (semiclassical spin dynamics) code, a C++ code developed for performing semiclassical spin dynamics simulations of radical pair systems.

This code uses the semiclassical spin dynamics method described in [1] (with additional developments described in [2]) to approximate the spin dynamics of recombining radical pairs. See references [1] and [2] for details of the method.

**Disclaimer:** _This code has been used as a testbed for a lot of ideas, so it's a little messy. If you want specific things tidied up or specfic features, or if you find bugs, added please let me know. If you try to use any of the functionality outside of the examples given, you do so at your own risk! Consider these untested!_

## Compilation

Compilation requires an up-to-date c++ compiler, I use the GNU compilers (version 9) with OpenMP, and cmake. From the scspin directory first run:
```
cmake -S src -B build
```
This creates the makes files in a folder "build". Then change into this directory and run make.
```
cd build
make
```
This should create the `scspin` program file in the `build` directory.

## Running scspin

The `scspin` program is run with an input file in the `.json` format, and the program is called with an input file by running the command:
```
scspin scspin_input_file.json
```
Several simulation types can be run, but the main one is a semiclassical simulation with anisotropic coupling tensors. This is specified as an "Anisotropic Simulation" in the json file.

Several examples are provided in the "examples" folder.

## Acknowledgements

I gratefully acknowledge the support of David Manolopoulos and Lachlan Lindoy in the development of this code.

## References
[1] Lewis, A. M., Manolopoulos, D. E. & Hore, P. J. Asymmetric recombination and electron spin relaxation in the semiclassical theory of radical pair reactions. J. Chem. Phys. 141, 044111 (2014).

[2] Lindoy, L. P., Fay, T. P. & Manolopoulos, D. E. Quantum mechanical spin dynamics of a molecular magnetoreceptor. J. Chem. Phys. 152, 164107 (2020).
