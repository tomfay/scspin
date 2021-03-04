# Example 1: An Anisotropic SC Simulation
Run this example with the command:
```
scspin input.json
```
This will produce a set of checkpoint files `checkpoint-n.dat` with n=1-10, the `anisosimout.dat` file and the `anisoyieldsout.dat` file.

## Input File Contents
Let us examine the input file `input.json`. This is in the `.json` file format, so I ad vise you look up the basics of this file format if anything doesn't make sense. [An intro to json can be found here.](https://www.tutorialspoint.com/json/json_quick_guide.htm) In my personal experience the biggest sources of errors in json inputs are spurious commas and brackets.
```
{
  "Anisotropic Simulation":{
    "Alt Norm": false,
    "Samples": 10000 ,
    "Step Size": 0.05 ,
    "Time Steps": 800 ,
    "System 1": {
      "Radicals": [
        {
          "Index": 0,
          "Omega": [0.0,0.0,1.0]
        },
        {
          "Index": 1,
          "Omega": [0.0,0.0,1.0],
          "Hyperfines": {
            "Number": 2,
            "Tensors": [
              [-0.5, 0.0, 0.0, 0.0, -0.5, 0.0, 0.0, 0.0, 1.0 ],
              [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 ]
            ] ,
            "Multiplicities": [3,2]
          }
        }
      ],
      "Rate Constants": {
        "Singlet Rate": 0.05,
        "Triplet Rate": 0.0,
        "Dephasing Rate": 0.0,
        "Spin Relaxation Rates": [0.0,0.0,0.0,0.0,0.0,0.0]
      },
      "Electron Couplings": {
        "Tensor": [0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 1.0 ]
      }
    }
  }
}
```

First the `{"Anisotropic Simulation": { ... } }` part specifies that the code will run a SC simulation of a radical pair with anisotropic coupling tensors. There should only be one "Simulation" block in each input.

Next we see some simulation parameters are specified. 
```
"Alt Norm": false,
"Samples": 10000 ,
"Step Size": 0.05 ,
"Time Steps": 800 ,
```
`"Alt Norm"` specifies how the electron spins are normalised in the nuclear spin evolution. "false" should generally be used, and this is the default.\
`"Samples"` specifies the number of Monte-Carlo samples used in the simulation. More samples will give more accurate, more converged results, but the simulation time scales linearly with the number of samples. Generally accurate results require at least 100000 samples, normally 1000000 samples.\
`"Step Size"` specifies the time step for the evolution with a Cayley type integrator for the semiclassical equations of motion. `scspin` works in arbitrary units. _You_ have to make sure the units are consistent and _you_ have to check that the time step is small enough that the equations are accurately integrated. A good rule-of-thumb is that the time step should be about 1/20 times the shortest time-scale in the problem.\
`"Time Steps"` specifies the number of integration steps performed, so the total simulation time is "Step Size" * "Time Steps".


Next we have the `"System 1": { ... }` block. Only "System 1" is needed and used for an "Anisotropic Simulation". This block specififies parameters for each of the radicals, and the electron spin dynamics. 

The `"Radicals:[...]"` specifies the parameters for each radical. Each radical must have an `"Index"` of 0 or 1. `"Omega"` is a three-element vector that specifies the x,y and z components of the Zeeman frequency for that radical electron spin. If the `"Hyperfines"` block is included, nuclear spins will also be included for that radical. The number of such couplings must be set by the `"Number"` block and the coupling tensors are specified in the `"Tensors"` block. The `"Tensors` block contains the hyperfine coupling tensors in a flattened form [_A<sub>xx</sub>, A<sub>xy</sub>, A<sub>xz</sub>, A<sub>yx</sub>, A<sub>yy</sub>, A<sub>yz</sub>, A<sub>zx</sub>, A<sub>zy</sub>, A<sub>zz</sub>_]. Finally the `"Multiplicities"` block specifies the spin multiplicities of the nuclear spins, i.e. each entry is 2 _I<sub>i,k</sub>_+1 for that spin.

The `"Rate Constants"` block specifies rate constants for the decay processes of the electron spins.\
`"Singlet Rate"` and `"Triplet Rate"` specify the singlet _k_<sub>S</sub> and triplet _k_<sub>T</sub>first order spin selective decay rates.\
`"Dephasing Rate"` specifies the electron spin singlet-triplet dephasing rate _k_<sub>STD</sub>.\
`"Spin Relaxation rates"` specifies six random fields type relaxation rates [_R_<sub>0,x</sub>,_R_<sub>0,y</sub>,_R_<sub>0,z</sub>,_R_<sub>1,x</sub>,_R_<sub>1,y</sub>,_R_<sub>1,z</sub>]. Each rate corresponds to a Random fields relaxation rate for a relaxation superoperator of the form L[_O_] = _R_<sub>i,a</sub> (_S_<sub>i,a</sub>._O_._S_<sub>i,a</sub> - (1/4)._O_). So for example the decay rate of _S_<sub>i,x</sub> is (_R_<sub>i,y</sub>+_R_<sub>i,z</sub>)/2. 

Finally the `"Electron Couplings"` block specifies the electron spin coupling. The `"Tensor"` block is all that is required, where the coupling tensor is specified such that the electron spin coupling term is **_S_**<sub>0</sub> .**T**.**_S_**<sub>1</sub>, where **T** is the tensor supplied.

## Output File Contents

The important output files are `anisosimout.dat` and `anisoyieldsout.dat`. 

`anisosimout.dat` contains five columns:
time, <P<sub>S</sub>(t)P<sub>S</sub>(0)>, <P<sub>T</sub>(t)P<sub>S</sub>(0)>, <P<sub>S</sub>(t)P<sub>T</sub>(0)>, <P<sub>T</sub>(t)P<sub>T</sub>(0)>

`anisoyieldsout.dat` contains quantum yields of the singlet and triplet reactions the first four numbers are: Y<sub>S</sub>(S), Y<sub>T</sub>(S)m, Y<sub>S</sub>(T) and Y<sub>T</sub>(T), where Y<sub>A</sub>(B) is the quantum yield of the A reaction for the B initial condition. The second four numbers are the <Y<sub>A</sub>(B)^2>, from which the standard error in each yield can be backed out. 
