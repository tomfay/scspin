# Example 1: An Anisotropic SC Simulation
Run this example with the command:
```
scspin input.json
```
This will produce a set of checkpoint files `checkpoint-n.dat` with n=1-10, the `anisosimout.dat` file and the `anisoyieldsout.dat` file.

## Input File Contents
Let us examine the input file `input.json`. This is int he `.json` file format, so I advise you look up the basics of this file format if anything doesn't make sense. [An intro to json can be found here.](https://www.tutorialspoint.com/json/json_quick_guide.htm)
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
"Alt Norm" specifies how the electron spins are normalised in the nuclear spin evolution. "false" should generally be used, and this is the default.
"Samples" specifies the number of Monte-Carlo samples used in the simulation. More samples will give more accurate, more converged results, but the simulation time scales linearly with the number of samples. Generally accurate results require at least 100000 samples, normally 1000000 samples.
"Step Size" specifies the time step for the evolution with a Cayley type integrator for the semiclassical equations of motion. `scspin` works in arbitrary units. _You_ have to make sure the units are consistent and _you_ have to check that the time step is small enough that the equations are accurately integrated. A good rule-of-thumb is that the time step should be about 1/20 times the shortest time-scale in the problem.
"Time Steps" specifies the number of integration steps performed, so the total simulation time is "Step Size" * "Time Steps".

Next we have the `"System 1": { ... }` block. Only "System 1" is needed and used for an "Anisotropic Simulation". This block specififies parameters for each of the radicals, and the electron spin dynamics.




## Output File Contents
