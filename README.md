# Builder

This software provides a modular statistical autodesigner that uses input parameters to calculate the optimal design of components from a library of modules.

## 📦 Installation

No special installation instructions, simply clone this repository to use it. Maybe I'll add a makefile at some point if we really need it.

## 🚀 How to use

Currently only the `nozzle` module is active. To use a module, first create a data file that specifies all the inputs the system expects. For example, in the case of the `nozzle` module, the system expects an input in the following format:

```
MODULE	NOZZLE

MFLOWR	100.0	5
CCTEMP	3600	100
CCPRES	5.00	1
AMPRES	0.05	0.01
EXPRES	0.05	0.01

K	1.20
M	24
```

where `MFLOWR` corresponds to the mass flow rate out of the nozzle and it's associated uncertainty (this part is optional, to ignore uncertainties simply enter it as 0), `CCTEMP` refers to the combustion chamber temperature, `CCPRES` refers to the combustion chamber pressure, `AMPRESS` refers to ambient pressure in the test environment, `EXPRESS` refers to the pressure of the exhaust at the point that it exits the nozzle, `k` is the specific heat ratio of the exhaust mixture, and `M` is the molar mass.

Note that the datafile's header, `MODULE`, tells the system to feed this data into the module `NOZZLE`. Datafile structure will vary from module to module.

In order to generate the output dxf, run the `builder.py` main file with the following arguments:

`builder.py path/to/datafile.dat path/to/outfile.dxf`

Builder will automatically identify the module to use and pass the necessary data.

## 💻 How to Contribute

Looking to help? Currently the main things that need to be done are developing new modules and implementing a system that allows us to dynamically load those new modules into the system.

I'm hoping to implement a module for designing various gear configurations based on gear ratios (e.g. PLANETARY, 50:1).

## 📣 Acknowledgements

Builder's DXF processer is pretty much ripped straight from the SDXF library by nycresistor and some guy named Stani. Kudos.

https://github.com/nycresistor/SDXF
