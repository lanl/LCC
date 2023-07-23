Kinetic Monte Carlo Crystal Growth Tool
===

About
-----
This folder contains a Python script ``growth_kmc.py`` which performs a Kinetic
Monte Carlo simulation of crystal growth from solution, starting with a spherical
seed. Currently, the code supports primative tetragonal lattices, with the intent
to generalize this. The core developer of this tool is Jacob Jeffries (jwjeffr@clemson.edu
or jwjeffr@lanl.gov).

Requirements
------------

-   The :code:`python` interpreter.

-   The external :code:`python` packages :code:`numpy` and code:`numba`.

The example input (:code:`example_input.json`) provided works for Python 3.9.12,
Numpy 1.21.6, and Numba 0.55.1. Other versions are not guaranteed to be functional.

Testing and running the code
----------------------------

The code can be tested with:

  python growth_kmc.py example_input.json

or:

  ./growth_kmc.py example_input.json

Two runs will be performed:

-   A short, small run which first compiles functions. This run data will be stored in
    :code:`small.dump` in the LAMMPS-style dump format.

-   A longer run with parameters provided in :code:`example_input.json`. The parameters
    are:

    Box dimensions = (30, 30, 70) (in lattice units)
    Number of steps = 100,000
    Dump every = 500 steps
    Dump file name = petn_growth.dump
    Initial seed radius = 75.0 angstroms
    Temperature = 300.0 kelvin
    First-neighbor cutoff distance = 7.0 angstroms
    Second-neighbor cutoff distance = 7.5 angstroms
    First-neighbor interaction energy = -0.291 electron volts
    Second-neighbor interaction energy = -0.186 electron volts

For new parameters, simply change the dictionary written in :code:`example_input.json` to
match your desired parameters.

Notes
-----

It's assumed that |a| = |b| = 9.088 angstroms and |c| = 6.737 angstroms, and
that the offset for the second molecule in the unit cell is (4.54348, 4.54346, 3.36908)
angstroms. These are specified in the :code:`initialize_simulation()` function starting at
line 134. These parameters are optimized for a PETN crystal, but are easy to change.
(TODO allow the user to specify the lattice in the input file. Create lattice classes).

This code is highly parallelized, and will use all available cores. If cores
are currently being used, your system might crash.
(TODO allow the user to specify the number of cores to use in the input file. Default
to all cores if not specified)

Adsorption and evaporation rates are both assumed to be Arrhenius with a prefactor of
1e+10 Hz. Adsorption rates are constant everywhere, with an activation energy of 0.9
electron volts, while evaporation rates are assumed to be governed by the energy change
when bonds are broken.
(TODO allow the user to modify the prefactors and the adsorption activation energy, or
come up with a better adsorption rate expression)
