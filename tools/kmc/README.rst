.. _PETN: https://en.wikipedia.org/wiki/Pentaerythritol_tetranitrate

.. role:: raw-math(raw)
    :format: latex html

.. _OVITO: https://www.ovito.org/

.. _rendering: https://ovito.org/manual/usage/rendering.html

.. _Kinetic Monte Carlo: https://en.wikipedia.org/wiki/Kinetic_Monte_Carlo

Kinetic Monte Carlo Crystal Growth Tool
----------------------------------------

About
######

This folder contains a Python script ``growth_kmc.py`` which performs a `Kinetic Monte Carlo`_ simulation of crystal growth from solution, starting with a spherical seed. The core developer of this tool is Jacob Jeffries (jwjeffr@clemson.edu or jwjeffr@lanl.gov).

Lattice Styles
##############

Two lattice styles are currently supported, which are TriclinicFourIndices and TriclinicThreeIndices.

The TriclinicFourIndices lattice style initializes sites at the positions:

:raw-math:`$$\mathbf{r}_{ijk\ell} = i\mathbf{a} + j\mathbf{b} + k\mathbf{c} + \ell\mathbf\Delta{r}$$`

where :raw-math:`$i$`, :raw-math:`$j$`, and :raw-math:`$k$` are any integers, and :raw-math:`$\ell$` is restricted to :raw-math:`$0$` and :raw-math:`$1$`. This style is for lattices whose unit cells contain two atoms or molecules, where :raw-math:`$\mathbf{r}_{ijk0}$` is the position of one atom/molecule in the unit cell, and :raw-math:`$\mathbf{r}_{ijk1}$` is the position of the other. The angles :raw-math:`$\alpha$`, :raw-math:`$\beta$`, and :raw-math:`$\gamma$` are also specifiable. See :code:`example_tri4.json` to see this lattice style used in an input file. Angles must be specified in degrees.

The OrthogonalFourIndices lattice style is a special case of the TriclinicFourIndices style, except with :raw-math:`$\alpha=\beta=\gamma=90^\circ$`. See :code:`example_orth4.json` to see this lattice style used in an input file.

The TriclinicThreeIndices lattice style initializes sites at the positions:

:raw-math:`$$\mathbf{r}_{ijk} = i\mathbf{a} + j\mathbf{b} + k\mathbf{c}$$`

where :raw-math:`$i$`, :raw-math:`$j$`, and :raw-math:`$k$` are any integers. This style is for lattices whose unit cells contain one atom or molecule. The angles :raw-math:`$\alpha$`, :raw-math:`$\beta$`, and :raw-math:`$\gamma$` are also specifiable. See :code:`example_tri3.json` to see this lattice style used in an input file. Angles must be specified in degrees.

The OrthogonalThreeIndices lattice style is a special case of the TriclinicThreeIndices style, except with :raw-math:`$\alpha=\beta=\gamma=90^\circ$`. See :code:`example_orth3.json` to see this lattice style used in an input file.

Energetics styles
#################

Three energetics styles are currently supported: IsotropicSecondNearest, AnisotropicThirdNearest, and AnisotropicThirdNearest. The anisotropic styles are both deprecated.

The IsotropicSecondNearest energetics style stores interaction energies between first and second nearest neighbors, specified by the first nearest cutoff and second nearest cutoff. In the input file, the cutoffs are specified as :code:`first_cutoff` and :code:`second_cutoff`, and the respective interaction energies are specified as :code:`first_energy` and :code:`second_energy`. See :code:`example_ortho4.json` to see this energetics style used in an input file.

The AnisotropicThirdNearest energetics style stores interaction energies between first, second, and third nearest neighbors. First nearest neighbor interactions depend on direction. In the input file, the first nearest neighbor interactions in the :code:`a`, :code:`b`, and :code:`c` directions are respectively specified by :code:`e_1a`, :code:`e_1b`, and :code:`e_1c`. The second nearest neighbor interactions in the :code:`b + c`, :code:`b - c`, :code:`c + a`, :code:`-c + a`, :code:`a + b`, and :code:`a - b` directions are respectively specified by :code:`e_2a`, :code:`e_2a_p`, :code:`e_2b`, :code:`e_2b_p`, :code:`e_2c`, and :code:`e_2c_p`. The third nearest neighbor interactions in the :code:`a + b + c`, :code:`a - b - c`, :code:`a - b + c`, and :code:`a + b - c` directions are respectively specified by :code:`e_31`, :code:`e_32`, :code:`e_33`, and :code:`e_34`.

The AnisotropicThirdNearestReconstruction energetics style is identical to the AnisotropicThirdNearest energetics style, except all second-nearest and third-nearest interactions are the same, respectively specified by :code:`second_nearest` and :code:`third_nearest`.

Requirements
##############

-   The :code:`python` interpreter.

-   The external :code:`python` packages :code:`numpy` and :code:`numba`.

The example input (:code:`example_ortho4.json`) provided works for Python 3.9.12, Numpy 1.21.6, and Numba 0.55.1. Other versions are not guaranteed to be functional.

Testing and running the code
#############################

The code can be tested with:

  python growth_kmc.py example_ortho4.json

or:

  ./growth_kmc.py example_ortho4.json

Two runs will be performed:

-   A short, small run which first compiles functions. This run data will be stored in :code:`small.dump` in the LAMMPS-style dump format.

-   A longer run with parameters provided in :code:`example_ortho4.json`. The parameters are:

    Box dimensions = (30, 30, 70) (in lattice units, so :raw-math:`$0 \leq i, j < 30$` and :raw-math:`$0 \leq k < 70$`)

    Number of steps = 300,000

    Dump every = 1,000 steps

    Dump file name = ortho_four.dump

    Initial configuration = spherical with a 75.0 angstrom radius

    Temperature = 300.0 kelvin

    Lattice type = OrthogonalFourIndices with a = b = 9.088 angstroms, c = 6.737 angstroms, and offset = (4.54348, 4.54346, 3.36908) angstroms

    Energetics type = IsotropicSecondNearest

    First-neighbor cutoff distance = 7.0 angstroms

    Second-neighbor cutoff distance = 7.5 angstroms

    First-neighbor interaction energy = -0.291 electron volts

    Second-neighbor interaction energy = -0.186 electron volts

    Adsorption prefactor = 1e+10 hertz

    Adsorption barrier = 0.9 electron volts

    Evaporation prefactor = 1e+10 hertz

    Number of cpus to use = all

    Log file name = ortho_four.log

For new parameters, simply change the dictionary written in :code:`ortho4.json` to match your desired parameters. Note that the energetics and geometric parameters specified in this file are optimized for a `PETN`_ crystal.

Notes
#####

This code is highly parallelized, and will use all available cores unless otherwise specified in the input file. If cores are currently being used, your system might crash. Specify a smaller number of cores with the :code:`num_cpus` input if necessary.
