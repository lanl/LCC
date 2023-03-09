
.. list-table:: 
  :header-rows: 1

  * - Issues
    - Pull Requests
    - CI
  * - .. image:: https://img.shields.io/github/issues/cnegre/LCC-1.svg
        :alt: GitHub issues
        :target: https://github.com/cnegre/LCC-1/issues
    - .. image:: https://img.shields.io/github/issues-pr/cnegre/LCC-1.svg
        :alt: GitHub pull requests
        :target: https://github.com/cnegre/LCC-1/pulls
    - .. image:: https://github.com/cnegre/LCC-1/actions/workflows/main.yml/badge.svg
        :alt: GitHub Actions
        :target: https://github.com/cnegre/LCC-1/actions


LCC
===

About
-----

Los Alamos Crystal Cut (LCC) is simple crystal builder. It is an easy-to-use 
and easy-to-develop code to make crystal solid/shape and slabs from a crystal lattice. 
Provided you have a ‘.pdb‘ file containing your lattice basis you can
create a solid or slab from command line. The core developer of this code is Christian Negre 
(cnegre@lanl.gov).



License
-------

© 2022. Triad National Security, LLC. All rights reserved. This program was produced under U.S. 
Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National 
Nuclear Security Administration. All rights in the program are reserved by Triad National Security, 
LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, 
irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute 
copies to the public, perform publicly and display publicly, and to permit others to do so.

This program is open source under the BSD-3 License.

Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice, this list of conditions and
the following disclaimer.
  2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
and the following disclaimer in the documentation and/or other materials provided with the
distribution.
  3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Requirements
------------

In order to follow this tutorial, we will assume that the reader have a
`LINUX` or `MAC` operative system with the following packages properly
installed:

-   The `git` program for cloning the codes.

-   A `C/C++` compiler (`gcc` and `g++` for example)

-   A `Fortran` compiler (`gfortran` for example)

-   The LAPACK and BLAS libraries (GNU `libblas` and `liblapack`
    for example)

-   The `python` interpreter (not essential).

-   The `pkgconfig` and `cmake` programs (not essential).

On an `x86_64` GNU/Linux Ubuntu 16.04 distribution the commands to be
typed are the following:

          $ sudo apt-get update
          $ sudo apt-get --yes --force-yes install gfortran gcc g++
          $ sudo apt-get --yes --force-yes install libblas-dev liblapack-dev
          $ sudo apt-get --yes --force-yes install cmake pkg-config cmake-data
          $ sudo apt-get --yes --force-yes install git python

**NOTE:** Through the course of this tutorial we will assume that the
follower will work and install the programs in the home directory
(`$HOME`).

Download and installation
---------------------------

We will need to clone the repository as follows:

          $ cd; git@github.com:lanl/LCC.git

Compiling PROGRESS and BML libraries
------------------------------------

The LCC code needs to be compiled with both
[PROGRESS](https://github.com/lanl/qmd-progress) and
[BML](https://github.com/lanl/bml) libraries. In this section we will
explain how to install both of these libraries and link the code against
them.

Scripts for quick installations can be found in the main folder.
In principle one should be able to install everything by typing:

        $ ./clone_libs.sh
        $ ./build_bml.sh
        $ ./build_progress.sh
        $ ./build.sh

Which will also build LCC with its binary file in `./src/lcc_main`.

Step-by-step install
--------------------

Clone the BML library (in your home directory) by doing[^1]:

        $ cd
        $ git clone git@github.com:lanl/bml.git

Take a loot at the `./scripts/example_build.sh` file which has a set of
instructions for configuring. Configure the installation by copying the
script into the main folder and run it:

        $ cp ./scripts/example_build.sh .
        $ sh example_build.sh

The `build.sh` script is called and the installation is configured by
creating the `build` directory. Go into the build directory and type:

        $ cd build
        $ make -j
        $ make install


To ensure bml is installed correctly type `$ make tests` or
`$ make test ARGS="-V"` to see details of the output. Series of tests
results should follow.

After BML is installed, return to you home folder and “clone” the
PROGRESS repository. To do this type:

        $ cd
        $ git clone git@github.com:lanl/qmd-progress.git

Once the folder is cloned, cd into that folder and use the
`example_build.sh` file to configure the installation by following the
same steps as for the bml library.

        $ sh example_build.sh
        $ cd build
        $ make; make install


You can test the installation by typing `$ make tests` in the same way
as it is done for BML.

LCC
---

Open the `Makefile` file in the `lcc/src` folder make sure the
path to both bml and progress libs are set correctly. NOTE: Sometimes,
depending on the architecture the libraries are installed in `/lib64`
instead of `/lib`. After the afforemention changes are done to the
`Makefile` file proceed compiling with the “make” command.

Contributors
------------

Christian Negre, email: cnegre@lanl.gov

Andrew Alvarado, email: aalvarado@lanl.gov


[^1]: In order to have access to the repository you should have a github
    account and make sure to add your public ssh key is added in the
    configuration windows of github account.

Contributing                                                                                                            
------------

Formally request to be added as a collaborator to the project by sending an email to cnegre@lanl.gov. 
After being added to the project do the followig:

  - Create a new branch with a proper name that can identify the new feature (git checkout -b "my_new_branch"
  - Make the changes or add your contributions to the new branch (git add newFile.F90 modifiedFile.F90)
  - Make sure the tests are passing (cd tests ; ./run_test.sh)
  - Commit the changes with proper commit messages (git commit -m "Adding a my new contribution")
  - Push the new branch to the repository (git push)
  - Go to repository on the github website and click on "create pull request"

SUGGESTION: Please, avoid commiting a large number of changes since it is difficult to review. Instead, 
add the changes gradually.


