
.. list-table:: 
  :header-rows: 1

  * - Issues
    - Pull Requests
    - CI
  * - .. image:: https://img.shields.io/github/issues/lanl/LCC.svg
        :alt: GitHub issues
        :target: https://github.com/cnegre/lanl/issues
    - .. image:: https://img.shields.io/github/issues-pr/lanl/LCC.svg
        :alt: GitHub pull requests
        :target: https://github.com/lanl/LCC/pulls
    - .. image:: https://github.com/lanl/LCC/actions/workflows/main.yml/badge.svg
        :alt: GitHub Actions
        :target: https://github.com/lanl/LCC/actions


LCC
===

About
-----

Los Alamos Crystal Cut (LCC) is simple crystal builder. It is an easy-to-use 
and easy-to-develop code to make crystal solid/shape and slabs from a crystal lattice. 
Provided you have a :code:`.pdb` file containing your lattice basis you can
create a solid or slab from command line. The core developer of this code is Christian Negre 
(cnegre@lanl.gov). This documentation has been approved for unlimited release with **LA-UR-23-28084**.


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

    1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

    3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

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
:code:`LINUX` or :code:`MAC` operative system with the following packages properly
installed:

-   The :code:`git` program for cloning the codes.

-   A :code:`C/C++` compiler (:code:`gcc` and :code:`g++` for example)

-   A :code:`Fortran` compiler (:code:`gfortran` for example)

-   The LAPACK and BLAS libraries (GNU :code:`libblas` and :code:`liblapack`
    for example)

-   The :code:`python` interpreter (not essential).

-   The :code:`pkgconfig` and :code:`cmake` programs (not essential).

On an :code:`x86_64` GNU/Linux Ubuntu 16.04 distribution the commands to be
typed are the following:::

         sudo apt-get update
         sudo apt-get --yes --force-yes install gfortran gcc g++
         sudo apt-get --yes --force-yes install libblas-dev liblapack-dev
         sudo apt-get --yes --force-yes install cmake pkg-config cmake-data
         sudo apt-get --yes --force-yes install git python

**NOTE:** Through the course of this tutorial we will assume that the
follower will work and install the programs in the home directory
(:code:`$HOME`).

Quick installation using `spack` 
-----------------------------------

Clone and setup the spack code::

       cd ~

       git clone git@github.com:spack/spack.git

       . spack/share/spack/setup-env.sh

Get info on the package::

       spack info lcc                    

Install the package, this will take a while because it'll install everything from scratch::

       spack install lcc

Load the lcc module::

      spack load lcc

Try lcc::

      cd tmp ; lcc_main

      spack install ovito

      spack load ovito

      cd /tmp

      echo "LCC{ ClusterType= Spheroid TypeOfLattice= FCC AAxis= 10.0 BAxis= 10.0 CAxis= 10.0 }"  | tee input.in  ; lcc_main input.in

      ovito coords.xyz

Download and installation
---------------------------

We will need to clone the repository as follows:::

      cd; git@github.com:lanl/LCC.git

Compiling PROGRESS and BML libraries
------------------------------------

The LCC code needs to be compiled with both
[PROGRESS](https://github.com/lanl/qmd-progress) and
[BML](https://github.com/lanl/bml) libraries. In this section we will
explain how to install both of these libraries and link the code against
them.

Scripts for quick installations can be found in the main folder.
In principle one should be able to install everything by typing:::

       ./clone_libs.sh
       ./build_bml.sh
       ./build_progress.sh
       ./build.sh

Which will also build LCC with its binary file in :code:`./src/lcc_main`.

Step-by-step installation
-------------------------

Clone the BML library (in your home directory) by doing [1]_::

        cd
        git clone git@github.com:lanl/bml.git

Take a loot at the :code:`./scripts/example_build.sh` file which has a set of
instructions for configuring. Configure the installation by copying the
script into the main folder and run it:::

        cp ./scripts/example_build.sh .
        sh example_build.sh

The :code:`build.sh` script is called and the installation is configured by
creating the :code:`build` directory. Go into the build directory and type::

       cd build
       make -j
       make install

To ensure bml is installed correctly type :code:`$ make tests` or
:code:`$ make test ARGS="-V"` to see details of the output. Series of tests
results should follow.

After BML is installed, return to you home folder and “clone” the
PROGRESS repository. To do this type::

        cd
        git clone git@github.com:lanl/qmd-progress.git

Once the folder is cloned, cd into that folder and use the
:code:`example_build.sh` file to configure the installation by following the
same steps as for the bml library.::

        sh example_build.sh
        cd build
        make; make install


You can test the installation by typing :code:`$ make tests` in the same way
as it is done for BML.

Open the :code:`Makefile` file in the :code:`lcc/src` folder make sure the
path to both bml and progress libs are set correctly. NOTE: Sometimes,
depending on the architecture the libraries are installed in :code:`/lib64`
instead of :code:`/lib`. After the aforementioned changes are done to the
:code:`Makefile` file proceed compiling with the “make” command.

Testing the code
-----------------

A test script can be run as follows::

  ./run_test

Quick example run 
-----------------
Assuming the code is installed in the :code:`$HOME` directory, we will run a simple example:::

        cd /tmp 
        echo "LCC{ ClusterType= Spheroid TypeOfLattice= FCC AAxis= 10.0 BAxis= 10.0 CAxis= 10.0 }"  | tee input.in  ; $HOME/LCC/build/lcc_main input.in

This will generate a spherical structure with an FCC lattice using default parameters.         
One can quickly get an input file sample by running the code without giving any input file. 
The available keywords can be listed by running :code:`lcc_main -h` 

Contributors
------------

Christian Negre, email: cnegre@lanl.gov

Andrew Alvarado, email: aalvarado@lanl.gov

Jacob Jeffries, email: jwjeffr@g.clemson.edu

.. [1] In order to have access to the repository you should have a github
    account and make sure to add your public ssh key is added in the
    configuration windows of github account.

Contributing                                                                                                            
------------

Formally request to be added as a collaborator to the project by sending an email to cnegre@lanl.gov. 
After being added to the project do the followig:

  - Create a new branch with a proper name that can identify the new feature (:code:`git checkout -b "my_new_branch"`
  - Make the changes or add your contributions to the new branch (:code:`git add newFile.F90 modifiedFile.F90`)
  - Make sure the tests are passing (:code:`cd tests ; ./run_test.sh`)
  - Commit the changes with proper commit messages (:code:`git commit -m "Adding a my new contribution"`)
  - Push the new branch to the repository (:code:`git push`)
  - Go to repository on the github website and click on "create pull request"

SUGGESTION: Please, avoid commiting a large number of changes since it is difficult to review. Instead, 
add the changes gradually.

Citing
------

If you find this code useful, we encourage you to cite us. Our project has a
citable DOI (`DOI:10.1088/1361-648X/acc294 <https://doi.org/10.1088/1361-648X/acc294>`_) 
with the following :code:`bibtex` snipped:

.. code-block:: bibtex

  @ARTICLE{lcc,
     title    = "A methodology to generate crystal-based molecular structures for
               atomistic simulations",
    author   = "Negre, Christian F A and Alvarado, Andrew and Singh, Himanshu and
               Finkelstein, Joshua and Martinez, Enrique and Perriot, Romain",
    abstract = "We propose a systematic method to construct crystal-based
               molecular structures often needed as input for computational
               chemistry studies. These structures include crystal 'slabs' with
               periodic boundary conditions (PBCs) and non-periodic solids such
               as Wulff structures. We also introduce a method to build crystal
               slabs with orthogonal PBC vectors. These methods are integrated
               into our code,Los Alamos Crystal Cut(LCC), which is open source
               and thus fully available to the community. Examples showing the
               use of these methods are given throughout the manuscript.",
    journal  = "J. Phys. Condens. Matter",
    volume   =  35,
    number   =  22,
    month    =  mar,
    year     =  2023,
    keywords = "crystal structures; extended structures; miller indices; quantum
               chemistry; unit cells",
    language = "en"
  }
