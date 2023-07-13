Input file choices
==================

In this secion we will describe the input file keywords. 
Every valid keword will use "cammel" syntax and will have 
and ``=`` sign right next to ie. For example, the following 
is a valid keyword syntax ``JobName= MyJob``. Coments need 
to have a ``#`` (hash) sign right nex to the phrase we want 
to comment. Example comment could be something like:
``#My comment``.  

`JobName=`
***********
This variable will indicate the name of the job we are sunning. 
It is just a tag to distinguish different outputs. 
As we mentioned before and example use sould be: ``JobName= MyJob``

`Verbose=`
*************
Controls the verbosity level of the output. If set to ``0`` no 
output is pronted out. If set to ``1``, only basic meesages of 
the current execution point of the code will be printed. 
If set to ``2``, information about basic quantities are also 
printed. If set to ``3``, all relevant possible info is printed.

`CoordsOutFile=`
****************
This will store the name of the output coordinates files. Basically 
if ``CoordsOutFile= coords`` two output file will be created: ``coords.xyz``
and ``coords.pdb``. 

  
`PrintCml=`
****************
By setiing ``PrintCml= T`` will also print 
create ``coords.cml`` which can be readed by avogadro. In order to have 
this option working one needs to install  `[openbabel] <https://openbabel.org/wiki/Main_Page>`_
In order to read a cml file one needs to have `[avogadro] <https://avogadro.cc/>`_
installed. On gnu linux:: 

    sudo apt-get avogadro
    sudo apt-get obabel

`ClusterType=` 
**************
This variable will define the type of shape/cluster/slab we 
want to construct. There are many options including 
``Bulk``, ``Planes``, ``Bravais`` and ``Spheroid``. We will explain 
all these in the following sention.

`ClusterType= Bulk`
*******************
This will just cut a "piece of bulk" by indicating how many lattice 
point we want. For example, the following will create a bulk/lattice with 
50 points on each a,b,c direction:: 

  ClusterType= Bulk 
  LatticePoints=  50                                                                                          

The following, instead, will create a bulk/lattice with 100 lattice points 
in the x direction and 50 on the rest::
                                                                                           
  LatticePointsX1=          1
  LatticePointsX2=          100
  LatticePointsY1=          1
  LatticePointsY2=          50
  LatticePointsZ1=          1
  LatticePointsZ2=          50

`ClusterType= Spheroid`
***********************
This will produce a "spheroid" center at the origin. 
And example follows:: 

  ClusterType= Spherid
  LatticePoints=  50   #This is necesary to construct the initial bulk
  AAxis=   1.0 #Radius in direction x
  BAxis=   2.0 #Radius in direction y
  CAxis=   2.0 #Radius in direction z

See section :ref:`regular` to see another example.

`ClusterType= Planes`
***********************
This will cut a shape using Miller indice. This is an important tool to 
construct a slab to study a surface. The cut does not gurantee periodicity.
In order to have a periodic structure different plane boudaries need to 
be tried and the structures needs to be checked using a molecular sivualizer. 
An example is given as follows:: 

   NumberOfPlanes=   6
   Planes[
    0  1  1  2.5
    0 -1 -1  1.5
    0 -1  1  4.5
    0  1 -1  3.5
    1  0  0  4.5
   -1  0  0  3.5
    ]

Three first number on each row indicate the Miller indices. The fourth number indicates how many 
Miller planes from the origin will be cut out. If the number of planes is 6, then the 
system tries to get the slab peridicity vectors since if the Miller planes are orthogonal 
to each other, the shape will be a "Parallelepiped". If instead the number diferent than 6, then 
the periodicity vectors are given by the "Boundaries" of the minimal box that contains the shape.

`CenterAtBox=` 
***********************
If set to ``T``, the shape will be centered at the box (the periodicity vectors 
of the shape/cluster)

`Reorient=` 
***********************
If set to ``T`` this, will reorient the shape, such that vector "a" will 
be aligned with the x dierction. This is important when making slabs 
needed to study a surface.

`AtomType=` 
***********************
This will sed the atom symbol if the lattice basis is not 
read from file.                                                                                                   

`TypeOfLattice=`
***********************
This will set the Lattice unit cell. if set to 
``SC`` or ``FCC`` either a simple cubic or face centered cubinc lattice is built provided 
we set ``LatticeConstanta=`` to the lattice constant 
value. For general unit cell we can set ``TypeOfLattice= Triclinic``, and provide 
the lattice parameters as in the following example::
 
  LatticeConstanta=   6.5329400000000000
  LatticeConstantb=   11.022100000000000
  LatticeConstantc=   7.3568800000000003
  LatticeAngleAlpha=   90.000000000000000
  LatticeAngleBeta=   102.65200000000000
  LatticeAngleGamma=   90.000000000000000

`RandomSeed=` 
***********************
To generate random positions in the lattice. This will 
need to be used in conjunction with ``RCoeff=`` which controll the degree 
of deviation from the lattice positions.

`PrimitiveFormat=`
***********************
This will indicate if the lattice needs to be constructed out 
of a,b,c and angle parameter or primitive lattice vectors. If 
``PrimitiveFormat= Angles`` (default), then the lattice parameters 
will need to be passed as in the following example:: 

  LatticeConstanta=   6.5329400000000000
  LatticeConstantb=   11.022100000000000
  LatticeConstantc=   7.3568800000000003
  LatticeAngleAlpha=   90.000000000000000
  LatticeAngleBeta=   102.65200000000000
  LatticeAngleGamma=   90.000000000000000

If instead, ``PrimitiveFormat= Vectors`` then the primitive vectors 
will need to be passed as in the following example:: 

  LatticeVectors[
    2.0 0   0      #First lattice vector         
    0.0 2.0 0 
    0.0 2.0 2.0 
  ]
                                                         
`UseLatticeBase=`
***********************
This is an important tool that allows us to "dress" every lattice point 
with a basis of choice. The basis is defined to be the minimal set of 
corrdinates and atom types needed to define a crystal system lattice point. 
The basis here will be red from file by providing the latticebase
``LatticeBaseFile=`` wich will contain our atom types and coordinates. 
If ``ReadLatticeFromFile=`` is set to ``T``, then, the lattice parameters will 
be read from the lattice basis file. If is set to ``F``, the the lattice 
parameters will need to be passed as explaines before. 
Another important keyword is the ``BaseFormat=``. If this is set to ``abc``, then 
the basis coordinates stored in the file are assumed to be given in fractional 
coordinates of the lattice parameters. If is set to ``xyz``, the it will be assumed 
to be given in catesian coordinates. 

`SymmetryOperations=`
*************************
If the basis needs to be constricted from symetry operation, 
then one needs to pass all these operation to the code 
as follows:: 

  SymmetryOperations= T
  NumberOfOperations= 4
  Translations[
    0 0 0  0.0
    1 1 1  0.5
    1 1 0  1.0
   -0.5 1.5 0.5 1.0
  ]
  Symmetries[
     0  0   0 
    -1  1  -1
    -1 -1  -1
     1 -1   1
  ]

The first block indictes the "translations" within the unit cell. The first three
rows indicating the directions of the translation and the fourth indicating the intensity. 
The second block indicates the symmetry of operations. For example, if an operation is indicated 
as :math:`(-x + 1/2, -y, -z)` then there will be a traslation ``0.5 0 0 1.0`` and a summetry ``-1 0 0``::

    NumberOfOperations=           0
    MaxCoordination=           1
    NumberOfIterations=           1
    Truncation=   1.0000000000000000E+040
    RCut=   20.000000000000000     
    RTol=   1.0000000000000000E-002
    CutAfterAddingBase=F                                                                                                   
    SeedFile=seed.pdb                                                                                            

