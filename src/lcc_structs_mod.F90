!> A module to handle the structures needed by the code.
!! \brief This module will be used to build and handle structures in the 
!! code.
!! @ingroup LCC
!!
module lcc_structs_mod

  use lcc_constants_mod
  use prg_system_mod

  implicit none

  private

  !> Build type
  type, public :: build_type  !< The type of build/run.

    !> Job name.
    character(len=20) :: job_name

    !> Output file name
    character(len=20) :: output_file_name

    !> Output file name for coordinates
    character(len=60), public :: coordsout_file

    !> Lattice base file name
    character(len=60), public :: latticebase_file

    !> Cut lattice using planes
    character(len=1) :: cut_by_planes

    !> Cut lattice after base is added 
    character(len=1) :: cut_with_base

    !> Read lattice from file
    character(len=1) :: read_lattice_from_file

    !> Use lattice base
    character(len=1) :: use_lattice_base

    !> Cluster (or solid) shape to be constructed
    character(len=60) :: cl_type

    !> Type of planes used for the cut
    character(len=60) :: planes_type

    !> File name for the seed used to grow a cluster
    character(len=60) :: seed_file

    !> Number of atoms
    integer :: N 
   
    !> Number of planes to use in the cut
    integer :: Nplanes

    !> Number of lattice points in +-(x, y, and z) directions 
    integer :: Nx1, Nx2, Ny1, Ny2, Nz1, Nz2

    !> Random seed
    integer :: seed

    !> Cluster number (if it is a solid with "magic" numbers) 
    integer :: cl_number
    
    !> Axis length if cluster is a spheroid
    real(dp) :: a_axis, b_axis, c_axis
    
    !> Coefficient used with random seed to create noise in coordinates
    real(dp) :: rcoeff

    !> Cutoff radius to build spheroids
    real(dp) :: r_cut
    
    !> Truncation for solids 
    real(dp) :: trunc

    !> Atom type (if specified on the input file)
    character(2) :: a_type

    !> Planes for the cut 
    real(dp), allocatable  ::  planes(:,:)
    
    !> Plenes weight factors
    real(dp), allocatable :: ploads(:)

    !> System seed to be grow on top
    type(system_type) :: syseed

    !> Number of atoms in cluster/slab
    integer :: Ncluster

    !> Atoms in the cluster/slab
    character(2), allocatable :: atom_in(:)
    character(2), allocatable :: atomname_in(:)
    integer, allocatable :: resindex_in(:)
    character(2), allocatable :: resname_in(:)

    !> Coordinates of the resulting cluster/slab
    real(dp), allocatable :: r_cluster(:,:)

    !> Max coordination number
    integer :: MaxCoordination
    
    !> Distance tolerance for distinguishing coordinates
    real(dp) :: rtol
    
    !> Number of iterations
    integer :: niter
    
    !> Verbose level
    integer :: verbose

    !> Center at box
    logical :: center
    
    !> Reorient first lattice vector toward x direction
    logical :: reorient
    
    !> Reorient first lattice vector toward x direction
    logical :: writecml

    !> To check periodicity
    logical :: checkperiod

    !> To compute RDFs
    character(5) :: rdfPair 

    !> Write LAMMPS input coordinates
    logical :: writelmp

    !> Use "number of interplanar distances" as unit of measurement for plane cut
    logical :: interPlanarDistances

    !> To build a slab out of regular vectors
    real(dp),allocatable :: slab(:,:)
    real(dp),allocatable :: sloads(:)
 
    !> To add randomness to coordinates
    logical :: randomCoordinates

    !> Width of a "chunk of lattice"
    integer :: width 

  end type build_type


  !> Lattice type to be read and extended
  type, public :: lattice_type  !< The type of lattice read from input.

    !> Lattice basis
    character(len=3) :: base_format

    !> The lattice primitive format (Angles of Vectors)
    character(len=60) ::  primitive_format 
    
    !> Type of lattice (sc, bcc, fcc, and triclinic)
    character(len=60) :: type_of_lattice

    !> Angles for triclinic lattice
    real(dp) :: angle_alpha, angle_beta, angle_gamma

    !> abc parameters for lattice
    real(dp) :: h_lattice_a, h_lattice_b, h_lattice_c
  
    !> Lattice vectors
    real(dp), allocatable :: lattice_vectors(:,:)

    !> Volume of the cell
    real(dp) :: volr

    !> Lattice reciprocal vectors
    real(dp), allocatable :: recip_vectors(:,:)

    !> Volume of the reciprocal cell
    real(dp) :: volk

    !> Number of atoms in the basis
    integer :: Nbase
    
    !> Basis atoms
    character(2), allocatable :: base_atom(:)

    !> Basis coordinates
    real(dp), allocatable :: r_base(:,:)
    
    !> System for the basis 
    type(system_type) :: sybase

    !> If there are symmetry operations to be performed
    logical :: bsopl

    !> Number of Symmetry operations
    integer :: Nop

    !> Translations to be performed 
    real(dp), allocatable :: bstr(:,:)

    !> Scaling factos (load) for the translation
    real(dp), allocatable :: bsopload(:)

    !> Symmetry operation (diagonal) 
    real(dp), allocatable :: bssym(:,:)

    !> Spicies index
    integer, allocatable :: spindex(:)
    
    !> System basis masses
    real(dp), allocatable :: base_mass(:)
    
    !> Residue index
    integer, allocatable :: resindex(:)

    !> To save the "bulk" positions
    real(dp), allocatable :: bulk(:,:)

    !> Check lattice 
    logical :: check

    !> Get optimal translations at symmetry operations
    logical :: getOptTrs

    !> To add randomness to each lattice position
    logical :: randomLattice

    !> To add random orientations to each lattice molecules/basis
    logical :: randrotations
    
    !> To set a particular density in [gr/cc]
    real(dp) :: setdensity

  end type lattice_type  

  
  !> Monte Carlo build type
  type, public :: mc_type  !< The type of lattice read from input.

    !> Job name
    character(len=20) :: job_name

    !> Job name
    character(len=20) :: run_type

    !> Number of cycles 
    integer :: num_cycles 

    !> Number of programs 
    integer :: num_prog 

    !> Number of initial cycles 
    integer :: num_init_cycles 

    !> Width of the crystal chunk
    integer :: width

    !> MC interation energies 
    !! hxmz means an interation with the neighbor in the a1, -a3 direction.
    integer :: hx, hy, hz, hxy, hxmy, hxz, hxmz, hyz, hymz, hxyz, hmxyz
    integer :: hxmyz, hxymz, h0 

  end type  mc_type

end module lcc_structs_mod

