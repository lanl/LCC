!> Module for printing through the code.
!!
module lcc_message_mod

  use lcc_constants_mod

  implicit none

  public :: lcc_print_ussage, lcc_print_message, lcc_print_warning
  public :: lcc_print_error, lcc_print_intVal, lcc_print_realVal
  public :: lcc_print_realVect

contains

  !> For printing the instructions on how to execute
  !! the code.
  !!
  subroutine lcc_print_ussage()
    implicit none
    write(*,*)""
    write(*,*)"Usage:"
    write(*,*)""
    write(*,*)"  $ clgen <inputfile>"
    write(*,*)""
    write(*,*)"<inputfile>: File containing the input parameters  "
    write(*,*)""
    write(*,*)" A sample_input.in file in now created."
    write(*,*)""
  end subroutine lcc_print_ussage

  !> Print a simple message.
  !! \param message Message to print.
  !! \param verbose Verbosity level.
  !!
  subroutine lcc_print_message(message,verbose)
    implicit none
    character(len=*), intent(in) ::  message
    integer, intent(in) :: verbose

    if(verbose > 0)then
      write(*,*)""
      write(*,*)">> ",trim(adjustl(message))
    endif

  end subroutine lcc_print_message

  !> Print a Warning (will not stop execution).
  !! \param at Name of the routine.
  !! \param message Message to print.
  !! \param verbose Verbosity level.
  !!
  subroutine lcc_print_warning(at,message,verbose)
    implicit none
    character(len=*), intent(in) :: at, message
    integer, intent(in) :: verbose

    if(verbose > 0)then
      write(*,*)""
      write(*,*)">> WARNING at subroutine ",trim(adjustl(at)),&
           &": ",trim(adjustl(message))
    endif

  end subroutine lcc_print_warning

  !> Print error (will stop execution).
  !! \param at Name of the routine.
  !! \param message Message to print.
  !!
  subroutine lcc_print_error(at,message)
    implicit none
    character(len=*), intent(in) ::  at, message

    write(*,*)""
    write(*,*)">> ERROR at subroutine ",trim(adjustl(at)),&
         &": ",trim(adjustl(message))
    stop

  end subroutine lcc_print_error

  !> Print integer magnitude.
  !! \param name Name of the magnitude.
  !! \param value Value to print.
  !! \param units Units of the magnitude.
  !!
  subroutine lcc_print_intVal(name,value,units,verbose)
    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in) :: value, verbose
    character(len=*), intent(in) :: units

    if(verbose > 1)then
      write(*,*)""
      write(*,*)">> ",name,"= ",value," [",units,"]"
    endif

  end subroutine lcc_print_intVal


  !> Print real magnitude.
  !! \param name Name of the magnitude.
  !! \param value Value to print.
  !! \param units Units of the magnitude.
  !!
  subroutine lcc_print_realVal(name,value,units,verbose)
    implicit none
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value
    integer, intent(in) :: verbose
    character(len=*), intent(in) :: units

    if(verbose > 1)then
      write(*,*)""
      write(*,*)">> ",name,"= ",value," [",units,"]"
    endif

  end subroutine lcc_print_realVal


  !> Print real vector.
  !! \param name Name of the quantities.
  !! \param vect Vector to print.
  !! \param units Units of the quantities.
  !! \param verbose Verbosity level.
  !!
  subroutine lcc_print_realVect(name,vect,units,verbose)
    implicit none
    character(len=*), intent(in) :: name
    real(dp), allocatable, intent(in) :: vect(:)
    integer, intent(in) :: verbose
    character(len=*), intent(in) :: units

    if(verbose > 2)then
      write(*,*)""
      write(*,*)">> ",name,units,"="
      write(*,*)vect(1),vect(2),vect(3)
    endif

  end subroutine lcc_print_realVect

  !> Print real vector.
  !! \param name Name of the quantities.
  !! \param mat Matrix to print.
  !! \param units Units of the quantities.
  !! \param verbose Verbosity level.
  !!
  subroutine lcc_print_realMat(name,mat,units,verbose)
    implicit none
    character(len=*), intent(in) :: name
    real(dp), allocatable, intent(in) :: mat(:,:)
    integer, intent(in) :: verbose
    character(len=*), intent(in) :: units
    character(40) :: f

    f = '(1X,1f5.3,1X,1f5.3,1X,1f5.3)'
    if(verbose > 2)then
      write(*,*)""
      write(*,*)">> ",name,units,"="
      write(*,f)mat(1,1:3)
      write(*,f)mat(2,1:3)
      write(*,f)mat(3,1:3)
    endif

  end subroutine lcc_print_realMat

  subroutine lcc_help()
    implicit none
    write(*,*)"LCC HELP"
    write(*,*)"========"
    write(*,*)""
    write(*,*)"ClusterType= Bulk "
    write(*,*)"    Will determine the type of cluster to build"
    write(*,*)"    Options are: Bulk, Planes, Slab, Spheroid, and BravaisGrowth"
    write(*,*)""
    write(*,*)"Verbose= 1 "
    write(*,*)"    Will control the verbosity of the output"
    write(*,*)"    0: No output printed execpt for input variables and errors"
    write(*,*)"    1: Status messages of the execution points will be printed"
    write(*,*)"    2: Scalars will be printed"
    write(*,*)"    3: Scalars and vectors will be printed"
    write(*,*)"    > 3: Everything is printed"
    write(*,*)""
    write(*,*)"ReadLatticeFromFile= F"
    write(*,*)"    Logical variable to indicate if the Lattice will be read from the "
    write(*,*)""
    write(*,*)"LatticeBaseFile= 'unit_cell.pdb' "
    write(*,*)"    Will tell the program the name of the file containing the lattice basis "
    write(*,*)""
    write(*,*)"LatticeConstanta= 1.0"
    write(*,*)"    Lattice parameter a to define the lattice"
    write(*,*)"    Only active when ReadLatticeFromFile= F "
    write(*,*)""
    write(*,*)"LatticeConstantab= 1.0"
    write(*,*)"    Lattice parameter b to define the lattice"
    write(*,*)"    Only active when ReadLatticeFromFile= F "
    write(*,*)""
    write(*,*)"LatticeConstantc= 1.0"
    write(*,*)"    Lattice parameter c to define the lattice"
    write(*,*)"    Only active when ReadLatticeFromFile= F "
    write(*,*)""
    write(*,*)"LatticeAngleAlpha= 90.0"
    write(*,*)"    Lattice parameter alpha to define the lattice"
    write(*,*)"    Only active when ReadLatticeFromFile= F "
    write(*,*)""
    write(*,*)"LatticeAngleBeta= 90.0"
    write(*,*)"    Lattice parameter beta to define the lattice"
    write(*,*)"    Only active when ReadLatticeFromFile= F "
    write(*,*)""
    write(*,*)"LatticeAngleGamma= 90.0"
    write(*,*)"    Lattice parameter gamma to define the lattice"
    write(*,*)"    Only active when ReadLatticeFromFile= F "
    write(*,*)""
    write(*,*)"SeedFile= seed.pdb" 
    write(*,*)"    Character variable to indicate the name of the seed structure to grow"
    write(*,*)"    from."
    write(*,*)""
    write(*,*)"Rdf= "
    write(*,*)"    Character variable to indicate if the Radial distribution functions needs"
    write(*,*)"    to be computed. Example: If Rdf= C-C, then the C-C radial distribition function will"
    write(*,*)"    be computed"
    write(*,*)""
    write(*,*)"RandomSeed= 123" 
    write(*,*)"    Integer variable to indicate the seed that will be used for random number generations."
    write(*,*)""
    write(*,*)"RandomCoordinates= F" 
    write(*,*)"    Logical variable to indicate if to add random displacements to the system's coordinate."
    write(*,*)"    It uses RCoeff= variable to set the max displacement."
    write(*,*)""
    write(*,*)"RandomLattice= F" 
    write(*,*)"    Logical variable to indicate if to add random displacements to the lattice points."
    write(*,*)"    It uses RCoeff= variable to set the max displacement."
    write(*,*)""
    write(*,*)"MaxCoordination= 0"
    write(*,*)"    Coordination to search for when Bravais method is used to grow from a seed coordinate."
    write(*,*)"    If MaxCoodination= 3 then a surface atom will be included in the growing shape only if" 
    write(*,*)"    it has a coordination hicher than 3."
    write(*,*)""
    write(*,*)"Reorient= F"
    write(*,*)"    If set to T it will reorient the shape so that the x axis will be alligned with the first"
    write(*,*)"    Miller plane"
    write(*,*)""
    write(*,*)"Slab[  ]"
    write(*,*)"    Active when TypeOfCluster= Slab"
    write(*,*)"    Contains the three PBC vectors and their lenghts needed to construct a slab. "
    write(*,*)"    An example follows:"
    write(*,*)"    Slab["
    write(*,*)"         1.0 0.0 0.0 10.0"
    write(*,*)"         0.0 1.0 1.0 10.0"
    write(*,*)"         0.0 0.0 1.0 10.0"
    write(*,*)"         ]"
    write(*,*)"    The latter means that the first PBC vector will have a lenght of 10.0 and a (1,0,0) direction."
    write(*,*)""
    write(*,*)"LatticePoints="
    write(*,*)"    Number of lattice points in x, y, and z diretion"
    write(*,*)""
    write(*,*)"LatticePointsX1="
    write(*,*)"    Number of point in -x direction"
    write(*,*)""
    write(*,*)"LatticePointsX2="
    write(*,*)"    Number of point in +x direction"
    write(*,*)""
    write(*,*)"LatticePointsY1="
    write(*,*)"    Number of point in -y direction"
    write(*,*)""
    write(*,*)"LatticePointsY2="
    write(*,*)"    Number of point in +y direction"
    write(*,*)""
    write(*,*)"LatticePointsZ1="
    write(*,*)"    Number of point in -z direction"
    write(*,*)""
    write(*,*)"LatticePointsZ2="
    write(*,*)"    Number of point in +z direction"
    write(*,*)""
    write(*,*)"NumberOfPlanes= 6"
    write(*,*)"    Total number of planes"
    write(*,*)"    Only active if ClusterType= Planes"
    write(*,*)""
    write(*,*)"Planes[ ]"
    write(*,*)"    List of planes and their position from the origin"
    write(*,*)"    An example follows:"
    write(*,*)"    Planes[ "
    write(*,*)"          1  0  1  10.0 "
    write(*,*)"         -1  0 -1  10.0 "
    write(*,*)"          0  1  0  10.0 "
    write(*,*)"          0 -1  0  10.0 "
    write(*,*)"          0  0  1  10.0 "
    write(*,*)"          0  0 -1  10.0 "
    write(*,*)"          ] "
    write(*,*)"    The latter means that the first plane will have a 1 0 1 orientation and will be placed at 10 from the origin"
    write(*,*)""
    write(*,*)"InterPlanarDistances= T"
    write(*,*)"    When set to T it will position the planes based on units of interplanar distances"
    write(*,*)""


  end subroutine lcc_help

end module lcc_message_mod
