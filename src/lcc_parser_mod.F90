!> This module controls the initialization of the variables.

module lcc_parser_mod

  use lcc_structs_mod
  use lcc_constants_mod
  use prg_kernelparser_mod
  use prg_system_mod
  use lcc_message_mod
  use lcc_aux_mod

  implicit none

  private

  public :: lcc_parse, lcc_make_sample_input, lcc_write_coords

contains

  !> Clustergen parser.
  !! \brief This module is used to parse all the input variables for this program.
  !! Adding a new input keyword to the parser:
  !! - If the variable is real, we have to increase nkey_re.
  !! - Add the keyword (character type) in the keyvector_re vector.
  !! - Add a default value (real type) in the valvector_re.
  !! - Define a new variable and pass the value through valvector_re(num)
  !! where num is the position of the new keyword in the vector.
  !! \param filename File name for the input.
  !! \param bld Build type. 
  !! \param ltt Lattice type. 
  !!
  subroutine lcc_parse(filename,bld,ltt)

    implicit none
    character(len=*), intent(in) :: filename
    type(build_type), intent(inout) :: bld
    type(lattice_type), intent(inout) :: ltt
    integer, parameter :: nkey_char = 14, nkey_int = 15, nkey_re = 13, nkey_log = 11
    integer :: i
    character(20) :: dummyc
    real(dp) :: angle_alpha_r, angle_beta_r, angle_gamma_r
    real(dp), allocatable :: params(:,:)

    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=100) :: &
         'JobName=', 'ClusterType=', 'TypeOfLattice=', 'PrimitiveFormat='&
         &, 'AtomType=', 'BaseFormat=', 'UseLatticeBase=', 'CutAfterAddingBase='&
         &, 'PlanesType=', 'LatticeBaseFile=', 'ReadLatticeFromFile=','CoordsOutFile='&
         &, 'SeedFile=', 'Rdf=' ]
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
         'MyJob', 'Bulk', 'Triclinic', 'Angles' &
         &, 'Au', 'abc', 'F', 'F'&
         &,'Miller', 'latticebase.pdb', 'F', 'coords', &
         &'seed.pdb', ' ']

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
         'ClusterNumber=', 'LatticePoints=', 'LatticePointsX1=' &
         &, 'LatticePointsX2=', 'LatticePointsY1=', 'LatticePointsY2='&
         &, 'LatticePointsZ1=', 'LatticePointsZ2=', 'NSpAtoms='&
         &, 'NumberOfPlanes=', 'RandomSeed=', "NumberOfOperations="&
         &, 'MaxCoordination=',"NumberOfIterations=","Verbose="]
    integer :: valvector_int(nkey_int) = (/ &
         1, -1, -10 &
         &, 10, -10, 10 &
         &, -10, 10, 2000&
         &, 6, 10, 0&
         &, 1, 1, 1/)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
         'Truncation=', 'LatticeConstanta=', 'LatticeConstantb=', 'LatticeConstantc='&
         &, 'LatticeAngleAlpha=', 'LatticeAngleBeta=', 'LatticeAngleGamma='&
         &, 'RCut=', 'AAxis=', 'BAxis=', 'CAxis=','RCoeff=','RTol=' ]
    real(dp) :: valvector_re(nkey_re) = (/&
         1.0d40, 4.08_dp, 4.08_dp, 4.08_dp&
         &, 90.0_dp, 90.0_dp, 90.0_dp&
         &, 20.0_dp, 1.0_dp, 1.0_dp, 1.0_dp,0.0_dp,0.01_dp/)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=100) :: &
         'SymmetryOperations=','CenterAtBox=','Reorient=','WriteCml=',&
         &'CheckLattice=','CheckPeriodicity=','OptimalTranslations=',&
         &'WriteLmp=','InterPlanarDistances=','RandomCoordinates=','RandomLattice=']
    logical :: valvector_log(nkey_log) = (/&
         .false.,.false.,.false.,.false.,&
         &.false.,.false.,.false.,.false.,.true.,.false.,.false./)

    !Start and stop characters
    character(len=50), parameter :: startstop(2) = [character(len=50) :: &
         'LCC{', '}']

    call prg_parsing_kernel(keyvector_char,valvector_char&
         ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
         keyvector_log,valvector_log,trim(filename),startstop)

    !Characters
    bld%job_name = valvector_char(1)
    bld%cl_type = valvector_char(2)
    ltt%type_of_lattice = valvector_char(3)
    ltt%primitive_format = valvector_char(4)
    bld%a_type = valvector_char(5)
    ltt%base_format = valvector_char(6)
    bld%use_lattice_base = valvector_char(7)
    bld%cut_with_base = valvector_char(8)
    bld%planes_type = valvector_char(9)
    bld%latticebase_file = valvector_char(10)
    bld%read_lattice_from_file = valvector_char(11)
    bld%coordsout_file = valvector_char(12)
    bld%seed_file = valvector_char(13)
    bld%rdfPair= valvector_char(14)

    !Integers
    bld%cl_number = valvector_int(1)
    bld%N = valvector_int(2)
    bld%Nx1 = valvector_int(3)
    bld%Nx2 = valvector_int(4)
    bld%Ny1 = valvector_int(5)
    bld%Ny2 = valvector_int(6)
    bld%Nz1 = valvector_int(7)
    bld%Nz2 = valvector_int(8)

    if(bld%N > 0)then
      bld%Nx1 = -int(bld%N/2); bld%Nx2 = int(bld%N/2)
      bld%Ny1 = -int(bld%N/2); bld%Ny2 = int(bld%N/2)
      bld%Nz1 = -int(bld%N/2); bld%Nz2 = int(bld%N/2)
    endif

    if(bld%Nx1 > bld%Nx2 .or. bld%Ny1 > bld%Ny2 .or. &
         & bld%Nz1 > bld%Nz2)then
      write(*,*)"ERROR: N(x/y/z)2 needs to be larger than N(x/y/z)1"
      stop
    endif

    !bld%n_spatoms = valvector_int(9)
    bld%Nplanes = valvector_int(10)
    bld%seed = valvector_int(11)
    ltt%Nop = valvector_int(12)
    bld%MaxCoordination = valvector_int(13)
    bld%niter = valvector_int(14)
    bld%verbose = valvector_int(15)

    !Real
    bld%trunc = valvector_re(1)
    ltt%h_lattice_a = valvector_re(2)
    ltt%h_lattice_b = valvector_re(3)
    ltt%h_lattice_c = valvector_re(4)
    ltt%angle_alpha = valvector_re(5)
    ltt%angle_beta = valvector_re(6)
    ltt%angle_gamma = valvector_re(7)
    bld%r_cut = valvector_re(8)
    bld%a_axis = valvector_re(9)
    bld%b_axis = valvector_re(10)
    bld%c_axis = valvector_re(11)
    bld%rcoeff = valvector_re(12)
    bld%rtol = valvector_re(13)

    !logical
    ltt%bsopl = valvector_log(1)
    bld%center = valvector_log(2)
    bld%reorient = valvector_log(3)
    bld%writecml = valvector_log(4)
    ltt%check = valvector_log(5)
    bld%checkperiod = valvector_log(6)
    ltt%getOptTrs = valvector_log(7)
    bld%writelmp = valvector_log(8)
    bld%interPlanarDistances = valvector_log(9)
    bld%randomCoordinates = valvector_log(10)
    ltt%randomLattice = valvector_log(11)

    if(bld%cl_type == 'Bulk') bld%checkperiod = .false.

    if(bld%cl_type == 'Planes' .or. bld%cl_type == 'PlanesMiller' )then
      write(*,*)"Reading planes ..."
      open(1, file=trim(filename))
      do i = 1,10000
        read(1,*) dummyc
        if(trim(dummyc) == "Planes[")then
          exit
        end if
        if(trim(dummyc) == "}")then
          write(*,*)'ERROR: No cutting planes defined'
          write(*,*)"Here is an example block you should add to the"
          write(*,*)"input file"
          write(*,*)""
          write(*,*)"Planes["
          write(*,*)"   0  1 0 10"
          write(*,*)"   0 -1 0 10"
          write(*,*)"]"
          stop
        end if
      end do
      allocate(bld%planes(3,bld%Nplanes))
      allocate(bld%ploads(bld%Nplanes))
      do i = 1,bld%Nplanes
        read(1,*) bld%planes(1,i),bld%planes(2,i),bld%planes(3,i),bld%ploads(i)
      end do
      write(*,*)""
      close(1)
    end if

    if(bld%cl_type == "Slab")then
      write(*,*)"Reading slab PBC vectors ..."
      open(1, file=trim(filename))
      do i = 1,10000
        read(1,*) dummyc
        if(adjustl(trim(dummyc)) == "Slab[")then
          exit
        endif
        if(adjustl(trim(dummyc)) == "}")then
          write(*,*)'ERROR: No Slab PBC vectors defined'
          write(*,*)"Here is an example block you should add to the"
          write(*,*)"input file"
          write(*,*)""
          write(*,*)"Slab["
          write(*,*)"   1.0 0.0 0.0 10.0"
          write(*,*)"   0.0 1.0 0.0 10.0"
          write(*,*)"   0.0 0.0 1.0 10.0"
          write(*,*)"]"
          stop
        end if
      enddo
      allocate(bld%slab(3,3))
      allocate(bld%sloads(3))
      do i = 1,3
        read(1,*) bld%slab(1,i),bld%slab(2,i),bld%slab(3,i),bld%sloads(i)
      end do
      write(*,*)""
      close(1)
    end if

    if(ltt%bsopl)then
      write(*,*)"Reading translation elements ..."
      open(1, file=trim(filename))
      do i = 1,10000
        read(1,*) dummyc
        if(trim(dummyc) == "Translations[")then
          exit
        end if
        if(trim(dummyc) == "}")then
          write(*,*)'ERROR: No Translation elements defined'
          write(*,*)"Here is an example block you should add to the"
          write(*,*)"input file"
          write(*,*)""
          write(*,*)"Translation["
          write(*,*)"   1  0 0 0.5"
          write(*,*)"   0 -1 0 0.5"
          write(*,*)"]"
          stop
        end if
      end do
      allocate(ltt%bstr(3,ltt%Nop))
      allocate(ltt%bsopload(ltt%Nop))
      do i = 1,ltt%Nop
        read(1,*) ltt%bstr(1,i),ltt%bstr(2,i),ltt%bstr(3,i),ltt%bsopload(i)
      end do
      write(*,*)""
      close(1)
    end if

    if(ltt%bsopl)then
      write(*,*)"Reading symmetry elements ..."
      open(1, file=trim(filename))
      do i = 1,10000
        read(1,*) dummyc
        if(trim(dummyc) == "Symmetries[")then
          exit
        end if
        if(trim(dummyc) == "}")then
          write(*,*)'ERROR: No Symmetry elements defined'
          write(*,*)"Here is an example block you should add to the"
          write(*,*)"input file"
          write(*,*)""
          write(*,*)"Symmetry["
          write(*,*)"   -1  1 -1 "
          write(*,*)"   1 -1 1"
          write(*,*)"]"
          stop
        end if
      end do
      allocate(ltt%bssym(3,ltt%Nop))
      do i = 1,ltt%Nop
        read(1,*) ltt%bssym(1,i),ltt%bssym(2,i),ltt%bssym(3,i)
      end do
      write(*,*)""
      close(1)
    end if

    !Asserts for the input file.
    if(ltt%type_of_lattice.eq.'Triclinic')then
      if(ltt%primitive_format.ne.'Angles')then
        if(ltt%primitive_format.ne.'Vectors')then
          stop 'WARNING: No PrimitiveFormat specification'
        end if
      end if
    end if

    if(bld%cut_with_base.eq.'T')then
      if(bld%use_lattice_base.eq.'F')then
        stop 'WARNING: UseLatticeBase must be consistent with CutAfterAddingBase option'
      end if
    end if

    allocate(ltt%lattice_vectors(3,3))

    if(ltt%primitive_format.eq.'Angles')then
      allocate(params(2,3))
      params(1,1) = ltt%h_lattice_a
      params(1,2) = ltt%h_lattice_b
      params(1,3) = ltt%h_lattice_c
      params(2,1) = ltt%angle_alpha
      params(2,2) = ltt%angle_beta
      params(2,3) = ltt%angle_gamma

      call lcc_parameters_to_vectors(params,ltt%lattice_vectors,bld%verbose)

      call lcc_vectors_to_parameters(ltt%lattice_vectors,params,bld%verbose)


    elseif(ltt%primitive_format.eq.'Vectors')then
      write(*,*)"Using primitive vectors ..."
      open(1, file=trim(filename))
      do i = 1,10000
        read(1,*) dummyc
        if(trim(dummyc) == "LatticeVectors[")then
          exit
        end if
        if(trim(dummyc) == "}")then
          write(*,*)'ERROR: No LatticeVectors are specified. These vectors are needed'
          write(*,*)'if PrimitiveFormat= Vectors'
          write(*,*)"Here is an example block you should add to the"
          write(*,*)"input file"
          write(*,*)""
          write(*,*)"LatticeVectors["
          write(*,*)"   1  0  1 "
          write(*,*)"   0  1  1 "
          write(*,*)"   0  0  1 "
          write(*,*)"]"
          stop
        end if
      end do
      do i = 1,3
        read(1,*)ltt%lattice_vectors(i,1), ltt%lattice_vectors(i,2), ltt%lattice_vectors(i,3)
      end do
      write(*,*)""
      close(1)
    else
      write(*,*)"ERROR PrimitiveFormat is not specified"
      stop
    endif

    if(trim(ltt%type_of_lattice) == "FCC" .or. trim(ltt%type_of_lattice) == "FCC") then
      if(abs(ltt%angle_alpha - 90.0_dp) > 0.000000000001_dp .or. &
           & abs(ltt%angle_beta - 90.0_dp) > 0.000000000001_dp .or. &
           & abs(ltt%angle_gamma - 90.0_dp) > 0.000000000001_dp)then
        write(*,*)"WARNING SC, and FCC will have Alpha = Beta = Gamma = 90"
      end if
    end if

    if(bld%cl_type == "Growth") then
      call prg_parse_system(bld%syseed,bld%seed_file)
    end if

    if(bld%cl_type == "Dress" .and. bld%use_lattice_base == "F")then
      STOP "UseLatticeBase must be set to T for ClusterType= Dress"
    end if

    call lcc_print_message("############### End of parameters read ################",&
         &bld%verbose)

  end subroutine lcc_parse


  !> Make a sample inputfile sample_input.in.
  !!
  subroutine lcc_make_sample_input()

    open(1,file="sample_input.in")
    write(1,'(A23)')"#Clustergen input file."
    write(1,*)""
    write(1,*)"CLGEN{                                                       "
    write(1,*)"                                                             "
    write(1,*)"  JobName=                 Clusters                          "
    write(1,*)"  ClusterType=             Planes            #or Spheroid/Bulk/Planes/PlanesMiller    "
    write(1,*)""
    write(1,*)'  LatticeBaseFile=        "latticebase.pdb"  #Lattice base unit cell   '
    write(1,*)"                                                              "
    write(1,*)"  TypeOfLattice=           Triclinic     #Or FCC/SC           "
    write(1,*)"  LatticePoints=           5             #Number of total lattice points in one dimention            "
    write(1,*)""
    write(1,*)"  #This section is only valid if ClusterType is set to Bulk                                          "
    write(1,*)"  LatticePointsX1=        -1             #Initial lattice  point in X direction                      "
    write(1,*)"  LatticePointsX2=         1             #Final lattice point in X direction                         "
    write(1,*)"  LatticePointsY1=        -1                                    "
    write(1,*)"  LatticePointsY2=         1                                    "
    write(1,*)"  LatticePointsZ1=        -1                                    "
    write(1,*)"  LatticePointsZ2=         1                                    "
    write(1,*)""
    write(1,*)"  PrimitiveFormat=         Angles        #Or Vectors            "
    write(1,*)"  AtomType=                X            #If the lattice basis is not provided                       "
    write(1,*)"  UseLatticeBase=          T             #Add a the basis to the lattice points                      "
    write(1,*)"  BaseFormat=              xyz           #Or abc                "
    write(1,*)"  CutAfterAddingBase=      F                                    "
    write(1,*)"                                                                "
    write(1,*)"  #If the lattice parameters are not specified in the LatticeBasisFile:                              "
    write(1,*)"  LatticeConstanta=        7.763                                "
    write(1,*)"  LatticeConstantb=        8.7109                               "
    write(1,*)"  LatticeConstantc=        10.8701                               "
    write(1,*)"  LatticeAngleAlpha=       90.0                                  "
    write(1,*)"  LatticeAngleBeta=        102.937                                "
    write(1,*)"  LatticeAngleGamma=       90.0                                  "
    write(1,*)"                                                                "
    write(1,*)"  #If PrimitiveFormat is set to Vectors:                        "
    write(1,*)"  Vectors[                                                      "
    write(1,*)"    1.0 0.0 0.0                                                 "
    write(1,*)"    0.0 1.0 0.0                                                 "
    write(1,*)"    0.0 0.0 1.0                                                 "
    write(1,*)"  ]                                                             "
    write(1,*)"                                                                "
    write(1,*)"  #If ClusterType is set to Planes:                             "
    write(1,*)"  PlanesType= Miller                                            "
    write(1,*)"  NumberOfPlanes= 6                                             "
    write(1,*)"  Planes[                                                       "
    write(1,*)"    -1.0   0.0   0.0  10.0   #x, y, and z plane direction; distance to origin  "
    write(1,*)"     1.0   0.0   0.0  10.0                                         "
    write(1,*)"     0.0  -1.0   0.0  10.0                                         "
    write(1,*)"     0.0   1.0   0.0  10.0                                         "
    write(1,*)"     0.0   0.0  -1.0  10.0                                         "
    write(1,*)"     0.0   0.0   1.0  10.0                                         "
    write(1,*)"  ]                                                             "
    write(1,*)"                                                                "
    write(1,*)"}                  "
    close(1)

  end subroutine lcc_make_sample_input

  !> Writes the coordinates to a file (coordsandbase.pdb)
  !! \param sy System type.
  !! \param bld Build type.
  !! \param coordsout_file File name to write the coordinates to.
  !! \param verbose Verbosity level.
  !!
  subroutine lcc_write_coords(sy,bld,coordsout_file,verbose)
    type(system_type) :: sy
    integer :: Ncluster
    character(len=*) :: coordsout_file
    character(2), allocatable :: atomname_in(:), atom_in(:), resname_in(:)
    real(dp), allocatable :: r_cluster(:,:)
    integer, allocatable :: resindex_in(:)
    type(build_type) :: bld
    integer, intent(in) :: verbose
    character(70) :: longChar

    if(.not.allocated(sy%coordinate))then
      write(*,*)"ERROR at lcc_write_coords: System coordinates are not allocated"
      stop
    endif

    if(sy%nats <= 0)sy%nats = size(sy%coordinate,dim=2)

    if(.not.allocated(sy%symbol))then
      call lcc_print_warning("lcc_write_coords",&
           &"Symbols are not provided. Using X instead",verbose)
      allocate(sy%symbol(sy%nats))
      sy%symbol = "X"
    endif

    if(.not. allocated(sy%resname))then
      call lcc_print_warning("lcc_write_coords",&
           &"Resnames are not provided. Using MOL instead",verbose)
      allocate(sy%resname(sy%nats))
      sy%resname = "M"
    endif

    if(.not.allocated(sy%resindex))then
      call lcc_print_warning("lcc_write_coords",&
           &"Resindex are not provided. Using 1 instead",verbose)
      allocate(sy%resindex(sy%nats))
      sy%resindex = 1
    endif

    if(.not.allocated(sy%atomname))then
      call lcc_print_warning("lcc_write_coords",&
           &"Atom names are not provided. Using X instead",verbose)
      allocate(sy%atomname(sy%nats))
      sy%atomname = "X"
    endif

    if(.not.allocated(sy%spindex))then
      call lcc_print_warning("lcc_write_coords",&
           &"Species indices are not provided. Using 1 instead",verbose)
      allocate(sy%spindex(sy%nats))
      sy%spindex = 1
    endif

    if(.not.allocated(sy%lattice_vector))then
      call lcc_print_warning("lcc_write_coords",&
           &"Lattice vectors are not provided. Using 10 ang instead",verbose)
      allocate(sy%lattice_vector(3,3))
      sy%lattice_vector = 10.0_dp
    endif

    call prg_write_system(sy,trim(coordsout_file),"pdb")
    call prg_write_system(sy,trim(coordsout_file),"xyz")

    if(bld%writecml)then
      longChar = "obabel -ipdb "//coordsout_file//".pdb -ocml -O"//coordsout_file//".cml "
      call system(longChar)
    endif

    if(bld%writelmp)then
      call prg_write_system(sy,trim(coordsout_file),"lmp")
    endif

  end subroutine lcc_write_coords

end module lcc_parser_mod
