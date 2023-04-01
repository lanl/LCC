!> Module to hold routines for handling the lattice and lattice base.
!!
module lcc_lattice_mod

  use bml
  use lcc_constants_mod
  use prg_system_mod
  use prg_extras_mod
  use prg_syrotation_mod
  use lcc_structs_mod
  use lcc_allocation_mod
  use lcc_message_mod

  implicit none

  public :: lcc_make_lattice, lcc_set_atom_type

contains

  !> Make a lattice depending on the input parameter.
  !! \brief This will make one of the following latices:
  !! SC: Simple cubic, FCC: Face center cubic, or Triclinic.
  !! \param bld Building structure (see lcc_structures_mod)
  !! \param ltt Lattice structure (see lcc_scturctures_mod)
  !! \param check If we want to check the basis for atom repetition.
  !! Note that checks can be expensive.
  !!
  subroutine lcc_make_lattice(bld,ltt,check,sy)
    integer :: mlsI
    type(lattice_type), intent (inout) :: ltt
    type(build_type), intent (inout) :: bld
    logical, intent(in) :: check
    type(system_type), intent(inout) :: sy
    real(dp), allocatable :: auxArr(:)

    call lcc_print_message("Making lattice ...",bld%verbose)

    if (bld%use_lattice_base == 'T') call lcc_read_base(bld,ltt,check,bld%verbose)

    if(ltt%type_of_lattice.eq.'SC')then
      call lcc_sc(bld%Nx1,bld%Nx2,bld%Ny1,bld%Ny2,bld%Nz1,bld%Nz2&
           &,ltt%h_lattice_a,sy%lattice_vector,sy%coordinate)
      sy%nats = size(sy%coordinate(1,:))
    elseif(ltt%type_of_lattice.eq.'FCC')then
      call lcc_fcc(bld%Nx1,bld%Nx2,bld%Ny1,bld%Ny2,bld%Nz1,bld%Nz2&
           &,ltt%h_lattice_a,sy%lattice_vector,sy%coordinate,bld%verbose)
      sy%nats = size(sy%coordinate(1,:))
    elseif(ltt%type_of_lattice.eq.'Triclinic')then
      call lcc_triclinic(bld%Nx1,bld%Nx2,bld%Ny1,bld%Ny2,bld%Nz1,bld%Nz2&
           &,ltt%lattice_vectors,sy%lattice_vector,sy%coordinate,bld%verbose)
      sy%nats = size(sy%coordinate(1,:))
    else
      call lcc_print_error("lcc_make_lattice","The requested TypeOfLattice is not implemented")
    endif

    if(ltt%randomLattice)then
      call lcc_add_randomness(sy%coordinate,ltt%lattice_vectors,bld%seed,bld%rcoeff)
    endif

    if(bld%checkperiod)then
      call lcc_reallocate_realMat(ltt%bulk,3,sy%nats)
      ltt%bulk = sy%coordinate
    endif
    
    allocate(auxArr(3))
    auxArr = ltt%lattice_vectors(1,:)
    call lcc_print_realVect("Lattice vector a1",auxArr," [Ang]",bld%verbose)
    auxArr = ltt%lattice_vectors(2,:)
    call lcc_print_realVect("Lattice vector a2",auxArr," [Ang]",bld%verbose)
    auxArr = ltt%lattice_vectors(3,:)
    call lcc_print_realVect("Lattice vector a3",auxArr," [Ang]",bld%verbose)

  end subroutine lcc_make_lattice


  !> Reading the basis from an input file.
  !! \brief This will read the coordinates for the basis from an input file
  !! If information about the lattice is contained, it will also be read.
  !! \param bld Building structure (see lcc_structures_mod).
  !! \param ltt Lattice structure (see lcc_scturctures_mod).
  !! \param check If we want to check the basis for atom repetition.
  !! \param verbose Verbose level.
  !! Note that checks can be expensive.
  !!
  subroutine lcc_read_base(bld,ltt,check,verbose)
    type(build_type), intent(inout) :: bld
    type(lattice_type), intent(inout) :: ltt
    type(system_type) :: sybase
    integer :: i,j,a,b,ai
    real(dp) :: massTot,volk,volr,newVol,scaleVect
    real(dp) :: orig(3), myTrs(3)
    logical, intent(in) :: check
    integer, intent(in) :: verbose

    !Parsing with progress library
    call prg_parse_system(ltt%sybase, bld%latticebase_file)

    orig = 0.5_dp
    !The basis in the unit cell may be different
    !than the one in the file due to symmetry operations
    ltt%nbase = ltt%sybase%nats
    allocate(ltt%base_atom(ltt%nbase))
    allocate(ltt%r_base(3,ltt%nbase))
    ltt%r_base = ltt%sybase%coordinate
    ltt%base_atom = ltt%sybase%symbol
    allocate(ltt%base_mass(ltt%nbase))
    ltt%base_mass = ltt%sybase%mass

    if (allocated(ltt%sybase%lattice_vector) .and. trim(bld%read_lattice_from_file) == "T") then
      call lcc_print_message("Picking up lattice parameters from file ...",verbose)
      ltt%lattice_vectors = ltt%sybase%lattice_vector
    else
      call lcc_print_message("Lattice parameters will be read from the main input file ...",verbose)
    end if

    !If there are NO symmetry operations
    if(.not.ltt%bsopl)then
      ltt%Nop = 1
      call lcc_reallocate_realVect(ltt%bsopload,ltt%Nop)
      call lcc_reallocate_realMat(ltt%bstr,3,ltt%Nop)
      call lcc_reallocate_realMat(ltt%bssym,3,ltt%Nop)
      ltt%bssym = 1.0_dp
      ltt%bstr = 0.0_dp
    endif

    call lcc_print_message("Applying translations ...",verbose)
    ltt%Nbase = ltt%Nop*ltt%Nbase
    if(allocated(ltt%base_atom)) deallocate(ltt%base_atom); allocate(ltt%base_atom(ltt%Nbase))
    call lcc_reallocate_realMat(ltt%r_base,3,ltt%Nbase)
    call lcc_reallocate_intVect(ltt%spindex,ltt%Nbase)
    call lcc_reallocate_intVect(ltt%resindex,ltt%Nbase)
    call lcc_reallocate_realVect(ltt%base_mass,ltt%Nbase)
    do i = 1,ltt%Nop
      a = (i-1)*ltt%sybase%nats + 1
      b = i*ltt%sybase%nats
      ltt%base_atom(a:b) = ltt%sybase%symbol(:)
      ltt%r_base(:,a:b) = ltt%sybase%coordinate(:,:)
      ltt%resindex(a:b) = i
      ltt%spindex(a:b) = ltt%sybase%spindex(:)
      ltt%base_mass(a:b) = ltt%sybase%mass(:)
    enddo

    ai = 0
    do i = 1,ltt%Nop
      do j = 1,ltt%sybase%nats
        ai = ai + 1
        if(ltt%base_format == "xyz")then
          !write(*,*)"Symmetry operations not implemented for xyz format"
          ltt%r_base(1,ai) =  ltt%bssym(1,i)*ltt%r_base(1,ai)
          ltt%r_base(1,ai) =  ltt%r_base(1,ai) + ltt%bsopload(i)*ltt%bstr(1,i)*ltt%lattice_vectors(1,1)
          ltt%r_base(1,ai) =  ltt%r_base(1,ai) + ltt%bsopload(i)*ltt%bstr(1,i)*ltt%lattice_vectors(2,1)
          ltt%r_base(1,ai) =  ltt%r_base(1,ai) + ltt%bsopload(i)*ltt%bstr(1,i)*ltt%lattice_vectors(3,1)

          ltt%r_base(2,ai) =  ltt%bssym(2,i)*ltt%r_base(2,ai)
          ltt%r_base(2,ai) =  ltt%r_base(2,ai) + ltt%bsopload(i)*ltt%bstr(2,i)*ltt%lattice_vectors(1,2)
          ltt%r_base(2,ai) =  ltt%r_base(2,ai) + ltt%bsopload(i)*ltt%bstr(2,i)*ltt%lattice_vectors(2,2)
          ltt%r_base(2,ai) =  ltt%r_base(2,ai) + ltt%bsopload(i)*ltt%bstr(2,i)*ltt%lattice_vectors(3,2)

          ltt%r_base(3,ai) =  ltt%bssym(3,i)*ltt%r_base(3,ai)
          ltt%r_base(3,ai) =  ltt%r_base(3,ai) + ltt%bsopload(i)*ltt%bstr(3,i)*ltt%lattice_vectors(1,3)
          ltt%r_base(3,ai) =  ltt%r_base(3,ai) + ltt%bsopload(i)*ltt%bstr(3,i)*ltt%lattice_vectors(2,3)
          ltt%r_base(3,ai) =  ltt%r_base(3,ai) + ltt%bsopload(i)*ltt%bstr(3,i)*ltt%lattice_vectors(3,3)
        else
          ltt%r_base(1,ai) =  ltt%bssym(1,i)*ltt%r_base(1,ai) + ltt%bsopload(i)*ltt%bstr(1,i)
          ltt%r_base(2,ai) =  ltt%bssym(2,i)*ltt%r_base(2,ai) + ltt%bsopload(i)*ltt%bstr(2,i)
          ltt%r_base(3,ai) =  ltt%bssym(3,i)*ltt%r_base(3,ai) + ltt%bsopload(i)*ltt%bstr(3,i)

          if(ltt%getOptTrs)then 
            if(i > 0 .and. j == 1)then
              call lcc_minimize_from(ltt%r_base(:,:),i,ai,ltt%sybase%nats,myTrs,verbose)
            endif
            ltt%r_base(:,ai) =  ltt%r_base(:,ai) + myTrs
          endif

        endif
      enddo
    enddo

    call prg_get_recip_vects(ltt%lattice_vectors,ltt%recip_vectors,ltt%volr,ltt%volk)

    massTot = sum(ltt%base_mass)

    call lcc_print_realVal("Volume of the cell",ltt%volr,"Ang^3",bld%verbose)
    call lcc_print_realVal("Total mass in the cell",massTot,"Au",bld%verbose)
    call lcc_print_realVal("Density",massTot/ltt%volr,"Au/Ang^3",bld%verbose)
    call lcc_print_realVal("Density",10.0_dp*massTot/(ltt%volr*6.022),"gr/cc",bld%verbose)

    if(check)then
      call lcc_print_message("Checking lattice structure...",bld%verbose)
      call lcc_check_basis(ltt%base_format,ltt%r_base,ltt%lattice_vectors,bld%verbose)
    endif

    if(ltt%setdensity > 0)then
      call lcc_print_message("Will scale lattive vetors to match input density ...",verbose)
      newVol = 1.0_dp/(6.022_dp*ltt%setdensity/(10.0_dp*massTot))
      scaleVect = (newVol/ltt%volr)**(1.0_dp/3.0_dp)
      ltt%lattice_vectors = scaleVect*ltt%lattice_vectors
    endif

  end subroutine lcc_read_base

  !> Routine to check for atom repetitions in basis
  !! \bnrief It will do all possible translations
  !! searching for atoms that could be repeated.
  !! \param base_format Basis format, if xyz of abc
  !! \param r_base Coordinates of the basis. r_base(1,7) means
  !! coordinate x of atom 7
  !! \param lattice_vectors Lattice vectors. WARNING, in this case
  !! lattice_vector(1,3) means the coordinate 3=z of vector 1.
  !!
  subroutine lcc_check_basis(base_format,r_base,lattice_vectors,verbose)
    implicit none
    integer :: Nbase,i,j,k,l,m
    real(dp), allocatable, intent(in) :: r_base(:,:),lattice_vectors(:,:)
    real(dp) :: r1(3),tr(3),distance
    character(len=*), intent(in) :: base_format
    integer, intent(in) :: verbose

    Nbase = size(r_base,dim=2)
    if(base_format == "abc")then
      !Check basis
      do i = 1,Nbase
        r1(:) = r_base(1,i)*lattice_vectors(1,:)
        r1(:) = r1 + r_base(2,i)*lattice_vectors(2,:)
        r1(:) = r1 + r_base(3,i)*lattice_vectors(3,:)
        do j = i+1, Nbase
          do k = -1,1
            do l = -1,1
              do m = -1,1
                tr(:) = r_base(1,j)*lattice_vectors(1,:)
                tr(:) = tr + r_base(2,j)*lattice_vectors(2,:)
                tr(:) = tr + r_base(3,j)*lattice_vectors(3,:)
                tr =  tr + k*lattice_vectors(1,:) + l*lattice_vectors(2,:) + m*lattice_vectors(3,:)
                distance = norm2(tr - r1)
                if(distance < 1.0d-1)then
                  write(*,*)"ERROR in the basis: An atom is equal to another by translation"
                  write(*,*)r_base(:,j),j
                  write(*,*)r_base(:,i),i
                  stop
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
    else
      write(*,*)"Check basis not yet implemented for xyz format"
    endif

  end subroutine lcc_check_basis

  !> Add a basis to the lattice.
  !! \brief This routine will add the basis to the system points previously
  !! cut from the lattice. This is the last step of the solid/shape/slab creation.
  !! \param ltt lattice_type See lcc_structs_mod
  !! \param sy system_type See progress library
  !!
  subroutine lcc_add_base_to_cluster(ltt,sy,verbose)
    integer, intent(in)    :: verbose
    integer                ::  cont, i, j
    integer                ::  natsCluster, natsBase
    integer, allocatable :: resindex_in(:)
    real(dp)               ::  x_component, y_component
    real(dp)               ::  z_component
    real(dp), allocatable  ::  coords_saved(:,:)
    real(dp), allocatable  ::  origin(:), r_cluster_in(:,:)
    real(dp) :: ran
    character(2), allocatable :: atom_in(:)
    character(3), allocatable :: resname_in(:)
    type(lattice_type) :: ltt
    type(system_type), intent(inout) :: sy
    type(rotation_type)               ::  rot
    integer, parameter :: seed=1234
    integer, allocatable :: seedin(:)
    integer :: ssize
    logical :: rotation

    call lcc_print_message("Adding basis to cluster ...",verbose)

    if(ltt%randrotations)then
      call lcc_print_message("Will apply random rotations to basis ...",verbose)
      call random_seed()
      call random_seed(size=ssize)
      allocate(seedin(ssize))
      seedin = seed
      call random_seed(PUT=seedin)
      allocate(coords_saved(3,ltt%nbase))
      coords_saved = ltt%r_base
      rot%patom1 = ltt%nbase
      rot%patom2 = 0
      rot%catom = 1
      rot%catom2 = 0
      rot%rotate_atoms(1) = 1; rot%rotate_atoms(2) = ltt%nbase
    endif

    natsBase = ltt%nbase
    natsCluster = ltt%nbase*sy%nats

    allocate(r_cluster_in(3,natsCluster))
    allocate(resindex_in(natsCluster))
    allocate(atom_in(natsCluster))
   
    if(ltt%base_format.eq.'abc')then
      cont = 0
      do i=1,sy%nats

        !Here we will add the posibility of a rotation
        if(ltt%randrotations)then
          call random_number(ran)
          ran = 2.0_dp*ran - 1.0_dp
          rot%v1(1)= ran
          call random_number(ran)
          ran = 2.0_dp*ran - 1.0_dp
          rot%v1(2)= ran
          call random_number(ran)
          ran = 2.0_dp*ran - 1.0_dp
          rot%v1(3)= ran
          !
          call random_number(ran)
          ran = 2.0_dp*ran - 1.0_dp
          rot%v2(1)= ran
          call random_number(ran)
          ran = 2.0_dp*ran - 1.0_dp
          rot%v2(2)= ran
          call random_number(ran)
          ran = 2.0_dp*ran - 1.0_dp
          rot%v2(3)= ran

          rot%vQ = 0.0_dp

          ltt%r_base = coords_saved
          call prg_rotate(rot,ltt%r_base,0)
        endif


        do j=1,natsBase
          cont = cont + 1
          r_cluster_in(1,cont) = sy%coordinate(1,i) + ltt%lattice_vectors(1,1)*ltt%r_base(1,j) + &
               &ltt%lattice_vectors(2,1)*ltt%r_base(2,j) + &
               &ltt%lattice_vectors(3,1)*ltt%r_base(3,j)

          r_cluster_in(2,cont) = sy%coordinate(2,i) + ltt%lattice_vectors(1,2)*ltt%r_base(1,j) + &
               &ltt%lattice_vectors(2,2)*ltt%r_base(2,j) + &
               &ltt%lattice_vectors(3,2)*ltt%r_base(3,j)

          r_cluster_in(3,cont) = sy%coordinate(3,i) + ltt%lattice_vectors(1,3)*ltt%r_base(1,j) + &
               &ltt%lattice_vectors(2,3)*ltt%r_base(2,j) + &
               &ltt%lattice_vectors(3,3)*ltt%r_base(3,j)

          atom_in(cont) = ltt%base_atom(j)
          resindex_in(cont) = ltt%resindex(j) + (i-1)*ltt%Nop !(i-1)*sy%nats
        enddo
      enddo

    elseif(ltt%base_format.eq.'xyz')then
      cont = 0
      do i=1,sy%nats

        if(ltt%randrotations)then
          call random_number(ran)
          ran = 2.0_dp*ran - 1.0_dp
          rot%v1(1)= ran
          call random_number(ran)
          ran = 2.0_dp*ran - 1.0_dp
          rot%v1(2)= ran
          call random_number(ran)
          ran = 2.0_dp*ran - 1.0_dp
          rot%v1(3)= ran
          !
          call random_number(ran)
          ran = 2.0_dp*ran - 1.0_dp
          rot%v2(1)= ran
          call random_number(ran)
          ran = 2.0_dp*ran - 1.0_dp
          rot%v2(2)= ran
          call random_number(ran)
          ran = 2.0_dp*ran - 1.0_dp
          rot%v2(3)= ran

          rot%vQ = 0.0_dp

          ltt%r_base = coords_saved
          call prg_rotate(rot,ltt%r_base,10)
        endif

        do j=1,natsBase
          cont = cont + 1
          r_cluster_in(1,cont) = sy%coordinate(1,i) + ltt%r_base(1,j)
          r_cluster_in(2,cont) = sy%coordinate(2,i) + ltt%r_base(2,j)
          r_cluster_in(3,cont) = sy%coordinate(3,i) + ltt%r_base(3,j)
          atom_in(cont) = ltt%base_atom(j)
          resindex_in(cont) = ltt%resindex(j) + (i-1)*ltt%Nop !+ (i-1)*sy%nats
        enddo
      enddo

    else
      write(*,*)'ERROR at lcc_add_base_to_cluster: Base Format not implemented'
      stop
    endif

    sy%nats = natsCluster
    call lcc_reallocate_realMat(sy%coordinate,3,natsCluster)
    sy%coordinate = r_cluster_in
    call lcc_reallocate_char2Vect(sy%symbol,natsCluster)
    sy%symbol = atom_in
    call lcc_reallocate_char3Vect(sy%atomname,natsCluster)
    sy%atomname = atom_in
    call lcc_reallocate_intVect(sy%resindex,natsCluster)
    sy%resindex = resindex_in

  end subroutine lcc_add_base_to_cluster

  !> Simple cubic (SC) lattice construction.
  !! \brief Constructs a "bulk" of Simple Cubic lattice.
  !! \param Nx1 Initial x lattice point.
  !! \param Nx2 Final x lattice point.
  !! \param Ny1 Initial y lattice point.
  !! \param Ny2 Final y lattice point.
  !! \param Nz1 Initial z lattice point.
  !! \param Nz2 Final z lattice point.
  !! \param h_lattice_a Lattice parameter.
  !! \param supra_lattice_vectors Lattice unit vectors of the resulting slab.
  !! \param r_sy Output system coordinates.
  !!
  subroutine lcc_sc(Nx1,Nx2,Ny1,Ny2,Nz1,Nz2,h_lattice_a,supra_lattice_vectors,r_sy)
    implicit none
    integer :: i,l,j,k,Ntop
    integer, intent(in) :: Nx1,Nx2,Ny1,Ny2,Nz1,Nz2
    real(dp), intent(in) :: h_lattice_a
    real(dp), allocatable, intent(inout) :: supra_lattice_vectors(:,:)
    real(dp), allocatable, intent(inout) :: r_sy(:,:)

    Ntop = (Nx2 - Nx1 + 1)*(Ny2 - Ny1 + 1)*(Nz2 - Nz1 + 1)

    allocate(r_sy(3,Ntop))

    l=0

    do i=Nx1,Nx2
      do j=Ny1,Ny2
        do k=Nz1,Nz2
          l=l+1
          r_sy(1,l)=(h_lattice_a)*i
          r_sy(2,l)=(h_lattice_a)*j
          r_sy(3,l)=(h_lattice_a)*k
        enddo
      enddo
    enddo

    call lcc_reallocate_realMat(supra_lattice_vectors,3,3)

    supra_lattice_vectors = 0.0_dp
    supra_lattice_vectors(1,1) = (Nx2 - Nx1 + 1)*h_lattice_a
    supra_lattice_vectors(2,2) = (Ny2 - Ny1 + 1)*h_lattice_a
    supra_lattice_vectors(3,3) = (Nz2 - Nz1 + 1)*h_lattice_a

  end subroutine lcc_sc

  !> Face center cubic (FCC) lattice construction.
  !! \brief Constructs a "bulk" of Face center cubic lattice.
  !! \param Nx1 Initial x lattice point.
  !! \param Nx2 Final x lattice point.
  !! \param Ny1 Initial y lattice point.
  !! \param Ny2 Final y lattice point.
  !! \param Nz1 Initial z lattice point.
  !! \param Nz2 Final z lattice point.
  !! \param h_lattice_a Lattice parameter.
  !! \param supra_lattice_vectors Lattice unit vectors of the resulting slab.
  !! \param r_sy Output system coordinates.
  !!
  subroutine lcc_fcc(Nx1,Nx2,Ny1,Ny2,Nz1,Nz2,h_lattice_a,supra_lattice_vectors,&
       &r_sy,verbose)
    implicit none
    integer :: i,l,j,k,Ntop
    integer, intent(in) :: Nx1,Nx2,Ny1,Ny2,Nz1,Nz2,verbose
    real(dp), intent(in) :: h_lattice_a
    real(dp), allocatable, intent(inout) :: supra_lattice_vectors(:,:)
    real(dp), allocatable, intent(inout) :: r_sy(:,:)

    call lcc_print_message("Making fcc lattice ...",verbose)

    Ntop = (Nx2 - Nx1 + 1)*(Ny2 - Ny1 + 1)*(Nz2 - Nz1 + 1)

    allocate(r_sy(3,Ntop))

    l=0

    do i=Nx1,Nx2
      do j=Ny1,Ny2
        do k=Nz1,Nz2
          l=l+1
          r_sy(1,l)=(h_lattice_a/2.0d0)*j + (h_lattice_a/2.0d0)*k
          r_sy(2,l)=(h_lattice_a/2.0d0)*i + (h_lattice_a/2.0d0)*k
          r_sy(3,l)=(h_lattice_a/2.0d0)*i + (h_lattice_a/2.0d0)*j
        enddo
      enddo
    enddo

    call lcc_reallocate_realMat(supra_lattice_vectors,3,3)

    supra_lattice_vectors = 0.0_dp
    supra_lattice_vectors(1,1) = (Nx2 - Nx1 + 1)*h_lattice_a
    supra_lattice_vectors(2,2) = (Ny2 - Ny1 + 1)*h_lattice_a
    supra_lattice_vectors(3,3) = (Nz2 - Nz1 + 1)*h_lattice_a

  end subroutine lcc_fcc


  !> Triclinic lattice construction.
  !! \brief Constructs a "bulk" of triclinic lattice.
  !! \param Nx1 Initial x lattice point.
  !! \param Nx2 Final x lattice point.
  !! \param Ny1 Initial y lattice point.
  !! \param Ny2 Final y lattice point.
  !! \param Nz1 Initial z lattice point.
  !! \param Nz2 Final z lattice point.
  !! \param lattice_vectors Lattice vectors.
  !! \param supra_lattice_vectors Lattice unit vectors of the resulting slab.
  !! \param r_sy Output system coordinates.
  !! Note: Unit cell representation has to be transformed from
  !! edges and angles to vetors before calling this routine.
  !! \todo A angles_to_vectors transformation will be available.
  !!
  subroutine lcc_triclinic(Nx1,Nx2,Ny1,Ny2,Nz1,Nz2,lattice_vectors,&
       &supra_lattice_vectors,r_sy,verbose)

    implicit none
    integer :: i,Ns,cont,l,j,k,Ntop
    integer, intent(in) :: Nx1,Nx2,Ny1,Ny2,Nz1,Nz2,verbose
    real(dp), allocatable, intent(inout) :: supra_lattice_vectors(:,:)
    real(dp), allocatable, intent(in) :: lattice_vectors(:,:)
    real(dp), allocatable, intent(inout) :: r_sy(:,:)

    l=0

    if(.not. allocated(supra_lattice_vectors)) then
      allocate(supra_lattice_vectors(3,3))
    endif

    Ntop = (Nx2 - Nx1 + 1)*(Ny2 - Ny1 + 1)*(Nz2 - Nz1 + 1)
    if(.not.(allocated(r_sy)))allocate(r_sy(3,Ntop))

    call lcc_print_message("Using (Nx1,Nx2,Ny1,Ny2,Nz1,Nz2)",verbose)

    do i=Nx1,Nx2
      do j=Ny1,Ny2
        do k=Nz1,Nz2
          l=l+1
          r_sy(1,l)=i*lattice_vectors(1,1) + j*lattice_vectors(2,1) + k*lattice_vectors(3,1)
          r_sy(2,l)=i*lattice_vectors(1,2) + j*lattice_vectors(2,2) + k*lattice_vectors(3,2)
          r_sy(3,l)=i*lattice_vectors(1,3) + j*lattice_vectors(2,3) + k*lattice_vectors(3,3)
        enddo
      enddo
    enddo

    call lcc_reallocate_realMat(supra_lattice_vectors,3,3)

    supra_lattice_vectors(1,:) = (Nx2-Nx1 + 1)*lattice_vectors(1,:)
    supra_lattice_vectors(2,:) = (Ny2-Ny1 + 1)*lattice_vectors(2,:)
    supra_lattice_vectors(3,:) = (Nz2-Nz1 + 1)*lattice_vectors(3,:)

  end subroutine lcc_triclinic


  !> Sets the atom type.
  !! \brief Sets the atom "symbol/type/name."
  !! \param a_type Atom symbol character.
  !! \param atom_symbol Atom symbols.
  !! \param atom_name Atom name.
  !! Note: Atom name is a tag that can distinguish
  !! atoms with same symbol.
  !!
  subroutine lcc_set_atom_type(a_type,atom_symbol,atom_name,nats)
    character(len=2), intent(in) :: a_type
    character(len=2), allocatable, intent(inout) :: atom_symbol(:)
    character(len=3), allocatable, intent(inout) :: atom_name(:)
    integer, intent(in) :: nats

    if(.not. allocated(atom_symbol))allocate(atom_symbol(nats))
    if(.not. allocated(atom_name))allocate(atom_name(nats))
    atom_symbol = a_type
    atom_name = a_type

  end subroutine lcc_set_atom_type


  !> Will add randomness to the system.
  !! \param r_inout System coordinates.
  !! \param lattice_vectors Lattice vectors.
  !! \param seed Random seed.
  !! rcoeff Coefficient for randomness.
  !!
  subroutine lcc_add_randomness(r_inout,lattice_vectors,seed,rcoeff)
    implicit none
    real(dp),allocatable,intent(inout) :: r_inout(:,:)
    real(dp), allocatable,intent(in) :: lattice_vectors(:,:)
    real(dp), intent(in) :: rcoeff
    real(dp) :: aux1,aux2,aux3,ran
    integer, intent(in) :: seed
    integer, allocatable :: seedin(:)
    integer :: nats, l,ssize

    call random_seed()
    call random_seed(size=ssize)
    allocate(seedin(ssize))
    seedin = seed
    call random_seed(PUT=seedin)

    nats = size(r_inout,dim=2)

    do l=1,nats
      call random_number(ran)
      aux1 = lattice_vectors(1,1)*rcoeff*(2.0_dp*ran - 1.0_dp)
      call random_number(ran)
      aux2 = lattice_vectors(2,1)*rcoeff*(2.0_dp*ran - 1.0_dp)
      call random_number(ran)
      aux3 = lattice_vectors(3,1)*rcoeff*(2.0_dp*ran - 1.0_dp)
      r_inout(1,l) = r_inout(1,l) + aux1 + aux2 + aux3

      call random_number(ran)
      aux1 = lattice_vectors(1,2)*rcoeff*(2.0_dp*ran - 1.0_dp)
      call random_number(ran)
      aux2 = lattice_vectors(2,2)*rcoeff*(2.0_dp*ran - 1.0_dp)
      call random_number(ran)
      aux3 = lattice_vectors(3,2)*rcoeff*(2.0_dp*ran - 1.0_dp)
      r_inout(2,l) = r_inout(2,l) + aux1 + aux2 + aux3

      call random_number(ran)
      aux1 = lattice_vectors(1,3)*rcoeff*(2.0_dp*ran - 1.0_dp)
      call random_number(ran)
      aux2 = lattice_vectors(2,3)*rcoeff*(2.0_dp*ran - 1.0_dp)
      call random_number(ran)
      aux3 = lattice_vectors(3,3)*rcoeff*(2.0_dp*ran - 1.0_dp)
      r_inout(3,l) = r_inout(3,l) + aux1 + aux2 + aux3
    enddo

  end subroutine lcc_add_randomness

  !> To get the best translation that minimizes the 
  !! distance to any previous fragment.
  !! \param xVar Coordinates of the full basis (including symmetry operations).
  !! \param i Fragmet being added at the "i" operation.
  !! \param ai Atom index to translate and get the optimal translation.
  !! \param nats Number of atoms in the fragment.
  !! \param trs Optimal translation.
  !!
  subroutine lcc_minimize_from(xVar,i,ai,nats,trs,verbose)

    implicit none
    real(dp) :: dist, minDist
    real(dp) ::  aux(3), myTrs(3)
    real(dp) :: xRef(3)
    real(dp), intent(inout) :: trs(3)
    real(dp), intent(in) :: xVar(:,:)
    integer, intent(in) :: i,ai,nats,verbose
    integer :: myK,k,l,m,ii

    call lcc_print_message("Getting optimal translation ...",verbose)

    minDist = 10000.0_dp
    do k = -1,1 !All possible translations
      do l = -1,1
        do m = -1,1
          do ii = 1,i !Previous equivalent atom
            xRef = xVar(:,ai - ii*nats)
            trs(1) = real(k,dp)
            trs(2) = real(l,dp)
            trs(3) = real(m,dp)
            aux =  xVar(:,ai) + trs(:)
            dist = norm2(aux - xRef)
            if(dist <= minDist)then
              myTrs = trs
              minDist = dist
            endif
          enddo
        enddo
      enddo
    enddo
    trs = myTrs
    if(verbose >= 1)then 
      call lcc_print_intVal("Optimal translations for operation",i,"",verbose)
      write(*,*)i,trs
    endif

  end subroutine lcc_minimize_from

end module lcc_lattice_mod
