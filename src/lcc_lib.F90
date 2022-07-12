!> Library module
!!
module lcc_lib

  use prg_timer_mod
  use prg_extras_mod
  use prg_system_mod
  use prg_syrotation_mod
  use lcc_message_mod
  use lcc_parser_mod
  use lcc_structs_mod
  use lcc_lattice_mod
  use lcc_aux_mod
  use lcc_build_mod
  use lcc_regular_mod
  use lcc_check_mod

  implicit none

  private 

  public :: lcc

  contains 

  subroutine lcc(readInputFile,inputFileName,syOut,writeOut,&
      &clType,planeIn)

  logical, intent(in) ::  readInputFile, writeOut
  character(len=60), intent(in) ::  inputFileName
  character(len=*), intent(in) ::  clType
  type(system_type),intent(inout) ::  syOut
  real(dp), allocatable, intent(in) :: planeIn(:,:)

  integer :: mlsI
  character(20) :: inputfile
  type(build_type) :: bld
  type(lattice_type) :: ltt
  type(system_type) :: sy
  real(dp), allocatable :: r(:,:),tmp_vectors(:,:)
  type(rotation_type) :: rot

  call getarg(1, inputfile)

  !If the argument passed is empty print ussage
  !and build a sample imput file.
  if(inputfile == "")then
    call lcc_print_ussage()
    call lcc_make_sample_input()
    stop
  end if

  !Read the build and lattice type
  mlsI = mls()
  call lcc_parse(inputfile,bld,ltt)
  call lcc_print_realVal("Time for lcc_parse",mls() - mlsI,"ms",bld%verbose)

  !Make lattice from unit cell.
  mlsI = mls()
  call lcc_make_lattice(bld,ltt,ltt%check,sy)

  call lcc_print_realVal("Time for make_lattice",mls() - mlsI,"ms",bld%verbose)

  !If the shape has to be cut from a "dressed lattice"
  if(bld%use_lattice_base.eq.'T'.and. bld%cut_with_base.eq.'T') then
    mlsI = mls()
    call lcc_add_base_to_cluster(ltt,sy)
    call lcc_print_realVal("Time for add_base_to_cluster",mls() - mlsI,"ms",bld%verbose)
  endif

  !Cout out the shape
  select case (bld%cl_type)
  case ('Bulk')
    call lcc_print_message("Shape will just be the bulk slab as constructed ...",bld%verbose)
  case ('BravaisGrowth')
    call lcc_bravais_growth(bld%niter,bld%rtol,bld%r_cut,bld%MaxCoordination,&
      &bld%seed_file,sy%coordinate)
    sy%nats = size(sy%coordinate,dim=2)
  case ('Spheroid')
    call lcc_spheroid(bld%a_axis,bld%b_axis,bld%c_axis,sy%coordinate)
    sy%nats = size(sy%coordinate,dim=2)
  case ("Planes")
    call lcc_plane_cut(bld%planes,bld%ploads,bld%interPlanarDistances,&
      &ltt%lattice_vectors,sy%lattice_vector,&
      &sy%resindex,sy%coordinate,bld%verbose)
    sy%nats = size(sy%coordinate,dim=2)
  case default
    call lcc_print_error("main","The requested TypeOfCluster is not implemented")
  end select
  
  if(bld%checkperiod)then
    call lcc_check_periodicity(sy%coordinate,sy%lattice_vector,ltt%bulk,0.0001_dp,bld%verbose)
  endif

  !If there is no lattice basis, we just set the atom
  !for every point in the lattice.
  if(bld%use_lattice_base.eq.'F') then
    call lcc_set_atom_type(bld%a_type,sy%symbol,sy%atomname,sy%nats)
  endif

  !If the shape has to be cut from the lattice but the basis has to be added a posteriori.
  if(bld%use_lattice_base.eq.'T'.and. bld%cut_with_base.eq.'F') then
    mlsI = mls()
    call lcc_add_base_to_cluster(ltt,sy)
    call lcc_print_realVal("Time for add_base_to_cluster",mls() - mlsI,"ms",bld%verbose)
  endif

  if(bld%reorient)then  
    call lcc_center_at_origin(sy%coordinate,bld%verbose)
    call lcc_canonical_basis(sy%lattice_vector,sy%coordinate,bld%verbose)
  endif

  if(bld%center)then
    call lcc_print_message("Centering at box ...",bld%verbose)
    call lcc_center_at_box(sy%lattice_vector,sy%coordinate,bld%verbose)
  endif

  !Write system
  call lcc_write_coords(sy,bld,"coords",bld%verbose)

  call lcc_print_message("CCL execution finilized!",1)

end subroutine lcc

end module lcc_lib
