!> Program for building periodic crystal slabs or shapes.
!!
program lcc_main

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
  use lcc_radial_mod

  implicit none
  integer :: mlsI
  character(50) :: inputfile
  type(build_type) :: bld
  type(lattice_type) :: ltt
  type(system_type) :: sy
  real(dp), allocatable :: r(:,:),tmp_vectors(:,:)
  type(rotation_type) :: rot

  call getarg(1, inputfile)

  !If the argument passed is empty, print ussage
  !and build a sample imput file.
  if(inputfile == "")then
    call lcc_print_ussage()
    call lcc_make_sample_input()
    stop
  end if

  !Print help info is argument is ~ help
  if(inputfile == "h" .or. inputfile == "-h" .or. &
       &inputfile == "--h" .or. inputfile == "help" .or. &
       &inputfile == "-help" .or. inputfile == "--help")then
    call lcc_help()
    stop
  endif

  !Read the build and lattice types
  mlsI = mls()
  call lcc_parse(inputfile,bld,ltt)
  call lcc_print_realVal("Time for lcc_parse",mls() - mlsI,"ms",bld%verbose)

  !Make lattice from unit cell
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
    if(bld%planes_type == "Miller")then
      call lcc_plane_cut(bld%planes,bld%ploads,bld%pgncheck,bld%interPlanarDistances,&
           &ltt%lattice_vectors,sy%lattice_vector,&
           &sy%resindex,sy%coordinate,bld%verbose)
      sy%nats = size(sy%coordinate,dim=2)
    elseif(bld%planes_type == "Regular")then
      call lcc_print_error("main","The requested PlanesType is not implemented yet")
    else
      call lcc_print_error("main","The requested PlanesType is not implemented. Available option are:&
           &Miller and Regular")
    endif
  case ("Slab")
    call lcc_build_slab(bld%slab,bld%sloads,ltt%lattice_vectors,&
       &sy%lattice_vector,sy%resindex,sy%coordinate,bld%verbose)
      sy%nats = size(sy%coordinate,dim=2)
  case default
    call lcc_print_error("main","The requested TypeOfCluster is not implemented. Available option are:&
         &Bulk, BravaisGrowth, Spheroid, Planes")
  end select

  if(bld%checkperiod)then
    call lcc_check_periodicity(sy%coordinate,sy%lattice_vector,ltt%bulk,0.0001_dp,bld%verbose)
  endif

  !If there is no lattice basis, we just set the atom
  !for every point in the lattice.
  if(bld%use_lattice_base.eq.'F') then
    call lcc_set_atom_type(bld%a_type,sy%symbol,sy%atomname,sy%nats)
  endif
  
  !If the shape has to be cut from the lattice but the basis has to be added a posteriori
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
    call lcc_center_at_box(sy%lattice_vector,sy%coordinate,bld%verbose)
  endif

  if(bld%randomCoordinates)then
    call lcc_add_randomness_to_coordinates(sy%coordinate,bld%seed,bld%rcoeff)
  endif

  !Write system
  call lcc_write_coords(sy,bld,"coords",bld%verbose)

  if(bld%rdfPair .ne. " ") call lcc_distance_matrix(sy%symbol,sy%coordinate,sy%lattice_vector,bld%rdfPair,bld%verbose)

  call lcc_print_message("CCL execution finalized!",bld%verbose)

end program lcc_main
