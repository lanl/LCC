!> Module for generating the shapes after lattice is constructed.
!!
module lcc_build_mod

  use lcc_constants_mod
  use lcc_lattice_mod
  use lcc_allocation_mod
  use prg_system_mod
  use lcc_aux_mod

  implicit none

  public :: lcc_bravais_growth, lcc_plane_cut, lcc_add_randomness_to_coordinates

contains

  !> For "growing" a crystal shape using Bravias
  !! type of growth teory.
  !! \param nCycles Number of shells to add.
  !! \param dTol Tolerance for distinguising the coordinates
  !! from the seed to the coodinates from the bulk.
  !! \param dTo Parameter to determine the coordination the incoming
  !! atom.
  !! \param tCoordination Target coordination. If coodination is larger
  !! than the target, the atom will be picked.
  !! \param seed_file Name of the file containing the seed.
  !! \param r_inout Input: Bulk lattice, Output: Crystal shape.
  !! \todo Optimize the routine.
  !!
  subroutine lcc_bravais_growth(nCycles,dTol,dTo,tCoordination,seed_file,r_inout)
    implicit none
    real(dp), allocatable, intent(inout) :: r_inout(:,:)
    real(dp), intent(in) :: dTol,dTo
    real(dp) :: d
    real(dp), allocatable :: r_tmp(:,:), r_seed(:,:)
    integer, intent(in) :: nCycles,tCoordination
    integer :: i,j,nBulk,nSeed,nTmp,cy,nCoordination
    logical, allocatable :: inSeed(:),inShell(:)
    character(len=*), intent(in) :: seed_file
    type(system_type) :: seed

    !Parsing with progress library
    call prg_parse_system(seed, trim(adjustl(seed_file)))

    write(*,*)"In bravais growth ..."

    nBulk = size(r_inout,dim=2)
    nSeed = seed%nats
    allocate(inSeed(nBulk))
    allocate(inShell(nBulk))
    allocate(r_tmp(3,nBulk))
    call lcc_reallocate_realMat(r_seed,3,nSeed)
    r_seed = seed%coordinate

    do cy=1,nCycles

      !Get the points that are not in (or very near) the points in r_in
      inSeed = .false.

      !$omp parallel do default (none) &
      !$omp shared(nBulk,nSeed,dTol,inSeed,r_seed,r_inout) &
      !$omp private(i,j,d)
      do i=1,nBulk
        do j=1,nSeed
          d = norm2(r_seed(:,j) - r_inout(:,i))
          if(d < dTol) inSeed(i) = .true.
        enddo
      enddo
      !$omp end parallel do

      r_tmp(:,1:nSeed) = r_seed(:,:)

      nTmp = nSeed

      !$omp parallel do default (none) &
      !$omp shared(nBulk,nSeed,dTo,inSeed,r_seed,r_inout) &
      !$omp shared(tCoordination,inShell) &
      !$omp private(i,j,d,nCoordination)
      do i=1,nBulk
        if(.not.inSeed(i))then
          nCoordination = 0
          do j=1,nSeed
            d = norm2(r_seed(:,j) - r_inout(:,i))
            if(d < dTo)then
              nCoordination = nCoordination + 1
            endif
          enddo
          if(nCoordination >= tCoordination)then
            inShell(i) = .true.
            inSeed(i) = .true.
          endif
        endif
      enddo
      !$omp end parallel do

      do i=1,nBulk
        if(inShell(i))then
          nTmp = nTmp + 1
          r_tmp(:,nTmp) = r_inout(:,i)
        endif
      enddo

      nSeed = nTmp
      call lcc_reallocate_realMat(r_seed,3,nSeed)
      r_seed = r_tmp(:,1:nTmp)
      write(*,*)"Cycle",cy,"Number of atoms",nTmp
    enddo

    call lcc_reallocate_realMat(r_inout,3,nTmp)
    r_inout(:,:) = r_tmp(:,1:nTmp)
    deallocate(inSeed)
    deallocate(inShell)
    deallocate(r_tmp)
    deallocate(r_seed)

  end subroutine lcc_bravais_growth

  !> Cutting a shape based on Miller planes.
  !! \brief A set of panes and distances is provided.
  !! \param planes List of planes to cut the shape with.
  !! \param ploads Distance from the origin to locate the plane.
  !! \param interPlanarDistance Use "interplanar distances" as measure for the cut.
  !! \param lattice_vectors Lattice vectors.
  !! \param cluster_lattice_vectors Lattice vectors of the shape.
  !! Note: this only makes sense if the planes make a parellelepiped.
  !! \param r_inout Coordinates in and out.
  !! \param verbose Verbosity level.
  !!
  subroutine lcc_plane_cut(planes,ploads,interPlanarDistances,lattice_vectors,&
       &cluster_lattice_vectors,resindex,r_inout,verbose)
    implicit none
    integer :: accept_index, j, l, nats
    real(dp) :: diff(3), dotprod, dr, modhkl
    real(dp) :: volk, volr, x_component
    real(dp) ::  y_component, z_component
    real(dp), allocatable  ::  abc_angles(:,:), plane_v(:,:), recip_vectors(:,:)
    real(dp) :: xBox(2),yBox(2),zBox(2)
    integer, allocatable :: resindex(:)
    integer :: i,Nplanes,Ntop,intAux,nearPlane
    real(dp), allocatable, intent(in) :: planes(:,:)
    real(dp), allocatable, intent(in) :: ploads(:)
    real(dp), allocatable, intent(inout) :: r_inout(:,:)
    real(dp), allocatable, intent(inout) :: lattice_vectors(:,:)
    real(dp), allocatable, intent(inout) :: cluster_lattice_vectors(:,:)
    real(dp), allocatable :: cluster_recip_vectors(:,:)
    real(dp), allocatable :: mod_cluster_lattice_vect(:)
    real(dp), allocatable :: myvect(:),vplanes(:,:),intplanes(:,:)
    real(dp), allocatable :: r_tmp(:,:),period(:)
    real(dp) :: PP,QP,t,density,area
    integer, intent(in) :: verbose
    logical, intent(in) :: interPlanarDistances
    logical :: inSurface

    call lcc_print_message("Cutting by Miller planes ...",verbose)

    allocate(recip_vectors(3,3))
    call prg_get_recip_vects(lattice_vectors,recip_vectors,volr,volk)

    call lcc_print_realMat("Reciprocal vector b1 ",recip_vectors,"[1/Ang] ",verbose)

    Nplanes = size(planes,dim=2)
    allocate(vplanes(3,Nplanes))
    allocate(intplanes(3,Nplanes))
    allocate(period(Nplanes))

    do j=1,Nplanes
      x_component =  planes(1,j)*recip_vectors(1,1) + planes(2,j)*recip_vectors(2,1) &
           &+ planes(3,j)*recip_vectors(3,1)
      y_component =  planes(1,j)*recip_vectors(1,2) + planes(2,j)*recip_vectors(2,2) &
           &+ planes(3,j)*recip_vectors(3,2)
      z_component =  planes(1,j)*recip_vectors(1,3) + planes(2,j)*recip_vectors(2,3) &
           &+ planes(3,j)*recip_vectors(3,3)

      modhkl = sqrt(x_component**2 + y_component**2 + z_component**2)
      period(j) =  2.0_dp*pi/modhkl

      if(interPlanarDistances)then
        call lcc_print_message("Using number of paralell planes as a distance measure ...",verbose)
        vplanes(1,j) =  period(j)*ploads(j)*x_component/modhkl
        vplanes(2,j) =  period(j)*ploads(j)*y_component/modhkl
        vplanes(3,j) =  period(j)*ploads(j)*z_component/modhkl
      else
        call lcc_print_message("Using euclidean distance as distance measure ...",verbose)
        vplanes(1,j) =  ploads(j)*x_component/modhkl
        vplanes(2,j) =  ploads(j)*y_component/modhkl
        vplanes(3,j) =  ploads(j)*z_component/modhkl
      endif

      intplanes(1,j) =  period(j)*floor(ploads(j))*x_component/modhkl
      intplanes(2,j) =  period(j)*floor(ploads(j))*y_component/modhkl
      intplanes(3,j) =  period(j)*floor(ploads(j))*z_component/modhkl

      write(*,*)"Plane:",planes(1,j),planes(2,j),planes(3,j)
      call lcc_print_realVal("Period in direction (Plane distance)", &
           &2.0_dp*pi/modhkl," Ang ",verbose)
      !call lcc_get_reticular_density(lattice_vectors,planes(:,j),density)
      area = volr/(2.0_dp*pi/modhkl)
      call lcc_print_realVal("Reticular area (V/d_hkl)", &
           &area,"Ang^2 ",verbose)
      call lcc_print_realVal("Reticular density", &
           &1.0_dp/area," points/Ang^2 ",verbose)

    enddo
    call lcc_reallocate_realMat(cluster_lattice_vectors,3,3)
    call lcc_reallocate_realMat(cluster_recip_vectors,3,3)
    call lcc_reallocate_realVect(mod_cluster_lattice_vect,3)


    if(mod(Nplanes,2) == 0 .and. Nplanes == 6)then
      call lcc_print_message("Computing boundaries assuming contiguous &
           &top and bottom planes (six planes forming a parallelepiped)",verbose)
      cluster_recip_vectors(1,:) = vplanes(:,1)
      cluster_recip_vectors(2,:) = vplanes(:,3)
      cluster_recip_vectors(3,:) = vplanes(:,5)

      cluster_recip_vectors(1,:)= cluster_recip_vectors(1,:)/norm2(cluster_recip_vectors(1,:))
      cluster_recip_vectors(2,:)= cluster_recip_vectors(2,:)/norm2(cluster_recip_vectors(2,:))
      cluster_recip_vectors(3,:)= cluster_recip_vectors(3,:)/norm2(cluster_recip_vectors(3,:))

      call prg_get_recip_vects(cluster_recip_vectors,cluster_lattice_vectors,volr,volk)
      PP = norm2(vplanes(:,1))**2
      QP = dot_product(cluster_lattice_vectors(1,:),vplanes(:,1))
      t = PP/QP
      cluster_lattice_vectors(1,:) = t*cluster_lattice_vectors(1,:)
      PP = norm2(vplanes(:,2))**2
      QP = dot_product(cluster_lattice_vectors(1,:),vplanes(:,2))
      t = PP/QP
      cluster_lattice_vectors(1,:) = cluster_lattice_vectors(1,:) - t*cluster_lattice_vectors(1,:)

      PP = norm2(vplanes(:,3))**2
      QP = dot_product(cluster_lattice_vectors(2,:),vplanes(:,3))
      t = PP/QP
      cluster_lattice_vectors(2,:) = t*cluster_lattice_vectors(2,:)
      PP = norm2(vplanes(:,4))**2
      QP = dot_product(cluster_lattice_vectors(2,:),vplanes(:,4))
      t = PP/QP
      cluster_lattice_vectors(2,:) = cluster_lattice_vectors(2,:) - t*cluster_lattice_vectors(2,:)

      PP = norm2(vplanes(:,5))**2
      QP = dot_product(cluster_lattice_vectors(3,:),vplanes(:,5))
      t = PP/QP
      cluster_lattice_vectors(3,:) = t*cluster_lattice_vectors(3,:)
      PP = norm2(vplanes(:,6))**2
      QP = dot_product(cluster_lattice_vectors(3,:),vplanes(:,6))
      t = PP/QP
      cluster_lattice_vectors(3,:) = cluster_lattice_vectors(3,:) - t*cluster_lattice_vectors(3,:)

    endif

    dr=0.1d0
    l=0
    Ntop = size(r_inout,dim=2)
    do i=1,Ntop  !Cut
      accept_index=0.0
      do j=1,Nplanes
        diff(:)=r_inout(:,i)-vplanes(:,j)
        dotprod = diff(1)*vplanes(1,j) + diff(2)*vplanes(2,j) + diff(3)*vplanes(3,j)
        if(dotprod.lt.0.0_dp)then
          accept_index=accept_index+1
        endif
      enddo
      if(accept_index.eq.Nplanes)l = l+1
    enddo

    nats = l

    call lcc_reallocate_realMat(r_tmp,3,Ntop)
    r_tmp = r_inout

    call lcc_reallocate_realMat(r_inout,3,nats)

    if(allocated(resindex)) deallocate(resindex)
    allocate(resindex(nats))

    l=0
    do i=1,Ntop  !Cut
      accept_index=0.0
      inSurface = .false.
      do j=1,Nplanes
        diff(:)=r_tmp(:,i)- vplanes(:,j)
        dotprod = diff(1)*vplanes(1,j) + diff(2)*vplanes(2,j) + diff(3)*vplanes(3,j)
        if(dotprod.lt.0.0_dp)then
          accept_index=accept_index+1
        endif
        dotprod = dotprod/(norm2(diff))
        if(dotprod.lt.0.0_dp .and. dotprod.gt.-5.01_dp) then
          inSurface = .true.
          nearPlane = j
        endif

      enddo
      if(accept_index.eq.Nplanes)then
        l = l+1
        r_inout(1,l) = r_tmp(1,i)
        r_inout(2,l) = r_tmp(2,i)
        r_inout(3,l) = r_tmp(3,i)
        xBox(1) = min(r_inout(1,l),xBox(1))
        xBox(2) = max(r_inout(1,l),xBox(2))
        yBox(1) = min(r_inout(2,l),yBox(1))
        yBox(2) = max(r_inout(2,l),yBox(2))
        zBox(1) = min(r_inout(3,l),zBox(1))
        zBox(2) = max(r_inout(3,l),zBox(2))
        if(inSurface) resindex(l) = nearPlane
      endif
    enddo

    if(.not.(mod(Nplanes,2) == 0 .and. Nplanes == 6))then
      call lcc_print_message("Computing boundaries assuming irregular &
           &shape",verbose)
      cluster_lattice_vectors = 0.0_dp
      cluster_lattice_vectors(1,1) = xBox(2) - xBox(1)
      cluster_lattice_vectors(2,2) = yBox(2) - yBox(1)
      cluster_lattice_vectors(3,3) = zBox(2) - zBox(1)
    endif

    deallocate(r_tmp)

  end subroutine lcc_plane_cut

  !> Cutting a shape based on PBC vectors.
  !! \brief A set of PBC vectors and distances is provided.
  !! \param planes List of planes to cut the shape with.
  !! \param ploads Distance from the origin to locate the plane.
  !! \param interPlanarDistance Use "interplanar distances" as measure for the cut.
  !! \param lattice_vectors Lattice vectors.
  !! \param cluster_lattice_vectors Lattice vectors of the shape.
  !! Note: this only makes sense if the planes make a parellelepiped.
  !! \param r_inout Coordinates in and out.
  !! \param verbose Verbosity level.
  !!
  subroutine lcc_build_slab(slab,sloads,lattice_vectors,&
       &cluster_lattice_vectors,resindex,r_inout,verbose)
    implicit none
    integer :: accept_index, j, l, nats
    real(dp) :: diff(3), dotprod, dr, normVect
    real(dp) :: volk, volr, x_component
    real(dp) ::  y_component, z_component
    real(dp), allocatable  ::  abc_angles(:,:), plane_v(:,:), recip_vectors(:,:)
    real(dp) :: xBox(2),yBox(2),zBox(2)
    integer, allocatable :: resindex(:)
    integer :: i,Nslab,Ntop,intAux,nearPlane
    real(dp), allocatable, intent(in) :: slab(:,:)
    real(dp), allocatable, intent(in) :: sloads(:)
    real(dp), allocatable, intent(inout) :: r_inout(:,:)
    real(dp), allocatable, intent(inout) :: lattice_vectors(:,:)
    real(dp), allocatable, intent(inout) :: cluster_lattice_vectors(:,:)
    real(dp), allocatable :: mod_cluster_lattice_vect(:)
    real(dp), allocatable :: r_tmp(:,:),pvect(:,:)
    real(dp) :: matM(3,3),invM(3,3),det
    integer, intent(in) :: verbose
    logical :: inSurface

    call lcc_print_message("Cutting a slab (In the direct lattice)...",verbose)

    Nslab = size(slab,dim=2)
    allocate(pvect(3,Nslab))

    do j=1,Nslab
      normVect = norm2(slab(:,j))
      x_component =  slab(1,j)*sloads(j)/normVect
      y_component =  slab(2,j)*sloads(j)/normVect
      z_component =  slab(3,j)*sloads(j)/normVect
      pvect(1,j) =  x_component
      pvect(2,j) =  y_component
      pvect(3,j) =  z_component
      matM(1,j) = x_component
      matM(2,j) = y_component
      matM(3,j) = z_component
    enddo

    det =     matM(1,1)*matM(2,2)*matM(3,3)  &
         - matM(1,1)*matM(2,3)*matM(3,2)  &
         - matM(1,2)*matM(2,1)*matM(3,3)  &
         + matM(1,2)*matM(2,3)*matM(3,1)  &
         + matM(1,3)*matM(2,1)*matM(3,2)  &
         - matM(1,3)*matM(2,2)*matM(3,1)

    invM(1,1) = +(matM(2,2)*matM(3,3)-matM(2,3)*matM(3,2))
    invM(1,2) = -(matM(2,1)*matM(3,3)-matM(2,3)*matM(3,1))
    invM(1,3) = +(matM(2,1)*matM(3,2)-matM(2,2)*matM(3,1))
    invM(2,1) = -(matM(1,2)*matM(3,3)-matM(1,3)*matM(3,2))
    invM(2,2) = +(matM(1,1)*matM(3,3)-matM(1,3)*matM(3,1))
    invM(2,3) = -(matM(1,1)*matM(3,2)-matM(1,2)*matM(3,1))
    invM(3,1) = +(matM(1,2)*matM(2,3)-matM(1,3)*matM(2,2))
    invM(3,2) = -(matM(1,1)*matM(2,3)-matM(1,3)*matM(2,1))
    invM(3,3) = +(matM(1,1)*matM(2,2)-matM(1,2)*matM(2,1))

    invM = transpose(invM) / DET

    dr=0.1d0
    l=0
    Ntop = size(r_inout,dim=2)

    do i=1,Ntop  !Cut
      diff(1)=invM(1,1)*r_inout(1,i) + invM(1,2)*r_inout(2,i) + invM(1,3)*r_inout(3,i)
      diff(2)=invM(2,1)*r_inout(1,i) + invM(2,2)*r_inout(2,i) + invM(2,3)*r_inout(3,i)
      diff(3)=invM(3,1)*r_inout(1,i) + invM(3,2)*r_inout(2,i) + invM(3,3)*r_inout(3,i)
      if(diff(1) .gt. 0.0_dp .and. diff(1).lt.1.0_dp)then
        if(diff(2) .gt. 0.0_dp .and. diff(2).lt.1.0_dp)then
          if(diff(3) .gt. 0.0_dp .and. diff(3).lt.1.0_dp)then
            l = l+1
          endif
        endif
      endif
    enddo

    nats = l

    call lcc_reallocate_realMat(r_tmp,3,Ntop)
    r_tmp = r_inout

    call lcc_reallocate_realMat(r_inout,3,nats)

    if(allocated(resindex)) deallocate(resindex)
    allocate(resindex(nats))

    l=0
    do i=1,Ntop  !Cut
      accept_index=0.0
      diff(1)=invM(1,1)*r_tmp(1,i) + invM(1,2)*r_tmp(2,i) + invM(1,3)*r_tmp(3,i)
      diff(2)=invM(2,1)*r_tmp(1,i) + invM(2,2)*r_tmp(2,i) + invM(2,3)*r_tmp(3,i)
      diff(3)=invM(3,1)*r_tmp(1,i) + invM(3,2)*r_tmp(2,i) + invM(3,3)*r_tmp(3,i)
      if(diff(1) .gt. 0.0_dp .and. diff(1).lt.1.0_dp)then
        if(diff(2) .gt. 0.0_dp .and. diff(2).lt.1.0_dp)then
          if(diff(3) .gt. 0.0_dp .and. diff(3).lt.1.0_dp)then
            l = l+1
            r_inout(1,l) = r_tmp(1,i)
            r_inout(2,l) = r_tmp(2,i)
            r_inout(3,l) = r_tmp(3,i)
          endif
        endif
      endif

    enddo

    resindex = 1
    cluster_lattice_vectors = 0.0_dp
    cluster_lattice_vectors(1,:) = pVect(:,1)
    cluster_lattice_vectors(2,:) = pVect(:,2)
    cluster_lattice_vectors(3,:) = pVect(:,3)

    deallocate(r_tmp)

  end subroutine lcc_build_slab

  !> Will add randomness to the system.
  !! \param r_inout System coordinates.
  !! \param lattice_vectors Lattice vectors.
  !! \param seed Random seed.
  !! rcoeff Coefficient for randomness.
  !!
  subroutine lcc_add_randomness_to_coordinates(r_inout,seed,rcoeff)
    implicit none
    real(dp),allocatable,intent(inout) :: r_inout(:,:)
    real(dp), intent(in) :: rcoeff
    real(dp) :: aux1,ran
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
      aux1 = rcoeff*(2.0_dp*ran - 1.0_dp)
      r_inout(1,l) = r_inout(1,l) + aux1 

      call random_number(ran)
      aux1 = rcoeff*(2.0_dp*ran - 1.0_dp)
      r_inout(2,l) = r_inout(2,l) + aux1 

      call random_number(ran)
      aux1 = rcoeff*(2.0_dp*ran - 1.0_dp)
      r_inout(3,l) = r_inout(3,l) + aux1 
    enddo

  end subroutine lcc_add_randomness_to_coordinates

end module lcc_build_mod
