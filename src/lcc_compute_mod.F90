!> Template module for contributing
!!
module lcc_compute_mod

  use bml
  use lcc_constants_mod
  use prg_system_mod
  use lcc_allocation_mod
  use lcc_message_mod
  use lcc_aux_mod

  implicit none

  public :: lcc_compute_roughness

contains

  !> Example subroutine.
  !! \param coords Coordinates.
  !! \param lattice_vectors Lattice vectors.
  !! \param isoval Parameter value to compute isosurface.
  !! \param rab Radius of the spherical probe.
  !! \param ni Number of discrete points on the a1 axis.
  !! \param nj Number of discrete points on the a2 axis.
  !! \param nk Number of discrete points on the a3 axis.
  !! \param verbose Verbosity level.
  !!
  subroutine lcc_compute_roughness(coords,lattice_vectors,isoval,rab,ni,nj,nk,verbose)
    implicit none
    integer, intent(in) :: verbose
    real(dp), allocatable, intent(in) :: lattice_vectors(:,:)
    real(dp), allocatable, intent(in) :: coords(:,:)
    integer, intent(in) :: ni, nj, nk
    real(dp), intent(in) :: isoval, rab
    real(dp) :: rijk(3)
    real(dp), allocatable :: a1(:), a2(:), a3(:)
    integer :: i, j, k, l, m, n, nats
    integer :: myindex
    real(dp) :: drc,volr,volk,dx,dy,dz,vol,rad
    real(dp) :: aPerp(3), S0, S1, myat(3)
    real(dp) :: gval, gvalTot
    real(8), allocatable :: gvalVect(:)
    real(dp) :: rough
    real(dp),allocatable :: x(:,:),y(:,:),z(:,:),rpar(:,:,:)
    real(dp) :: dS
    real(dp) :: ru(3), rv(3), ruXrv(3)

    nats = size(coords,dim=2)
    call lcc_print_message("Computing roughness ...",verbose)
    
    allocate(a1(3))
    allocate(a2(3))
    allocate(a3(3))
    a1 = lattice_vectors(1,:)
    a2 = lattice_vectors(2,:)
    a3 = lattice_vectors(3,:)
    
    !Discretization
    !ni = 40 ; nj = 80 ; nk = 80
    dx = 1.0_dp ; dy = 1.0_dp ; dz = 1.0_dp
    dx = dx/real(ni)
    dy = dy/real(nj)
    dz = dz/real(nk)

    myindex = 0
    !Rad probe
    !rad = 1.0_dp
    allocate(gvalVect(ni*nj*nk))

    !$omp parallel do default (none) &
    !$omp shared(ni,nj,nk,dx,dy,dz) &
    !$omp shared(a1,a2,a3,coords,nats) &
    !$omp shared(rad,gvalVect) &
    !$omp private(i,j,k,l,rijk,myindex) &
    !$omp private(drc,myat,gval,gvalTot) 
    do i=1,ni
      do j=1,nj
        do k=1,nk
          myindex = (i-1)*nj*nk + (j-1)*nk + k
          rijk = a1*i*dx + a2*j*dy + a3*k*dz + a1
          gvalTot = 0.0_dp
          do l = 1,nats
            do m = -1,1
              do n = -1,1
                myat = (coords(:,l)) + m*a2 + n*a3
                drc = norm2(rijk - myat)
                gval = 10*exp(-0.1*drc**2)/exp(0.1_dp)
                gvalTot = gvalTot + gval
              enddo
            enddo
          enddo
          gvalVect(myindex) = gvalTot
        enddo
      enddo
    enddo
    !$omp end parallel do

    allocate(rpar(ni,nj,nk))
    do i=1,ni
      do j=1,nj
        do k=1,nk
          myindex = (i-1)*nj*nk + (j-1)*nk + k
          rijk = a1*i*dx + a2*j*dy + a3*k*dz + a1
          rpar(i,j,k) = gvalVect(myindex)
        enddo
      enddo
    enddo

    !isoval = 50.0_dp
    allocate(x(nj,nk))
    allocate(y(nj,nk))
    allocate(z(nj,nk))
    do k = 1,nk
      do j = 1,nj
        do i = 1,ni
          !Get highest point after isoval
          if(rpar(i,j,k) < isoval)then
            x(j,k) = a1(1)*i*dx + a2(1)*j*dy + a3(1)*k*dz + a1(1)
            y(j,k) = a1(2)*i*dx + a2(2)*j*dy + a3(2)*k*dz + a1(2)
            z(j,k) = a1(3)*i*dx + a2(3)*j*dy + a3(3)*k*dz + a1(3)
            exit
          endif
        enddo
      enddo
    enddo

    !A mask to plot the effective area
    open(101,file="mask.xyz",status='UNKNOWN')
    write(101,*)nj*nk
    write(101,*)""
    do k = 1,nk
      do j = 1,nj
        write(101,*)"H",x(j,k),y(j,k),z(j,k)
      enddo
    enddo
    close(101)

    !Integral
    S1 = 0.0_dp
    do k = 2,nk
      do j = 2,nj
        ru(1) = (y(j,k) - y(j-1,k))/dy
        ru(2) = (z(j,k) - z(j-1,k))/dy 
        ru(3) = (x(j,k) - x(j-1,k))/dy
        
        rv(1) = (y(j,k) - y(j,k-1))/dz
        rv(2) = (z(j,k) - z(j,k-1))/dz 
        rv(3) = (x(j,k) - x(j,k-1))/dz

        !Cross prod 
        ruXrv = crossProd(ru,rv)
        !Surf element
        dS = norm2(ruXrv)*dy*dz
        S1 = S1 + dS
      enddo
    enddo

    aPerp = crossProd(a2,a3)
    S0 = norm2(aPerp)

    call lcc_print_realVal("Flat area (S0)",S0,"Ang^2",verbose)
    call lcc_print_realVal("Effective area (S1)",S1,"Ang^2",verbose)
    rough = S1/S0
    call lcc_print_realVal("Roughness parameter (S1/S0)",rough,"",verbose)

  end subroutine lcc_compute_roughness

end module lcc_compute_mod
