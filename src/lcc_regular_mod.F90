!> Module for generating regular shapes after lattice is constructed.
!!
module lcc_regular_mod

  use lcc_constants_mod
  use lcc_lattice_mod
  use lcc_allocation_mod
  use prg_system_mod

  implicit none

  public :: lcc_spheroid

contains

  !> For building spheroidal shapes out of a bulk 
  !! lattice.
  !! \param a_axis Lenght in the x direction.
  !! \param b_axis Lenght in the y direction.
  !! \param c_axis Lenght in the z direction.
  !! \param r_inout Input and output coordinates.
  !!
  subroutine lcc_spheroid(a_axis,b_axis,c_axis,r_inout)
    implicit none
    real(dp), allocatable :: r_inout(:,:)
    real(dp) :: a_axis,b_axis,c_axis
    real(dp) :: x_param,y_param,z_param,r_param
    integer :: i,nTmp,nats
    real(dp), allocatable :: r_tmp(:,:)

    write(*,*)'In spheroid ...'

    nTmp = size(r_inout,dim=2)
    call lcc_reallocate_realMat(r_tmp,3,nTmp)
    nats = 0
    do i=1,nTmp
      x_param = r_inout(1,i)/a_axis
      y_param = r_inout(2,i)/b_axis
      z_param = r_inout(3,i)/c_axis
      r_param = sqrt(x_param**2 + y_param**2 + z_param**2)
      if(r_param <= 1)then 
        nats = nats + 1
        r_tmp(:,nats) = r_inout(:,i)
      endif 
    enddo 

    call lcc_reallocate_realMat(r_inout,3,nats)
    r_inout(:,:) = r_tmp(:,1:nats)
    deallocate(r_tmp)

  end subroutine lcc_spheroid

end module lcc_regular_mod 
