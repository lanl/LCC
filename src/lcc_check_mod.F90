!> Module for checking operations routines
!!
module lcc_check_mod

  use lcc_message_mod
  use lcc_constants_mod

  implicit none 

  public :: lcc_check_periodicity

  contains

  !> Check the periodicity.
  !! \brief Will use a "brute force" approach to check periodidity.
  !! \param r_in Input coordinates.
  !! \param lattice_vectors Translation vectors for the slab.
  !! \param r_ref Reference or "bulk structure from where the 
  !! shape was cut. 
  !! \param verbose Verbosity level.
  !!
  subroutine lcc_check_periodicity(r_in,lattice_vectors,r_ref,tol,verbose)
    implicit none 
    real(dp),allocatable,intent(in) :: r_in(:,:) 
    real(dp),allocatable,intent(in) :: lattice_vectors(:,:)
    real(dp),allocatable,intent(in) :: r_ref(:,:)
    integer, intent(in) :: verbose
    integer :: i,j,nats,natsRef
    real(dp), intent(in) :: tol
    real(dp) :: distance, jump(3)
    logical :: inBulk

    call lcc_print_message("Checking periodicity of the slab ...",verbose)

    nats = size(r_in,dim=2)
    natsRef = size(r_ref,dim=2)

    inBulk = .false.

    !This will check if points correspond to a bulk point 
    !by translation.
    do i = 1,nats
      !Make a jump on x y and z
      jump = r_in(:,i) + lattice_vectors(1,:) + &
        & lattice_vectors(2,:) + lattice_vectors(3,:)
      
      !Check if the jump is on the reference structure.
      do j = 1,natsRef
        distance = norm2(jump(:) - r_ref(:,j))
        if (distance < tol)then 
          inBulk = .true.
          exit
        endif
      enddo

      if(.not.inBulk)then 
          call lcc_print_error("lcc_check_periodicity","Periodicity condition &
          &is not fullfilled for the shape")
      endif
      inBulk = .false.

    enddo

    call lcc_print_message("Periodicity condition is OK ...",verbose)

  end subroutine lcc_check_periodicity

end module lcc_check_mod


