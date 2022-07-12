!> Module for allocation operations. 
!!
module lcc_allocation_mod

  use bml
  use lcc_constants_mod

  implicit none

  public :: lcc_reallocate_realVect, lcc_reallocate_realMat
  public :: lcc_reallocate_intVect, lcc_reallocate_intMat
  public :: lcc_reallocate_char2Vect, lcc_reallocate_char3Vect

contains

  !> To reallocate a real vector.
  !! \brief This will reallocate a vector
  !! If it is already allocated, a deallocation
  !! will first happen.
  !! \param vect Real 1D array.
  !! \param ndim Dimension to reallocate the vector to.
  !!
  subroutine lcc_reallocate_realVect(vect,ndim)
    integer, intent(in) :: ndim
    real(dp), allocatable, intent(inout) :: vect(:)

    if(allocated(vect))then 
     deallocate(vect) 
    endif
    allocate(vect(ndim))

  end subroutine lcc_reallocate_realVect

  !> To reallocate a real mxn matrix.
  !! \brief This will reallocate a matrix
  !! If it is already allocated, a deallocation
  !! will first happen.
  !! \param mat Real 2D array.
  !! \param mnim First dimension to realocate the matrix to.
  !! \param ndim Second dimension to reallocate the matrix to.
  !!
  subroutine lcc_reallocate_realMat(mat,mdim,ndim)
    integer, intent(in) :: ndim, mdim
    real(dp), allocatable, intent(inout) :: mat(:,:)

    if(allocated(mat))then
     deallocate(mat)
    endif
    allocate(mat(mdim,ndim))
  
  end subroutine lcc_reallocate_realMat

  !> To reallocate a real vector.
  !! \brief This will reallocate a vector
  !! If it is already allocated, a deallocation
  !! will first happen.
  !! \param vect Integer 1D array.
  !! \param ndim Dimension to reallocate the vector to.
  !!
  subroutine lcc_reallocate_intVect(vect,ndim)
    integer, intent(in) :: ndim
    integer, allocatable, intent(inout) :: vect(:)

    if(allocated(vect))then
     deallocate(vect)
    endif
    allocate(vect(ndim))

  end subroutine lcc_reallocate_intVect

  !> To reallocate an integer mxn matrix.
  !! \brief This will reallocate a matrix.
  !! If it is already allocated, a deallocation
  !! will first happen.
  !! \param mat Integer 2D array.
  !! \param mnim First dimension to realocate the matrix to.
  !! \param ndim Second dimension to reallocate the matrix to.
  !!
  subroutine lcc_reallocate_intMat(mat,mdim,ndim)
    integer, intent(in) :: ndim, mdim
    integer, allocatable, intent(inout) :: mat(:,:)

    if(allocated(mat))then
     deallocate(mat)
    endif
    allocate(mat(mdim,ndim))

  end subroutine lcc_reallocate_intMat

  !> To reallocate a character vector.
  !! \brief This will reallocate a character len=2 vector
  !! If it is already allocated, a deallocation
  !! will first happen.
  !! \param vect Character(2) 1D array.
  !! \param ndim Dimension to reallocate the vector to.
  !!
  subroutine lcc_reallocate_char2Vect(vect,ndim)
    integer, intent(in) :: ndim
    character(2), allocatable, intent(inout) :: vect(:)

    if(allocated(vect))then
     deallocate(vect)
    endif
    allocate(vect(ndim))

  end subroutine lcc_reallocate_char2Vect

  !> To reallocate a character vector
  !! \brief This will reallocate a character len=3 vector.
  !! If it is already allocated, a deallocation
  !! will first happen.
  !! \param vect Character(3) 1D array.
  !! \param ndim Dimension to reallocate the vector to.
  !!
  subroutine lcc_reallocate_char3Vect(vect,ndim)
    integer, intent(in) :: ndim
    character(3), allocatable, intent(inout) :: vect(:)

    if(allocated(vect))then
     deallocate(vect)
    endif
    allocate(vect(ndim))

  end subroutine lcc_reallocate_char3Vect

end module lcc_allocation_mod
