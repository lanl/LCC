!> A module to handle the constants needed by the code.
!! \brief This module will be used to store the constants needed
!! in the code
!!
!! @ingroup LCC
!!
module lcc_constants_mod

  implicit none

  !> Precision used troughout the code
  integer, parameter, public :: dp = kind(1.0d0)

  !> Pi number
  real(dp),parameter :: pi = 3.14159265358979323846264338327950_dp

end module lcc_constants_mod

