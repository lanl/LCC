!> Template module for contributing
!!
module lcc_template_mod

  use bml
  use lcc_constants_mod
  use prg_system_mod
  use lcc_allocation_mod
  use lcc_message_mod

  implicit none

  public :: lcc_template_subroutine

contains

  !> Example subroutine.
  !! \param dummy Example variable
  !! \param verbose Verbosity level.
  !!
  subroutine lcc_template_subroutine(dummy, verbose)
    implicit none
    character(*), intent(in)   ::  dummy
    integer, intent(in) :: verbose

    call lcc_print_message("Centering at box ...",verbose)
    stop

  end subroutine lcc_template_subroutine

end module lcc_template_mod
