!> Module for Monte Carlo related routines 
!!
module lcc_mc_mod
  use lcc_constants_mod
  use lcc_allocation_mod
  use lcc_aux_mod

!> Maximize: This checks the acceptance
!!
    subroutine lcc_check_system(r,iter,temp,cost,cost0)
      implicit none 
      integer :: i
      integer, intent(in) :: iter
      real(dp) :: DeltV, DeltvB, ran
      real(dp), intent(in) :: temp
      real(dp), allocatable, intent(inout) :: r(:,:)

      if(.not.allocated(r))then 
        stop "Coordinates not allocated"
      endif

      Beta =1.0_dp/(temp*kb) !Inverse temperature
      deltV=cost-cost0  !Checking acceptance
      deltvB=beta*deltV

      call RANDOM_NUMBER(ran)

      if(deltvB.LT.75.0_dp)then
        if(deltv.LE.0.0_dp)then !If its favorable we accept the move
          r(m,1)=rxinew
          r(m,2)=ryinew
          r(m,3)=rzinew
          accept = accept + 1 !We count the acceptance rate
          cost0 = cost
        elseif(exp(-deltvB).GT.ran)then
          r(m,1)=rxinew
          r(m,2)=ryinew
          r(m,3)=rzinew
          accept = accept + 1
          cost0 = cost
        else
          r(m,1)=rxiold
          r(m,2)=ryiold
          r(m,3)=rziold
        endif
      else
          r(m,1)=rxiold
          r(m,2)=ryiold
          r(m,3)=rziold
      endif

      write(*,*)iter,temp,real(accept)/real(iter),cost0

    end subroutine lcc_check_system

! END OF THE KERNEL ***********************************************     

end module lcc_mc_mod
