!> Program to move and wrapp the strcture
!!
!! \ingroup PROGRAMS
!!
!! Example using this program:
!!
!!     \verbatim  moveandwrapp coordinate.xyz x 1.0 \endverbatim
!!
!!
program moveandwrapp

  !PROGRESS lib modes.
  use prg_system_mod
  integer, parameter                :: dp = kind(1.0d0)
  integer                           :: i
  type(system_type)                 :: sy
  character(30)                     :: filein, namein
  character(3)                      :: extin
  character(2)                      :: axis
  character(4)                      :: lenghtC 
  real(dp)                          :: lenght, projb1, projb2, projb3
  real(dp)                          :: matM(3,3),invM(3,3), DET 
  character(1), allocatable         :: tempc(:)
  character(len=30)                 :: tempcflex
  integer                           :: lenc

  call getarg(1, filein)
  call getarg(2, axis)
  call getarg(3, lenghtC)

  if(filein == "")then
    write(*,*)""
    write(*,*)"Usage:"
    write(*,*)""
    write(*,*)"  $ moveandwrapp <filein> <axis> <lenght> "
    write(*,*)""
    write(*,*)"<filein>:  Input coordinates file "
    write(*,*)"<axis>: x, y, or z (slab axis)"
    write(*,*)"<lenght>: Length of the move"
    write(*,*)""
    stop
  endif

  read(lenghtC,*) lenght

  lenc=len(adjustl(trim(filein)))
  if(.not.allocated(tempc))allocate(tempc(lenc))
  tempcflex = adjustl(trim(filein))
  namein = adjustl(trim(tempcflex(1:lenc-4)))
  extin = adjustl(trim(tempcflex(lenc-2:lenc+1)))

  call prg_parse_system(sy,adjustl(trim(namein)),extin)

  !Moving
  do i = 1, sy%nats
    if(trim(adjustl(axis)) == "x")then 
      sy%coordinate(:,i) = sy%coordinate(:,i) + lenght*sy%lattice_vector(1,:)/norm2(sy%lattice_vector(1,:)) 
    elseif(trim(adjustl(axis)) == "y")then
      sy%coordinate(:,i) = sy%coordinate(:,i) + lenght*sy%lattice_vector(2,:)/norm2(sy%lattice_vector(2,:)) 
    elseif(trim(adjustl(axis)) == "z")then
      sy%coordinate(:,i) = sy%coordinate(:,i) + lenght*sy%lattice_vector(3,:)/norm2(sy%lattice_vector(3,:)) 
    else
      STOP "Axis should be either x,y, or z."
    endif
  enddo

    do j=1,3
      matM(1,j) = sy%lattice_vector(j,1)
      matM(2,j) = sy%lattice_vector(j,2)
      matM(3,j) = sy%lattice_vector(j,3)
    enddo

    DET =     matM(1,1)*matM(2,2)*matM(3,3)  &
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



  !Wrapping
  do i = 1, sy%nats 

    !To b1,b2,b3 basis

    projb1 = invM(1,1)*sy%coordinate(1,i) + invM(1,2)*sy%coordinate(2,i) + invM(1,3)*sy%coordinate(3,i)
    projb2 = invM(2,1)*sy%coordinate(1,i) + invM(2,2)*sy%coordinate(2,i) + invM(2,3)*sy%coordinate(3,i)
    projb3 = invM(3,1)*sy%coordinate(1,i) + invM(3,2)*sy%coordinate(2,i) + invM(3,3)*sy%coordinate(3,i)

    !Wrapping
    if(projb1 > 1.0_dp)then
      projb1 = projb1 - 1.0_dp
    elseif(projb1 < 0.0_dp)then
      projb1 = projb1 + 1.0_dp
    endif

    if(projb2 > 1.0_dp)then
      projb2 = projb2 - 1.0_dp
    elseif(projb2 < 0.0_dp)then
      projb2 = projb2 + 1.0_dp
    endif

    if(projb3 > 1.0_dp)then
      projb3 = projb3 - 1.0_dp
    elseif(projb3 < 0.0_dp)then
      projb3 = projb3 + 1.0_dp
    endif

    !Back to canonical
    sy%coordinate(:,i) = projb1*sy%lattice_vector(1,:) + projb2*sy%lattice_vector(2,:) + &
     & projb3*sy%lattice_vector(3,:)

  enddo

  call prg_write_system(sy,adjustl(trim(namein))//"_moved",extin)

end program moveandwrapp

