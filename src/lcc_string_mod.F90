!> Module for manipulating strings
!!
module lcc_string_mod

  use lcc_constants_mod

  implicit none

  public :: lcc_get_word, lcc_split_string

contains

  !> Cut a word from string.
  !! \param string Full string.
  !! \param posh Cut from position.
  !! \param post Cut to position.
  !! \param word Extracted word.
  !!
  subroutine lcc_get_word(string,posh,post,word)
    implicit none
    character(len=*), intent(in)      ::  string
    character(20), intent(inout)      ::  word
    character(1), allocatable         ::  tempc(:)
    character(len=30)                 ::  tempcflex
    integer, intent(in)               ::  posh, post
    integer                           ::  lenc

    lenc=len(adjustl(trim(string)))
    if(.not.allocated(tempc))allocate(tempc(lenc))
    tempcflex = adjustl(trim(string))
    word = adjustl(trim(tempcflex(posh:post)))

  end subroutine lcc_get_word

  !> Split a string in two words uning a delimiter.
  !! \param string Full string.
  !! \param delimit Delimiter.
  !! \param head First word.
  !! \param tail Last word.
  !!
  subroutine lcc_split_string(string,delimit,head,tail)
    implicit none
    character(len=*), intent(in)      ::  string
    character(20), intent(inout)      ::  head, tail
    character(1), intent(in)          ::  delimit
    character(1), allocatable         ::  tempc(:)
    character(len=30)                 ::  tempcflex
    character(1)                      ::  myChar
    integer                           ::  lenc,i,pos

    lenc=len(adjustl(trim(string)))
    if(.not.allocated(tempc))allocate(tempc(lenc))
    tempcflex = adjustl(trim(string))
    do i=1,lenc
      myChar = adjustl(trim(tempcflex(i:i)))
      if(myChar == delimit)then 
        pos = i
        exit
      endif
    enddo

    head = adjustl(trim(tempcflex(1:pos-1)))
    tail = adjustl(trim(tempcflex(pos+1:lenc)))

  end subroutine lcc_split_string


end module lcc_string_mod
