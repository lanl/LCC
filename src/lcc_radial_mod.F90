module lcc_radial_mod

  use lcc_constants_mod
  use prg_dos_mod
  use lcc_aux_mod
  use prg_openfiles_mod
  use lcc_string_mod

  implicit none

  private

  public :: lcc_distance_matrix

contains

  !> Radial distribution function 
  !! \brief computed the Rdf of a pair of atomic species.
  !! \param atom_symbols Atomic symbols from the system 
  !! \param r_in coordinates of the system
  !! \param lattice_vectors Lattice vectors to determin the points density
  !! \param atI Symbol I 
  !! \param atJ Symbol J
  !!
  subroutine lcc_distance_matrix(atom_symbols,r_in,lattice_vectors,atIJ,verbose)
    implicit none
    character(2), allocatable, intent(in) :: atom_symbols(:)
    real(dp), allocatable, intent(in) :: r_in(:,:)
    character(5), intent(in) :: atIJ
    character(20) :: atI, atJ
    character(20) :: word
    integer :: i,j,nats,cont,npts,io
    real(dp), allocatable :: dij(:),loads(:),lattice_vectors(:,:),recip_vectors(:,:)
    character(20) :: filename
    real(dp) :: dr,dmin,dmax,volr,volk,density
    integer, optional, intent(in) :: verbose

    if(present(verbose)) call lcc_print_message("Computing Rdf ...",verbose)

    call lcc_split_string(atIJ,"-",atI,atJ)

    nats = size(atom_symbols,dim=1)
    allocate(dij(nats*nats))
    allocate(loads(size(dij)))

    !Compute density
    call prg_get_recip_vects(lattice_vectors,recip_vectors,volr,volk)
    density = real(nats,dp)/volr
    cont = 0
    loads = 1.0_dp
    dij = 1.0D10
    do i=1,nats
      if(atom_symbols(i) == trim(adjustl(atI)))then
        do j=1,nats
          if(atom_symbols(j) == trim(adjustl(atJ)) .and. (i .ne. j) )then
            cont = cont + 1
            dij(cont) = norm2(r_in(:,i) - r_in(:,j))
            loads(cont) = 1.0_dp/(nats*density*4*pi*dij(cont)**2)
          endif
        enddo
      endif
    enddo

    filename = trim(adjustl(atI))//"-"//trim(adjustl(atJ))//".rdf"

    call prg_open_file(io,filename)

    dmax = 20.0_dp
    dmin = 0.2_dp

    npts = 1000
    dr = (dmax-dmin)/real(npts)

    write(io,*)"#  r    d"
    do i = 1, npts
      write(io,*) dmin + dr*i,lorentz(dmin + dr*i, dij, loads, 0.2_dp)
    end do

    close(io)

  end subroutine lcc_distance_matrix


  !> Lorentzian Function
  !! \brief Computes:
  !! \f$ L(\epsilon) = \sum_{k} \frac{\omega(k)\Gamma}{2 \pi}\frac{1}{(\epsilon - \epsilon_k)^2 + (\Gamma/2)^2} \f$
  !! \param energy Energy point.
  !! \param eigenvals Eigenvalues of the system.
  !! \param Gamma Lorentz function broadening.
  !!
  real(dp) function lorentz(energy, eigenvals, loads, Gamma)
    implicit none
    integer               ::  Nstates, k
    real(dp)              ::  auxfactor, auxterm, pi
    real(dp), intent(in)  ::  Gamma, eigenvals(:), energy, loads(:)

    Nstates = size(eigenvals,dim=1)
    pi = 3.14159265358979323846264338327950_dp

    !Lorentz parameters
    auxfactor = Gamma/(2.0_dp*pi)
    auxterm = (Gamma/2.0_dp)**2
    lorentz = 0.0_dp

    do k = 1, Nstates
      lorentz = lorentz + loads(k)/((energy-eigenvals(k))**2 + auxterm)
    end do

    lorentz = auxfactor*lorentz

  end function lorentz

end module lcc_radial_mod
