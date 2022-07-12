!> Module for auxiliary operations routines
!!
module lcc_aux_mod

  use bml
  use lcc_constants_mod
  use prg_system_mod
  use lcc_allocation_mod
  use lcc_message_mod

  implicit none

  public :: lcc_get_coordination, lcc_vectors_to_parameters
  public :: lcc_parameters_to_vectors, lcc_center_at_box
  public :: lcc_canonical_basis, lcc_center_at_origin
  public :: lcc_get_reticular_density

contains

  !> Transforms the lattice vectors into lattice parameters.
  !! \param lattice_vector 3x3 array containing the lattice vectors.
  !! lattice_vector(1,:) = \f$ \overrightarrow{a} \f$
  !! \param abc_angles 2x3 array containing the lattice parameters.
  !! abc_angles(1,1) = a, abc_angles(1,2) = b and abc_angles(1,3) = c
  !! abc_angles(2,1) = \f$ \alpha \f$, abc_angles(2,2) = \f$ \beta \f$, and abc_angles(2,3) = \f$ \gamma \f$.
  !! \param verbose Verbosity level.
  !!
  subroutine lcc_vectors_to_parameters(lattice_vector,abc_angles,verbose)
    implicit none
    real(dp)               ::  angle_alpha, angle_beta, angle_gamma, lattice_a
    real(dp)               ::  lattice_b, lattice_c, pi
    real(dp), intent(in)   ::  lattice_vector(3,3)
    real(dp), intent(out)  ::  abc_angles(2,3)
    integer, intent(in) :: verbose

    pi = 3.14159265358979323846264338327950_dp

    lattice_a = sqrt(lattice_vector(1,1)**2+lattice_vector(1,2)**2+lattice_vector(1,3)**2)
    lattice_b = sqrt(lattice_vector(2,1)**2+lattice_vector(2,2)**2+lattice_vector(2,3)**2)
    lattice_c = sqrt(lattice_vector(3,1)**2+lattice_vector(3,2)**2+lattice_vector(3,3)**2)

    angle_gamma = dot_product(lattice_vector(1,:),lattice_vector(2,:))/(lattice_a*lattice_b)
    angle_beta = dot_product(lattice_vector(1,:),lattice_vector(3,:))/(lattice_a*lattice_c)
    angle_alpha = dot_product(lattice_vector(2,:),lattice_vector(3,:))/(lattice_b*lattice_c)

    angle_alpha = 360.0_dp*acos(angle_alpha)/(2.0_dp*pi)
    angle_beta = 360.0_dp*acos(angle_beta)/(2.0_dp*pi)
    angle_gamma = 360.0_dp*acos(angle_gamma)/(2.0_dp*pi)

    abc_angles(1,1) = lattice_a
    abc_angles(1,2) = lattice_b
    abc_angles(1,3) = lattice_c

    abc_angles(2,1) = angle_alpha
    abc_angles(2,2) = angle_beta
    abc_angles(2,3) = angle_gamma

  end subroutine lcc_vectors_to_parameters

  !> Transforms the lattice parameters into lattice vectors.
  !! \param abc_angles 2x3 array containing the lattice parameters.
  !! abc_angles(1,1) = a, abc_angles(1,2) = b, and abc_angles(1,3) = c
  !! abc_angles(2,1) = \f$ \alpha \f$ , abc_angles(2,2) = \f$ \beta \f$ and abc_angles(2,3) = \f$ \gamma \f$
  !! \param lattice_vector 3x3 array containing the lattice vectors.
  !! lattice_vector(1,:) = \f$ \overrightarrow{a} \f$
  !! \param verbose Verbosity level.
  !!
  subroutine lcc_parameters_to_vectors(abc_angles,lattice_vector,verbose)
    implicit none
    real(dp)               ::  angle_alpha, angle_beta, angle_gamma, lattice_a
    real(dp)               ::  lattice_b, lattice_c, pi
    real(dp), intent(in)   ::  abc_angles(2,3)
    real(dp), intent(out)  ::  lattice_vector(3,3)
    integer, intent(in)    ::  verbose

    pi = 3.14159265358979323846264338327950_dp

    lattice_a = abc_angles(1,1)
    lattice_b = abc_angles(1,2)
    lattice_c = abc_angles(1,3)
    angle_alpha = abc_angles(2,1)
    angle_beta = abc_angles(2,2)
    angle_gamma = abc_angles(2,3)


    angle_alpha = 2.0_dp*pi*angle_alpha/360.0_dp
    angle_beta = 2.0_dp*pi*angle_beta/360.0_dp
    angle_gamma = 2.0_dp*pi*angle_gamma/360.0_dp

    lattice_vector(1,1) = lattice_a
    lattice_vector(1,2)=0
    lattice_vector(1,3)=0

    lattice_vector(2,1)=lattice_b*cos(angle_gamma)
    lattice_vector(2,2)=lattice_b*sin(angle_gamma)
    lattice_vector(2,3)=0

    lattice_vector(3,1)=lattice_c*cos(angle_beta)
    lattice_vector(3,2)=lattice_c*( cos(angle_alpha)-cos(angle_gamma)* &
         cos(angle_beta) )/sin(angle_gamma)
    lattice_vector(3,3)=sqrt(lattice_c**2 - lattice_vector(3,1)**2 - lattice_vector(3,2)**2)

  end subroutine lcc_parameters_to_vectors

  !> Get the coordination of an atom.
  !! \brief Will count how many atoms are around a particular
  !! atom (coordination number) given a set radius.
  !! \param r_at Coodinates of the atom for which we need the coordination.
  !! \param r_env Coordinated of the environment sorounding atom at r_at.
  !! \param thres Threshod distance to find coordinations.
  !! \param cnum Coordination number (output).
  !!
  subroutine lcc_get_coordination(r_at,r_env,thresh,cnum)
    implicit none
    real(dp), intent(in) :: r_at(3)
    real(dp), allocatable, intent(in) :: r_env(:,:)
    real(dp), intent(in) :: thresh
    integer, intent(inout) :: cnum
    integer :: nenv,i
    real(dp) :: distance

    cnum = 0
    nenv = size(r_env,dim=2)

    do i=1,nenv
      distance = norm2(r_env(:,i) - r_at)
      if(distance <= thresh)then
        nenv = nenv + 1
      endif
    enddo

  end subroutine lcc_get_coordination

  !> To "canonical base" transformation.
  !! \brief This will reorient the shape/slab so that the
  !! first translation vector is alligned with x.
  !! \param lattice_vectors Translation vectors for the shape/slab.
  !! \param r_inout Coordinates to be transformed.
  !! \param verbose Verbosity level.
  !!
  subroutine lcc_canonical_basis(lattice_vectors,r_inout,verbose)
    implicit none
    real(dp), allocatable, intent(inout) :: lattice_vectors(:,:)
    real(dp), allocatable, intent(inout) :: r_inout(:,:)
    real(dp) :: lattice_vectors_ref(3,3),lattice_vectors_ref_norm(3,3)
    real(dp) :: lattice_vectors_inv(3,3),lattice_vectors_norm(3,3)
    real(dp) :: matrixTransf(3,3)
    real(dp) :: abc_angles(2,3),rx,ry,rz
    integer :: i,nats
    integer, intent(in) :: verbose
    type(system_type) :: sy

    !To parameters
    call lcc_vectors_to_parameters(lattice_vectors,abc_angles,verbose)
    allocate(sy%lattice_vector(3,3))
    allocate(sy%coordinate(3,4))
    sy%lattice_vector(1,1:3) = lattice_vectors(1,1:3)
    sy%lattice_vector(2,1:3) = lattice_vectors(2,1:3)
    sy%lattice_vector(3,1:3) = lattice_vectors(3,1:3)
    sy%coordinate(1:3,1) = lattice_vectors(1,1:3)
    sy%coordinate(1:3,2) = lattice_vectors(2,1:3)
    sy%coordinate(1:3,3) = lattice_vectors(3,1:3)
    sy%coordinate(1:3,4) = 0.0
    allocate(sy%symbol(4))
    sy%symbol = "I"
    sy%symbol(4) = "F"
    sy%nats = 4

    !Back to vector (to have an aligned basis set)
    call lcc_parameters_to_vectors(abc_angles,lattice_vectors_ref,verbose)

    !Normalize vectors
    lattice_vectors_ref_norm(1,:) = lattice_vectors_ref(1,:)/norm2(lattice_vectors_ref(1,:))
    lattice_vectors_ref_norm(2,:) = lattice_vectors_ref(2,:)/norm2(lattice_vectors_ref(2,:))
    lattice_vectors_ref_norm(3,:) = lattice_vectors_ref(3,:)/norm2(lattice_vectors_ref(3,:))
    lattice_vectors_norm(1,:) = lattice_vectors(1,:)/norm2(lattice_vectors(1,:))
    lattice_vectors_norm(2,:) = lattice_vectors(2,:)/norm2(lattice_vectors(2,:))
    lattice_vectors_norm(3,:) = lattice_vectors(3,:)/norm2(lattice_vectors(3,:))

    !V_ref^T = M * V_orig^T => M = V_ref^T *V_orig^T^-1
    lattice_vectors_inv = inv(transpose(lattice_vectors_norm))

    matrixTransf = matmul(transpose(lattice_vectors_ref_norm),lattice_vectors_inv)

    nats = size(r_inout,dim=2)
    do i=1,nats
      rx = dot_product(matrixTransf(1,:),r_inout(:,i))
      ry = dot_product(matrixTransf(2,:),r_inout(:,i))
      rz = dot_product(matrixTransf(3,:),r_inout(:,i))
      r_inout(1,i) = rx; r_inout(2,i) = ry; r_inout(3,i) = rz
    enddo

    do i=1,sy%nats
      rx = dot_product(matrixTransf(1,:),sy%coordinate(:,i))
      ry = dot_product(matrixTransf(2,:),sy%coordinate(:,i))
      rz = dot_product(matrixTransf(3,:),sy%coordinate(:,i))
      sy%coordinate(1,i) = rx;
      sy%coordinate(2,i) = ry;
      sy%coordinate(3,i) = rz
    enddo

    sy%lattice_vector = lattice_vectors_ref

    lattice_vectors = lattice_vectors_ref

  end subroutine lcc_canonical_basis

  !> Cetering the system inside the lattice box.
  !! \brief This will move the coordinates so that the geometric
  !! center of the system is at the center of the box.
  !! \param lattice_vectors Translation vectors for the shape/slab.
  !! \param r_inout Coordinates to be transform.
  !! \param verbose Verbosity level.
  !!
  subroutine lcc_center_at_box(lattice_vectors,r_inout,verbose)
    implicit none
    real(dp), allocatable, intent(in) :: lattice_vectors(:,:)
    real(dp), allocatable, intent(inout) :: r_inout(:,:)
    integer :: nats, i
    integer, intent(in) :: verbose
    real(dp) :: geomCent(3), boxGeomCent(3)

    call lcc_print_message("Centering at box ...",verbose)

    nats = size(r_inout,dim=2)

    geomCent = 0.0_dp
    do i = 1,nats
      geomCent = geomCent + r_inout(:,i)
    enddo
    geomCent = geomCent/real(nats,dp)

    boxGeomCent = (lattice_vectors(1,:) + lattice_vectors(2,:) + lattice_vectors(3,:))/2.0_dp

    do i = 1,nats
      r_inout(:,i) = r_inout(:,i) - geomCent +  boxGeomCent
    enddo

    geomCent = 0.0_dp
    do i = 1,nats
      geomCent = geomCent + r_inout(:,i)
    enddo
    geomCent = geomCent/real(nats,dp)

  end subroutine lcc_center_at_box

  !> Cetering the system at the origin.
  !! \brief This will move the coordinates so that the geometric
  !! center of the system is at (0,0,0).
  !! \param r_inout Coordinates to be transform.
  !! \param verbose Verbosity level.
  !!
  subroutine lcc_center_at_origin(r_inout,verbose)
    implicit none
    real(dp), allocatable :: r_inout(:,:)
    integer :: nats,i
    real(dp) :: geomCent(3)
    integer, intent(in) :: verbose

    nats = size(r_inout,dim=2)

    geomCent = 0.0_dp
    do i = 1,nats
      geomCent = geomCent + r_inout(:,i)
    enddo
    geomCent = geomCent/real(nats,dp)

    do i = 1,nats
      r_inout(:,i) = r_inout(:,i) - geomCent
    enddo

  end subroutine lcc_center_at_origin

  !> Computes the inverse of a matrix using an LU decomposition.
  !! \param A nxn Matrix to be inverted.
  !! \param Ainv Inverse of matrix A 
  !!
  function inv(A) result(invOfA)
    external dgetrf 
    external dgetri 
    real(dp), intent(in) :: A(:,:)
    real(dp), allocatable :: invOfA(:,:)
    real(dp), allocatable :: work(:)
    integer, allocatable :: ipiv(:)
    integer :: n, info

    n = size(A,dim=1)
    allocate(work(n))
    allocate(ipiv(n))
    invOfA = A

    call dgetrf(n, n, invOfA, n, ipiv, info)
    call dgetri(n, invOfA, n, ipiv, work, n, info)

    deallocate(work)
    deallocate(ipiv)

  end function inv

  !> Get the reticular density of a particular hkl face:
  !! This soubroutine computes:
  !! \param lattice_vectors Lattice vectors for the system.
  !! \param hkl_in Vector containing h, k, and l.
  !! \param density Reticular density.
  !! 
  subroutine lcc_get_reticular_density(lattice_vectors,hkl_in,density)
    implicit none
    real(dp), intent(in)                 ::  lattice_vectors(:,:)
    real(dp), intent(in)                 ::  hkl_in(:)
    real(dp), allocatable                ::  hkl(:)
    real(dp), intent(out)                ::  density
    real(dp) :: tol, area
    real(dp) :: v1(3)

    allocate(hkl(3)) ; hkl = hkl_in

    tol = 0.0001_dp
    !  (*00) case
    if((abs(hkl(2)) < tol) .and. (abs(hkl(3)) < tol))then
      hkl = hkl/norm2(hkl)
      v1 = crossProd(lattice_vectors(2,:),lattice_vectors(3,:))
      area = norm2(v1)
      !  (0*0) case
    elseif((abs(hkl(1)) < tol) .and. (abs(hkl(3)) < tol))then
      hkl = hkl/norm2(hkl)
      v1 = crossProd(lattice_vectors(1,:),lattice_vectors(3,:))
      write(*,*)"LATTICE",lattice_vectors(1,:)
      write(*,*)"LATTICE",lattice_vectors(3,:)
      area = norm2(v1)
      write(*,*)"AREA",area
      !  (00*) case
    elseif((abs(hkl(1)) < tol) .and. (abs(hkl(2)) < tol))then
      hkl = hkl/norm2(hkl)
      v1 = crossProd(lattice_vectors(1,:),lattice_vectors(2,:))
      area = norm2(v1)
      !  (**0) case
    elseif((abs(hkl(1)) > tol) .and. (abs(hkl(2)) > tol) .and. (abs(hkl(3)) < tol))then
      if(abs(hkl(1)) <= abs(hkl(2)))then
        v1 = (-hkl(2))*lattice_vectors(1,:)
        v1 = v1 + (1.0_dp)*lattice_vectors(2,:)
        v1 = crossProd(v1,lattice_vectors(3,:))
        area = norm2(v1)
      else
        v1 = (-hkl(1))*lattice_vectors(2,:)
        v1 = v1 + (1.0_dp)*lattice_vectors(1,:)
        v1 = crossProd(v1,lattice_vectors(3,:))
        area = norm2(v1)
      endif
    elseif((abs(hkl(1)) > tol) .and. (abs(hkl(3)) > tol) .and. (abs(hkl(2)) < tol))then
      if(abs(hkl(1)) <= abs(hkl(3)))then
        v1 = (-hkl(3))*lattice_vectors(1,:)
        v1 = v1 + (1.0_dp)*lattice_vectors(3,:)
        v1 = crossProd(v1,lattice_vectors(2,:))
        area = norm2(v1)
      else
        v1 = (-hkl(1))*lattice_vectors(3,:)
        v1 = v1 + (1.0_dp)*lattice_vectors(1,:)
        v1 = crossProd(v1,lattice_vectors(2,:))
        area = norm2(v1)
      endif
    elseif((abs(hkl(2)) > tol) .and. (abs(hkl(3)) > tol) .and. (abs(hkl(1)) < tol))then
      if(abs(hkl(2)) <= abs(hkl(3)))then
        v1 = (-hkl(3))*lattice_vectors(2,:)
        v1 = v1 + (1.0_dp)*lattice_vectors(3,:)
        v1 = crossProd(v1,lattice_vectors(1,:))
        area = norm2(v1)
      else
        v1 = (-hkl(2))*lattice_vectors(3,:)
        v1 = v1 + (1.0_dp)*lattice_vectors(2,:)
        v1 = crossProd(v1,lattice_vectors(1,:))
        area = norm2(v1)
      endif

    endif

    density = 1.0_dp/area

  end subroutine lcc_get_reticular_density

  function crossProd(r1,r2) result(r3)
    implicit none
    real(dp), intent(in)  :: r1(:), r2(:)
    real(dp), allocatable  :: r3(:)

    !| i       j     k   |
    !| r1(1) r1(2) r1(3) |
    !| r2(1) r2(2) r2(3) |

    if(.not. allocated(r3)) allocate(r3(3))
    r3(1) =r1(2)*r2(3) - r2(2)*r1(3)
    r3(2) = -(r1(1)*r2(3) - r1(3)*r2(1))
    r3(3) = r1(1)*r2(2) - r1(2)*r2(1)

  end function crossProd

end module lcc_aux_mod
