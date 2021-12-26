module pde

  ! This module contains the kinetic equation PDE, aka the Boltzmann equation,
  ! or whatever other equation (Vlasov / BGK / ...) according to the collision
  ! operator that we employ.
  ! Hyper2D is still two-dimensional, meaning that the kinetic equation will 
  ! be 1D1V: 1D in space and 1V in the velocity. 
  ! Dimension "x" refers to the space coordinate, while dimension "y" is a 
  ! velocity "v".
   
  use global_module

  implicit none

  real(kind=8), parameter :: M = 9.109d-31     ! [kg] particle mass
  real(kind=8), parameter :: q = -1.602176d-19 ! [C] particle charge

  real(kind=8), parameter :: kB  = 1.3806503d-23 ! [J/K] Boltzmann constant
  real(kind=8), parameter :: PI  = 3.1415926535 
 
  integer, parameter :: Neq = 1 ! Number of equations: (x, v, f)

  ! Name of primitive variables, used ONLY for exporting the solution to VTK file
  ! NOTE: To initialize it here, all entries need to have the same length!!!
  !       I'm using three letters for simplicity.
  character(len=20), dimension(Neq) :: prim_names = (/'f'/)

  contains

  ! ============================================================

  subroutine initialize_solution(U)
 
    implicit none

    real(kind=8), dimension(:,:,:), intent(inout) :: U

    integer      :: i, j
    real(kind=8) :: x, p, gam, v

    ! Initialize all cells with a single-velocity (1V) Maxwellian VDF
    ! do j = 1, Ny
    !   v = compute_v_from_j(j)
    !   do i = 1, Nx
    !     ! x = x_min + dx/2.0 + real(i-1)*dx
    !     U(1,i,j) = n0*sqrt(M/(2.0*PI*kB*T0))*exp(-M/(2.0*kB*T0)*(v - u0)*(v - u0))
    !   end do
    ! end do

    ! Init with a box
    U(1,:,:) = 1.0d-10
    do j = 1, Ny ! floor(Ny/2.0 - Ny/20.0), floor(Ny/2.0 + Ny/20.0)

      p   = compute_p_from_j(j)
      gam = sqrt(1 + p*p/M/c)
      v   = p/gam/M

      do i = 1, Nx ! floor(Nx/2.0 - Nx/20.0), floor(Nx/2.0 + Nx/20.0)

        x = compute_x_from_i(i)

        if ((x .ge. 0.0) .and. (x.lt.c/10.0) .and. (v .ge. 0.0) .and. (v .le. c/2.5) ) then
          U(1,i,j) = 1.0
        end if
      end do
    end do

  end subroutine 

  ! ============================================================

  function compute_x_from_i(i)

    implicit none

    integer, intent(in) :: i
    real(kind=8)        :: compute_x_from_i

    compute_x_from_i = x_min + dx/2.0 + real(i-1)*dx

  end function 

  ! ============================================================

  function compute_p_from_j(j)

    implicit none

    integer, intent(in) :: j
    real(kind=8)        :: compute_p_from_j

    compute_p_from_j = y_min + dy/2.0 + real(j-1)*dy

  end function 

  ! ============================================================

  subroutine assign_BCs(U)
 
    implicit none

    real(kind=8), dimension(:,:,:), intent(inout) :: U

    integer :: i, j

    ! Periodic BCs 
    if (X_PERIODIC_BOOL) then
      U(:,1,:)    = U(:,Nx-3,:)
      U(:,2,:)    = U(:,Nx-2,:)
      U(:,Nx-1,:) = U(:,3,:)
      U(:,Nx,:)   = U(:,4,:)
    end if

  end subroutine 

  ! =========================================================================

  function F_field(x)

    implicit none

    real(kind=8), intent(in) :: x
    real(kind=8) :: F_field

    real(kind=8) :: k

    ! Cosinusoidal electric field profile
    F_field = 2.0d-21 ! [N]

  end function

end module
