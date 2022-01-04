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

  real(kind=8), parameter :: M = 6.6337d-26   ! [kg] particle mass
  real(kind=8), parameter :: q = 1.602176d-19 ! [C] particle charge

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
    real(kind=8) :: x, v

    ! Initialize all cells with a single-velocity (1V) Maxwellian VDF
    do j = 1, Ny

      v = compute_v_from_j(j)

      do i = 1, Nx

        ! x = x_min + dx/2.0 + real(i-1)*dx

        U(1,i,j) = n0*sqrt(M/(2.0*PI*kB*T0))*exp(-M/(2.0*kB*T0)*(v - u0)*(v - u0))

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

  function compute_v_from_j(j)

    implicit none

    integer, intent(in) :: j
    real(kind=8)        :: compute_v_from_j

    compute_v_from_j = y_min + dy/2.0 + real(j-1)*dy

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

  function E_field(x)

    implicit none

    real(kind=8), intent(in) :: x
    real(kind=8) :: E_field

    real(kind=8) :: k

    ! Cosinusoidal electric field profile
    k       = 2.0*PI*4.0/(x_max-x_min)
    E_field = 10.0*k*cos(k*x)

  end function

!!   ! ==========================================================================
!! 
!!   subroutine compute_moments(U, i, n, u_ave, T)
!! 
!!     implicit none
!! 
!!     real(kind=8), dimension(Neq,Nx,Ny), intent(in)  :: U
!!     integer,      intent(in)  :: i
!!     real(kind=8), intent(out) :: n, u_ave, T
!! 
!!     integer :: j
!! 
!!     real(kind=8) :: n_u, rhoE, v
!! 
!!     ! For every position along x, integrate in velocity
!!     n     = 0.0
!!     n_u   = 0.0
!!     u_ave = 0.0
!!     T     = 0.0
!!     rhoE  = 0.0
!! 
!!     do j = 1, Ny
!! 
!!       v    = compute_v_from_j(j)
!! 
!!       n     = n     + U(1,i,j)*dy
!!       n_u   = n_u   + v*U(1,i,j)*dy
!!       rhoE  = rhoE  + M*v*v/2.0*U(1,i,j)*dy
!! 
!!     end do
!! 
!!     u_ave = n_u/n ! Average velocity
!!     T     = (rhoE - n*M*u_ave*u_ave/2.0)*2.0/(n*kB)
!! 
!!   end subroutine

end module
