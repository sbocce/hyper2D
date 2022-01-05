module integration

  use pde
  use global_module

  implicit none
  
  real(kind=8), dimension(Neq,Nx) :: U_new ! Working variable

  contains

  ! ======================================================================== 
  
  subroutine forward_Euler_step(U, dt)

    ! This function performs one step of the Forward Euler explicit time integrator 

    implicit none
 
    real(kind=8), dimension(Neq,Nx), intent(inout)  :: U

    real(kind=8), intent(in) :: dt
    integer :: i, j, eqID
    
    real(kind=8), dimension(Neq) :: F_W, F_E

    ! First, assign BCs. Wall BCs need to be updated in time, taking the value
    ! from neighboring cells.
    call assign_BCs(U) ! See the pde.f03 module

    ! Update solution using the forward Euler integrator - INTERNAL CELLS only
    ! (I am using 4 ghost cells per side. This is useful for further generalizations
    !  to higher order methods)
    do i = 5, Nx-4

      ! Choose one of these fluxes
      call compute_fluxes_HLL(U, i, F_W, F_E)
      ! call compute_fluxes_Rusanov(U, i, F_W, F_E)

      ! Advect remaining stuff
      U_new(:,i) = U(:,i) - dt/dx*(F_E - F_W)  ! No source term

      ! Check that the solution did not diverge
      do eqID = 1, Neq
        if (isnan(U_new(eqID,i))) then 
          print*, 'Solution diverged, try with a smaller time step! Aborting.'
          print*, 'Solution that diverged: ', U_new(:,i)
          print*, 'in cell i = ', i
          stop
        end if
      end do

    end do

    ! Save solution (internal cells ONLY! Do not overwrite ghost cells!)
    U(:,5:Nx-4) = U_new(:,5:Nx-4)

  end subroutine

  ! ======================================================================== 
 
  subroutine compute_fluxes_HLL(U, i, F_W, F_E)

    ! Computes HLL numerical fluxes for the cell "i"
    ! 'F_W' and 'F_E' are the complete fluxes at the cell interfaces, to be then
    ! added in the time-discretized PDE.
    ! The variables 'F_L' and 'F_R' are temporary variables used for the 
    ! fluxes at the left and right of an interface.

    implicit none

    real(kind=8), dimension(Neq,Nx), intent(in)  :: U
    real(kind=8), dimension(Neq),    intent(out) :: F_W, F_E
    integer, intent(in) :: i

    integer :: eqID
    real(kind=8), dimension(Neq) :: U_im3, U_im2, U_im1, U_i, U_ip1, U_ip2, U_ip3
    real(kind=8), dimension(Neq) :: U_L, U_R, F_L, F_R

    real(kind=8) :: ws_min_L, ws_max_L, ws_min_R, ws_max_R, ws_min, ws_max ! wave speeds

    ! ---- West interface -----

    U_L = U(:,i-1)
    U_R = U(:,i)

    call compute_flux_ws_x(U_L, F_L, ws_max_L, ws_min_L)
    call compute_flux_ws_x(U_R, F_R, ws_max_R, ws_min_R)

    ws_min = MIN(ws_min_L, ws_min_R)
    ws_max = MAX(ws_max_L, ws_max_R)

    ! HLL fluxes
    if (ws_min .ge. 0.0) then
      F_W = F_L
    else if (ws_max .lt. 0.0) then
      F_W = F_R
    else
      F_W = (ws_min*ws_max*(U_R - U_L) + ws_max*F_L - ws_min*F_R)/(ws_max - ws_min);
    end if

    ! ---- East interface -----

    U_L = U(:,i)
    U_R = U(:,i+1)

    call compute_flux_ws_x(U_L, F_L, ws_max_L, ws_min_L)
    call compute_flux_ws_x(U_R, F_R, ws_max_R, ws_min_R)

    ws_min = MIN(ws_min_L, ws_min_R)
    ws_max = MAX(ws_max_L, ws_max_R)

    ! HLL fluxes
    if (ws_min .ge. 0.0) then
      F_E = F_L
    else if (ws_max .lt. 0.0) then
      F_E = F_R
    else
      F_E = (ws_min*ws_max*(U_R - U_L) + ws_max*F_L - ws_min*F_R)/(ws_max - ws_min);
    end if

  end subroutine

  ! ======================================================================== 

  subroutine compute_fluxes_Rusanov(U, i, F_W, F_E)

    ! Computes HLL numerical fluxes for the cell "i"
    ! 'F_W' and 'F_E' are the complete fluxes at the cell interfaces, to be then
    ! added in the time-discretized PDE.
    ! The variables 'F_L' and 'F_R' are temporary variables used for the 
    ! fluxes at the left and right of an interface.

    implicit none

    real(kind=8), dimension(Neq,Nx), intent(in)  :: U
    real(kind=8), dimension(Neq),    intent(out) :: F_W, F_E
    integer, intent(in) :: i

    integer :: eqID
    real(kind=8), dimension(Neq) :: U_L, U_R, F_L, F_R
    real(kind=8) :: ws_min_L, ws_max_L, ws_min_R, ws_max_R, ws_min, ws_max ! wave speeds

    ! ---- West interface -----

    U_L = U(:,i-1)
    U_R = U(:,i)

    call compute_flux_ws_x(U_L, F_L, ws_max_L, ws_min_L)
    call compute_flux_ws_x(U_R, F_R, ws_max_R, ws_min_R)

    ws_max = MAX(ABS(ws_max_L), ABS(ws_min_L), ABS(ws_max_R), ABS(ws_min_R))

    ! West interface
    F_W    = 0.5*(F_R + F_L) - ws_max/2.0*(U_R-U_L) ! Rusanov flux

    ! ---- East interface -----

    U_L = U(:,i)
    U_R = U(:,i+1)

    call compute_flux_ws_x(U_L, F_L, ws_max_L, ws_min_L)
    call compute_flux_ws_x(U_R, F_R, ws_max_R, ws_min_R)

    ws_max = MAX(ABS(ws_max_L), ABS(ws_min_L), ABS(ws_max_R), ABS(ws_min_R))

    ! West interface
    F_E    = 0.5*(F_R + F_L) - ws_max/2.0*(U_R-U_L) ! Rusanov flux

  end subroutine

end module
