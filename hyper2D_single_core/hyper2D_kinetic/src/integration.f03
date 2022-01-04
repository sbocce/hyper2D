module integration

  use pde
  use global_module

  implicit none
  
  real(kind=8), dimension(Neq,Nx,Ny) :: U_new ! Working variable

  contains

  ! ======================================================================== 
  
  subroutine forward_Euler_step(U, dt)

    ! This function performs one step of the Forward Euler explicit time integrator 

    implicit none
 
    real(kind=8), dimension(Neq,Nx,Ny), intent(inout)  :: U

    real(kind=8), intent(in) :: dt
    integer :: i, j, eqID
    
    real(kind=8), dimension(Neq) :: F_N, F_S, F_W, F_E 

    ! First, assign BCs. Wall BCs need to be updated in time, taking the value
    ! from neighboring cells.
    call assign_BCs(U) ! See the pde.f03 module

    ! Update solution using the forward Euler integrator - INTERNAL CELLS only
    do j = 3, Ny-2
      do i = 3, Nx-2

        call compute_fluxes_Mieussens(U, i, j, F_N, F_S, F_W, F_E)

        U_new(:,i,j) = U(:,i,j) - dt/dx*(F_E - F_W) - dt/dy*(F_N - F_S)
        
        ! Check that the solution did not diverge
        do eqID = 1, Neq
          if (isnan(U_new(eqID,i,j))) then 
            print*, 'Solution diverged, try with a smaller time step! Aborting.'
            stop
          end if
        end do

      end do
    end do

    ! Save solution (internal cells ONLY! Do not overwrite ghost cells!)
    U(:,3:Nx-2,3:Ny-2) = U_new(:,3:Nx-2,3:Ny-2)

  end subroutine

  ! ======================================================================== 
  
  subroutine midpoint_Euler_step(U, dt)

    ! This function performs one step of the Forward Euler explicit time integrator 

    implicit none
 
    real(kind=8), dimension(Neq,Nx,Ny), intent(inout)  :: U

    real(kind=8), intent(in) :: dt
    integer :: i, j, eqID
    
    real(kind=8), dimension(Neq) :: F_N, F_S, F_W, F_E 

    ! First, assign BCs. Wall BCs need to be updated in time, taking the value
    ! from neighboring cells.
    call assign_BCs(U) ! See the pde.f03 module

    ! Compute first half step - Create U_new that is actually U_half
    do j = 3, Ny-2
      do i = 3, Nx-2
        call compute_fluxes_Mieussens(U, i, j, F_N, F_S, F_W, F_E)
        U_new(:,i,j) = U(:,i,j) - dt/2.0/dx*(F_E - F_W) - dt/2.0/dy*(F_N - F_S)
      end do
    end do

    call assign_BCs(U_new) ! See the pde.f03 module

    ! Compute second half step - Use U_new (not U!!!! Be careful!) to compute the fluxes, 
    ! and update U
    do j = 3, Ny-2
      do i = 3, Nx-2

        ! Compute fluxes from U_new (= U_half). 
        ! If you compute fluxes directly from U, you make a mistake, because you are
        ! updating U dynamically! You would compute weird fluxes!
        call compute_fluxes_Mieussens(U_new, i, j, F_N, F_S, F_W, F_E)
        U(:,i,j) = U(:,i,j) - dt/dx*(F_E - F_W) - dt/dy*(F_N - F_S)
        
        ! Check that the solution did not diverge
        do eqID = 1, Neq
          if (isnan(U(eqID,i,j))) then 
            print*, 'Solution diverged, try with a smaller time step! Aborting.'
            stop
          end if
        end do

      end do
    end do

  end subroutine



  ! ======================================================================== 

  subroutine compute_fluxes_Mieussens(U, i, j, F_N, F_S, F_W, F_E)

    ! Computes Rusanov numerical fluxes for the cell (i,j)

    implicit none

    real(kind=8), dimension(Neq,Nx,Ny), intent(in)  :: U
    real(kind=8), dimension(Neq),       intent(out) :: F_N, F_S, F_W, F_E 
    integer, intent(in) :: i, j

    integer :: eqID
    real(kind=8), dimension(Neq) :: U_L, U_R, F_L, F_R
    real(kind=8) :: x_now, v_now, E_now, Psi

    x_now = compute_x_from_i(i)
    v_now = compute_v_from_j(j)
    E_now = E_field(compute_x_from_i(i))

    ! North flux
    call slope_fun(U(1,i,j-1),U(1,i,j),U(1,i,j+1),U(1,i,j+2), Psi)
    F_N = 0.5*(q*E_now/M*U(1,i,j) + q*E_now/M*U(1,i,j+1) - abs(q*E_now/M)*(U(1,i,j+1) - U(1,i,j) - Psi))

    ! South flux
    call slope_fun(U(1,i,j-2),U(1,i,j-1),U(1,i,j),U(1,i,j+1), Psi)
    F_S = 0.5*(q*E_now/M*U(1,i,j-1) + q*E_now/M*U(1,i,j) - abs(q*E_now/M)*(U(1,i,j) - U(1,i,j-1) - Psi))

    ! East flux
    call slope_fun(U(1,i-1,j),U(1,i,j),U(1,i+1,j),U(1,i+2,j), Psi)
    F_E = 0.5*(v_now*U(1,i,j) + v_now*U(1,i+1,j) - abs(v_now)*(U(1,i+1,j) - U(1,i,j) - Psi))

    ! West flux
    call slope_fun(U(1,i-2,j),U(1,i-1,j),U(1,i,j),U(1,i+1,j), Psi)
    F_W = 0.5*(v_now*U(1,i-1,j) + v_now*U(1,i,j) - abs(v_now)*(U(1,i,j) - U(1,i-1,j) - Psi))

  end subroutine

  ! =======================================================================

  subroutine slope_fun(U_m2, U_m1, U_p1, U_p2, Psi)

    implicit none
   
    real(kind=8), intent(in)  :: U_m2, U_m1, U_p1, U_p2
    real(kind=8), intent(out) :: Psi

    Psi = 0.0 ! Initialize as first order

    if (second_order_bool) then
      Psi = minmod((U_m1-U_m2), (U_p1-U_m1), (U_p2-U_p1))
    end if

  end subroutine

  ! =======================================================================

  function minmod(a,b,c)

    implicit none

    real(kind=8), intent(in) :: a, b, c 
 
    real(kind=8) :: minmod

    if ( (a .ge. 0.0) .and. (b .ge. 0.0) .and. (c .ge. 0.0) ) then
      minmod = min(a,b,c)
    else if ( (a .le. 0.0) .and. (b .le. 0.0) .and. (c .le. 0.0) ) then
      minmod = -min(abs(a), abs(b), abs(c))
    else
      minmod = 0.0
    end if

  end function

end module
