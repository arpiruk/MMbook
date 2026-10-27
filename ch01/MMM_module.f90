module MMM_module
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private

   integer, parameter, public :: rk = real64

   public :: get_erf_cellvalue
   public :: update_periodic_ghosts_pv, update_periodic_ghosts_cv
   public :: get_spatial_derivative_inplace

contains

   elemental function get_erf_cellvalue(a, b, rho, c0, c1) result(value)
      real(rk), intent(in) :: a, b, rho, c0, c1
      real(rk) :: value
      real(rk) :: a_minus, a_plus, b_minus, b_plus
      real(rk) :: erf_group, exp_group, denominator
      real(rk), parameter :: pi = acos(-1.0_rk)

      a_minus = rho * (a - c0)
      a_plus  = rho * (a + c1)
      b_minus = rho * (b - c0)
      b_plus  = rho * (b + c1)

      erf_group = sqrt(pi) * ( &
           b_minus * erf(b_minus) &
         - b_plus  * erf(b_plus)  &
         - a_minus * erf(a_minus) &
         + a_plus  * erf(a_plus) )

      exp_group = &
         - exp(-a_minus*a_minus) &
         + exp(-b_minus*b_minus) &
         + exp(-a_plus*a_plus)   &
         - exp(-b_plus*b_plus)

      denominator = 2.0_rk * rho * sqrt(pi) * (a - b)
      value = (erf_group + exp_group) / denominator
   end function get_erf_cellvalue

   subroutine update_periodic_ghosts_pv(phi_arr, n, ng)
      real(rk), intent(inout) :: phi_arr(:)
      integer, intent(in) :: n, ng

      phi_arr(1:ng) = phi_arr(n+1:n+ng)
      phi_arr(n+ng+1:n+2*ng+1) = phi_arr(ng+1:2*ng+1)
   end subroutine update_periodic_ghosts_pv

   subroutine update_periodic_ghosts_cv(phi_arr, n, ng)
      real(rk), intent(inout) :: phi_arr(:)
      integer, intent(in) :: n, ng

      phi_arr(1:ng) = phi_arr(n+1:n+ng)
      phi_arr(n+ng+1:n+2*ng) = phi_arr(ng+1:2*ng)
   end subroutine update_periodic_ghosts_cv

   subroutine get_spatial_derivative_inplace(phi_pv, phi_cv, u, dx, n, ng, rhs_pv, rhs_cv, scheme)
      real(rk), intent(inout) :: phi_pv(:), phi_cv(:)
      real(rk), intent(in) :: u, dx
      integer, intent(in) :: n, ng
      real(rk), intent(inout) :: rhs_pv(:), rhs_cv(:)
      character(len=*), intent(in) :: scheme
      real(rk) :: updx, c1coef, c2coef
      integer :: i
      
       real(rk) :: cv1, cv2, cv3, cv4
       real(rk) :: pv1, pv2, pv3

      call update_periodic_ghosts_pv(phi_pv, n, ng)
      call update_periodic_ghosts_cv(phi_cv, n, ng)

      rhs_pv = 0.0_rk
      rhs_cv = 0.0_rk

      updx = u / dx

      do i = ng + 1, n + ng
         rhs_cv(i) = updx * (phi_pv(i+1) - phi_pv(i))
      end do

      c1coef = 2.0_rk * u / dx
      c2coef = 0.5_rk * u / dx

      select case (trim(adjustl(scheme)))
      case ('central4')
         do i = ng + 1, n + ng
            rhs_pv(i) = c1coef * (phi_cv(i) - phi_cv(i-1)) &
                      - c2coef * (phi_pv(i+1) - phi_pv(i-1))
         end do

      case ('central6')
         do i = ng + 1, n + ng
            rhs_pv(i) = updx * ( &
                 (9.0_rk / 4.0_rk)   * (phi_cv(i)   - phi_cv(i-1)) &
               + (1.0_rk / 36.0_rk)  * (phi_cv(i+1) - phi_cv(i-2)) &
               - (2.0_rk / 3.0_rk)   * (phi_pv(i+1) - phi_pv(i-1)) )
         end do

      case ('central8')
         do i = ng + 1, n + ng
            rhs_pv(i) = updx * ( &
                 (5.0_rk / 2.0_rk)   * (phi_cv(i)   - phi_cv(i-1)) &
               + (7.0_rk / 54.0_rk)  * (phi_cv(i+1) - phi_cv(i-2)) &
               - (8.0_rk / 9.0_rk)   * (phi_pv(i+1) - phi_pv(i-1)) &
               - (1.0_rk / 36.0_rk)  * (phi_pv(i+2) - phi_pv(i-2)) )
         end do
         
      case ('central10')
         do i = ng + 1, n + ng
            rhs_pv(i) = updx * ( &
                 (1.0_rk / 900.0_rk) * (phi_cv(i+2) - phi_cv(i-3)) &
               + (7.0_rk / 36.0_rk)  * (phi_cv(i+1) - phi_cv(i-2)) &
               + (47.0_rk / 18.0_rk) * (phi_cv(i)   - phi_cv(i-1)) &
               - 1.0_rk              * (phi_pv(i+1) - phi_pv(i-1)) &
               - (1.0_rk / 20.0_rk)  * (phi_pv(i+2) - phi_pv(i-2)) )
         end do        
         
   case ('central12')
    updx = u / dx
         cv1 = updx*49.0_rk /18.0_rk    ; pv1 = -updx*9.0_rk/8.0_rk
         cv2 = updx*41.0_rk/144.0_rk    ; pv2 = -updx*9.0_rk/100.0_rk
         cv3 = updx*157.0_rk/ 18000.0_rk; pv3 = -updx*1.0_rk/600.0_rk
         
          do i = ng + 1, n + ng
            rhs_pv(i) = &
              ! --- Beta Terms (Cell Values) ---
                  cv1 * (phi_cv(i  ) - phi_cv(i-1)) & ! \tau = 1 (i+1/2 and i-1/2)
                + cv2 * (phi_cv(i+1) - phi_cv(i-2)) & ! \tau = 2 (i+3/2 and i-3/2)
                + cv3 * (phi_cv(i+2) - phi_cv(i-3)) & ! \tau = 3 (i+5/2 and i-5/2)
              ! --- Gamma Terms (Point Values) ---
                + pv1 * (phi_pv(i+1) - phi_pv(i-1)) & ! \tau = 1
                + pv2 * (phi_pv(i+2) - phi_pv(i-2)) & ! \tau = 2
                + pv3 * (phi_pv(i+3) - phi_pv(i-3))   ! \tau = 3
         end do    
         case ('central14')
         updx = u / dx
         cv1 = updx * 401.0_rk  / 144.0_rk   ; pv1 = -updx * 6.0_rk / 5.0_rk
         cv2 = updx * 1241.0_rk / 3600.0_rk  ; pv2 = -updx * 3.0_rk / 25.0_rk
         cv3 = updx * 301.0_rk  / 18000.0_rk ; pv3 = -updx * 2.0_rk / 525.0_rk
         cv4 = updx * 1.0_rk    / 19600.0_rk
         
         do i = ng + 1, n + ng
            rhs_pv(i) = &
                ! --- Beta Terms (Cell Values) ---
                  cv1 * (phi_cv(i  ) - phi_cv(i-1)) & ! \tau = 1
                + cv2 * (phi_cv(i+1) - phi_cv(i-2)) & ! \tau = 2
                + cv3 * (phi_cv(i+2) - phi_cv(i-3)) & ! \tau = 3
                + cv4 * (phi_cv(i+3) - phi_cv(i-4)) & ! \tau = 4
                ! --- Gamma Terms (Point Values) ---
                + pv1 * (phi_pv(i+1) - phi_pv(i-1)) & ! \tau = 1
                + pv2 * (phi_pv(i+2) - phi_pv(i-2)) & ! \tau = 2
                + pv3 * (phi_pv(i+3) - phi_pv(i-3))   ! \tau = 3
         end do               
    

      case ('upwind1')
         do i = ng + 1, n + ng
            rhs_pv(i) = updx * (2.0_rk * (phi_pv(i) - phi_cv(i-1)))
         end do

      case ('upwind3')
         do i = ng + 1, n + ng
            rhs_pv(i) = u * ( &
                 -3.5_rk * phi_cv(i-1) &
                 +0.5_rk * phi_cv(i)   &
                 +1.0_rk * phi_pv(i-1) &
                 +2.0_rk * phi_pv(i) ) / dx
         end do

      case ('upwind5')
         do i = ng + 1, n + ng
            rhs_pv(i) = updx * ( &
                 -(1.0_rk/18.0_rk)  * phi_cv(i-2) &
                 -(55.0_rk/18.0_rk) * phi_cv(i-1) &
                 +(13.0_rk/9.0_rk)  * phi_cv(i)   &
                 + phi_pv(i-1) + phi_pv(i)         &
                 -(1.0_rk/3.0_rk)   * phi_pv(i+1) )
         end do

      case ('upwind7')
         do i = ng + 1, n + ng
            rhs_pv(i) = updx * ( &
                 -(53.0_rk/216.0_rk)  * phi_cv(i-2) &
                 -(725.0_rk/216.0_rk) * phi_cv(i-1) &
                 +(355.0_rk/216.0_rk) * phi_cv(i)   &
                 +(1.0_rk/72.0_rk)    * phi_cv(i+1) &
                 +(1.0_rk/18.0_rk)    * phi_pv(i-2) &
                 +(4.0_rk/3.0_rk)     * phi_pv(i-1) &
                 + phi_pv(i)                           &
                 -(4.0_rk/9.0_rk)     * phi_pv(i+1) )
         end do

      case default
         error stop 'Unknown spatial scheme.'
      end select
   end subroutine get_spatial_derivative_inplace

end module MMM_module





