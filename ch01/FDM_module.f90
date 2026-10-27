module FDM_module
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none

   integer, parameter :: rk = real64

type :: work_arrays
      real(rk), allocatable :: tmp_pv(:)
      real(rk), allocatable :: k1_pv(:), k2_pv(:), k3_pv(:), k4_pv(:)
      real(rk), allocatable :: k5_pv(:), k6_pv(:), k7_pv(:), k8_pv(:) ! Added for RK5/7
! New registers for RK8
      real(rk), allocatable :: k9_pv(:) 
      real(rk), allocatable :: k10_pv(:) 
      real(rk), allocatable :: k11_pv(:)   ,  k12_pv(:)  , k13_pv(:)  , k14_pv(:)  , k15_pv(:)     , k16_pv(:)
      real(rk), allocatable :: rhs_inner(:)
      real(rk), allocatable :: mat_a(:,:), mat_b(:,:)
   end type work_arrays

contains

   subroutine allocate_work(work, npx, ncell)
      type(work_arrays), intent(inout) :: work
      integer, intent(in) :: npx, ncell

      allocate(work%tmp_pv(npx))
      allocate(work%k1_pv(npx), work%k2_pv(npx), work%k3_pv(npx), work%k4_pv(npx))
      allocate(work%k5_pv(npx), work%k6_pv(npx), work%k7_pv(npx), work%k8_pv(npx))  

 
      allocate(work%k9_pv(npx) )
       allocate(work%k10_pv(npx)) 
      allocate(work%k11_pv(npx) )
      allocate(work%k12_pv(npx), work%k13_pv(npx), work%k14_pv(npx), work%k15_pv(npx), work%k16_pv(npx)) 
      
      allocate(work%rhs_inner(ncell))
      allocate(work%mat_a(ncell,ncell), work%mat_b(ncell,ncell))



      work%tmp_pv = 0.0_rk
      work%k1_pv = 0.0_rk
      work%k2_pv = 0.0_rk
      work%k3_pv = 0.0_rk
      work%k4_pv = 0.0_rk
      work%k5_pv = 0.0_rk; work%k6_pv = 0.0_rk
      work%k7_pv = 0.0_rk; work%k8_pv = 0.0_rk   
      work%k9_pv = 0.0_rk; work%k10_pv = 0.0_rk 
      work%k11_pv = 0.0_rk; work%k12_pv = 0.0_rk     
            work%k13_pv = 0.0_rk; work%k14_pv = 0.0_rk  
                  work%k15_pv = 0.0_rk; work%k16_pv = 0.0_rk  
      work%rhs_inner = 0.0_rk
      work%mat_a = 0.0_rk
      work%mat_b = 0.0_rk
   end subroutine allocate_work

   subroutine update_periodic_ghosts_pv(phi_arr, n, ng)
      real(rk), intent(inout) :: phi_arr(:)
      integer, intent(in) :: n, ng

      ! Python equivalent:
      ! phi_arr[0:ng] = phi_arr[N:N+ng]
      ! phi_arr[N+ng:N+2*ng+1] = phi_arr[ng:2*ng+1]
      
      phi_arr(1:ng) = phi_arr(n+1:n+ng)
      phi_arr(n+ng+1:n+2*ng+1) = phi_arr(ng+1:2*ng+1)
   end subroutine update_periodic_ghosts_pv

   subroutine get_spatial_derivative_inplace(phi_pv, u, dx, n, ng, rhs_pv, scheme)
      real(rk), intent(inout) :: phi_pv(:)
      real(rk), intent(in) :: u, dx
      integer, intent(in) :: n, ng
      real(rk), intent(inout) :: rhs_pv(:)
      character(len=*), intent(in) :: scheme
      integer :: j
      real(rk) :: c1
      real(rk) :: pv1,pv2,pv3,pv4,pv5,pv6,pv7

      call update_periodic_ghosts_pv(phi_pv, n, ng)
      rhs_pv = 0.0_rk

      select case (trim(scheme))
      case ('central2')
         do j = ng + 1, n + ng
            rhs_pv(j) = u * (phi_pv(j+1) - phi_pv(j-1)) / (2.0_rk * dx)
         end do

      case ('central4')
         c1 = u / (12.0_rk * dx)
         do j = ng + 1, n + ng
            rhs_pv(j) = c1 * (-phi_pv(j+2) + 8.0_rk*phi_pv(j+1) &
                              -8.0_rk*phi_pv(j-1) + phi_pv(j-2))
         end do

      case ('central6')
         c1 = u / (60.0_rk * dx)
         do j = ng + 1, n + ng
            rhs_pv(j) = c1 * (           -phi_pv(j-3) +  9.0_rk * phi_pv(j-2) - 45.0_rk * phi_pv(j-1) &
                                         +phi_pv(j+3) -  9.0_rk * phi_pv(j+2) + 45.0_rk * phi_pv(j+1)  )
         end do
         
!         soldfEX6 := Vector[row](6, {(1) = -(1/60)/h, (2) = (3/20)/h, (3) = -(3/4)/h, (4) = (3/4)/h, (5) = -(3/20)/h, (6) = (1/60)/h})
         

      case ('central8')
         c1 = u / (840.0_rk * dx)
         do j = ng + 1, n + ng
            rhs_pv(j) = c1 * (+3.0_rk * phi_pv(j-4) -32.0_rk * phi_pv(j-3)  +168.0_rk * phi_pv(j-2) -672.0_rk * phi_pv(j-1) &  
                              -3.0_rk * phi_pv(j+4) +32.0_rk * phi_pv(j+3)  -168.0_rk * phi_pv(j+2) +672.0_rk * phi_pv(j+1) &
                                )
         end do
         
      case ('central10')
         c1 =  u / (2520.0_rk * dx)
         do j = ng + 1, n + ng
            rhs_pv(j) = c1 * ( -2.0_rk * phi_pv(j-5) +  25.0_rk * phi_pv(j-4) - 150.0_rk * phi_pv(j-3) + 600.0_rk * phi_pv(j-2) - 2100.0_rk * phi_pv(j-1) &
                               +2.0_rk * phi_pv(j+5) -  25.0_rk * phi_pv(j+4) + 150.0_rk * phi_pv(j+3) - 600.0_rk * phi_pv(j+2) + 2100.0_rk * phi_pv(j+1) )
         end do         

	case ('central12')
         c1 =  u / (27720.0_rk * dx)
         do j = ng + 1, n + ng
            rhs_pv(j) = c1 * (                                         &
               +     5.0_rk * phi_pv(j-6)  - 5.0_rk   * phi_pv(j+6)    &
               -    72.0_rk * phi_pv(j-5)  + 72.0_rk  * phi_pv(j+5)    &
               +   495.0_rk * phi_pv(j-4)  - 495.0_rk * phi_pv(j+4)    &
               -  2200.0_rk * phi_pv(j-3) + 2200.0_rk * phi_pv(j+3)    &
               +  7425.0_rk * phi_pv(j-2) - 7425.0_rk * phi_pv(j+2)    &
               - 23760.0_rk * phi_pv(j-1)+ 23760.0_rk * phi_pv(j+1) )
               
        
               
 !!!!5/h,-72/h, 495/h, -2200/h, 7425/h, -23760/h, 23760/h,-7425/h, 2200/h, -495/h, 72/h, -5h    
               
               
         end do
 
 	case ('central14')
             c1 =  u /  dx 
	     pv7=   1.0_rk*c1 / 24024.0_rk
	     pv6=  -7.0_rk*c1  / 10296.0_rk
	     pv5=   7.0_rk*c1  / 1320.0_rk
	     pv4=  -7.0_rk*c1  / 264.0_rk
	     pv3=   7.0_rk*c1  / 72.0_rk
	     pv2=  -7.0_rk*c1  / 24.0_rk
	     pv1=   7.0_rk*c1  / 8.0_rk
         do j = ng + 1, n + ng
            rhs_pv(j) =                                  &
               +  pv7 * (phi_pv(j+7)  - phi_pv(j-7) )    &            
               +  pv6 * (phi_pv(j+6)  - phi_pv(j-6) )    &
               +  pv5 * (phi_pv(j+5)  - phi_pv(j-5) )    &
               +  pv4 * (phi_pv(j+4)  - phi_pv(j-4) )    &
               +  pv3 * (phi_pv(j+3)  - phi_pv(j-3) )    &
               +  pv2 * (phi_pv(j+2)  - phi_pv(j-2) )    &
               +  pv1 * (phi_pv(j+1)  - phi_pv(j-1) )
               

        
!  -1/(24024*h), 7/(10296*h), -7/(1320*h), 7/(264*h), -7/(72*h), 7/(24*h), -7/(8*h), 7/(8*h), -7/(24*h), 7/(72*h), -7/(264*h), 7/(1320*h), -7/(10296*h), 1/(24024*h)
         end do
 
 
    
 
	 case ('central101')
	   c1 = u / (2520.0_rk * dx)
	   do j = ng + 1, n + ng
	      rhs_pv(j) = c1 * ( -2.0_rk * phi_pv(j-5)  &
		                + 25.0_rk * phi_pv(j-4)  &
		                -150.0_rk * phi_pv(j-3)  &
		                +600.0_rk * phi_pv(j-2)  &
		                -2100.0_rk * phi_pv(j-1) &
		                +2100.0_rk * phi_pv(j+1) &
		                -600.0_rk * phi_pv(j+2)  &
		                +150.0_rk * phi_pv(j+3)  &
		                - 25.0_rk * phi_pv(j+4)  &
		                +  2.0_rk * phi_pv(j+5) )
	   end do     

      case ('upwind1')
         do j = ng + 1, n + ng
            rhs_pv(j) = u * (phi_pv(j) - phi_pv(j-1)) / dx
         end do

      case ('upwind3')
         c1 = u / (6.0_rk * dx)
         do j = ng + 1, n + ng
            rhs_pv(j) = c1 * (2.0_rk*phi_pv(j+1) + 3.0_rk*phi_pv(j) &
                              -6.0_rk*phi_pv(j-1) + phi_pv(j-2))
         end do

      case ('bl34')
         do j = ng + 1, n + ng
            rhs_pv(j) = (u/(6.0_rk*dx)) * (2.0_rk*phi_pv(j+1) + 3.0_rk*phi_pv(j) &
                        -6.0_rk*phi_pv(j-1) + phi_pv(j-2)) &
                      + (u/(12.0_rk*dx)) * (-phi_pv(j+2) + 8.0_rk*phi_pv(j+1) &
                        -8.0_rk*phi_pv(j-1) + phi_pv(j-2))
         end do

      case default
         error stop 'Unknown spatial scheme.'
      end select

      call update_periodic_ghosts_pv(rhs_pv, n, ng)
   end subroutine get_spatial_derivative_inplace

   subroutine advance_one_step_inplace(phi_pv, u, dx, dt, n, ng, work, method, scheme)
      real(rk), intent(inout) :: phi_pv(:)
      real(rk), intent(in) :: u, dx, dt
      integer, intent(in) :: n, ng
      type(work_arrays), intent(inout) :: work
      character(len=*), intent(in) :: method, scheme
      real(rk) :: factor
      integer :: i
      real(rk) :: s

      select case (trim(method))
      case ('euler')
         call get_spatial_derivative_inplace(phi_pv, u, dx, n, ng, work%k1_pv, scheme)
         phi_pv = phi_pv - dt * work%k1_pv
         call update_periodic_ghosts_pv(phi_pv, n, ng)

      case ('rk4')
         call get_spatial_derivative_inplace(phi_pv, u, dx, n, ng, work%k1_pv, scheme)

         work%tmp_pv = phi_pv - (0.5_rk * dt) * work%k1_pv
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k2_pv, scheme)

         work%tmp_pv = phi_pv - (0.5_rk * dt) * work%k2_pv
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k3_pv, scheme)

         work%tmp_pv = phi_pv - dt * work%k3_pv
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k4_pv, scheme)

         factor = dt / 6.0_rk
         phi_pv = phi_pv - factor*work%k1_pv - 2.0_rk*factor*work%k2_pv &
                         - 2.0_rk*factor*work%k3_pv - factor*work%k4_pv
         call update_periodic_ghosts_pv(phi_pv, n, ng)

case ('rk5')
         ! 6-Stage, 5th-Order Explicit Runge-Kutta Scheme (Butcher)
         call get_spatial_derivative_inplace(phi_pv, u, dx, n, ng, work%k1_pv, scheme)

         work%tmp_pv = phi_pv - dt * (1.0_rk/4.0_rk) * work%k1_pv
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k2_pv, scheme)

         work%tmp_pv = phi_pv - dt * ((1.0_rk/8.0_rk)*work%k1_pv + (1.0_rk/8.0_rk)*work%k2_pv)
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k3_pv, scheme)

         work%tmp_pv = phi_pv - dt * (-(1.0_rk/2.0_rk)*work%k2_pv + 1.0_rk*work%k3_pv)
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k4_pv, scheme)

         work%tmp_pv = phi_pv - dt * ((3.0_rk/16.0_rk)*work%k1_pv + (9.0_rk/16.0_rk)*work%k4_pv)
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k5_pv, scheme)

         work%tmp_pv = phi_pv - dt * (-(3.0_rk/7.0_rk)*work%k1_pv + (2.0_rk/7.0_rk)*work%k2_pv &
                                    + (12.0_rk/7.0_rk)*work%k3_pv - (12.0_rk/7.0_rk)*work%k4_pv &
                                    + (8.0_rk/7.0_rk)*work%k5_pv)
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k6_pv, scheme)

         ! Final step combination
         phi_pv = phi_pv - dt * ((7.0_rk/90.0_rk)*work%k1_pv + (32.0_rk/90.0_rk)*work%k3_pv &
                               + (12.0_rk/90.0_rk)*work%k4_pv + (32.0_rk/90.0_rk)*work%k5_pv &
                               + (7.0_rk/90.0_rk)*work%k6_pv)
         call update_periodic_ghosts_pv(phi_pv, n, ng)

case ('rk7')
         ! 9-Stage, 7th-Order Explicit Runge-Kutta Scheme (Butcher, p. 196)
         
         ! Stage 1
         call get_spatial_derivative_inplace(phi_pv, u, dx, n, ng, work%k1_pv, scheme)

         ! Stage 2: c2 = 1/6
         work%tmp_pv = phi_pv - dt * (1.0_rk / 6.0_rk) * work%k1_pv
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k2_pv, scheme)

         ! Stage 3: c3 = 1/3
         work%tmp_pv = phi_pv - dt * (1.0_rk / 3.0_rk) * work%k2_pv
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k3_pv, scheme)

         ! Stage 4: c4 = 1/2
         work%tmp_pv = phi_pv - dt * ((1.0_rk / 8.0_rk) * work%k1_pv + (3.0_rk / 8.0_rk) * work%k3_pv)
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k4_pv, scheme)

         ! Stage 5: c5 = 2/11
         work%tmp_pv = phi_pv - dt * ((148.0_rk / 1331.0_rk) * work%k1_pv &
                                    + (150.0_rk / 1331.0_rk) * work%k3_pv &
                                    - (56.0_rk / 1331.0_rk)  * work%k4_pv)
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k5_pv, scheme)

         ! Stage 6: c6 = 2/3
         work%tmp_pv = phi_pv - dt * (-(404.0_rk / 243.0_rk)   * work%k1_pv &
                                    - (170.0_rk / 27.0_rk)    * work%k3_pv &
                                    + (4024.0_rk / 1701.0_rk) * work%k4_pv &
                                    + (10648.0_rk / 1701.0_rk)* work%k5_pv)
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k6_pv, scheme)

         ! Stage 7: c7 = 6/7
         work%tmp_pv = phi_pv - dt * ((2466.0_rk / 2401.0_rk)    * work%k1_pv &
                                    + (1242.0_rk / 343.0_rk)     * work%k3_pv &
                                    - (19176.0_rk / 16807.0_rk)  * work%k4_pv &
                                    - (51909.0_rk / 16807.0_rk)  * work%k5_pv &
                                    + (1053.0_rk / 2401.0_rk)    * work%k6_pv)
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k7_pv, scheme)

         ! Stage 8: c8 = 1/7
         work%tmp_pv = phi_pv - dt * ((5.0_rk / 154.0_rk)       * work%k1_pv &
                                    + (96.0_rk / 539.0_rk)      * work%k4_pv &
                                    - (1815.0_rk / 20384.0_rk)  * work%k5_pv &
                                    - (405.0_rk / 2464.0_rk)    * work%k6_pv &
                                    + (49.0_rk / 1144.0_rk)     * work%k7_pv)
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k8_pv, scheme)

         ! Stage 9: c9 = 1.0 (Note: Requires temporary k9 storage or reusing rhs_inner)
         ! We allocate work%rhs_inner(ncell) in your module, but since we need an array of size npx, 
         ! make sure work%rhs_inner is resized to npx or add an explicit work%k9_pv buffer.
         ! Assuming work%rhs_inner is updated to size npx, we use it as a scratchpad here:
         work%tmp_pv = phi_pv - dt * (-(113.0_rk / 32.0_rk)    * work%k1_pv &
                                    - (195.0_rk / 22.0_rk)    * work%k3_pv &
                                    + (32.0_rk / 7.0_rk)      * work%k4_pv &
                                    + (29403.0_rk / 3584.0_rk)* work%k5_pv &
                                    - (729.0_rk / 512.0_rk)   * work%k6_pv &
                                    + (1029.0_rk / 1408.0_rk) * work%k7_pv &
                                    + (21.0_rk / 16.0_rk)     * work%k8_pv)
         ! We will hijack work%rhs_inner here for the 9th stage derivative computation
         ! Alternatively, to keep it clean, name it as a separate collection step:
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k3_pv, scheme) ! Overwriting k3 since it's no longer needed in weights

         ! Final step combination (using b-weights vector)
         phi_pv = phi_pv - dt * ((32.0_rk / 105.0_rk)             * work%k4_pv &
                               + (1771561.0_rk / 6289920.0_rk)    * work%k5_pv &
                               + (243.0_rk / 2560.0_rk)           * work%k6_pv &
                               + (16807.0_rk / 74880.0_rk)        * work%k7_pv &
                               + (77.0_rk / 1440.0_rk)            * work%k8_pv & 
                               + (11.0_rk / 270.0_rk)             * work%k3_pv ) ! work%k3_pv contains stage 9 evaluation
         call update_periodic_ghosts_pv(phi_pv, n, ng)
 
case ('rk8')
         ! Cooper-Verner 11-stage, 8th-order explicit RK adapted for pure FDM
         s = sqrt(21.0_rk)

         ! Stage 1
         call get_spatial_derivative_inplace(phi_pv, u, dx, n, ng, work%k1_pv, scheme)

         ! Stage 2
         work%tmp_pv = phi_pv - dt * ((1.0_rk / 2.0_rk) * work%k1_pv)
         call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k2_pv, scheme)

         ! Stage 3
         work%tmp_pv = phi_pv - dt * ((1.0_rk / 4.0_rk) * work%k1_pv + (1.0_rk / 4.0_rk) * work%k2_pv)
         call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k3_pv, scheme)

         ! Stage 4, c = 1/2 - sqrt(21)/14
         work%tmp_pv = phi_pv - dt * ((1.0_rk / 7.0_rk) * work%k1_pv &
            + (-1.0_rk / 14.0_rk + 3.0_rk * s / 98.0_rk) * work%k2_pv &
            + (3.0_rk / 7.0_rk - 5.0_rk * s / 49.0_rk) * work%k3_pv)
         call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k4_pv, scheme)

         ! Stage 5
         work%tmp_pv = phi_pv - dt * ((11.0_rk / 84.0_rk - s / 84.0_rk) * work%k1_pv &
            + (2.0_rk / 7.0_rk - 4.0_rk * s / 63.0_rk) * work%k3_pv &
            + (1.0_rk / 12.0_rk + s / 252.0_rk) * work%k4_pv)
         call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k5_pv, scheme)

         ! Stage 6
         work%tmp_pv = phi_pv - dt * ((5.0_rk / 48.0_rk - s / 48.0_rk) * work%k1_pv &
            + (1.0_rk / 4.0_rk - s / 36.0_rk) * work%k3_pv &
            + (-77.0_rk / 120.0_rk - 7.0_rk * s / 180.0_rk) * work%k4_pv &
            + (63.0_rk / 80.0_rk + 7.0_rk * s / 80.0_rk) * work%k5_pv)
         call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k6_pv, scheme)

         ! Stage 7
         work%tmp_pv = phi_pv - dt * ((5.0_rk / 21.0_rk + s / 42.0_rk) * work%k1_pv &
            + (-48.0_rk / 35.0_rk - 92.0_rk * s / 315.0_rk) * work%k3_pv &
            + (211.0_rk / 30.0_rk + 29.0_rk * s / 18.0_rk) * work%k4_pv &
            + (-36.0_rk / 5.0_rk - 23.0_rk * s / 14.0_rk) * work%k5_pv &
            + (9.0_rk / 5.0_rk + 13.0_rk * s / 35.0_rk) * work%k6_pv)
         call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k7_pv, scheme)

         ! Stage 8
         work%tmp_pv = phi_pv - dt * ((1.0_rk / 14.0_rk) * work%k1_pv &
            + (1.0_rk / 9.0_rk + s / 42.0_rk) * work%k5_pv &
            + (13.0_rk / 63.0_rk + s / 21.0_rk) * work%k6_pv &
            + (1.0_rk / 9.0_rk) * work%k7_pv)
         call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k8_pv, scheme)

         ! Stage 9
         work%tmp_pv = phi_pv - dt * ((1.0_rk / 32.0_rk) * work%k1_pv &
            + (91.0_rk / 576.0_rk + 7.0_rk * s / 192.0_rk) * work%k5_pv &
            + (11.0_rk / 72.0_rk) * work%k6_pv &
            + (-385.0_rk / 1152.0_rk + 25.0_rk * s / 384.0_rk) * work%k7_pv &
            + (63.0_rk / 128.0_rk - 13.0_rk * s / 128.0_rk) * work%k8_pv)
         call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k9_pv, scheme)

         ! Stage 10
         work%tmp_pv = phi_pv - dt * ((1.0_rk / 14.0_rk) * work%k1_pv &
            + (1.0_rk / 9.0_rk) * work%k5_pv &
            + (-733.0_rk / 2205.0_rk + s / 15.0_rk) * work%k6_pv &
            + (515.0_rk / 504.0_rk - 37.0_rk * s / 168.0_rk) * work%k7_pv &
            + (-51.0_rk / 56.0_rk + 11.0_rk * s / 56.0_rk) * work%k8_pv &
            + (132.0_rk / 245.0_rk - 4.0_rk * s / 35.0_rk) * work%k9_pv)
         call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k10_pv, scheme)

         ! Stage 11
         work%tmp_pv = phi_pv - dt * ((-7.0_rk / 3.0_rk - 7.0_rk * s / 18.0_rk) * work%k5_pv &
            + (-2.0_rk / 5.0_rk - 28.0_rk * s / 45.0_rk) * work%k6_pv &
            + (-91.0_rk / 24.0_rk + 53.0_rk * s / 72.0_rk) * work%k7_pv &
            + (301.0_rk / 72.0_rk - 53.0_rk * s / 72.0_rk) * work%k8_pv &
            + (28.0_rk / 45.0_rk + 28.0_rk * s / 45.0_rk) * work%k9_pv &
            + (49.0_rk / 18.0_rk + 7.0_rk * s / 18.0_rk) * work%k10_pv)
         call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
         call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k11_pv, scheme)

         ! Final 8th-order assembly
         phi_pv = phi_pv - dt * ((1.0_rk / 20.0_rk) * work%k1_pv &
            + (49.0_rk / 180.0_rk) * work%k8_pv &
            + (16.0_rk / 45.0_rk) * work%k9_pv &
            + (49.0_rk / 180.0_rk) * work%k10_pv &
            + (1.0_rk / 20.0_rk) * work%k11_pv)

         call update_periodic_ghosts_pv(phi_pv, n, ng) 
case ('rk9')

   ! ---- Stage 1 -------------------------------------------------
   call get_spatial_derivative_inplace(phi_pv, u, dx, n, ng, work%k1_pv, scheme)

   ! ---- Stage 2 -------------------------------------------------
   ! c2 ~ 0.097024
   work%tmp_pv = phi_pv - dt*(0.03462_rk * work%k1_pv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k2_pv, scheme)

   ! ---- Stage 3 -------------------------------------------------
   ! c3 ~ 0.145537
   work%tmp_pv = phi_pv - dt*( &
       (-0.03893354388572875_rk)*work%k1_pv + &
       ( 0.13595789452450918_rk)*work%k2_pv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k3_pv, scheme)

   ! ---- Stage 4 -------------------------------------------------
   ! c4 = 0.561
   work%tmp_pv = phi_pv - dt*( &
       (0.03638413148954267_rk)*work%k1_pv + &
       (0.10915239446862801_rk)*work%k3_pv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k4_pv, scheme)

   ! ---- Stage 5 -------------------------------------------------
   ! c5 ~ 0.229008
   work%tmp_pv = phi_pv - dt*( &
       ( 2.0257639143939694_rk)*work%k1_pv + &
       (-7.638023836496291_rk )*work%k3_pv + &
       ( 6.173259922102322_rk )*work%k4_pv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k5_pv, scheme)

   ! ---- Stage 6 -------------------------------------------------
   ! c6 ~ 0.544992
   work%tmp_pv = phi_pv - dt*( &
       (0.05112275589406061_rk )*work%k1_pv + &
       (0.17708237945550218_rk )*work%k4_pv + &
       (0.0008027762409222536_rk)*work%k5_pv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k6_pv, scheme)

   ! ---- Stage 7 -------------------------------------------------
   ! c7 = 0.645
   work%tmp_pv = phi_pv - dt*( &
       ( 0.13160063579752163_rk)*work%k1_pv + &
       (-0.2957276252669636_rk )*work%k4_pv + &
       ( 0.08781378035642955_rk)*work%k5_pv + &
       ( 0.6213052975225274_rk )*work%k6_pv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k7_pv, scheme)

   ! ---- Stage 8 -------------------------------------------------
   ! c8 = 0.48375
   work%tmp_pv = phi_pv - dt*( &
       (0.07166666666666667_rk)*work%k1_pv + &
       (0.33055335789153195_rk)*work%k6_pv + &
       (0.2427799754418014_rk )*work%k7_pv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k8_pv, scheme)

   ! ---- Stage 9 -------------------------------------------------
   ! c9 = 0.06757
   work%tmp_pv = phi_pv - dt*( &
       ( 0.071806640625_rk    )*work%k1_pv + &
       ( 0.3294380283228177_rk)*work%k6_pv + &
       ( 0.1165190029271823_rk)*work%k7_pv + &
       (-0.034013671875_rk    )*work%k8_pv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k9_pv, scheme)

   ! ---- Stage 10 ------------------------------------------------
   ! c10 = 0.25
   work%tmp_pv = phi_pv - dt*( &
       ( 0.04836757646340646_rk )*work%k1_pv  + &
       ( 0.03928989925676164_rk )*work%k6_pv  + &
       ( 0.10547409458903446_rk )*work%k7_pv  + &
       (-0.021438652846483126_rk)*work%k8_pv  + &
       (-0.10412291746271944_rk )*work%k9_pv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k10_pv, scheme)

   ! ---- Stage 11 ------------------------------------------------
   ! c11 ~ 0.659065
   work%tmp_pv = phi_pv - dt*( &
       (-0.026645614872014785_rk)*work%k1_pv  + &
       ( 0.03333333333333333_rk )*work%k6_pv  + &
       (-0.1631072244872467_rk  )*work%k7_pv  + &
       ( 0.03396081684127761_rk )*work%k8_pv  + &
       ( 0.1572319413814626_rk  )*work%k9_pv  + &
       ( 0.21522674780318796_rk )*work%k10_pv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k11_pv, scheme)

   ! ---- Stage 12 ------------------------------------------------
   ! c12 = 0.8206
   work%tmp_pv = phi_pv - dt*( &
       ( 0.03689009248708622_rk )*work%k1_pv  + &
       (-0.1465181576725543_rk  )*work%k6_pv  + &
       ( 0.2242577768172024_rk  )*work%k7_pv  + &
       ( 0.02294405717066073_rk )*work%k8_pv  + &
       (-0.0035850052905728597_rk)*work%k9_pv + &
       ( 0.08669223316444385_rk )*work%k10_pv + &
       ( 0.43838406519683376_rk )*work%k11_pv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k12_pv, scheme)

   ! ---- Stage 13 ------------------------------------------------
   ! c13 = 0.9012
   work%tmp_pv = phi_pv - dt*( &
       (-0.4866012215113341_rk)*work%k1_pv  + &
       (-6.304602650282853_rk )*work%k6_pv  + &
       (-0.2812456182894729_rk)*work%k7_pv  + &
       (-2.679019236219849_rk )*work%k8_pv  + &
       ( 0.5188156639241577_rk)*work%k9_pv  + &
       ( 1.3653531876033418_rk)*work%k10_pv + &
       ( 5.8850910885039465_rk)*work%k11_pv + &
       ( 2.8028087862720628_rk)*work%k12_pv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng) 
   call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k13_pv, scheme)

   ! ---- Stage 14 ------------------------------------------------
   ! c14 ~ c13  (correction stage)
   work%tmp_pv = phi_pv - dt*( &
       ( 0.4185367457753472_rk)*work%k1_pv  + &
       ( 6.724547581906459_rk )*work%k6_pv  + &
       (-0.42544428016461133_rk)*work%k7_pv + &
       ( 3.3432791530012653_rk)*work%k8_pv  + &
       ( 0.6170816631175374_rk)*work%k9_pv  + &
       (-0.9299661239399329_rk)*work%k10_pv + &
       (-6.099948804751011_rk )*work%k11_pv + &
       (-3.002206187889399_rk )*work%k12_pv + &
       ( 0.2553202529443446_rk)*work%k13_pv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k14_pv, scheme)

   ! ---- Stage 15 ------------------------------------------------
   ! c15 ~ c10
   work%tmp_pv = phi_pv - dt*( &
        (-0.7793740861228848_rk )*work%k1_pv  + &
        (-13.937342538107776_rk )*work%k6_pv  + &
        ( 1.2520488533793563_rk )*work%k7_pv  + &
        (-14.691500408016868_rk )*work%k8_pv  + &
        (-0.494705058533141_rk  )*work%k9_pv  + &
        ( 2.2429749091462368_rk )*work%k10_pv + &
        ( 13.367893803828643_rk )*work%k11_pv + &
        ( 14.396650486650687_rk )*work%k12_pv + &
        (-0.79758133317768_rk   )*work%k13_pv + &
        ( 0.4409353709534278_rk )*work%k14_pv)

   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k15_pv, scheme)

   ! ---- Stage 16 ------------------------------------------------
   ! c16 = 1  (embedded stage; computed but not used in final update)
   !call get_spatial_derivative_inplace(work%tmp_pv, u, dx, n, ng, work%k16_pv, scheme)


   work%tmp_pv = phi_pv - dt*( &
        ( 2.0580513374668867_rk )*work%k1_pv  + &
        ( 22.357937727968032_rk )*work%k6_pv  + &
        ( 0.9094981099755646_rk )*work%k7_pv  + &
        ( 35.89110098240264_rk  )*work%k8_pv  + &
        (-3.442515027624454_rk  )*work%k9_pv  + &
        (-4.865481358036369_rk  )*work%k10_pv + &
        (-18.909803813543427_rk )*work%k11_pv + &
        (-34.26354448030452_rk  )*work%k12_pv + &
        ( 1.2647565216956427_rk )*work%k13_pv)
 
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
 


   ! ---- Final 9th-order assembly --------------------------------
   phi_pv = phi_pv - dt*( &
       ( 0.014611976858423152_rk)*work%k1_pv  + &
       (-0.3915211862331339_rk  )*work%k8_pv  + &
       ( 0.23109325002895065_rk )*work%k9_pv  + &
       ( 0.12747667699928525_rk )*work%k10_pv + &
       ( 0.2246434176204158_rk  )*work%k11_pv + &
       ( 0.5684352689748513_rk  )*work%k12_pv + &
       ( 0.058258715572158275_rk)*work%k13_pv + &
       ( 0.13643174034822156_rk )*work%k14_pv + &
       ( 0.030570139830827976_rk)*work%k15_pv)

   call update_periodic_ghosts_pv(phi_pv, n, ng) 
 

      case ('implicit_euler')
         call build_implicit_matrices(u, dx, dt, n, scheme, method, work%mat_a, work%mat_b)
         work%rhs_inner = phi_pv(ng+1:n+ng)
         call solve_dense(work%mat_a, work%rhs_inner, n)
         do i = 1, n
            phi_pv(ng+i) = work%rhs_inner(i)
         end do
         call update_periodic_ghosts_pv(phi_pv, n, ng)

      case ('cn')
         call build_implicit_matrices(u, dx, dt, n, scheme, method, work%mat_a, work%mat_b)
         work%rhs_inner = matmul(work%mat_b, phi_pv(ng+1:n+ng))
         call solve_dense(work%mat_a, work%rhs_inner, n)
         do i = 1, n
            phi_pv(ng+i) = work%rhs_inner(i)
         end do
         call update_periodic_ghosts_pv(phi_pv, n, ng)

      case default
         error stop 'Unknown integration method.'
      end select
   end subroutine advance_one_step_inplace

   subroutine build_implicit_matrices(u, dx, dt, n, scheme, method, a, b)
      real(rk), intent(in) :: u, dx, dt
      integer, intent(in) :: n
      character(len=*), intent(in) :: scheme, method
      real(rk), intent(out) :: a(n,n), b(n,n)
      real(rk), allocatable :: m(:,:)
      integer :: i, jp1, jm1, jp2, jm2
      real(rk) :: c

      allocate(m(n,n))
      m = 0.0_rk
      a = 0.0_rk
      b = 0.0_rk

      do i = 1, n
         jp1 = periodic_index(i+1, n)
         jm1 = periodic_index(i-1, n)
         jp2 = periodic_index(i+2, n)
         jm2 = periodic_index(i-2, n)

         select case (trim(scheme))
         case ('central4')
            c = u / (12.0_rk * dx)
            m(i,jp1) = m(i,jp1) + 8.0_rk*c
            m(i,jm1) = m(i,jm1) - 8.0_rk*c
            m(i,jp2) = m(i,jp2) - c
            m(i,jm2) = m(i,jm2) + c

         case ('upwind1')
            c = u / dx
            m(i,i)   = m(i,i)   + c
            m(i,jm1) = m(i,jm1) - c

         case ('upwind3')
            c = u / (6.0_rk * dx)
            m(i,jp1) = m(i,jp1) + 2.0_rk*c
            m(i,i)   = m(i,i)   + 3.0_rk*c
            m(i,jm1) = m(i,jm1) - 6.0_rk*c
            m(i,jm2) = m(i,jm2) + c

         case default
            error stop 'Implicit methods currently support central4, upwind1, and upwind3.'
         end select
      end do

      call identity_matrix(a, n)
      call identity_matrix(b, n)

      if (trim(method) == 'implicit_euler') then
         a = a + dt * m
      else if (trim(method) == 'cn') then
         a = a + 0.5_rk * dt * m
         b = b - 0.5_rk * dt * m
      else
         error stop 'Method does not require implicit matrices.'
      end if

      deallocate(m)
   end subroutine build_implicit_matrices

   integer function periodic_index(i, n) result(ip)
      integer, intent(in) :: i, n
      ip = modulo(i - 1, n) + 1
   end function periodic_index

   subroutine identity_matrix(mat, n)
      real(rk), intent(out) :: mat(n,n)
      integer, intent(in) :: n
      integer :: i
      mat = 0.0_rk
      do i = 1, n
         mat(i,i) = 1.0_rk
      end do
   end subroutine identity_matrix

   subroutine solve_dense(a_in, b, n)
      ! Simple Gaussian elimination with partial pivoting.
      real(rk), intent(in) :: a_in(n,n)
      real(rk), intent(inout) :: b(n)
      integer, intent(in) :: n
      real(rk), allocatable :: a(:,:), rowtmp(:)
      real(rk) :: factor, tmp, pivot_abs
      integer :: i, j, k, pivot

      allocate(a(n,n), rowtmp(n))
      a = a_in

      do k = 1, n - 1
         pivot = k
         pivot_abs = abs(a(k,k))
         do i = k + 1, n
            if (abs(a(i,k)) > pivot_abs) then
               pivot = i
               pivot_abs = abs(a(i,k))
            end if
         end do

         if (pivot_abs <= epsilon(1.0_rk)) error stop 'Singular matrix in solve_dense.'

         if (pivot /= k) then
            rowtmp = a(k,:)
            a(k,:) = a(pivot,:)
            a(pivot,:) = rowtmp
            tmp = b(k)
            b(k) = b(pivot)
            b(pivot) = tmp
         end if

         do i = k + 1, n
            factor = a(i,k) / a(k,k)
            a(i,k:n) = a(i,k:n) - factor * a(k,k:n)
            b(i) = b(i) - factor * b(k)
         end do
      end do

      if (abs(a(n,n)) <= epsilon(1.0_rk)) error stop 'Singular matrix in solve_dense.'

      do i = n, 1, -1
         if (i < n) b(i) = b(i) - sum(a(i,i+1:n) * b(i+1:n))
         b(i) = b(i) / a(i,i)
      end do

      deallocate(a, rowtmp)
   end subroutine solve_dense

end module FDM_module
