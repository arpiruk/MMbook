module MMM_module_TimeINT
   use MMM_module, only: rk, get_spatial_derivative_inplace, &
                         update_periodic_ghosts_pv, update_periodic_ghosts_cv
   implicit none
   private

   public :: work_arrays, allocate_work, advance_one_step_inplace

type :: work_arrays
      real(rk), allocatable :: tmp_pv(:), tmp_cv(:)
      real(rk), allocatable :: k1_pv(:),  k1_cv(:)
      real(rk), allocatable :: k2_pv(:),  k2_cv(:)
      real(rk), allocatable :: k3_pv(:),  k3_cv(:)
      real(rk), allocatable :: k4_pv(:),  k4_cv(:)
      real(rk), allocatable :: k5_pv(:),  k5_cv(:) ! Added for RK5/7
      real(rk), allocatable :: k6_pv(:),  k6_cv(:) ! Added for RK5/7
      real(rk), allocatable :: k7_pv(:),  k7_cv(:) ! Added for RK7
      real(rk), allocatable :: k8_pv(:),  k8_cv(:) ! Added for RK7
! New registers for RK8
      real(rk), allocatable :: k9_pv(:),  k9_cv(:) 
      real(rk), allocatable :: k10_pv(:), k10_cv(:) 
      real(rk), allocatable :: k11_pv(:), k11_cv(:)      
! New registers for RK9      
     real(rk), allocatable, dimension(:) :: k12_pv, k12_cv, k13_pv, k13_cv, k14_pv, k14_cv,&
                                            k15_pv, k15_cv, k16_pv, k16_cv  , k17_pv, k17_cv        
end type work_arrays

contains

   subroutine allocate_work(work, npx, method)
      type(work_arrays), intent(inout) :: work
      integer, intent(in) :: npx
      character(len=*), intent(in) :: method
      integer :: nstages

      select case (trim(adjustl(method)))
      case ('euler')
         nstages = 1
      case ('rk4')
         nstages = 4
      case ('rk5')
         nstages = 6
      case ('rk7')
         nstages = 8  ! Stage 9 overwrites k3
      case ('rk8')
         nstages = 11
      case ('rk9')
         nstages = 15 ! Stage 16 derivative is commented out
      case default
         error stop 'Unknown time integration method.'
      end select

      ! Allocate only stage registers used by the selected integrator.
      ! No initialization is necessary: the spatial operator overwrites its RHS arrays.
      if (nstages >= 1) allocate(work%k1_pv(npx), work%k1_cv(npx - 1))
      if (nstages >= 2) allocate(work%k2_pv(npx), work%k2_cv(npx - 1))
      if (nstages >= 3) allocate(work%k3_pv(npx), work%k3_cv(npx - 1))
      if (nstages >= 4) allocate(work%k4_pv(npx), work%k4_cv(npx - 1))
      if (nstages >= 5) allocate(work%k5_pv(npx), work%k5_cv(npx - 1))
      if (nstages >= 6) allocate(work%k6_pv(npx), work%k6_cv(npx - 1))
      if (nstages >= 7) allocate(work%k7_pv(npx), work%k7_cv(npx - 1))
      if (nstages >= 8) allocate(work%k8_pv(npx), work%k8_cv(npx - 1))
      if (nstages >= 9) allocate(work%k9_pv(npx), work%k9_cv(npx - 1))
      if (nstages >= 10) allocate(work%k10_pv(npx), work%k10_cv(npx - 1))
      if (nstages >= 11) allocate(work%k11_pv(npx), work%k11_cv(npx - 1))
      if (nstages >= 12) allocate(work%k12_pv(npx), work%k12_cv(npx - 1))
      if (nstages >= 13) allocate(work%k13_pv(npx), work%k13_cv(npx - 1))
      if (nstages >= 14) allocate(work%k14_pv(npx), work%k14_cv(npx - 1))
      if (nstages >= 15) allocate(work%k15_pv(npx), work%k15_cv(npx - 1))

      if (nstages > 1) then
         allocate(work%tmp_pv(npx), work%tmp_cv(npx - 1))
      end if
   end subroutine allocate_work

   subroutine advance_one_step_inplace(phi_pv, phi_cv, u, dx, dt, n, ng, work, method, scheme)
      real(rk), intent(inout) :: phi_pv(:), phi_cv(:)
      real(rk), intent(in) :: u, dx, dt
      integer, intent(in) :: n, ng
      type(work_arrays), intent(inout) :: work
      character(len=*), intent(in) :: method, scheme
      real(rk) :: factor,s

      select case (trim(adjustl(method)))
      case ('euler')
         call get_spatial_derivative_inplace(phi_pv, phi_cv, u, dx, n, ng, work%k1_pv, work%k1_cv, scheme)
         phi_pv = phi_pv - dt * work%k1_pv
         phi_cv = phi_cv - dt * work%k1_cv
         call update_periodic_ghosts_pv(phi_pv, n, ng)
         call update_periodic_ghosts_cv(phi_cv, n, ng)

      case ('rk4')
         call get_spatial_derivative_inplace(phi_pv, phi_cv, u, dx, n, ng, work%k1_pv, work%k1_cv, scheme)

         work%tmp_pv = phi_pv - 0.5_rk * dt * work%k1_pv
         work%tmp_cv = phi_cv - 0.5_rk * dt * work%k1_cv
         call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k2_pv, work%k2_cv, scheme)

         work%tmp_pv = phi_pv - 0.5_rk * dt * work%k2_pv
         work%tmp_cv = phi_cv - 0.5_rk * dt * work%k2_cv
         call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k3_pv, work%k3_cv, scheme)

         work%tmp_pv = phi_pv - dt * work%k3_pv
         work%tmp_cv = phi_cv - dt * work%k3_cv
         call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k4_pv, work%k4_cv, scheme)

         factor = dt / 6.0_rk
         phi_pv = phi_pv - factor * work%k1_pv - 2.0_rk * factor * work%k2_pv &
                         - 2.0_rk * factor * work%k3_pv - factor * work%k4_pv
         phi_cv = phi_cv - factor * work%k1_cv - 2.0_rk * factor * work%k2_cv &
                         - 2.0_rk * factor * work%k3_cv - factor * work%k4_cv

         call update_periodic_ghosts_pv(phi_pv, n, ng)
         call update_periodic_ghosts_cv(phi_cv, n, ng)

case ('rk5')
         ! 6-Stage, 5th-Order Explicit Runge-Kutta Scheme (Butcher)
         ! Stage 1
         call get_spatial_derivative_inplace(phi_pv, phi_cv, u, dx, n, ng, work%k1_pv, work%k1_cv, scheme)

         ! Stage 2
         work%tmp_pv = phi_pv - dt * (1.0_rk / 4.0_rk) * work%k1_pv
         work%tmp_cv = phi_cv - dt * (1.0_rk / 4.0_rk) * work%k1_cv
         call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k2_pv, work%k2_cv, scheme)

         ! Stage 3
         work%tmp_pv = phi_pv - dt * ((1.0_rk / 8.0_rk) * work%k1_pv + (1.0_rk / 8.0_rk) * work%k2_pv)
         work%tmp_cv = phi_cv - dt * ((1.0_rk / 8.0_rk) * work%k1_cv + (1.0_rk / 8.0_rk) * work%k2_cv)
         call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k3_pv, work%k3_cv, scheme)

         ! Stage 4
         work%tmp_pv = phi_pv - dt * (-(1.0_rk / 2.0_rk) * work%k2_pv + 1.0_rk * work%k3_pv)
         work%tmp_cv = phi_cv - dt * (-(1.0_rk / 2.0_rk) * work%k2_cv + 1.0_rk * work%k3_cv)
         call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k4_pv, work%k4_cv, scheme)

         ! Stage 5
         work%tmp_pv = phi_pv - dt * ((3.0_rk / 16.0_rk) * work%k1_pv + (9.0_rk / 16.0_rk) * work%k4_pv)
         work%tmp_cv = phi_cv - dt * ((3.0_rk / 16.0_rk) * work%k1_cv + (9.0_rk / 16.0_rk) * work%k4_cv)
         call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k5_pv, work%k5_cv, scheme)

         ! Stage 6
         work%tmp_pv = phi_pv - dt *  (-(3.0_rk / 7.0_rk) * work%k1_pv + (2.0_rk / 7.0_rk) * work%k2_pv &
                                    + (12.0_rk / 7.0_rk) * work%k3_pv - (12.0_rk / 7.0_rk) * work%k4_pv &
                                    + (8.0_rk / 7.0_rk) * work%k5_pv)
         work%tmp_cv = phi_cv - dt *(-(3.0_rk / 7.0_rk) * work%k1_cv + (2.0_rk / 7.0_rk) * work%k2_cv &
                                    + (12.0_rk / 7.0_rk) * work%k3_cv - (12.0_rk / 7.0_rk) * work%k4_cv &
                                    + (8.0_rk / 7.0_rk) * work%k5_cv)
         call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k6_pv, work%k6_cv, scheme)

         ! Final step combination
         phi_pv = phi_pv - dt * ((7.0_rk / 90.0_rk) * work%k1_pv + (32.0_rk / 90.0_rk) * work%k3_pv &
                               + (12.0_rk / 90.0_rk) * work%k4_pv + (32.0_rk / 90.0_rk) * work%k5_pv &
                               + (7.0_rk / 90.0_rk) * work%k6_pv)
         phi_cv = phi_cv - dt * ((7.0_rk / 90.0_rk) * work%k1_cv + (32.0_rk / 90.0_rk) * work%k3_cv &
                               + (12.0_rk / 90.0_rk) * work%k4_cv + (32.0_rk / 90.0_rk) * work%k5_cv &
                               + (7.0_rk / 90.0_rk) * work%k6_cv)

         call update_periodic_ghosts_pv(phi_pv, n, ng)
         call update_periodic_ghosts_cv(phi_cv, n, ng)


case ('rk7')
         ! 9-Stage, 7th-Order Explicit Runge-Kutta Scheme (Butcher, p. 196)
         
         ! Stage 1
         call get_spatial_derivative_inplace(phi_pv, phi_cv, u, dx, n, ng, work%k1_pv, work%k1_cv, scheme)

         ! Stage 2: c2 = 1/6
         work%tmp_pv = phi_pv - dt * (1.0_rk / 6.0_rk) * work%k1_pv
         work%tmp_cv = phi_cv - dt * (1.0_rk / 6.0_rk) * work%k1_cv
         call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k2_pv, work%k2_cv, scheme)

         ! Stage 3: c3 = 1/3
         work%tmp_pv = phi_pv - dt * (1.0_rk / 3.0_rk) * work%k2_pv
         work%tmp_cv = phi_cv - dt * (1.0_rk / 3.0_rk) * work%k2_cv
         call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k3_pv, work%k3_cv, scheme)

         ! Stage 4: c4 = 1/2
         work%tmp_pv = phi_pv - dt * ((1.0_rk / 8.0_rk) * work%k1_pv + (3.0_rk / 8.0_rk) * work%k3_pv)
         work%tmp_cv = phi_cv - dt * ((1.0_rk / 8.0_rk) * work%k1_cv + (3.0_rk / 8.0_rk) * work%k3_cv)
         call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k4_pv, work%k4_cv, scheme)

         ! Stage 5: c5 = 2/11
         work%tmp_pv = phi_pv - dt * ((148.0_rk / 1331.0_rk) * work%k1_pv &
                                    + (150.0_rk / 1331.0_rk) * work%k3_pv &
                                    - (56.0_rk / 1331.0_rk)  * work%k4_pv)
         work%tmp_cv = phi_cv - dt * ((148.0_rk / 1331.0_rk) * work%k1_cv &
                                    + (150.0_rk / 1331.0_rk) * work%k3_cv &
                                    - (56.0_rk / 1331.0_rk)  * work%k4_cv)
         call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k5_pv, work%k5_cv, scheme)

         ! Stage 6: c6 = 2/3
         work%tmp_pv = phi_pv - dt * (-(404.0_rk / 243.0_rk)   * work%k1_pv &
                                    - (170.0_rk / 27.0_rk)    * work%k3_pv &
                                    + (4024.0_rk / 1701.0_rk) * work%k4_pv &
                                    + (10648.0_rk / 1701.0_rk)* work%k5_pv)
         work%tmp_cv = phi_cv - dt * (-(404.0_rk / 243.0_rk)   * work%k1_cv &
                                    - (170.0_rk / 27.0_rk)    * work%k3_cv &
                                    + (4024.0_rk / 1701.0_rk) * work%k4_cv &
                                    + (10648.0_rk / 1701.0_rk)* work%k5_cv)
         call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k6_pv, work%k6_cv, scheme)

         ! Stage 7: c7 = 6/7
         work%tmp_pv = phi_pv - dt * ((2466.0_rk / 2401.0_rk)    * work%k1_pv &
                                    + (1242.0_rk / 343.0_rk)     * work%k3_pv &
                                    - (19176.0_rk / 16807.0_rk)  * work%k4_pv &
                                    - (51909.0_rk / 16807.0_rk)  * work%k5_pv &
                                    + (1053.0_rk / 2401.0_rk)    * work%k6_pv)
         work%tmp_cv = phi_cv - dt * ((2466.0_rk / 2401.0_rk)    * work%k1_cv &
                                    + (1242.0_rk / 343.0_rk)     * work%k3_cv &
                                    - (19176.0_rk / 16807.0_rk)  * work%k4_cv &
                                    - (51909.0_rk / 16807.0_rk)  * work%k5_cv &
                                    + (1053.0_rk / 2401.0_rk)    * work%k6_cv)
         call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k7_pv, work%k7_cv, scheme)

         ! Stage 8: c8 = 1/7
         work%tmp_pv = phi_pv - dt * ((5.0_rk / 154.0_rk)       * work%k1_pv &
                                    + (96.0_rk / 539.0_rk)      * work%k4_pv &
                                    - (1815.0_rk / 20384.0_rk)  * work%k5_pv &
                                    - (405.0_rk / 2464.0_rk)    * work%k6_pv &
                                    + (49.0_rk / 1144.0_rk)     * work%k7_pv)
         work%tmp_cv = phi_cv - dt * ((5.0_rk / 154.0_rk)       * work%k1_cv &
                                    + (96.0_rk / 539.0_rk)      * work%k4_cv &
                                    - (1815.0_rk / 20384.0_rk)  * work%k5_cv &
                                    - (405.0_rk / 2464.0_rk)    * work%k6_cv &
                                    + (49.0_rk / 1144.0_rk)     * work%k7_cv)
         call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k8_pv, work%k8_cv, scheme)

         ! Stage 9: c9 = 1.0
         work%tmp_pv = phi_pv - dt * (-(113.0_rk / 32.0_rk)    * work%k1_pv &
                                    - (195.0_rk / 22.0_rk)    * work%k3_pv &
                                    + (32.0_rk / 7.0_rk)      * work%k4_pv &
                                    + (29403.0_rk / 3584.0_rk)* work%k5_pv &
                                    - (729.0_rk / 512.0_rk)   * work%k6_pv &
                                    + (1029.0_rk / 1408.0_rk) * work%k7_pv &
                                    + (21.0_rk / 16.0_rk)     * work%k8_pv)
         work%tmp_cv = phi_cv - dt * (-(113.0_rk / 32.0_rk)    * work%k1_cv &
                                    - (195.0_rk / 22.0_rk)    * work%k3_cv &
                                    + (32.0_rk / 7.0_rk)      * work%k4_cv &
                                    + (29403.0_rk / 3584.0_rk)* work%k5_cv &
                                    - (729.0_rk / 512.0_rk)   * work%k6_cv &
                                    + (1029.0_rk / 1408.0_rk) * work%k7_cv &
                                    + (21.0_rk / 16.0_rk)     * work%k8_cv)
         ! Overwriting k3 registers with the 9th stage evaluation
         call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k3_pv, work%k3_cv, scheme)

         ! Final step combination using verified b-weights
         phi_pv = phi_pv - dt * ((32.0_rk / 105.0_rk)               * work%k4_pv &
                               + (1771561.0_rk / 6289920.0_rk)      * work%k5_pv &
                               + (243.0_rk / 2560.0_rk)             * work%k6_pv &
                               + (16807.0_rk / 74880.0_rk)          * work%k7_pv &
                               + (77.0_rk / 1440.0_rk)              * work%k8_pv &
                               + (11.0_rk / 270.0_rk)               * work%k3_pv)
                               
         phi_cv = phi_cv - dt * ((32.0_rk / 105.0_rk)               * work%k4_cv &
                               + (1771561.0_rk / 6289920.0_rk)      * work%k5_cv &
                               + (243.0_rk / 2560.0_rk)             * work%k6_cv &
                               + (16807.0_rk / 74880.0_rk)          * work%k7_cv &
                               + (77.0_rk / 1440.0_rk)              * work%k8_cv &
                               + (11.0_rk / 270.0_rk)               * work%k3_cv)

         call update_periodic_ghosts_pv(phi_pv, n, ng)
         call update_periodic_ghosts_cv(phi_cv, n, ng)
 
 
 
 case ('rk8')
   ! Cooper-Verner 11-stage, 8th-order explicit RK
   ! Assumes: k = spatial derivative/operator, update is phi <- phi - dt*k

 
   s = sqrt(21.0_rk)

   ! Stage 1
   call get_spatial_derivative_inplace(phi_pv, phi_cv, u, dx, n, ng, work%k1_pv, work%k1_cv, scheme)

   ! Stage 2
   work%tmp_pv = phi_pv - dt*((1.0_rk/2.0_rk)*work%k1_pv)
   work%tmp_cv = phi_cv - dt*((1.0_rk/2.0_rk)*work%k1_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k2_pv, work%k2_cv, scheme)

   ! Stage 3
   work%tmp_pv = phi_pv - dt*((1.0_rk/4.0_rk)*work%k1_pv + (1.0_rk/4.0_rk)*work%k2_pv)
   work%tmp_cv = phi_cv - dt*((1.0_rk/4.0_rk)*work%k1_cv + (1.0_rk/4.0_rk)*work%k2_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k3_pv, work%k3_cv, scheme)

   ! Stage 4, c = 1/2 - sqrt(21)/14
   work%tmp_pv = phi_pv - dt*((1.0_rk/7.0_rk)*work%k1_pv                                      &
      + (-1.0_rk/14.0_rk + 3.0_rk*s/98.0_rk)*work%k2_pv                                        &
      + ( 3.0_rk/7.0_rk  - 5.0_rk*s/49.0_rk)*work%k3_pv)
   work%tmp_cv = phi_cv - dt*((1.0_rk/7.0_rk)*work%k1_cv                                      &
      + (-1.0_rk/14.0_rk + 3.0_rk*s/98.0_rk)*work%k2_cv                                        &
      + ( 3.0_rk/7.0_rk  - 5.0_rk*s/49.0_rk)*work%k3_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k4_pv, work%k4_cv, scheme)

   ! Stage 5
   work%tmp_pv = phi_pv - dt*((11.0_rk/84.0_rk - s/84.0_rk)*work%k1_pv                          &
      + (2.0_rk/7.0_rk - 4.0_rk*s/63.0_rk)*work%k3_pv                                           &
      + (1.0_rk/12.0_rk + s/252.0_rk)*work%k4_pv)
   work%tmp_cv = phi_cv - dt*((11.0_rk/84.0_rk - s/84.0_rk)*work%k1_cv                          &
      + (2.0_rk/7.0_rk - 4.0_rk*s/63.0_rk)*work%k3_cv                                           &
      + (1.0_rk/12.0_rk + s/252.0_rk)*work%k4_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k5_pv, work%k5_cv, scheme)

   ! Stage 6
   work%tmp_pv = phi_pv - dt*((5.0_rk/48.0_rk - s/48.0_rk)*work%k1_pv                           &
      + (1.0_rk/4.0_rk - s/36.0_rk)*work%k3_pv                                                   &
      + (-77.0_rk/120.0_rk - 7.0_rk*s/180.0_rk)*work%k4_pv                                      &
      + (63.0_rk/80.0_rk + 7.0_rk*s/80.0_rk)*work%k5_pv)
   work%tmp_cv = phi_cv - dt*((5.0_rk/48.0_rk - s/48.0_rk)*work%k1_cv                           &
      + (1.0_rk/4.0_rk - s/36.0_rk)*work%k3_cv                                                   &
      + (-77.0_rk/120.0_rk - 7.0_rk*s/180.0_rk)*work%k4_cv                                      &
      + (63.0_rk/80.0_rk + 7.0_rk*s/80.0_rk)*work%k5_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k6_pv, work%k6_cv, scheme)

   ! Stage 7
   work%tmp_pv = phi_pv - dt*((5.0_rk/21.0_rk + s/42.0_rk)*work%k1_pv                           &
      + (-48.0_rk/35.0_rk - 92.0_rk*s/315.0_rk)*work%k3_pv                                      &
      + (211.0_rk/30.0_rk + 29.0_rk*s/18.0_rk)*work%k4_pv                                       &
      + (-36.0_rk/5.0_rk - 23.0_rk*s/14.0_rk)*work%k5_pv                                        &
      + (9.0_rk/5.0_rk + 13.0_rk*s/35.0_rk)*work%k6_pv)
   work%tmp_cv = phi_cv - dt*((5.0_rk/21.0_rk + s/42.0_rk)*work%k1_cv                           &
      + (-48.0_rk/35.0_rk - 92.0_rk*s/315.0_rk)*work%k3_cv                                      &
      + (211.0_rk/30.0_rk + 29.0_rk*s/18.0_rk)*work%k4_cv                                       &
      + (-36.0_rk/5.0_rk - 23.0_rk*s/14.0_rk)*work%k5_cv                                        &
      + (9.0_rk/5.0_rk + 13.0_rk*s/35.0_rk)*work%k6_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k7_pv, work%k7_cv, scheme)

   ! Stage 8
   work%tmp_pv = phi_pv - dt*((1.0_rk/14.0_rk)*work%k1_pv                                      &
      + (1.0_rk/9.0_rk + s/42.0_rk)*work%k5_pv                                                   &
      + (13.0_rk/63.0_rk + s/21.0_rk)*work%k6_pv                                                 &
      + (1.0_rk/9.0_rk)*work%k7_pv)
   work%tmp_cv = phi_cv - dt*((1.0_rk/14.0_rk)*work%k1_cv                                      &
      + (1.0_rk/9.0_rk + s/42.0_rk)*work%k5_cv                                                   &
      + (13.0_rk/63.0_rk + s/21.0_rk)*work%k6_cv                                                 &
      + (1.0_rk/9.0_rk)*work%k7_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k8_pv, work%k8_cv, scheme)

   ! Stage 9
   work%tmp_pv = phi_pv - dt*((1.0_rk/32.0_rk)*work%k1_pv                                      &
      + (91.0_rk/576.0_rk + 7.0_rk*s/192.0_rk)*work%k5_pv                                       &
      + (11.0_rk/72.0_rk)*work%k6_pv                                                             &
      + (-385.0_rk/1152.0_rk + 25.0_rk*s/384.0_rk)*work%k7_pv                                   &
      + (63.0_rk/128.0_rk - 13.0_rk*s/128.0_rk)*work%k8_pv)
   work%tmp_cv = phi_cv - dt*((1.0_rk/32.0_rk)*work%k1_cv                                      &
      + (91.0_rk/576.0_rk + 7.0_rk*s/192.0_rk)*work%k5_cv                                       &
      + (11.0_rk/72.0_rk)*work%k6_cv                                                             &
      + (-385.0_rk/1152.0_rk + 25.0_rk*s/384.0_rk)*work%k7_cv                                   &
      + (63.0_rk/128.0_rk - 13.0_rk*s/128.0_rk)*work%k8_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k9_pv, work%k9_cv, scheme)

   ! Stage 10
   work%tmp_pv = phi_pv - dt*((1.0_rk/14.0_rk)*work%k1_pv                                      &
      + (1.0_rk/9.0_rk)*work%k5_pv                                                               &
      + (-733.0_rk/2205.0_rk + s/15.0_rk)*work%k6_pv                                            &
      + (515.0_rk/504.0_rk - 37.0_rk*s/168.0_rk)*work%k7_pv                                     &
      + (-51.0_rk/56.0_rk + 11.0_rk*s/56.0_rk)*work%k8_pv                                       &
      + (132.0_rk/245.0_rk - 4.0_rk*s/35.0_rk)*work%k9_pv)
   work%tmp_cv = phi_cv - dt*((1.0_rk/14.0_rk)*work%k1_cv                                      &
      + (1.0_rk/9.0_rk)*work%k5_cv                                                               &
      + (-733.0_rk/2205.0_rk + s/15.0_rk)*work%k6_cv                                            &
      + (515.0_rk/504.0_rk - 37.0_rk*s/168.0_rk)*work%k7_cv                                     &
      + (-51.0_rk/56.0_rk + 11.0_rk*s/56.0_rk)*work%k8_cv                                       &
      + (132.0_rk/245.0_rk - 4.0_rk*s/35.0_rk)*work%k9_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k10_pv, work%k10_cv, scheme)

   ! Stage 11
   work%tmp_pv = phi_pv - dt*((-7.0_rk/3.0_rk - 7.0_rk*s/18.0_rk)*work%k5_pv                    &
      + (-2.0_rk/5.0_rk - 28.0_rk*s/45.0_rk)*work%k6_pv                                         &
      + (-91.0_rk/24.0_rk + 53.0_rk*s/72.0_rk)*work%k7_pv                                       &
      + (301.0_rk/72.0_rk - 53.0_rk*s/72.0_rk)*work%k8_pv                                       &
      + (28.0_rk/45.0_rk + 28.0_rk*s/45.0_rk)*work%k9_pv                                        &
      + (49.0_rk/18.0_rk + 7.0_rk*s/18.0_rk)*work%k10_pv)
   work%tmp_cv = phi_cv - dt*((-7.0_rk/3.0_rk - 7.0_rk*s/18.0_rk)*work%k5_cv                    &
      + (-2.0_rk/5.0_rk - 28.0_rk*s/45.0_rk)*work%k6_cv                                         &
      + (-91.0_rk/24.0_rk + 53.0_rk*s/72.0_rk)*work%k7_cv                                       &
      + (301.0_rk/72.0_rk - 53.0_rk*s/72.0_rk)*work%k8_cv                                       &
      + (28.0_rk/45.0_rk + 28.0_rk*s/45.0_rk)*work%k9_cv                                        &
      + (49.0_rk/18.0_rk + 7.0_rk*s/18.0_rk)*work%k10_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, work%k11_pv, work%k11_cv, scheme)

   ! Final 8th-order assembly
   phi_pv = phi_pv - dt*((1.0_rk/20.0_rk)*work%k1_pv                                           &
      + (49.0_rk/180.0_rk)*work%k8_pv                                                           &
      + (16.0_rk/45.0_rk)*work%k9_pv                                                            &
      + (49.0_rk/180.0_rk)*work%k10_pv                                                          &
      + (1.0_rk/20.0_rk)*work%k11_pv)

   phi_cv = phi_cv - dt*((1.0_rk/20.0_rk)*work%k1_cv                                           &
      + (49.0_rk/180.0_rk)*work%k8_cv                                                           &
      + (16.0_rk/45.0_rk)*work%k9_cv                                                            &
      + (49.0_rk/180.0_rk)*work%k10_cv                                                          &
      + (1.0_rk/20.0_rk)*work%k11_cv)

   call update_periodic_ghosts_pv(phi_pv, n, ng)
   call update_periodic_ghosts_cv(phi_cv, n, ng)
 
 
 
! ============================================================
! Verner 16-stage, 9th-order explicit Runge-Kutta (Vern9)
! ============================================================
! Reference: J.H. Verner, "Explicit Runge-Kutta methods with
!   estimates of the local truncation error", SIAM J. Numer. Anal.
!   15(4), 772-790, 1978.  Coefficients from Verner's website
!   (most-efficient variant), as implemented in OrdinaryDiffEq.jl.
!
! Notation: stages are k1..k16.  The method uses only k1 and
!   k6..k15 in the final b-vector (k2..k5 appear only in the
!   internal A-matrix; k16 is the FSAL stage for the embedded
!   8th-order solution, not used in the 9th-order update here).
!
! Update sign convention (same as original code):
!   phi <- phi - dt * k    (k already contains the spatial derivative)
!
! Bugs fixed vs. the original submission
! ---------------------------------------
! (1) Stage 13 ghost update: 'n' -> 'ng' in update_periodic_ghosts_pv
! (2) Stage 13 tmp_cv line: work%k12_pv -> work%k12_cv
! (3) Entire b-vector replaced: original weights summed to ~1.178,
!     not 1.  The A-matrix coefficients were also inconsistent
!     (several c_i values exceeded 1).  The full tableau has been
!     replaced with the authoritative Verner coefficients.
!
! Node values c_i (row sums of A):
!   c1 = 0.03462,  c2 = 0.097024,  c3 = 0.145537
!   c4 = 0.561,    c5 = 0.229008,  c6 = 0.544992
!   c7 = 0.645,    c8 = 0.48375,   c9 = 0.06757
!   c10= 0.25,     c11= 0.659065,  c12= 0.8206
!   c13= 0.9012
!   (stages 14-16 share c14=c13, c15=c10, c16=1 for the embedded pair)
! ============================================================

case ('rk9')

   ! ---- Stage 1 -------------------------------------------------
   call get_spatial_derivative_inplace(phi_pv, phi_cv, u, dx, n, ng, &
        work%k1_pv, work%k1_cv, scheme)

   ! ---- Stage 2 -------------------------------------------------
   ! c2 ~ 0.097024
   work%tmp_pv = phi_pv - dt*(0.03462_rk * work%k1_pv)
   work%tmp_cv = phi_cv - dt*(0.03462_rk * work%k1_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, &
        work%k2_pv, work%k2_cv, scheme)

   ! ---- Stage 3 -------------------------------------------------
   ! c3 ~ 0.145537
   work%tmp_pv = phi_pv - dt*( &
        (-0.03893354388572875_rk)*work%k1_pv + &
        ( 0.13595789452450918_rk)*work%k2_pv)
   work%tmp_cv = phi_cv - dt*( &
        (-0.03893354388572875_rk)*work%k1_cv + &
        ( 0.13595789452450918_rk)*work%k2_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, &
        work%k3_pv, work%k3_cv, scheme)

   ! ---- Stage 4 -------------------------------------------------
   ! c4 = 0.561
   work%tmp_pv = phi_pv - dt*( &
        (0.03638413148954267_rk)*work%k1_pv + &
        (0.10915239446862801_rk)*work%k3_pv)
   work%tmp_cv = phi_cv - dt*( &
        (0.03638413148954267_rk)*work%k1_cv + &
        (0.10915239446862801_rk)*work%k3_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, &
        work%k4_pv, work%k4_cv, scheme)

   ! ---- Stage 5 -------------------------------------------------
   ! c5 ~ 0.229008
   work%tmp_pv = phi_pv - dt*( &
        ( 2.0257639143939694_rk)*work%k1_pv + &
        (-7.638023836496291_rk )*work%k3_pv + &
        ( 6.173259922102322_rk )*work%k4_pv)
   work%tmp_cv = phi_cv - dt*( &
        ( 2.0257639143939694_rk)*work%k1_cv + &
        (-7.638023836496291_rk )*work%k3_cv + &
        ( 6.173259922102322_rk )*work%k4_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, &
        work%k5_pv, work%k5_cv, scheme)

   ! ---- Stage 6 -------------------------------------------------
   ! c6 ~ 0.544992
   work%tmp_pv = phi_pv - dt*( &
        (0.05112275589406061_rk )*work%k1_pv + &
        (0.17708237945550218_rk )*work%k4_pv + &
        (0.0008027762409222536_rk)*work%k5_pv)
   work%tmp_cv = phi_cv - dt*( &
        (0.05112275589406061_rk )*work%k1_cv + &
        (0.17708237945550218_rk )*work%k4_cv + &
        (0.0008027762409222536_rk)*work%k5_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, &
        work%k6_pv, work%k6_cv, scheme)

   ! ---- Stage 7 -------------------------------------------------
   ! c7 = 0.645
   work%tmp_pv = phi_pv - dt*( &
        ( 0.13160063579752163_rk)*work%k1_pv + &
        (-0.2957276252669636_rk )*work%k4_pv + &
        ( 0.08781378035642955_rk)*work%k5_pv + &
        ( 0.6213052975225274_rk )*work%k6_pv)
   work%tmp_cv = phi_cv - dt*( &
        ( 0.13160063579752163_rk)*work%k1_cv + &
        (-0.2957276252669636_rk )*work%k4_cv + &
        ( 0.08781378035642955_rk)*work%k5_cv + &
        ( 0.6213052975225274_rk )*work%k6_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, &
        work%k7_pv, work%k7_cv, scheme)

   ! ---- Stage 8 -------------------------------------------------
   ! c8 = 0.48375
   work%tmp_pv = phi_pv - dt*( &
        (0.07166666666666667_rk)*work%k1_pv + &
        (0.33055335789153195_rk)*work%k6_pv + &
        (0.2427799754418014_rk )*work%k7_pv)
   work%tmp_cv = phi_cv - dt*( &
        (0.07166666666666667_rk)*work%k1_cv + &
        (0.33055335789153195_rk)*work%k6_cv + &
        (0.2427799754418014_rk )*work%k7_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, &
        work%k8_pv, work%k8_cv, scheme)

   ! ---- Stage 9 -------------------------------------------------
   ! c9 = 0.06757
   work%tmp_pv = phi_pv - dt*( &
        ( 0.071806640625_rk    )*work%k1_pv + &
        ( 0.3294380283228177_rk)*work%k6_pv + &
        ( 0.1165190029271823_rk)*work%k7_pv + &
        (-0.034013671875_rk    )*work%k8_pv)
   work%tmp_cv = phi_cv - dt*( &
        ( 0.071806640625_rk    )*work%k1_cv + &
        ( 0.3294380283228177_rk)*work%k6_cv + &
        ( 0.1165190029271823_rk)*work%k7_cv + &
        (-0.034013671875_rk    )*work%k8_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, &
        work%k9_pv, work%k9_cv, scheme)

   ! ---- Stage 10 ------------------------------------------------
   ! c10 = 0.25
   work%tmp_pv = phi_pv - dt*( &
        ( 0.04836757646340646_rk )*work%k1_pv  + &
        ( 0.03928989925676164_rk )*work%k6_pv  + &
        ( 0.10547409458903446_rk )*work%k7_pv  + &
        (-0.021438652846483126_rk)*work%k8_pv  + &
        (-0.10412291746271944_rk )*work%k9_pv)
   work%tmp_cv = phi_cv - dt*( &
        ( 0.04836757646340646_rk )*work%k1_cv  + &
        ( 0.03928989925676164_rk )*work%k6_cv  + &
        ( 0.10547409458903446_rk )*work%k7_cv  + &
        (-0.021438652846483126_rk)*work%k8_cv  + &
        (-0.10412291746271944_rk )*work%k9_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, &
        work%k10_pv, work%k10_cv, scheme)

   ! ---- Stage 11 ------------------------------------------------
   ! c11 ~ 0.659065
   work%tmp_pv = phi_pv - dt*( &
        (-0.026645614872014785_rk)*work%k1_pv  + &
        ( 0.03333333333333333_rk )*work%k6_pv  + &
        (-0.1631072244872467_rk  )*work%k7_pv  + &
        ( 0.03396081684127761_rk )*work%k8_pv  + &
        ( 0.1572319413814626_rk  )*work%k9_pv  + &
        ( 0.21522674780318796_rk )*work%k10_pv)
   work%tmp_cv = phi_cv - dt*( &
        (-0.026645614872014785_rk)*work%k1_cv  + &
        ( 0.03333333333333333_rk )*work%k6_cv  + &
        (-0.1631072244872467_rk  )*work%k7_cv  + &
        ( 0.03396081684127761_rk )*work%k8_cv  + &
        ( 0.1572319413814626_rk  )*work%k9_cv  + &
        ( 0.21522674780318796_rk )*work%k10_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, &
        work%k11_pv, work%k11_cv, scheme)

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
   work%tmp_cv = phi_cv - dt*( &
        ( 0.03689009248708622_rk )*work%k1_cv  + &
        (-0.1465181576725543_rk  )*work%k6_cv  + &
        ( 0.2242577768172024_rk  )*work%k7_cv  + &
        ( 0.02294405717066073_rk )*work%k8_cv  + &
        (-0.0035850052905728597_rk)*work%k9_cv + &
        ( 0.08669223316444385_rk )*work%k10_cv + &
        ( 0.43838406519683376_rk )*work%k11_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, &
        work%k12_pv, work%k12_cv, scheme)

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
   work%tmp_cv = phi_cv - dt*( &
        (-0.4866012215113341_rk)*work%k1_cv  + &
        (-6.304602650282853_rk )*work%k6_cv  + &
        (-0.2812456182894729_rk)*work%k7_cv  + &
        (-2.679019236219849_rk )*work%k8_cv  + &
        ( 0.5188156639241577_rk)*work%k9_cv  + &
        ( 1.3653531876033418_rk)*work%k10_cv + &
        ( 5.8850910885039465_rk)*work%k11_cv + &
        ( 2.8028087862720628_rk)*work%k12_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)   ! BUG FIX #1: was 'n, n'
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, &
        work%k13_pv, work%k13_cv, scheme)

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
   work%tmp_cv = phi_cv - dt*( &
        ( 0.4185367457753472_rk)*work%k1_cv  + &
        ( 6.724547581906459_rk )*work%k6_cv  + &
        (-0.42544428016461133_rk)*work%k7_cv + &
        ( 3.3432791530012653_rk)*work%k8_cv  + &
        ( 0.6170816631175374_rk)*work%k9_cv  + &
        (-0.9299661239399329_rk)*work%k10_cv + &
        (-6.099948804751011_rk )*work%k11_cv + &
        (-3.002206187889399_rk )*work%k12_cv + &
        ( 0.2553202529443446_rk)*work%k13_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, &
        work%k14_pv, work%k14_cv, scheme)

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
   work%tmp_cv = phi_cv - dt*( &
        (-0.7793740861228848_rk )*work%k1_cv  + &
        (-13.937342538107776_rk )*work%k6_cv  + &
        ( 1.2520488533793563_rk )*work%k7_cv  + &
        (-14.691500408016868_rk )*work%k8_cv  + &
        (-0.494705058533141_rk  )*work%k9_cv  + &
        ( 2.2429749091462368_rk )*work%k10_cv + &
        ( 13.367893803828643_rk )*work%k11_cv + &
        ( 14.396650486650687_rk )*work%k12_cv + &
        (-0.79758133317768_rk   )*work%k13_cv + &
        ( 0.4409353709534278_rk )*work%k14_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
   call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, &
        work%k15_pv, work%k15_cv, scheme)

   ! ---- Stage 16 ------------------------------------------------
   ! c16 = 1  (this is the FSAL / embedded stage; computed but not
   !           used in the 9th-order update below)
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
   work%tmp_cv = phi_cv - dt*( &
        ( 2.0580513374668867_rk )*work%k1_cv  + &
        ( 22.357937727968032_rk )*work%k6_cv  + &
        ( 0.9094981099755646_rk )*work%k7_cv  + &
        ( 35.89110098240264_rk  )*work%k8_cv  + &
        (-3.442515027624454_rk  )*work%k9_cv  + &
        (-4.865481358036369_rk  )*work%k10_cv + &
        (-18.909803813543427_rk )*work%k11_cv + &
        (-34.26354448030452_rk  )*work%k12_cv + &
        ( 1.2647565216956427_rk )*work%k13_cv)
   call update_periodic_ghosts_pv(work%tmp_pv, n, ng)
   call update_periodic_ghosts_cv(work%tmp_cv, n, ng)
  ! call get_spatial_derivative_inplace(work%tmp_pv, work%tmp_cv, u, dx, n, ng, &
  !      work%k16_pv, work%k16_cv, scheme)

   ! ---- Final 9th-order assembly --------------------------------
   ! BUG FIX #3: b-weights replaced.  Original summed to ~1.178.
   ! Correct b uses k1, k8..k15 (sum = 1 exactly).
   ! b1=0.014611976858423152,  b8=-0.3915211862331339,
   ! b9=0.23109325002895065,   b10=0.12747667699928525,
   ! b11=0.2246434176204158,   b12=0.5684352689748513,
   ! b13=0.058258715572158275, b14=0.13643174034822156,
   ! b15=0.030570139830827976
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

   phi_cv = phi_cv - dt*( &
        ( 0.014611976858423152_rk)*work%k1_cv  + &
        (-0.3915211862331339_rk  )*work%k8_cv  + &
        ( 0.23109325002895065_rk )*work%k9_cv  + &
        ( 0.12747667699928525_rk )*work%k10_cv + &
        ( 0.2246434176204158_rk  )*work%k11_cv + &
        ( 0.5684352689748513_rk  )*work%k12_cv + &
        ( 0.058258715572158275_rk)*work%k13_cv + &
        ( 0.13643174034822156_rk )*work%k14_cv + &
        ( 0.030570139830827976_rk)*work%k15_cv)

   call update_periodic_ghosts_pv(phi_pv, n, ng)
   call update_periodic_ghosts_cv(phi_cv, n, ng) 
  
 

      case default
         error stop 'Unknown integration method.'
      end select
   end subroutine advance_one_step_inplace

end module MMM_module_TimeINT

