program MMM_1Dconvect

   use MMM_module, only: rk,    get_erf_cellvalue
   use MMM_module_TimeINT, only:  allocate_work, work_arrays , advance_one_step_inplace                    
                         

   implicit none

   !====================================================================
   ! 1. Variable declarations
   !====================================================================

   integer :: ng, ncell, n_steps, npx
   integer :: i, unit, ios

   real(rk), parameter :: pi = acos(-1.0_rk)

   ! Domain and simulation parameters
   real(rk) :: l, l0, l1, dx, u
   real(rk) :: dt, cfl, cfl0, t_final
   real(rk) :: rho, c0, c1

   ! Performance and error measurements
   real(rk) :: start_time, end_time, solve_time
   real(rk) :: pv_linf, cv_linf
   real(rk) :: l2_ratioPV, l2_ratioCV

   ! Simulation configuration
   character(len=16)  :: scase, method, scheme
   character(len=256) :: filename

   ! Solution arrays
   real(rk), allocatable :: x(:), xc(:)
   real(rk), allocatable :: phi_pv(:), phi_cv(:)
   real(rk), allocatable :: phi_pv_initial(:), phi_cv_initial(:)

   type(work_arrays) :: work

   !====================================================================
   ! 2. Read simulation setup
   !====================================================================

   call read_setup('setup.txt', ng, ncell, cfl, u, t_final, &
                   scase, method, scheme)

   ! Multi-moment method uses half the nominal number of cells.
   ncell = ncell / 2

   !====================================================================
   ! 3. Define computational domain and initial conditions
   !====================================================================

   select case (trim(scase))

   case ('sine')

      l       = 6.0_rk * pi
      dx      = l / real(ncell, rk)
      u       = 1.0_rk
      t_final = 60.0_rk * pi

      l0 = -real(ng, rk) * dx
      l1 = l + real(ng, rk) * dx

   case ('erf')

      l  = 40.0_rk
      dx = l / real(ncell, rk)

      l0 = -20.0_rk - real(ng, rk) * dx
      l1 = l0 + l + 2.0_rk * real(ng, rk) * dx

      rho = 10.0_rk
      c0  = 2.5_rk
      c1  = 2.5_rk

   case default

      error stop 'Unknown simulation case.'

   end select

   !====================================================================
   ! 4. Determine time-step size
   !====================================================================

   dt      = cfl * dx / u
   n_steps = int(t_final / dt)

   ! Adjust dt so that the final time is reached exactly.
   dt   = t_final / real(n_steps, rk)
   cfl0 = cfl
   cfl  = u * dt / dx

   ! Total number of grid points, including ghost cells.
   npx = ncell + 2 * ng + 1

   !====================================================================
   ! 5. Allocate arrays and construct computational grid
   !====================================================================

   allocate(x(npx), xc(npx))
   allocate(phi_pv(npx), phi_cv(npx - 1))
   allocate(phi_pv_initial(npx), phi_cv_initial(npx - 1))

   call allocate_work(work, npx, method)

   do i = 1, npx
      x(i)  = l0 + real(i - 1, rk) * (l1 - l0) / real(npx - 1, rk)
      xc(i) = x(i) + 0.5_rk * dx
   end do

   !====================================================================
   ! 6. Initialize point values and cell averages
   !====================================================================

   select case (trim(scase))

   case ('sine')

      do i = 1, npx
         phi_pv(i) = sin(x(i))
      end do

      do i = 1, npx - 1
         phi_cv(i) = (cos(x(i)) - cos(x(i) + dx)) / dx
      end do

   case ('erf')

      do i = 1, npx
         phi_pv(i) = 0.5_rk * (erf(rho * (x(i) + c1)) - &
                             erf(rho * (x(i) - c0)))
      end do

      do i = 1, npx - 1
         phi_cv(i) = get_erf_cellvalue(x(i), x(i + 1), rho, c0, c1)
      end do

   end select

   ! Preserve the initial solution for error evaluation.
   phi_pv_initial = phi_pv
   phi_cv_initial = phi_cv

   !====================================================================
   ! 7. Advance solution in time
   !====================================================================

   print *, 'Starting simulation...'

   write(*, '(A,A,A,A,A,I0,A,I0,A,F8.4)') &
      'Scheme=MMM ', trim(uppercase(scheme)), &
      ', Method=', trim(uppercase(method)), &
      ', Ncell=', ncell, &
      ', steps=', n_steps, &
      ', CFL=', cfl

   call cpu_time(start_time)

!************************* MAIN COMPUTING LOOP ************************
   do i = 1, n_steps
      call advance_one_step_inplace( &
         phi_pv, phi_cv, u, dx, dt, ncell, ng, work, method, scheme)
   end do
!**********************************************************************

   call cpu_time(end_time)

   solve_time = end_time - start_time

   !====================================================================
   ! 8. Compute error and solution norms
   !====================================================================

   ! Point-value maximum error
   pv_linf = maxval(abs( &
      phi_pv(ng + 1:ncell + ng + 1) - &
      phi_pv_initial(ng + 1:ncell + ng + 1)))

   ! Cell-average maximum error
   cv_linf = maxval(abs( &
      phi_cv(ng + 1:ncell + ng) - &
      phi_cv_initial(ng + 1:ncell + ng)))

   ! Point-value L2 norm ratio
   l2_ratioPV = sqrt(sum( &
      phi_pv(ng + 1:ncell + ng + 1)**2)) / &
      sqrt(sum(phi_pv_initial(ng + 1:ncell + ng + 1)**2))

   ! Cell-average L2 norm ratio
   l2_ratioCV = sqrt(sum( &
      phi_cv(ng + 1:ncell + ng)**2)) / &
      sqrt(sum(phi_cv_initial(ng + 1:ncell + ng)**2))

   !====================================================================
   ! 9. Display simulation results
   !====================================================================

   write(*, '(A,ES12.4,A)') &
      'L-inf Error at Time=', t_final, ':'

   write(*, '(A,2ES12.4)') &
      '# PV L-inf Error, L2-ERR : ', pv_linf, 1.0_rk - l2_ratioPV

   write(*, '(A,2ES12.4)') &
      '# CV L-inf Error, L2-ERR : ', cv_linf, 1.0_rk - l2_ratioCV

   write(*, '(A,2ES24.16)') &
      '# L2 RatioPV-CV       : ', l2_ratioPV, l2_ratioCV

   write(*, '(A,F10.4,A)') &
      '# Solve time: ', solve_time, ' seconds'

   !====================================================================
   ! 10. Save results to file
   !====================================================================

   write(filename, '(A,A,A,A,A,F7.5,A,I0,A)') &
      'advection_MMM_', trim(scheme), '_', trim(method), &
      '_CFL_', cfl0, '_', ncell, 'cells.dat'

   open(newunit=unit, file=trim(filename), &
        status='replace', action='write', iostat=ios)

   if (ios /= 0) error stop 'Could not open output file.'

   ! File header
   write(unit, '(A,A,A,A,A,I0)') &
      '# Simulation Profile: MMM ', trim(uppercase(scheme)), &
      ', Method: ', trim(uppercase(method)), &
      ', Ncell: ', ncell

   write(unit, '(A,F10.4)') &
      '# Solve Time (s) : ', solve_time

   write(unit, '(A,2ES12.4)') &
      '# PV L-inf Error, L2-ERR : ', pv_linf, 1.0_rk - l2_ratioPV

   write(unit, '(A,2ES12.4)') &
      '# CV L-inf Error, L2-ERR : ', cv_linf, 1.0_rk - l2_ratioCV

   write(unit, '(A,2ES24.16)') &
      '# L2 RatioPV-CV       : ', l2_ratioPV, l2_ratioCV

   ! Summary record
   write(unit, '(A)') &
      '#      ncell,        solve_time,        pv_linf,          cv_linf,   1.0-l2_ratioPV,   1.0-l2_ratioCV'

   write(unit, '(A,I8,5ES19.10)') &
      '# REPORT: ', ncell, solve_time, pv_linf, cv_linf, &
      1.0_rk - l2_ratioPV, 1.0_rk - l2_ratioCV

   ! Solution profile
   write(unit, '(A)') &
      '#x phi_pv phi_pv_exact xc phi_cv phi_cv_exact'

   do i = ng + 1, ncell + ng + 1

      if (i <= ncell + ng) then

         write(unit, '(6ES18.8)') &
            x(i), phi_pv(i), phi_pv_initial(i), &
            xc(i), phi_cv(i), phi_cv_initial(i)

      else

         write(unit, '(3ES18.8,3A18)') &
            x(i), phi_pv(i), phi_pv_initial(i), &
            'NaN', 'NaN', 'NaN'

      end if

   end do

   close(unit)

   print *, 'Data successfully saved to ', trim(filename)

contains

   !====================================================================
   ! Convert a character string to uppercase
   !====================================================================

   pure function uppercase(s) result(out)

      character(len=*), intent(in) :: s
      character(len=len(s))     :: out

      integer :: j, ia

      do j = 1, len(s)

         ia = iachar(s(j:j))

         if (ia >= iachar('a') .and. ia <= iachar('z')) then
            out(j:j) = achar(ia - 32)
         else
            out(j:j) = s(j:j)
         end if

      end do

   end function uppercase

   !====================================================================
   ! Read simulation configuration from file
   !====================================================================

   subroutine read_setup(fname, ng, ncell, cfl, u, t_final, &
                         scase, method, scheme)

      implicit none

      character(len=*), intent(in)  :: fname

      integer,          intent(out) :: ng, ncell
      real(rk),         intent(out) :: cfl, u, t_final

      character(len=*), intent(out) :: scase, method, scheme

      integer :: iu

      open(newunit=iu, file=fname, status='old', action='read')

      read(iu, *) ng
      read(iu, *) ncell
      read(iu, *) cfl
      read(iu, *) u
      read(iu, *) t_final

      read(iu, '(A)') scase
      read(iu, '(A)') method
      read(iu, '(A)') scheme

      close(iu)

      scase  = adjustl(trim(scase))
      method = adjustl(trim(method))
      scheme = adjustl(trim(scheme))

      print *, 'ng, ncell, cfl : ', ng, ncell, cfl

   end subroutine read_setup

end program MMM_1Dconvect
