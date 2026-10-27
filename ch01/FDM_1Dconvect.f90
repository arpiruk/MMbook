program FDM_1Dconvect
   use, intrinsic :: iso_fortran_env, only: real64
   use FDM_module, only: rk, work_arrays, allocate_work, advance_one_step_inplace
   implicit none

   integer :: ng, ncell, n_steps, npx, i, unit, ios
   real(rk) :: l, dx, u, t_final, cfl, cfl0, l0, l1, dt, rho, c0, c1
   real(rk) :: start_time, end_time, solve_time, pv_linf, l2_ratio
   real(rk), parameter :: pi = acos(-1.0_rk)
   character(len=16) :: scase, method, scheme
   character(len=256) :: filename
   real(rk), allocatable :: x(:), phi_pv(:), phi_pv_initial(:)
   type(work_arrays) :: work


   !scase = 'erf'
   !method = 'rk4'
   !scheme = 'central4'
   call read_setup('setup.txt', ng, ncell, cfl, u, t_final, scase, method, scheme)    


   write(6,*) ng
   write(6,*) ncell
   write(6,*) cfl
   write(6,*) u
   write(6,*) t_final
   write(6,'(A)') scase
   write(6,'(A)') method
   write(6,'(A)') scheme
 
   select case (trim(scase))
   case ('sine')
      l = 2.0_rk * pi
      dx = l / real(ncell, rk)
      l0 = 0.0_rk - real(ng, rk) * dx
      l1 = l + real(ng, rk) * dx
   case ('erf')
      l = 40.0_rk
      dx = l / real(ncell, rk)
      l0 = -20.0_rk - real(ng, rk) * dx
      l1 = l0 + l + 2.0_rk * real(ng, rk) * dx
      rho = 10.0_rk
      c0 = 2.5_rk
      c1 = 2.5_rk
   case default
      error stop 'Unknown Scase.'
   end select
    
   dt = cfl * (dx / u)
   n_steps = int(t_final / dt)
   write(6,*) "n_steps " , n_steps
   dt = t_final / real(n_steps, rk)
   cfl0 = cfl
   cfl = u * dt / dx
   npx = ncell + 2 * ng + 1

   allocate(x(npx), phi_pv(npx), phi_pv_initial(npx))
   call allocate_work(work, npx, ncell)

   do i = 1, npx
      x(i) = l0 + real(i-1, rk) * (l1 - l0) / real(npx - 1, rk)
   end do

   select case (trim(scase))
   case ('sine')
      do i = 1, npx
         phi_pv(i) = sin(x(i))
      end do
   case ('erf')
      do i = 1, npx
         phi_pv(i) = 0.5_rk * (erf(rho * (x(i) + c1)) - erf(rho * (x(i) - c0)))
      end do
   end select

   phi_pv_initial = phi_pv

   print *, 'Starting simulation... w nTimeINT = ', n_steps
   write(*,'(A,A,A,A,A,I0,A,I0,A,F8.4)') 'Scheme=', trim(uppercase(scheme)), ', Method=', &
      trim(uppercase(method)), ', Ncell=', ncell, ', steps=', n_steps, ', CFL=', cfl

   call cpu_time(start_time)
   do i = 1, n_steps
      call advance_one_step_inplace(phi_pv, u, dx, dt, ncell, ng, work, method, scheme)
   end do
   call cpu_time(end_time)
   solve_time = end_time - start_time

   pv_linf = maxval(abs(phi_pv(ng+1:ncell+ng+1) - phi_pv_initial(ng+1:ncell+ng+1)))
   l2_ratio = sqrt(sum(phi_pv(ng+1:ncell+ng)**2)) / sqrt(sum(phi_pv_initial(ng+1:ncell+ng)**2))

   write(*,'(A,ES12.4,A)')'#L-inf Error at Time=', t_final, ':'
   write(*,'(A,ES12.4)')  '#  -> PV Error: ', pv_linf
   write(*,'(A,2ES24.16)') '# L2 Ratio       : ', l2_ratio, 1.0_rk-l2_ratio
   write(*,'(A,F10.4,A)') '#Solve time: ', solve_time, ' seconds'

   write(filename,'(A,A,A,A,A,F7.5,A,I0,A)') 'advection_FDM_', trim(scheme), '_', trim(method), '_CFL_', cfl0, '_',ncell, 'cells.dat'
   open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
   if (ios /= 0) error stop 'Could not open output file.'

   write(unit,'(A,A,A,A,A,I0)') '# Simulation Profile: FDM ', trim(uppercase(scheme)), &
      ', Method: ', trim(uppercase(method)), ', Ncell: ', ncell
   write(unit,'(A,F10.4)')   '# Solve Time (s) : ', solve_time
   write(unit,'(A,ES12.4)')  '# PV L-inf Error : ', pv_linf
   write(unit,'(A,2ES24.16)') '# L2 Ratio       : ', l2_ratio, 1.0_rk-l2_ratio
   write(unit,'(A)')             '#             ncell,        solve_time,       pv_linf,     1.0-l2_ratio'   
   write(unit,'(A,I8,3ES19.10)') '# REPORT: ',  ncell,  solve_time, pv_linf  ,1.0-l2_ratio  
   write(unit,'(A)') '#x     phi_num      phi_exact'

   do i = ng + 1, ncell + ng + 1
      write(unit,'(3ES18.8)') x(i), phi_pv(i), phi_pv_initial(i)
   end do
   close(unit)
   print *, 'Data successfully saved to ', trim(filename)

contains
 subroutine read_setup(fname, ng, ncell, cfl, u, t_final, scase, method, scheme)

   implicit none

   character(len=*), intent(in)  :: fname
   integer,          intent(out) :: ng   
   integer,          intent(out) :: ncell
   real(rk),         intent(out) :: cfl
   real(rk),         intent(out) :: u
   real(rk),         intent(out) :: t_final
   character(len=*), intent(out) :: scase
   character(len=*), intent(out) :: method
   character(len=*), intent(out) :: scheme

   integer :: iu

   open(newunit=iu, file=fname, status='old', action='read')

   read(iu,*) ng
   read(iu,*) ncell
   read(iu,*) cfl
   read(iu,*) u
   read(iu,*) t_final
   read(iu,'(A)') scase
   read(iu,'(A)') method
   read(iu,'(A)') scheme

   close(iu)

   scase  = adjustl(trim(scase))
   method = adjustl(trim(method))
   scheme = adjustl(trim(scheme))

end subroutine read_setup 

   pure function uppercase(s) result(out)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: out
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

end program FDM_1Dconvect
