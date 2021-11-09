module params
  ! Variables to be shared by all program units
  real :: pad_finest
  

  
  integer :: levelmin, levelmax, nlevel, npad
  ! integer :: ix, iy, iz
  ! integer :: nx, ny, nz

  real :: omega_m, omega_v, boxlen, H0, astart
  real :: omega_m_w, omega_v_w, H0_w, astart_w
end module params

program ref_mask

  ! NOTES
  !
  ! off_abs is the absolute offset of the grids at ilevel from the
  ! origin, measured in ilevel grids
  !
  ! off_rel is the relative offset between the grids at ilevel and the
  ! grids at ilevel-1, measured in ilevel-1 grids

  ! TODO
  !
  ! perhaps bulding the refinement hierarchy is necessary all the
  ! way down past levelmin, in case the zoom region isn't consistent
  ! below there?

  use params

  implicit none
  
  type :: masks
     integer, pointer :: data(:, :, :)
  end type masks

  type(masks), allocatable :: rm(:)

  real, allocatable :: slab(:, :)
  real, allocatable :: data_w(:, :, :), data_coarse(:, :, :)

  character(len=800):: nml_file
  character(len=100) :: g9p_path, mus_path
  character(len=4):: lnum
  character(len=6):: ldim
  character(len=3):: xyz = 'xyz'
  
  ! integer, allocatable :: n(:, :)
  ! integer, allocatable :: off_abs(:, :), off_rel(:, :)

  integer :: mv
  real :: dx
  real, dimension(3) :: x = 0.0
  integer, dimension(3) :: i_fine = 0
  
  integer :: ilevel
  integer :: ii, jj, kk, i

  ! Convenience and header variables for writing
  real :: dx_w, tol
  real :: pvar_val
  real, dimension(3) :: off_w, xo_mus
  integer, dimension(3) :: n_w, n_coarse, off_abs
  integer, dimension(3) :: n_mus
  ! g9p header variables
  integer :: n1_g9p
  real :: dx_g9p, x1o_g9p, x2o_g9p, x3o_g9p, dx_mus

  real :: ainv, adot, tmp1, growth, tmp, fupper, flower, growthVel, vFact
  integer :: ixc, iyc, izc, ixv, iyv, izv, k, kv, iax, nf
  integer :: ix, iy, iz, nx, ny, nz

  real, external :: func_dtdacube, rombint

  namelist /mask_params/ levelmin, levelmax, g9p_path, mus_path, pvar_val

  ! Read the namelist
  if(command_argument_count() .lt. 1) then
     print*,'Usage: ./write_g9p namelist'
     stop 1
  endif
  
  call get_command_argument(1, nml_file)

  open(11,file=trim(nml_file),status='old')
  read(11,nml=mask_params)
  close(11)
  print *,'Read the namelist:'
  print *,'---- levelmin', levelmin, 'levelmax', levelmax
  ! print *,'---- ix', ix, 'iy', iy, 'iz', iz
  ! print *,'---- nx', nx, 'ny', ny, 'nz', nz
  ! print *,'---- x', xi, 'y',  yi, 'z',  zi, 'r', r
  ! print *,'---- npad', npad
  print *,'---- g9p path ', trim(g9p_path)
  print *,'---- music path ', trim(mus_path)
  print *,''

  nlevel = levelmax - levelmin
  pvar_val = 0.5
  
  ! Read in some parameters from the g9p file
  write(ldim, '(I5)') 2**levelmin
  open(13,file=trim(g9p_path)//'g9p'//trim(adjustl(ldim))//'_delta',form='unformatted',status='old')
  read(13) n1_g9p, n1_g9p, n1_g9p, dx_g9p, x1o_g9p, x2o_g9p, x3o_g9p, astart, omega_m, omega_v, H0
  print *, 'Read in parameters from ', 'g9p'//trim(adjustl(ldim))//'_delta'
  print *, '---- astart', astart, 'zstart', (1./astart) - 1.
  print *, '---- H0', H0, 'omega_m', omega_m, 'omega_v', omega_v
  print *, ''
  close(13)
  ! Calculate box length from the cell spacing on the coarsest level
  boxlen = (2**levelmin) * dx_g9p

  
  ! allocate(n(nlevel+1, 3))
  ! allocate(off_abs(nlevel+1, 3))
  ! allocate(off_rel(nlevel+1, 3))

  ! allocate(n_w(3), off_w(3), off_abs(3))

  ! n_w = (/ 0., 0., 0. /)
  ! off_w = (/ 0., 0., 0. /)
  ! off_abs = (/ 0., 0., 0. /)
  
  ! Now write the rest of the IC files
  do ilevel = nlevel+1, 1, -1
     print *, 'Writing IC files for level ', ilevel + levelmin -1
     ! Utility variables for loading files
     ! Maximum possible number of cells in this level along each dimension
     write(lnum, '(I4)') ilevel + levelmin - 1 + 1000
     write(ldim, '(I5)') 2**(ilevel + levelmin - 1)
     
     open(10,file=trim(g9p_path)//'g9p'//trim(adjustl(ldim))//'_delta',form='unformatted',status='old')
     read(10) n1_g9p, n1_g9p, n1_g9p, dx_g9p, x1o_g9p, x2o_g9p, x3o_g9p, astart_w, omega_m_w, omega_v_w, H0_w
     close(10)

     open(10,file=trim(mus_path)//'level_'//lnum(2:4)//'/ic_refmap', form='unformatted')
     read(10) n_mus(1), n_mus(2), n_mus(3), dx_mus, xo_mus(1), xo_mus(2), xo_mus(3), astart, omega_m, omega_v, H0

     dx_w = dx_g9p ! offset and cell spacing comes from g9p
     
     do i=1,3
        off_abs(i) = nint(xo_mus(i) / dx_mus)  ! ilevel grid cells
        ! print*, off_abs(i)
        ! print*, xo_mus(i)
        ! print*, dx_mus
        print*, '-- off_abs diff', real(off_abs(i)) - (xo_mus(i) / dx_mus)
        off_w(i) = off_abs(i) * dx_w  ! offset and cell spacing comes from g9p
        n_w(i) = n_mus(i)  ! number of cells from music
     end do

     ! Write out the new refmaps with updated parameters
     open(11,file=trim(mus_path)//'level_'//lnum(2:4)//'/ic_refmap_new', form='unformatted')
     open(12,file=trim(mus_path)//'level_'//lnum(2:4)//'/ic_pvar_00001_new', form='unformatted')
     write(11) n_w(1), n_w(2), n_w(3), dx_w, off_w(1), off_w(2), off_w(3), astart_w, omega_m_w, omega_v_w, H0_w
     write(12) n_w(1), n_w(2), n_w(3), dx_w, off_w(1), off_w(2), off_w(3), astart_w, omega_m_w, omega_v_w, H0_w

     allocate(slab(n_w(1), n_w(2)))
     do k=1,n_w(3)
        read(10) slab
        write(11) slab
        write(12) slab * pvar_val
     end do

     close(10)
     close(11)
     close(12)
     deallocate(slab)
     
     
     ! Bits of the following do loop is (almost) directly pinched from
     ! Sergey Pilipenko's cubic_mask3.f90
     do iax = 0,3
        if(iax.eq.0) then
           open(13,file=trim(g9p_path)//'g9p'//trim(adjustl(ldim))//'_delta',form='unformatted',status='old')
           open(12,file=trim(mus_path)//'level_'//lnum(2:4)//'/ic_deltab',form='unformatted')
        else 
           open(13,file=trim(g9p_path)//'g9p'//trim(adjustl(ldim))//'_vel'//xyz(iax:iax),form='unformatted',status='old')
           open(12,file=trim(mus_path)//'level_'//lnum(2:4)//'/ic_velc'//xyz(iax:iax),form='unformatted')
           open(14,file=trim(mus_path)//'level_'//lnum(2:4)//'/ic_posc'//xyz(iax:iax),form='unformatted')
        endif
        read(13) n1_g9p, n1_g9p, n1_g9p, dx_g9p, x1o_g9p, x2o_g9p, x3o_g9p, astart, omega_m, omega_v, H0

        ! Calculate grid offsets in comoving Mpc and extent in grid cells
        ! do ii = 1, 3
        !    n_w(ii) = n(ilevel, ii)
        !    off_w(ii) = off_abs(ilevel, ii) * dx_g9p
        ! end do

        ! print*, 'dx_g9p', dx_g9p, 'dx_mus', dx_mus

        tol = 1e-5
        
        if (abs(dx_w - dx_g9p) .gt. tol) then
           print*, 'dx different', dx_w, dx_g9p
           call exit(1)
        end if
        if (abs(omega_m_w - omega_m) .gt. tol) then
           print*, 'omega_m different', omega_m_w, omega_m
           call exit(1)
        end if
        if (abs(omega_v_w - omega_v) .gt. tol) then
           print*, 'omega_v different', omega_v_w, omega_v
           call exit(1)
        end if
        if (abs(H0_w - H0) .gt. tol) then
           call exit(1)
           print*, 'H0 different', H0_w, H0
        end if
        if (abs(astart_w - astart) .gt. tol) then
           print*, 'astart different', astart_w, astart
           call exit(1)
        end if

        write(12) n_w(1), n_w(2), n_w(3), dx_w, off_w(1), off_w(2), off_w(3), astart_w, omega_m_w, omega_v_w, H0_w
        if(iax.gt.0) write(14) n_w(1), n_w(2), n_w(3), dx_w, off_w(1), off_w(2), off_w(3), astart_w, omega_m_w, omega_v_w, H0_w

        allocate(slab(n1_g9p, n1_g9p))

        ainv = 1. / astart
        adot = sqrt(omega_m*(ainv-1) + omega_v*(astart*astart-1) + 1)
        tmp1 = rombint(func_dtdacube, 1e-6, astart, 1e-5)
        growth = 2.5 * omega_m * adot * tmp1 / astart
        tmp  = 1. - omega_m - omega_v
        fupper = 2.5 * omega_m / growth - 1.5 * omega_m * ainv - tmp
        flower = omega_m * ainv + omega_v * astart * astart + tmp
        growthVel = fupper/flower

        vFact = 1.  / ( (adot) * 100. * (growthVel) )

        ixv = nint(x1o_g9p/dx_g9p)
        iyv = nint(x2o_g9p/dx_g9p)
        izv = nint(x3o_g9p/dx_g9p)

        ! Calculate i*c, since on levelmin we start from 1
        if (ilevel .ne. 1) then
           ixc = off_abs(1)
           iyc = off_abs(2)
           izc = off_abs(3)
        else
           ixc = 1
           iyc = 1
           izc = 1
        end if

        do k = 1, n1_g9p
           read(13) slab
           kv = k + izv
           ! if (ilevel .eq. 1) kv
           ! LC - swapped .ge. for .gt. and .lt. for .le. to account for 1-indexing and added +1
           ! LC 7/5 - swapped back
           if((kv .ge. izc) .and. (kv .lt. izc + n_w(3))) then
              write(12) slab(ixc-ixv:ixc+n_w(1)-1-ixv, &
                   iyc-iyv:iyc+n_w(2)-1-iyv)
              if(iax.gt.0) then
                 write(14) slab(ixc-ixv:ixc+n_w(1)-1-ixv, &
                      iyc-iyv:iyc+n_w(2)-1-iyv)  * vFact
              endif
           endif
        enddo
        deallocate(slab)

        close(13)
        close(12)
        if(iax.gt.0) close(14)
     enddo
  end do

  ! Finally, write the namelist
  call write_namelist()
  
  ! Using internal subroutines since the size of n, off_abs, and off_rel
  ! will be dynamically allocated and hence we need expicit interfaces
contains
  
  subroutine write_namelist
    use params

    integer :: ilev
    character(len=64) :: nexpand, str_levelmin
    character(len=3), dimension(10) :: i_nexpand

    open(20, file='nml.nml')!, status='new')

    ! RUN_PARAMS
    write(20, '(a)') '&RUN_PARAMS'
    write(20, '(a)') 'cosmo=.true.'
    write(20, '(a)') 'pic=.true.'
    write(20, '(a)') 'poisson=.true.'
    write(20, '(a)') 'hydro=.true.'
    write(20, '(a)') 'nrestart=0'
    write(20, '(a)') 'nremap=10'
    write(20, '(a)') 'ncontrol=1'
    write(20, '(a)') 'verbose=.false.'
    write(20, '(a)') '/'
    write(20, '(a)') ''

    ! INIT_PARAMS
    write(20, '(a)') '&INIT_PARAMS'
    write(20, '(a)') "filetype='grafic'"
    do ilev = levelmin, levelmax
       write(20, fmt='(a, i1, a, i3.3, a)') "initfile(",ilev-levelmin+1,")='level_",ilev,"'"
    end do
    write(20, '(a)') '/'
    write(20, '(a)') ''

    ! AMR_PARAMS
    write(20, '(a)') '&AMR_PARAMS'
    write(str_levelmin, '(i2)') levelmin
    str_levelmin = trim(adjustl(str_levelmin))
    write(20, fmt='(a, a)') 'levelmin=', str_levelmin
    write(20, fmt='(a, i2)') 'levelmax=', levelmax + 9
    write(20, '(a)') 'ngridmax=100000'
    write(20, '(a)') 'npartmax=200000'
    if (levelmin .eq. levelmax) then
       write(20, '(a)') 'nexpand=1'
    else
       nexpand='nexpand='
       do ilev=levelmin, levelmax-2
          write(i_nexpand(ilev-levelmin+1), '(i2, a)') (npad-1)/2, ','
          nexpand = trim(adjustl(nexpand))//trim(adjustl(i_nexpand(ilev-levelmin+1)))
       end do
       nexpand = trim(nexpand)//'1,1'
       write(20, '(a)') trim(nexpand)
    end if

    write(20, '(a)') '/'
    write(20, '(a)') ''

    ! REFINE_PARAMS
    write(20, '(a)') '&REFINE_PARAMS'
    write(20, fmt='(a, i2, a)') 'm_refine=', levelmax - levelmin + 10,'*8.,'
    write(20, fmt='(a, i1)') 'ivar_refine=',6  ! ic_pvar_00001
    write(20, fmt='(a, f5.3)') 'var_cut_refine=', 0.01
    write(20, fmt='(a, e11.5)') 'mass_cut_refine=', 2.0/(2.0 ** (levelmax * 3.0))
    write(20, '(a)') 'interpol_var=1'
    write(20, '(a)') 'interpol_type=0'
    write(20, '(a)') '/'
    write(20, '(a)') ''

    ! OUTPUT_PARAMS
    write(20, '(a)') '&OUTPUT_PARAMS'
    write(20, '(a)') '/'
    write(20, '(a)') ''

    ! HYDRO_PARAMS
    write(20, '(a)') '&HYDRO_PARAMS'
    write(20, '(a)') '/'
    write(20, '(a)') ''

    close(20)

  end subroutine write_namelist
 
end program ref_mask

! The following two functions are also directly pinched from Sergey
! Pilipenko's cubic_mask3.f90
real function func_dtdacube(a)
  !  real,intent(in)::a
  use params
  
  func_dtdacube = 1. / (omega_m * (1./a - 1.0) + omega_v * (a * a - 1.0) + 1.0)**1.5
end function func_dtdacube

real function rombint(f,a,b,tol)
        parameter (MAXITER=30,MAXJ=5)
        implicit real (a-h,o-z)
        dimension g(MAXJ+1)
        real a, b, tol
        external f
!        interface
!          real function f(x)
!            real,intent(in)::x
!          end function f
!        end interface
        h=0.5d0*(b-a)
        gmax=h*(f(a)+f(b))
        g(1)=gmax
        nint1=1
        error=1.0d20
        i=0
10        i=i+1
          if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) go to 40
!  Calculate next trapezoidal rule approximation to integral.
          g0=0.0d0
            do 20 k=1,nint1
            g0=g0+f(a+(k+k-1)*h)
20        continue
          g0=0.5d0*g(1)+h*g0
          h=0.5d0*h
          nint1=nint1+nint1
          jmax=min(i,MAXJ)
          fourj=1.0d0
            do 30 j=1,jmax
!  Use Richardson extrapolation.
            fourj=4.0d0*fourj
            g1=g0+(g0-g(j))/(fourj-1.0d0)
            g(j)=g0
            g0=g1
30        continue
          if (abs(g0).gt.tol) then
            error=1.0d0-gmax/g0
          else
            error=gmax
          end if
          gmax=g0
          g(jmax+1)=g0
        go to 10
40      rombint=g0
        if (i.gt.MAXITER.and.abs(error).gt.tol) write(*,*) 'Rombint failed to converge; integral, error=',rombint,error
        return
end function rombint
