! converts a restart file from CICE 4.1 to 5.0 format
! single processor (not parallelized)
! core restarts only (not including tracers except Tsfc)

! search for "modify as needed" to find potential configuration changes.
! NB: restarts will likely work only for exactly the same configuration
!     (BL99 thermo in particular).

! intel compiler:
! ifort convert_restarts.f90 -o convert_restarts -i4 -r8 -convert big_endian -assume byterecl -mcmodel=medium -shared-intel

!  ifort convert_restarts.f90 -o convert_restarts -i4 -r8 -convert big_endian -assume byterecl -mcmodel=medium -shared-intel -O0 -fpe0 -debug -ftrapuv `nf-config --fflags --flibs`
!
      program convert_restarts_new

      implicit none

      integer, parameter :: char_len_long  = 256, &
                            log_kind  = kind(.true.), &
                            int_kind  = selected_int_kind(6), &
                            dbl_kind  = selected_real_kind(13)

!!!!!!! modify as needed (begin)
      ! these values are for the standard gx1 configuration
      integer (kind=int_kind), parameter :: &
                            ncat = 5, &       ! number of thickness cats
                            nilyr = 4, &      ! number of ice layers
                            nslyr = 1, &      ! number of snow layers
                            nx_block = 4500, & ! global grid size
                            ny_block = 3297

      ! new layer number and layer variables
      integer (kind=int_kind), parameter                                :: nilyrnew = 7
      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyrnew,ncat) :: qicen, sicen
      real (kind=dbl_kind), dimension(nilyr)                            :: ain
      real (kind=dbl_kind), dimension(nilyrnew)                         :: aout

      ! these values are for the standard gx3 configuration
!      integer (kind=int_kind), parameter :: &
!                            ncat = 5, &       ! number of thickness cats
!                            nilyr = 4, &      ! number of ice layers
!                            nslyr = 1, &      ! number of snow layers
!                            nx_block = 100, & ! global grid size
!                            ny_block = 116

      ! flags
      logical (kind=log_kind), parameter :: &
         oceanmixed_ice = .true., & ! if true, read/write ocean mixed layer fields
         heat_capacity  = .true., & ! if true, ice has nonzero heat capacity
!#ifdef UNDEPRECATE_0LAYER
                                    ! if false, use zero-layer thermodynamics
!#else
                                    ! heat_capacity = .false. (zero-layer thermodynamics)
                                    ! has been deprecated in CICE and Icepack
!#endif
         diag = .true.              ! write min/max diagnostics for fields

      ! file names
      character (len=char_len_long), parameter :: &
         iced_4_1 = 'rtofs_glo.t00z.n-06.restart_cice', &  ! gx1
         iced_5_0 = 'rtofs_glo.t00z.n-06.restart_cice_converted'
!         iced_4_1 = '/scratch/eclare/tmp/restarts/iced_gx1_v4.0_kcatbound0', &  ! gx1
!         iced_5_0 = '/scratch/eclare/tmp/restarts/iced_gx1_v4.0_kcatbound0_converted'
!         iced_4_1 = 'iced_gx3_v4.0_kcatbound0', &  ! gx3
!         iced_5_0 = 'iced_gx3_v4.0_kcatbound0_converted''

!!!!!!! modify as needed (end)

      ! array sizes
      integer (kind=int_kind), parameter :: &
         ntilyr    = ncat*nilyr, & ! number of ice layers in all categories
         ntslyr    = ncat*nslyr, & ! number of snow layers in all categories
         max_ntrcr =   1         & ! 1 = surface temperature
                   + nilyr       & ! ice salinity
                   + nilyr       & ! ice enthalpy
                   + nslyr         ! snow enthalpy

      integer (kind=int_kind), dimension(ncat) :: &
         ilyr1          , & ! starting ice layer number for each category
         slyr1              ! starting snow layer number for each category

      integer (kind=int_kind) :: &
         ntrcr          , & ! number of required tracers
         i,j,n,k,m          ! indices

      ! tracer indices
      integer (kind=int_kind) :: &
         nt_Tsfc  , & ! ice/snow temperature
         nt_qice  , & ! volume-weighted ice enthalpy (in layers)
         nt_qsno  , & ! volume-weighted snow enthalpy (in layers)
         nt_sice      ! volume-weighted ice bulk salinity (CICE grid layers)

      ! time info
      integer (kind=int_kind) :: &
         istep1    ! counter, number of steps at current timestep

      real (kind=dbl_kind) :: &
         time           , & ! total elapsed time (s)
         time_forc          ! time of last forcing update (s)

      ! restart fields
      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ncat) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,max_ntrcr,ncat) :: &
         trcrn     ! tracers
                   ! 1: surface temperature of ice/snow (C)

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ntilyr) :: &
         eicen     ! energy of melting for each ice layer  (J/m^2)

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ntslyr) :: &
         esnon     ! energy of melting for each snow layer (J/m^2)

      real (kind=dbl_kind), dimension(nx_block,ny_block) :: &
         uvel    , & ! x-component of velocity (m/s)
         vvel    , & ! y-component of velocity (m/s)
         swvdr   , & ! sw down, visible, direct  (W/m^2)
         swvdf   , & ! sw down, visible, diffuse (W/m^2)
         swidr   , & ! sw down, near IR, direct  (W/m^2)
         swidf   , & ! sw down, near IR, diffuse (W/m^2)
         scale_factor, &! scaling factor for shortwave components
         strocnxT, & ! ice-ocean stress, x-direction
         strocnyT, & ! ice-ocean stress, y-direction
         stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4, & ! sigma12
         sst     , & ! sea surface temperature (C)
         frzmlt  , & ! freezing/melting potential (W/m^2)
         iceumask    ! ice extent mask (U-cell)

      ! flags
      logical (kind=log_kind) :: &
         l_brine     ! if true, treat brine pocket effects

      ! numbers
      real (kind=dbl_kind), parameter :: &
         spval_dbl = 1.0e30_dbl_kind, & ! special value (double precision)
         puny = 1.0e-11_dbl_kind, &
         pi   = 3.14159265358979323846_dbl_kind,&! pi
         c0   = 0.0_dbl_kind, &
         c1   = 1.0_dbl_kind, &
         c2   = 2.0_dbl_kind, &
         p5   = 0.5_dbl_kind

      ! physical parameters
      real (kind=dbl_kind), parameter :: &
         Lsub      = 2.835e6_dbl_kind ,&! latent heat, sublimation freshwater (J/kg)
         Lvap      = 2.501e6_dbl_kind ,&! latent heat, vaporization freshwater (J/kg)
         Lfresh    = Lsub-Lvap        ,&! latent heat of melting of fresh ice (J/kg)
         cp_ice    = 2106._dbl_kind   ,&! specific heat of fresh ice (J/kg/K)
         rhos      = 330.0_dbl_kind   ,&! density of snow (kg/m^3)
         hs_min  = 1.e-4_dbl_kind, &  ! min snow thickness for computing Tsno (m)
         nsal    = 0.407_dbl_kind, &
         msal    = 0.573_dbl_kind, &
         min_salin = 0.1_dbl_kind, &  ! threshold for brine pocket treatment
         saltmax = 3.2_dbl_kind       ! max salinity at ice base (ppt)

      ! useful temps
      real (kind=dbl_kind), dimension(nilyr+1) :: &
         salin       ! initial salinity  profile (ppt)

      real (kind=dbl_kind)    :: &
         Tmin, Tmax,    & ! min and max snow temperature
         zTsn,          & ! snow temperature
         zn,            & ! thickness
         rnslyr,        & ! real(nslyr)
         rnilyr           ! real(nilyr)

      ! count tracers
         nt_Tsfc = 1           ! index tracers, starting with Tsfc = 1
         ntrcr = 1             ! count tracers, starting with Tsfc = 1
         nt_qice = ntrcr + 1
         ntrcr = ntrcr + nilyr ! qice in nilyr layers
         nt_qsno = ntrcr + 1
         ntrcr = ntrcr + nslyr ! qsno in nslyr layers
         nt_sice = ntrcr + 1
         ntrcr = ntrcr + nilyr ! sice in nilyr layers
         if (ntrcr /= max_ntrcr) write (*,*) 'Tracer number mismatch'

      ! ice salinity profile
      if (saltmax > min_salin .and. heat_capacity) then
         l_brine = .true.
      else
         l_brine = .false.
      endif

      if (l_brine) then
         do k = 1, nilyr
            zn = (real(k,kind=dbl_kind)-p5) /  &
                  real(nilyr,kind=dbl_kind)
            salin(k)=(saltmax/c2)*(c1-cos(pi*zn**(nsal/(msal+zn))))
!            salin(k)=saltmax ! for isosaline ice ! modify as needed
         enddo
         salin(nilyr+1) = saltmax
      else
         do k = 1, nilyr+1
            salin(k) = c0
         enddo
      endif
      do k = 1, nilyr
         trcrn(:,:,nt_sice+k-1,:) = salin(k)
      enddo

      ! read cice 4.1 restart data
      call restartfile_4_1 (iced_4_1)

      ! convert eicen, esnon to qicen, qsnon
      ilyr1(1) = 1
      slyr1(1) = 1
      do n = 2, ncat
         ilyr1(n) = ilyr1(n-1) + nilyr
         slyr1(n) = slyr1(n-1) + nslyr
      enddo

      rnslyr = real(nslyr, kind=dbl_kind)
      rnilyr = real(nilyr, kind=dbl_kind)
      Tmin = -100.  ! minimum allowed snow temperature
      do j = 1, ny_block
      do i = 1, nx_block
      do n = 1, ncat
      if (aicen(i,j,n) > puny) then
      if (vsnon(i,j,n)/aicen(i,j,n) > hs_min) then
         do k = 1, nslyr
            ! qsn, esnon < 0
            trcrn(i,j,nt_qsno+k-1,n) = &
               min(esnon(i,j,slyr1(n)+k-1)*rnslyr/vsnon(i,j,n),-rhos*Lfresh)
            Tmax = -trcrn(i,j,nt_qsno+k-1,n)*puny*rnslyr/(rhos*cp_ice*vsnon(i,j,n))
            if (.not. heat_capacity) then
               trcrn(i,j,nt_qsno+k-1,n) = -rhos * Lfresh
               Tmax = puny
            endif
!!!!!!! modify as needed (begin)
! if your restarts do not work, try uncommenting this section
!            ! snow temperature
!            zTsn = (Lfresh + trcrn(i,j,nt_qsno+k-1,n)/rhos)/cp_ice
!            ! zap entire snow volume if temperature is out of bounds
!            if (zTsn < Tmin .or. zTsn > Tmax) then
!               print*, 'zapping snow volume ', i,j,n,vsnon(i,j,n)
!               vsnon(i,j,n) = c0
!               do m = 1, nslyr
!                  trcrn(i,j,nt_qsno+m-1,n) = c0
!               enddo
!            endif
!!!!!!! modify as needed (begin)
         enddo
      else
         vsnon(i,j,n) = c0
         do k = 1, nslyr
            trcrn(i,j,nt_qsno+k-1,n) = c0
         enddo
      endif
      endif
      do k = 1, nilyr
!         if (vicen(i,j,n) > puny) then
         if (aicen(i,j,n) > puny) then  ! matches v4.1
            trcrn(i,j,nt_qice+k-1,n) = eicen(i,j,ilyr1(n)+k-1)*rnilyr/vicen(i,j,n)
         else
            trcrn(i,j,nt_qice+k-1,n) = c0
         endif
      enddo
      enddo
      enddo
      enddo

      ! write cice 5.0 restart data
      !call dumpfile_5_0 (iced_5_0)

      ! write cice 5.0 netcdf restart file, original layers
      !call write_netcdf(trim(iced_5_0)//'.nc', .false.)

      ! interpolate from nilyr to nilyrnew
      call interpolate_newlayers(trcrn,nt_sice,sicen)
      call interpolate_newlayers(trcrn,nt_qice,qicen)

      ! write cice 5.0 netcdf restart file, new layers
      call write_netcdf(trim(iced_5_0)//'.intp.nc', .true.)

!=======================================================================

      contains

!=======================================================================
!---! these subroutines write/read Fortran unformatted data files
!=======================================================================

      subroutine restartfile_4_1 (ice_ic)

      character (*) :: ice_ic

      integer (kind=int_kind) :: i, j, k, n

      character(len=char_len_long) :: &
         filename

      integer (kind=int_kind), parameter :: nu_restart = 20

         filename = ice_ic
         open(nu_restart,file=filename,form='unformatted')

         write(*,*) 'Using restart dump=', trim(filename)
         read (nu_restart) istep1,time,time_forc
         write(*,*) 'Restart read at istep=',istep1,time,time_forc

      !-----------------------------------------------------------------
      ! state variables
      ! Tsfc is the only tracer read in this file.  All other
      ! tracers are in their own dump/restart files.
      !-----------------------------------------------------------------
      do n=1,ncat
              write(*,*) 'cat ',n, &
                               ' min/max area, vol ice, vol snow, Tsfc'

         call ice_read(nu_restart,aicen(:,:,n),diag)
         call ice_read(nu_restart,vicen(:,:,n),diag)
         call ice_read(nu_restart,vsnon(:,:,n),diag)
         call ice_read(nu_restart,trcrn(:,:,nt_Tsfc,n),diag)
         !print *,n,trcrn(370:391,3000,1,n)
      enddo

           write(*,*) 'min/max eicen for each layer'
      do k=1,ntilyr
         call ice_read(nu_restart,eicen(:,:,k),diag)
      enddo

           write(*,*) 'min/max esnon for each layer'
      do k=1,ntslyr
         call ice_read(nu_restart,esnon(:,:,k),diag)
      enddo

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
           write(*,*) 'min/max velocity components'

      call ice_read(nu_restart,uvel,diag)
      call ice_read(nu_restart,vvel,diag)

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------
         write(*,*) 'radiation fields'

      call ice_read(nu_restart,scale_factor,diag)
      call ice_read(nu_restart,swvdr,diag)
      call ice_read(nu_restart,swvdf,diag)
      call ice_read(nu_restart,swidr,diag)
      call ice_read(nu_restart,swidf,diag)

      !-----------------------------------------------------------------
      ! ocean stress
      !-----------------------------------------------------------------
           write(*,*) 'min/max ocean stress components'

      call ice_read(nu_restart,strocnxT,diag)
      call ice_read(nu_restart,strocnyT,diag)

      !-----------------------------------------------------------------
      ! internal stress
      !-----------------------------------------------------------------
           write(*,*) 'internal stress components'

      call ice_read(nu_restart,stressp_1,diag) ! stressp_1
      call ice_read(nu_restart,stressp_3,diag) ! stressp_3

      call ice_read(nu_restart,stressp_2,diag) ! stressp_2
      call ice_read(nu_restart,stressp_4,diag) ! stressp_4

      call ice_read(nu_restart,stressm_1,diag) ! stressm_1
      call ice_read(nu_restart,stressm_3,diag) ! stressm_3

      call ice_read(nu_restart,stressm_2,diag) ! stressm_2
      call ice_read(nu_restart,stressm_4,diag) ! stressm_4

      call ice_read(nu_restart,stress12_1,diag) ! stress12_1
      call ice_read(nu_restart,stress12_3,diag) ! stress12_3

      call ice_read(nu_restart,stress12_2,diag) ! stress12_2
      call ice_read(nu_restart,stress12_4,diag) ! stress12_4

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
           write(*,*) 'ice mask for dynamics'

      call ice_read(nu_restart,iceumask,diag)

      ! for mixed layer model
      if (oceanmixed_ice) then
              write(*,*) 'min/max sst, frzmlt'

         call ice_read(nu_restart,sst,diag)
         call ice_read(nu_restart,frzmlt,diag)
      endif

      !-----------------------------------------------------------------
      ! mask out all-land blocks
      !-----------------------------------------------------------------
      where (aicen        > 0.5*spval_dbl) aicen        = c0
      where (vicen        > 0.5*spval_dbl) vicen        = c0
      where (vsnon        > 0.5*spval_dbl) vsnon        = c0
      where (trcrn        > 0.5*spval_dbl) trcrn        = c0
      where (eicen        > 0.5*spval_dbl) eicen        = c0
      where (esnon        > 0.5*spval_dbl) esnon        = c0
      where (uvel         > 0.5*spval_dbl) uvel         = c0
      where (vvel         > 0.5*spval_dbl) vvel         = c0
      where (scale_factor > 0.5*spval_dbl) scale_factor = c0
      where (swvdr        > 0.5*spval_dbl) swvdr        = c0
      where (swvdf        > 0.5*spval_dbl) swvdf        = c0
      where (swidr        > 0.5*spval_dbl) swidr        = c0
      where (swidf        > 0.5*spval_dbl) swidf        = c0
      where (strocnxT     > 0.5*spval_dbl) strocnxT     = c0
      where (strocnyT     > 0.5*spval_dbl) strocnyT     = c0
      where (stressp_1    > 0.5*spval_dbl) stressp_1    = c0
      where (stressp_2    > 0.5*spval_dbl) stressp_2    = c0
      where (stressp_3    > 0.5*spval_dbl) stressp_3    = c0
      where (stressp_4    > 0.5*spval_dbl) stressp_4    = c0
      where (stressm_1    > 0.5*spval_dbl) stressm_1    = c0
      where (stressm_2    > 0.5*spval_dbl) stressm_2    = c0
      where (stressm_3    > 0.5*spval_dbl) stressm_3    = c0
      where (stressm_4    > 0.5*spval_dbl) stressm_4    = c0
      where (stress12_1   > 0.5*spval_dbl) stress12_1   = c0
      where (stress12_2   > 0.5*spval_dbl) stress12_2   = c0
      where (stress12_3   > 0.5*spval_dbl) stress12_3   = c0
      where (stress12_4   > 0.5*spval_dbl) stress12_4   = c0
      where (iceumask     > 0.5*spval_dbl) iceumask     = c0
      if (oceanmixed_ice) then
         where (sst       > 0.5*spval_dbl) sst          = c0
         where (frzmlt    > 0.5*spval_dbl) frzmlt       = c0
      endif

      close(nu_restart)

      end subroutine restartfile_4_1

!=======================================================================

      subroutine dumpfile_5_0(filename_spec)

      character(len=char_len_long), intent(in) :: filename_spec

      integer (kind=int_kind) :: i, j, k, n

      character(len=char_len_long) :: filename

      integer (kind=int_kind), parameter :: nu_dump = 21

        filename = trim(filename_spec)
        open(nu_dump,file=filename,form='unformatted')

        write(nu_dump) istep1,time,time_forc
        write(*,*) 'Writing ',trim(filename)
        write(*,*) 'Restart written ',istep1,time,time_forc

      !-----------------------------------------------------------------
      ! state variables
      ! Tsfc is the only tracer written to this file.  All other
      ! tracers are written to their own dump/restart files.
      !-----------------------------------------------------------------

      do n=1,ncat
              write(*,*) 'cat ',n, &
                               ' min/max area, vol ice, vol snow, Tsfc'

         call ice_write(nu_dump,aicen(:,:,n),diag)
         call ice_write(nu_dump,vicen(:,:,n),diag)
         call ice_write(nu_dump,vsnon(:,:,n),diag)
         call ice_write(nu_dump,trcrn(:,:,nt_Tsfc,n),diag)

           write(*,*) 'min/max sicen for each layer'
         do k=1,nilyr
            call ice_write(nu_dump,trcrn(:,:,nt_sice+k-1,n),diag)
         enddo

           write(*,*) 'min/max qicen for each layer'
         do k=1,nilyr
            call ice_write(nu_dump,trcrn(:,:,nt_qice+k-1,n),diag)
         enddo

           write(*,*) 'min/max qsnon for each layer'
         do k=1,nslyr
            call ice_write(nu_dump,trcrn(:,:,nt_qsno+k-1,n),diag)
         enddo
      enddo ! ncat

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
           write(*,*) 'min/max velocity components'

      call ice_write(nu_dump,uvel,diag)
      call ice_write(nu_dump,vvel,diag)

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------
         write(*,*) 'radiation fields'

      call ice_write(nu_dump,scale_factor,diag)
      call ice_write(nu_dump,swvdr,diag)
      call ice_write(nu_dump,swvdf,diag)
      call ice_write(nu_dump,swidr,diag)
      call ice_write(nu_dump,swidf,diag)

      !-----------------------------------------------------------------
      ! ocean stress (for bottom heat flux in thermo)
      !-----------------------------------------------------------------
           write(*,*) 'min/max ocean stress components'

      call ice_write(nu_dump,strocnxT,diag)
      call ice_write(nu_dump,strocnyT,diag)

      !-----------------------------------------------------------------
      ! internal stress
      !-----------------------------------------------------------------
           write(*,*) 'internal stress components'

      call ice_write(nu_dump,stressp_1,diag)
      call ice_write(nu_dump,stressp_3,diag)
      call ice_write(nu_dump,stressp_2,diag)
      call ice_write(nu_dump,stressp_4,diag)

      call ice_write(nu_dump,stressm_1,diag)
      call ice_write(nu_dump,stressm_3,diag)
      call ice_write(nu_dump,stressm_2,diag)
      call ice_write(nu_dump,stressm_4,diag)

      call ice_write(nu_dump,stress12_1,diag)
      call ice_write(nu_dump,stress12_3,diag)
      call ice_write(nu_dump,stress12_2,diag)
      call ice_write(nu_dump,stress12_4,diag)

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
           write(*,*) 'ice mask for dynamics'

      call ice_write(nu_dump,iceumask,diag)

      ! for mixed layer model
      if (oceanmixed_ice) then
              write(*,*) 'min/max sst, frzmlt'

         call ice_write(nu_dump,sst,diag)
         call ice_write(nu_dump,frzmlt,diag)
      endif

      close(nu_dump)

      end subroutine dumpfile_5_0

!=======================================================================

      subroutine ice_read(nu, work_g1, diag)

      integer (kind=int_kind), intent(in) :: &
           nu                ! unit number

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
           intent(out) :: &
           work_g1           ! output array (real, 8-byte)

      logical (kind=log_kind), intent(in) :: &
           diag              ! if true, write diagnostic output

      integer (kind=int_kind) :: i, j

      real (kind=dbl_kind) :: &
         amin, amax         ! min and max values of input array

               read(nu) ((work_g1(i,j),i=1,nx_block),j=1,ny_block)

      if (diag) then
         amin = minval(work_g1)
         amax = maxval(work_g1, mask = work_g1 /= spval_dbl)
         write(*,*) ' read_global ',nu, amin, amax
      endif

      end subroutine ice_read

!=======================================================================

      subroutine ice_write(nu, work_g1, diag)

      integer (kind=int_kind), intent(in) :: &
           nu                ! unit number

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
           intent(in) :: &
           work_g1           ! input array (real, 8-byte)

      logical (kind=log_kind), intent(in) :: &
           diag              ! if true, write diagnostic output
!
      integer (kind=int_kind) :: i, j

      real (kind=dbl_kind) :: &
         amin, amax     ! min and max values of ouput array

            write(nu) ((work_g1(i,j),i=1,nx_block),j=1,ny_block)

         if (diag) then
            amin = minval(work_g1)
            amax = maxval(work_g1)
            write(*,*) ' write_global ', nu, amin, amax
         endif

      end subroutine ice_write

      subroutine write_netcdf(fname, new_layers)
        use netcdf

        character(len=*), intent(in) :: fname
        logical         , intent(in) :: new_layers

        ! local variables
        integer, parameter :: ncatvars = 4, nlyrvars = 2, nsnwvars = 1, n2dvars = 23
        integer            :: n, nilyr_out
        integer            :: ncid, rc, idimid, jdimid, kdimid, id
        integer            :: dim3(3), dim2(2)
        character(len=20)  :: vname
        character(len= 3)  :: nchar

        real (kind=dbl_kind), dimension (nx_block,ny_block,ncat) :: tmpvar
        real (kind=dbl_kind), dimension (nx_block,ny_block)      :: coszen

        character(len=20), dimension(ncatvars) ::  catvars = (/"aicen", "vicen", "vsnon", "Tsfcn"/)
        character(len=20), dimension(nlyrvars) ::  lyrvars = (/"sice", "qice"/)
        character(len=20), dimension(nsnwvars) ::  snwvars = (/"qsno"/)
        character(len=20), dimension(n2dvars)  :: varnames = &
             (/ "uvel        ", "vvel        ", "scale_factor", "coszen      ",  &
                "swvdr       ", "swvdf       ", "swidr       ", "swidf       ",  &
                "stressp_1   ", "stressp_2   ", "stressp_3   ", "stressp_4   ",  &
                "stressm_1   ", "stressm_2   ", "stressm_3   ", "stressm_4   ",  &
                "stress12_1  ", "stress12_2  ", "stress12_3  ", "stress12_4  ",  &
                "strocnxT    ", "strocnyT    ", "iceumask    "/)

        coszen = c0
        if (new_layers) then
           nilyr_out = nilyrnew
        else
           nilyr_out = nilyr
        end if
        print *,'creating netcdf restart ',trim(fname)

        rc = nf90_create(trim(fname), nf90_64bit_offset, ncid)
        rc = nf90_def_dim(ncid,   'ni',  nx_block, idimid)
        rc = nf90_def_dim(ncid,   'nj',  ny_block, jdimid)
        rc = nf90_def_dim(ncid, 'ncat',  ncat,     kdimid)

        ! define category variables
        dim3(:) = (/idimid, jdimid, kdimid/)
        do n = 1,ncatvars
           vname = trim(catvars(n))
           rc = nf90_def_var(ncid, vname, nf90_double, dim3, id)
        end do
        ! define ice layer variables
        do n = 1,nlyrvars
           do k = 1,nilyr_out
              write(nchar,'(i3.3)') k
              vname = trim(lyrvars(n))//trim(nchar)
              rc = nf90_def_var(ncid, vname, nf90_double, dim3, id)
           end do
        end do
        ! define snow layer variables
        do n = 1,nsnwvars
           do k = 1,nslyr
              write(nchar,'(i3.3)') k
              vname = trim(snwvars(n))//trim(nchar)
              rc = nf90_def_var(ncid, vname, nf90_double, dim3, id)
           end do
        end do
        !define 2d variables
        dim2(:) = (/idimid, jdimid/)
        do n = 1,n2dvars
           vname = trim(varnames(n))
           rc = nf90_def_var(ncid, vname, nf90_double, dim2, id)
        end do
        rc = nf90_enddef(ncid)

        ! write variables
        rc = nf90_inq_varid(ncid,  'aicen', id)
        rc = nf90_put_var(ncid, id, aicen)
        print *,'aicen ',minval(aicen), maxval(aicen)

        rc = nf90_inq_varid(ncid,  'vicen', id)
        rc = nf90_put_var(ncid, id, vicen)
        print *,'vicen ',minval(vicen), maxval(vicen)

        rc = nf90_inq_varid(ncid,  'vsnon', id)
        rc = nf90_put_var(ncid, id, vsnon)
        print *,'vsnon ',minval(vsnon), maxval(vsnon)

        tmpvar(:,:,:) = trcrn(:,:,nt_Tsfc,:)
        rc = nf90_inq_varid(ncid,  'Tsfcn', id)
        rc = nf90_put_var(ncid, id, tmpvar)
        print *,'Tsfcn ',minval(tmpvar), maxval(tmpvar)

        n = 1   ! sice
        do k = 1,nilyr_out
           write(nchar,'(i3.3)') k
           vname = trim(lyrvars(n))//trim(nchar)
           if (new_layers) then
              tmpvar(:,:,:) = sicen(:,:,k,:)
           else
              tmpvar(:,:,:) = trcrn(:,:,nt_sice+k-1,:)
           end if
           print '(a,i5,2g15.7)',trim(vname)//' layer ',k,minval(tmpvar),maxval(tmpvar)
           rc = nf90_inq_varid(ncid, trim(vname), id)
           rc = nf90_put_var(ncid, id, tmpvar)
        end do
        n = 2   ! qice
        do k = 1,nilyr_out
           write(nchar,'(i3.3)') k
           vname = trim(lyrvars(n))//trim(nchar)
           if (new_layers) then
              tmpvar(:,:,:) = qicen(:,:,k,:)
           else
              tmpvar(:,:,:) = trcrn(:,:,nt_qice+k-1,:)
           endif
           print '(a,i5,2g15.7)',trim(vname)//' layer ',k,minval(tmpvar),maxval(tmpvar)
           rc = nf90_inq_varid(ncid, trim(vname), id)
           rc = nf90_put_var(ncid, id, tmpvar)
        end do

        n = 1  ! qsno
        do k = 1,nslyr
           write(nchar,'(i3.3)') k
           vname = trim(snwvars(n))//trim(nchar)
           tmpvar(:,:,:) = trcrn(:,:,nt_qsno+k-1,:)
           print '(a,i5,2g15.7)',trim(vname)//' tracer index ',nt_qsno+k-1,minval(tmpvar),maxval(tmpvar)
           rc = nf90_inq_varid(ncid, trim(vname), id)
           rc = nf90_put_var(ncid, id, tmpvar)
        end do

        rc = nf90_inq_varid(ncid, 'uvel', id)
        rc = nf90_put_var(ncid, id, uvel)
        print '(a,2g15.7)', 'uvel', minval(uvel), maxval(uvel)
        rc = nf90_inq_varid(ncid, 'vvel', id)
        rc = nf90_put_var(ncid, id, vvel)
        print '(a,2g15.7)', 'vvel', minval(vvel), maxval(vvel)

        rc = nf90_inq_varid(ncid, 'scale_factor', id)
        rc = nf90_put_var(ncid, id, scale_factor)
        print '(a,2g15.7)', 'scale_factor', minval(scale_factor), maxval(scale_factor)
        rc = nf90_inq_varid(ncid, 'coszen', id)
        rc = nf90_put_var(ncid, id, coszen)
        print '(a,2g15.7)', 'coszen', minval(coszen), maxval(coszen)

        rc = nf90_inq_varid(ncid, 'swvdr', id)
        rc = nf90_put_var(ncid, id, swvdr)
        print '(a,2g15.7)', 'swvdr', minval(swvdr), maxval(swvdr)
        rc = nf90_inq_varid(ncid, 'swvdf', id)
        rc = nf90_put_var(ncid, id, swvdf)
        print '(a,2g15.7)', 'swvdf', minval(swvdf), maxval(swvdf)
        rc = nf90_inq_varid(ncid, 'swidr', id)
        rc = nf90_put_var(ncid, id, swidr)
        print '(a,2g15.7)', 'swidr', minval(swidr), maxval(swidr)
        rc = nf90_inq_varid(ncid, 'swidf', id)
        rc = nf90_put_var(ncid, id, swidf)
        print '(a,2g15.7)', 'swidf', minval(swidf), maxval(swidf)

        rc = nf90_inq_varid(ncid, 'strocnxT', id)
        rc = nf90_put_var(ncid, id, strocnxT)
        print '(a,2g15.7)', 'strocnxT', minval(strocnxT), maxval(strocnxT)
        rc = nf90_inq_varid(ncid, 'strocnyT', id)
        rc = nf90_put_var(ncid, id, strocnyT)
        print '(a,2g15.7)', 'strocnyT', minval(strocnyT), maxval(strocnyT)

        rc = nf90_inq_varid(ncid, 'stressp_1', id)
        rc = nf90_put_var(ncid, id, stressp_1)
        print '(a,2g15.7)', 'stressp_1', minval(stressp_1), maxval(stressp_1)
        rc = nf90_inq_varid(ncid, 'stressp_2', id)
        rc = nf90_put_var(ncid, id, stressp_2)
        print '(a,2g15.7)', 'stressp_2', minval(stressp_2), maxval(stressp_2)
        rc = nf90_inq_varid(ncid, 'stressp_3', id)
        rc = nf90_put_var(ncid, id, stressp_3)
        print '(a,2g15.7)', 'stressp_3', minval(stressp_3), maxval(stressp_3)
        rc = nf90_inq_varid(ncid, 'stressp_4', id)
        rc = nf90_put_var(ncid, id, stressp_4)
        print '(a,2g15.7)', 'stressp_4', minval(stressp_4), maxval(stressp_4)

        rc = nf90_inq_varid(ncid, 'stressm_1', id)
        rc = nf90_put_var(ncid, id, stressm_1)
        print '(a,2g15.7)', 'stressm_1', minval(stressm_1), maxval(stressm_1)
        rc = nf90_inq_varid(ncid, 'stressm_2', id)
        rc = nf90_put_var(ncid, id, stressm_2)
        print '(a,2g15.7)', 'stressm_2', minval(stressm_2), maxval(stressm_2)
        rc = nf90_inq_varid(ncid, 'stressm_3', id)
        rc = nf90_put_var(ncid, id, stressm_3)
        print '(a,2g15.7)', 'stressm_3', minval(stressm_3), maxval(stressm_3)
        rc = nf90_inq_varid(ncid, 'stressm_4', id)
        rc = nf90_put_var(ncid, id, stressm_4)
        print '(a,2g15.7)', 'stressm_4', minval(stressm_4), maxval(stressm_4)

        rc = nf90_inq_varid(ncid, 'stress12_1', id)
        rc = nf90_put_var(ncid, id, stress12_1)
        print '(a,2g15.7)', 'stress12_1', minval(stress12_1), maxval(stress12_1)
        rc = nf90_inq_varid(ncid, 'stress12_2', id)
        rc = nf90_put_var(ncid, id, stress12_2)
        print '(a,2g15.7)', 'stress12_2', minval(stress12_2), maxval(stress12_2)
        rc = nf90_inq_varid(ncid, 'stress12_3', id)
        rc = nf90_put_var(ncid, id, stress12_3)
        print '(a,2g15.7)', 'stress12_3', minval(stress12_3), maxval(stress12_3)
        rc = nf90_inq_varid(ncid, 'stress12_4', id)
        rc = nf90_put_var(ncid, id, stress12_4)
        print '(a,2g15.7)', 'stress12_4', minval(stress12_4), maxval(stress12_4)

        rc = nf90_inq_varid(ncid, 'iceumask', id)
        rc = nf90_put_var(ncid, id, iceumask)
        print '(a,2g15.7)', 'iceumask', minval(iceumask), maxval(iceumask)

        rc = nf90_close(ncid)

      end subroutine write_netcdf

      subroutine interpolate_newlayers(ain,ntstart,aout)

        real (kind = dbl_kind),  intent(in)  :: ain(nx_block,ny_block,max_ntrcr,ncat)
        integer (kind=int_kind), intent(in)  :: ntstart
        real (kind = dbl_kind),  intent(out) :: aout(nx_block,ny_block,nilyrnew,ncat)

        ! local variables
        integer (kind=int_kind) :: kk
        real (kind = dbl_kind)  :: delh,delhnew,h,slope
        real (kind = dbl_kind)  :: xi(nilyr), fi(nilyr)         ! input layer axis and values
        real (kind = dbl_kind)  :: xo(nilyrnew), fo(nilyrnew)   ! output layer axis and values

        aout = 0.0d0
        ! define data axis for old and new layers
        h = 1.0d0
        delh = h/float(nilyr)
        delhnew = h/float(nilyrnew)

        xi(1) = delh/2.0d0
        do k = 2,nilyr
           xi(k) = xi(k-1) + delh
        end do
        !print *,xi
        xo(1) = delhnew/2.0d0
        do k = 2,nilyrnew
           xo(k) = xo(k-1) + delhnew
        end do
        !print *,xo

        ! i = 3000; j = 3000 ;
        ! do k = 1,nilyr
        !    print '(i4,6e14.5)',k,xi(k),(ain(i,j,ntstart+k-1,n),n=1,5)
        ! end do

        do n = 1,ncat
           do j = 1,ny_block
              do i = 1,nx_block
                 ! nilyr values
                 do k = 1,nilyr
                    fi(k) = ain(i,j,ntstart+k-1,n)
                 end do
                 ! nilyrnew values; piecewise linear interpolation
                 ! set new top and bottom values
                 fo(1)        = fi(1)
                 fo(nilyrnew) = fi(nilyr)
                 do k = 2,nilyrnew-1
                    do kk = 1,nilyr-1
                       if (xo(k) .gt. xi(kk) .and. xo(k) .lt. xi(kk+1) ) then
                          slope = (fi(kk+1) - fi(kk))/(xi(kk+1) - xi(kk))
                          fo(k) = fi(kk) + slope * (xo(k) - xi(kk))
                       end if
                    end do
                 end do
                 aout(i,j,:,n) = fo(:)
              end do
           end do
        end do

        ! i = 3000; j = 3000 ; n = 5
        ! do k = 1,nilyrnew
        !    print '(i4,6e14.5)',k,xo(k),(aout(i,j,k,n),n=1,5)
        ! end do
        ! print *

      end subroutine interpolate_newlayers

!=======================================================================

      end program convert_restarts_new

!=======================================================================
