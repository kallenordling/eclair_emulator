!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_emulator.f90
!!
!! \brief
!! Eclair emulator module
!! 
!!
!! \author Kalle Nordling (FMI)
!!
!! \responsible_coder
!! Kalle Nordling, kalle.nordling@fmi.fi
!! Jukka-Pekka Keskinen, jukka-pekka.keskinen@fmi.fi
!!
!! \revision_history
!!   - K. Nordling  (FMI) - Original code (2018-2019)
!!   - K. Nordling  (FMI) - Added cdnc and cos_mu variabiables to emulator inputs
!!   - J-P Keskinen (FMI) - Train and standardisation subroutines are now 
!!                          more general. (6/2019)
!!   - J-P Keskinen (FMI) - Inclusion of the precipitation emulator. (10/2019)
!!   - K. Nordling (FMI)  - New cloud top/base subrutine from P.Räisänen added for testing(01/2020)
!!   - K. Nordling (FMI)  - Only values inside training data are accpeted to emulator (02/2020)
!!   - K. Nordling (FMI)  - included 3D mask where emulator values are used (02/2020)
!!   - K. Nordling (FMI)  - CALL standardize_emulator_inputs - modified so that input and output vectors are seperate (02/2020)
!! \limitations
!!
!! \details
!!
!! \belongs_to
!!
!! \copyright
!! Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!! licencing agreement to be found at:
!! https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!! The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef __xlC__
@PROCESS HOT
@PROCESS XLF90(NOSIGNEDZERO)
#else
#define FSEL(a,b,c) MERGE(b,c,(a) >= 0._wp)
#define SWDIV_NOCHK(a,b) ((a)/(b))
#endif

MODULE mo_emulator

  USE mo_kind,                 ONLY : wp
  USE mo_physical_constants,   ONLY : vtmpc1, cpd, grav
  USE m_gp
  USE m_gp_example
  USE mo_mpi,                  ONLY: p_parallel,p_parallel_io,p_io,p_bcast
  USE mo_netcdf,               ONLY: nf_max_name, nf_open, NF_NOERR, nf_nowrite, nf_inq_dimid, nf_inq_dimlen, &
                              nf_get_var_double, nf_get_vara_double, nf_close, nf_inq_varid
  USE mo_decomposition,        ONLY: ldc=>local_decomposition
  USE mo_exception,            ONLY: finish
  USE mo_tracer,               ONLY: get_tracer
  USE mo_physical_constants,   ONLY : cpd, grav, rgrav, rd, alv, als, rv
  USE mo_memory_g3b,           ONLY : cos_mu
  USE mo_time_control,         ONLY : lstart, lresume
  USE mo_control,              ONLY: emul_cloud  ! Switch for cloud top/base calculation
  !USE m_gp_dence,	       ONLY: initDenseGP
  IMPLICIT NONE
 

  PRIVATE
  !Define emulators, global
  CLASS(BaseGP), allocatable :: gp_emu_w_n, gp_emu_w_d, gp_emu_precip_n, gp_emu_precip_d
  ! Define arrays needed by the standardisation of emulator input and output
  REAL(wp), DIMENSION(500,6) :: ex1_w_n, ex1_precip_n
  REAL(wp), DIMENSION(497,7) :: ex1_precip_d
  REAL(wp), DIMENSION(497,7) :: ex1_w_d
  REAL(wp), DIMENSION(497)   :: t1_w_d
  REAL(wp), DIMENSION(500)   :: t1_w_n, t1_precip_n, t1_precip_d

  ! type definition
  TYPE t_emulator
    REAL(wp),  POINTER             :: emu_mask(:,:)      ! emulator mask 
    REAL(wp),  POINTER             :: emu_mask_3d(:,:,:) ! emulator mask 3d 
    REAL(wp),  POINTER             :: autoconv_int(:,:)  ! Autoconversion [m]
    REAL(wp),  POINTER             :: accretion1_int(:,:)! accretion term 1 [m]
    REAL(wp),  POINTER             :: accretion2_int(:,:)! accretion term 2 [m]
    REAL(wp),  POINTER             :: rainevap_int(:,:)  ! rain evaporation [m]
    REAL(wp),  POINTER             :: cloud_top(:,:)     ! emulated cloud top [m]
    REAL(wp),  POINTER             :: cloud_base(:,:)    ! emulated cloud base [m]
    REAL(wp),  POINTER             :: hpalim700(:,:)     ! 700hpa limit debug[m]
    REAL(wp),  POINTER             :: clw_profile(:,:,:) ! clw profile debug[m]
    REAL(wp),  POINTER             :: tpot_inv(:,:)      ! potential temperature inversion strengt [m]
    REAL(wp),  POINTER             :: tpot_pbl(:,:)      ! potential temperature at pbl [m]
    REAL(wp),  POINTER             :: h2o_inv(:,:)       ! h2o inversion
    REAL(wp),  POINTER             :: h2o_pbl(:,:)       ! h2o mixin ration at pbl[kg]
    REAL(wp),  POINTER             :: pbl_num(:,:)       ! particle number at pbl[m3]
    REAL(wp),  POINTER             :: pbl_h(:,:)         ! pbl height [hPa] 
    REAL(wp),  POINTER             :: emu_w(:,:)         ! emulated updraft velocity [m/s]
    REAL(wp),  POINTER             :: tpot(:,:,:)        ! potential temperature[k]
    REAL(wp),  POINTER             :: q_total(:,:,:)     ! total water content[kg/kg]
    REAL(wp),  POINTER             :: emu_lwp(:,:)       ! lwp input for emulator
    REAL(wp),  POINTER             :: emu_cdnc(:,:)      ! cdnc input for emulator
    REAL(wp),  POINTER             :: ice(:,:,:)         ! ice content
    REAL(wp),  POINTER             :: emu_precip(:,:)    ! emulated rain production [kg/kg]
    REAL(wp),  POINTER             :: pxlb_unemul(:,:) ! Water no emulated [kg?]
    REAL(wp),  POINTER             :: emu_diag(:,:,:)    ! diagnostigs where emulator inputs are out of trainin data
    REAL(wp),  POINTER             :: emu_diag2(:,:,:)    ! diagnostigs values where emulator inputs are out of trainin data
    REAL(wp),  POINTER             :: emu_raut(:,:,:)     ! autoconversion
    REAL(wp),  POINTER             :: emu_rac1(:,:,:)     ! accretion 1
    REAL(wp),  POINTER             :: emu_rac2(:,:,:)     ! accretion 2
    REAL(wp),  POINTER             :: emu_precip3d(:,:,:) ! emulated precipitation in 3d
  END TYPE t_emulator
  ! emulator_stream=23
  ! emulatorvars is used to check namelist input for valid names
  INTEGER, PARAMETER           :: nemulatorvars=31
  CHARACTER(LEN=32)            :: emulatorvars(1:nemulatorvars)= &
                                (/'emu_mask'      ,&
				  'emu_mask_3d'   ,&
				  'autoconv_int'  ,&
				  'accretion1_int',&
				  'accretion2_int',&
				  'rainevap_int'  ,&
				  'cloud_top'     ,&
				  'hpalim700'     ,&
				  'clw_profile'   ,&
				  'cloud_base'    ,&
				  'tpot_inv'      ,&
                                  'tpot_pbl'      ,&
                                  'h2o_inv'       ,&
                                  'h2o_pbl'       ,&  
                                  'pbl_num'       ,& 
                                  'pbl_h'         ,&  
                                  'emu_w'         ,&
                                  'tpot'          ,&
                                  'q_total'       ,&
				  'emu_lwp'       ,&
				  'emu_cdnc'      ,&
                                  'ice'           ,&
                                  'cos_mu'        ,&
                                  'emu_precip'    ,&
                                  'pxlb_unemul'   ,&
				  'emu_diag'      ,&
				  'emu_diag2'     ,&
                                  'emu_raut'      ,&
                                  'emu_rac1'      ,& 
                                  'emu_rac2'      ,&
                                  'emu_precip3d'  /)


  ! variable pointers
  TYPE(t_emulator)               :: emulators
  PUBLIC :: emulator,initEmulatorStream,train,emulators,initEmulator
CONTAINS

  SUBROUTINE train(gpemu,ind)
    ! This subroutine reads in a pre-trained emulator and the used training data.

    INTEGER,INTENT(IN)                                  :: ind
    CLASS(BaseGP), ALLOCATABLE, INTENT(OUT)             :: gpemu


    iF (.NOT. ALLOCATED(gpemu)) allocate(gpemu, source = DenseGP(ind))

  END SUBROUTINE train


  SUBROUTINE  readTraining(train_filename,train_n,train_dim,train_in,train_out)

    CHARACTER(len=*), INTENT(IN)                        :: train_filename
    INTEGER                                             :: uu, train_n, i, train_dim
    INTEGER, DIMENSION(train_n)                         :: a
    REAL(wp), DIMENSION(train_n,train_dim), INTENT(OUT) :: train_in
    REAL(wp), DIMENSION(train_n), INTENT(OUT)           :: train_out
    ! Read in emulator training data. These are used only in the
    ! standardisation and only their means are used.
    open(newunit=uu, file=train_filename)
    read (uu,*) (train_in(i,1:train_dim), a(i), train_out(i), i=1,train_n)
    close(uu)

  END SUBROUTINE readTraining

  SUBROUTINE initEmulator()
    !This routine initialize emulators, reads emulator files
       CALL initDenseGP('gp_w_night.out',1)
       CALL initDenseGP('gp_w_day.out',2)
       CALL initDenseGP('gp_precip_night_v2.out',3)
       CALL initDenseGP('gp_precip_day_v2.out',4) 
    IF (p_parallel_io) THEN   
       CALL readTraining('DATA_w_night',500,6,ex1_w_n,t1_w_n)
       CALL readTraining('DATA_w_day',497,7,ex1_w_d,t1_w_d)
       CALL readTraining('DATA_precip_night_v2',500,6,ex1_precip_n,t1_precip_n)
       CALL readTraining('DATA_precip_day_v2',500,7,ex1_precip_d,t1_precip_d)
    END IF
    IF (p_parallel) THEN
       CALL p_bcast (ex1_w_n, p_io)
       CALL p_bcast (t1_w_n, p_io)
       CALL p_bcast (ex1_w_d, p_io)
       CALL p_bcast (t1_w_d, p_io)
       CALL p_bcast (ex1_precip_n, p_io)
       CALL p_bcast (t1_precip_n, p_io)
       CALL p_bcast (ex1_precip_d, p_io)
       CALL p_bcast (t1_precip_d, p_io)
    END IF 

  END SUBROUTINE initEmulator

  SUBROUTINE initEmulatorStream
    
    USE mo_control,             ONLY: nlev
    USE mo_submodel_streams,    ONLY: emulator_lpost, emulator_tinterval, emulatornam
    USE mo_string_utls,         ONLY: st1_in_st2_proof
    USE mo_util_string,         ONLY: tolower
    USE mo_exception,           ONLY: finish
    USE mo_memory_base,         ONLY: t_stream, new_stream, &
                                      default_stream_setting, &
                                      add_stream_element, &
                                      AUTO, BELOWSUR
    ! local variables
    INTEGER, PARAMETER             :: ndefault = 31
    CHARACTER(LEN=32)            :: defnam(1:nemulatorvars)= &
                                (/'emu_mask'      ,&
				  'emu_mask_3d'   ,&
				  'autoconv_int'  ,&
				  'accretion1_int',&
				  'accretion2_int',&
				  'rainevap_int'  ,&
				  'cloud_top'     ,&
				  'hpalim700'     ,&
				  'clw_profile'   ,&
				  'cloud_base'    ,&
				  'tpot_inv'      ,&
                                  'tpot_pbl'      ,&
                                  'h2o_inv'       ,&
                                  'h2o_pbl'       ,&  
                                  'pbl_num'       ,& 
                                  'pbl_h'         ,&  
                                  'emu_w'         ,&
                                  'tpot'          ,&
                                  'q_total'       ,&
				  'emu_lwp'       ,&
				  'emu_cdnc'      ,&
                                  'ice'           ,&
                                  'cos_mu'        ,&
                                  'emu_precip'    ,&
                                  'pxlb_unemul'   ,&
				  'emu_diag'      ,&
				  'emu_diag2'     ,&
                                  'emu_raut'      ,&
                                  'emu_rac1'      ,&
                                  'emu_rac2'      ,&
                                  'emu_precip3d'  /)
 
    TYPE (t_stream), POINTER       :: emulator_stream
    INTEGER                        :: ierr
    LOGICAL                        :: lpost
    CHARACTER(LEN=32) :: file
    INTEGER :: iret, ncid, DimId, VarId, xdmy
    INTEGER :: inml, iunit
    !-- handle ALL and DEFAULT options
    IF (TRIM(tolower(emulatornam(1))) == 'all')     emulatornam(1:nemulatorvars) = emulatorvars(:)
    IF (TRIM(tolower(emulatornam(1))) == 'default') emulatornam(1:ndefault) = defnam(:)

    !-- check that all variable names from namelist are valid
    IF (.NOT. st1_in_st2_proof( emulatornam, emulatorvars, ierr=ierr) ) THEN
      IF (ierr > 0) CALL finish ( 'ini_emulator_stream', 'variable '// &
                                  emulatornam(ierr)//' does not exist in emulator stream' )
    END IF

    !-- open new stream
    !SF #383: registering emulator to rerun storage
    CALL new_stream (emulator_stream,'emulator',lpost=emulator_lpost,lrerun=.TRUE., &
         interval=emulator_tinterval)
    CALL default_stream_setting (emulator_stream, lrerun = .TRUE., &
         contnorest = .TRUE., table = 199, &
         laccu = .false., code = AUTO)

    !-- add individual variables to stream
    lpost = st1_in_st2_proof( 'emu_mask', emulatornam)
    CALL add_stream_element (emulator_stream, 'emu_mask', emulators%emu_mask, &
         longname = 'Mask where emulator is used', &
         units = '', lpost = lpost)

    lpost = st1_in_st2_proof( 'emu_mask_3d', emulatornam)
    CALL add_stream_element (emulator_stream, 'emu_mask_3d', emulators%emu_mask_3d, &
         longname = 'Mask where emulator is used vertical mask included', &
         units = '', lpost = lpost)


    lpost = st1_in_st2_proof( 'emu_diag', emulatornam)
    CALL add_stream_element (emulator_stream, 'emu_diag1', emulators%emu_diag, &
         longname = 'diagnosis where input variaibles is out of trainign set', &
         units = '', lpost = lpost)

    lpost = st1_in_st2_proof( 'emu_diag2', emulatornam)
    CALL add_stream_element (emulator_stream, 'emu_diag2', emulators%emu_diag2, &
         longname = 'diagnosis, values of diag1 points', &
         units = '', lpost = lpost)


    lpost = st1_in_st2_proof( 'hpalim700', emulatornam)
    CALL add_stream_element (emulator_stream, 'hpalim700', emulators%hpalim700, &
         longname = 'hpalim700', &
         units = '', lpost = .TRUE.)

    lpost = st1_in_st2_proof( 'emu_lwp', emulatornam)
    CALL add_stream_element (emulator_stream, 'emu_lwp', emulators%emu_lwp, &
         longname = 'emu_lwp', &
         units = '', lpost = .TRUE.)

    lpost = st1_in_st2_proof( 'emu_cdnc', emulatornam)
    CALL add_stream_element (emulator_stream, 'emu_cdnc', emulators%emu_cdnc, &
         longname = 'emu_cdnc', &
         units = '', lpost = .TRUE.)

    lpost = st1_in_st2_proof( 'clw_profile', emulatornam)
    CALL add_stream_element (emulator_stream, 'clw_profile', emulators%clw_profile, &
         longname = 'clw_profile', &
         units = 'K', lpost = .TRUE.,LACCU=.false. )


    CALL add_stream_element (emulator_stream, 'autoconv_int', emulators%autoconv_int, &
         longname = 'autoconversion', &
         units = '', lpost = .TRUE.,laccu=.TRUE.)

    CALL add_stream_element (emulator_stream, 'accretion1_int', emulators%accretion1_int, &
         longname = 'accreation term 1', &
         units = '', lpost = .TRUE.,laccu=.TRUE.)

    CALL add_stream_element (emulator_stream, 'accretion2_int', emulators%accretion2_int, &
         longname = 'accreation term 2', &
         units = '', lpost = .TRUE.,laccu=.TRUE.)

    CALL add_stream_element (emulator_stream, 'rainevap_int', emulators%rainevap_int, &
         longname = 'Rain evaporations', &
         units = '', lpost = .TRUE.,laccu=.TRUE.)

    lpost = st1_in_st2_proof( 'tpot_inv', emulatornam)
    CALL add_stream_element (emulator_stream, 'tpot_inv', emulators%tpot_inv, &
         longname = 'Potential temperature inversion', &
         units = 'K', lpost = lpost)
    
    lpost = st1_in_st2_proof( 'tpot_pbl', emulatornam)
    CALL add_stream_element (emulator_stream, 'tpot_pbl', emulators%tpot_pbl, &   
         longname = 'potential temperatuer at boundary layer', &
         units = 'K', lpost = lpost)
    
    lpost = st1_in_st2_proof( 'h2o_inv', emulatornam)
    CALL add_stream_element (emulator_stream, 'h2o_inv', emulators%h2o_inv, &
         longname = 'h2o inversion', &
         units = 'kg', lpost = lpost)
    
    lpost = st1_in_st2_proof( 'h2o_pbl', emulatornam)
    CALL add_stream_element (emulator_stream, 'h2o_pbl', emulators%h2o_pbl, &
         longname = 'h2o in plb', &
         units = 'm3', lpost = lpost)

    lpost = st1_in_st2_proof( 'cloud_top', emulatornam)
    CALL add_stream_element (emulator_stream, 'cloud_top', emulators%cloud_top, &
         longname = 'Cloud top level', &
         units = 'level', lpost = lpost)

    lpost = st1_in_st2_proof( 'cloud_base', emulatornam)
    CALL add_stream_element (emulator_stream, 'cloud_base', emulators%cloud_base, &
         longname = 'Cloud base level', &
         units = 'level', lpost = lpost)

    lpost = st1_in_st2_proof( 'pbl_num', emulatornam)
    CALL add_stream_element (emulator_stream, 'pbl_num', emulators%pbl_num, &
         longname = 'number of particles in pbl', &
         units = 'kg kg-1', lpost=lpost)

    lpost = st1_in_st2_proof( 'pbl_h', emulatornam)
    CALL add_stream_element (emulator_stream, 'pbl_h', emulators%pbl_h, &
         longname = 'pbl height', &
         units = 'hPa', lpost=lpost)

    lpost = st1_in_st2_proof( 'emu_w', emulatornam)
    CALL add_stream_element (emulator_stream, 'emu_w', emulators%emu_w, &
         longname = 'updraft velocity', &
         units = 'm/s', lpost=lpost)

    CALL add_stream_element (emulator_stream, 'tpot', emulators%tpot, &
         longname = 'potential temperature', &
         units = 'K', lpost=lpost)
    
    CALL add_stream_element (emulator_stream, 'q_total', emulators%q_total, &
         longname = 'total water content', &
         units = 'kg/kg', lpost=lpost)

    CALL add_stream_element (emulator_stream, 'ice', emulators%ice, &
         longname = 'ice content', &
         units = 'kg/kg', lpost=lpost)

    CALL add_stream_element (emulator_stream, 'cos_mu', cos_mu, &
         longname = 'Cosine of solar zenith angle', &
         units = ' ', lpost=lpost)

    CALL add_stream_element (emulator_stream, 'emu_precip', emulators%emu_precip, &
         longname = 'Emulated precipitation production (autoconversion+accretion)', &
         units = 'kg/kg', lpost=lpost)

    CALL add_stream_element (emulator_stream, 'pxlb_unemul', emulators%pxlb_unemul, &
         longname = 'Amount of unemulated water', &
         units = 'kg', lpost=lpost)

    CALL add_stream_element (emulator_stream, 'emu_raut', emulators%emu_raut, &
         longname = 'Autoconversion after emulation', &
         units = 'kg/kg', lpost=lpost)

    CALL add_stream_element (emulator_stream, 'emu_rac1', emulators%emu_rac1, &
         longname = 'First accretion after emulation', &
         units = 'kg/kg', lpost=lpost)

    CALL add_stream_element (emulator_stream, 'emu_rac2', emulators%emu_rac2, &
         longname = 'Second accretion after emulation', &
         units = 'kg/kg', lpost=lpost)

    CALL add_stream_element (emulator_stream, 'emu_precip3d', emulators%emu_precip3d, &
         longname = 'Emulated precipitation production (autoconversion+accretion), now in 3D!', &
         units = 'kg/kg', lpost=lpost)


  END SUBROUTINE initEmulatorStream

!>
!! This subroutine calculates vertical interacls from ground to certain level
!! used to calculate lwp inside boundary layer
  SUBROUTINE calculate_lwp(x,klev,cloud_top,zdpg,integral) 
       INTEGER, INTENT(IN)    :: cloud_top, klev
       REAL(wp),INTENT(IN)   ::  x(klev),zdpg(klev)
       REAL(wp),INTENT(OUT)   :: integral
       REAL(wp)    :: zintegral
       INTEGER                :: jk
       zintegral = 0
       DO 100 jk = cloud_top,klev
          zintegral  = zintegral  + x(jk) *zdpg(jk)
       100 END DO
       integral=zintegral
  END SUBROUTINE calculate_lwp



  SUBROUTINE find_cloud2(lev,indhi,x,cloud_top,cloud_base)
! Input: lev   = number of layers
!        indhi = index of the highest layer considered 
!                (e.g., the 700 hPa level)
!        x(lev) = cloud liquid water content (g/kg)
! Output: cloud_top  = full-level index for the lowermost layer of the lowest cloud 
!         cloud_base = full-level index for the uppermost layer of the lowest cloud 
   
!    INTEGER, PARAMETER :: dp = selected_real_kind(15, 307)
! This would be better:
    USE mo_kind, ONLY : dp
    INTEGER, INTENT(IN)   :: lev, indhi
    REAL(dp), INTENT(in)  :: x(lev)
    INTEGER, INTENT(out)  :: cloud_top, cloud_base

    INTEGER :: ilev

! Initialize to undefined
    cloud_top=-999
    cloud_base=-999
    DO ilev=lev,indhi,-1
      IF (x(ilev) >= 0.01) THEN
        cloud_base=ilev
        EXIT
      END IF
    ENDDO 
    IF (cloud_base >= indhi) THEN
      DO ilev=cloud_base,indhi,-1
        IF (x(ilev) < 0.01) EXIT
      ENDDO 
      cloud_top=MAX(indhi,ilev+1)
    END IF

    RETURN
  END SUBROUTINE find_cloud2

  SUBROUTINE find_cloud(x,lev,cloud_top,cloud_base)
        integer, parameter :: dp = selected_real_kind(15, 307)
	INTEGER, INTENT(IN)   :: lev
        REAL(dp),INTENT(IN)   :: x(lev)
        REAL(dp)     :: x1(lev)
	INTEGER, INTENT(OUT)  :: cloud_base,cloud_top
        INTEGER               :: i,ind
        cloud_top=0
        cloud_base=0
	!write(*,*) x,'find cloud profile'
	x1=x
	x1 = x1(size(x1):1:-1) !reverse order
	WHERE(x1 >= 0.01) x1 = 999
	WHERE(x1 < 0.01) x1=-999
	cloud_base =  SIZE(x1,1)-MAXLOC(x1,1)
        cloud_top =  SIZE(x1,1)-(MINLOC(x1(MAXLOC(x1,1):SIZE(x1,1)),1)-1+MAXLOC(x1,1)-1)-1 !change to calculate right could top
  END SUBROUTINE find_cloud

!>
!! This subroutine computes mask where emulator is applied for each timestep.
!! If region is inside training space emulator is also applied
!! retruns maks and emulated value of what ever we are emulating.
  SUBROUTINE emulator(kproma, kbdim, ktdia, klev, klevp1, krow,ktrac, &
                      pfull, phalf, &
		      pxlm1,pxim1,pqm1,pcdnc)
  USE mo_memory_g3b,    ONLY: aps, slm, seaice, tpot
! Input for the emulator 

    INTEGER, INTENT(IN)    :: kbdim, klevp1, klev, kproma, ktdia,krow,ktrac
    !INTEGER                :: id_cdnc

    REAL(wp), INTENT(in) :: &
      pfull(kbdim,klev),    & ! Full-level pressure (Pa)  
      phalf(kbdim,klevp1),  & ! Half-level pressure (Pa)  
      pxlm1(kbdim,klev),    & ! cloud liquid water content 
      pxim1(kbdim,klev),    & ! cloud ice content 
      pqm1(kbdim,klev),     & ! specific humidity
      pcdnc(kbdim,klev)       !< CDNC

    !Local variables
    INTEGER                :: zfull700(kbdim),zhalf700(kbdim),id_cdnc ! index for nearest level at 700hpa level. 

    INTEGER                :: jk,        &
                              cloud_top, &                      ! index for top of the lowest cloud layer
                              cloud_base,inv_ind_hi,inv_ind_lo                      ! index for basa of lowest cloud layer,index's for inversion calculations
    REAL(wp)                   :: zlwp700(kbdim),zlwp,ziwp,ziwp700, &            !LWP inside 700hpa layer,total lwp and iwp
                              pdp(kbdim,klev), pdpg(kbdim,klev), &    !pressure difference
                              g_rcp,zx(klev),zx1(klevp1),ztotalq(klev)           !zx is just temporary vector,ztotalq is total water

    REAL(wp) :: meanlt, stdlt
    REAL(wp), DIMENSION(:), ALLOCATABLE :: emuInputVec,emuInputVecStandard
!    REAL(dp), DIMENSION(:), ALLOCATABLE :: lt
    INTEGER :: n

    ALLOCATE(REAL(wp) ::emuInputVec(7))
    ALLOCATE(REAL(wp) ::emuInputVecStandard(7))

   g_rcp = 1._wp / grav
   pdp(1:kproma,ktdia:klev)           = phalf(1:kproma,ktdia+1:klevp1) - phalf(1:kproma,ktdia:klev)
   pdpg(1:kproma,ktdia:klev)          = g_rcp * pdp(1:kproma,ktdia:klev)!*100
   emulators%tpot(:,:,krow) = tpot(:,:,krow)
   emulators%q_total(:,:,krow) = pxlm1+pqm1
   emulators%ice(:,:,krow) = pxim1
   emulators%clw_profile(:,:,krow) = pxlm1(:,:)

    !>>jpk
    ! Train emulator
    IF (lstart .OR. lresume) THEN
         CALL train(gp_emu_w_n,1)
         CALL train(gp_emu_w_d,2)
         CALL train(gp_emu_precip_n,3)
         CALL train(gp_emu_precip_d,4)
     END IF
!    allocate(lt(501))
    !<<jpk


    !emulator is only trained low level clouds (clouds below 700hpa level). 1. step find index for level 700hpa both for half levels and full levels
    zfull700 = -1
    zhalf700 = -1
    DO 101 jk = 1,kproma
	   zx =pfull(jk,:)
     	   WHERE(zx .GE. 70000) zx = 999
	   zx1 =phalf(jk,:)
     	   WHERE(zx1 .GE. 70000) zx1 = 999
           zfull700(jk) = MINLOC(MERGE(0,1,zx == 999),DIM=1)
           zhalf700(jk) = MINLOC(MERGE(0,1,zx1 == 999),DIM=1)
    101 END DO

    !loop over grid points
    !emulators%emu_mask_3d(:,:,1:klev) = 0.0_wp
    emulators%emu_diag(1:kproma,1:klev,krow) = 0.0_wp
    emulators%emu_diag2(1:kproma,1:klev,krow) = 0.0_wp
    emulators%emu_precip3d(1:kproma,1:klev,krow) = 0.0_wp
    DO 102 jk =1,kproma
	    emulators%emu_mask(jk,krow) = 0.0_wp
	    emulators%emu_mask_3d(jk,1:klev,krow) = 0.0_wp
	    emulators%emu_lwp(jk,krow) = -999
            emulators%pxlb_unemul(jk,krow)=0.0_wp
            emulators%emu_w(jk,krow) = -999
            emulators%emu_precip(jk,krow) = -999
            emulators%emu_raut(jk,:,krow) = -999
            emulators%emu_rac1(jk,:,krow) = -999
            emulators%emu_rac2(jk,:,krow) = -999
	    !2.eliminaite points over ice and land, 
	    IF((slm(jk,krow)  .GE. 1.0_wp) .OR. seaice(jk,krow) .GE. 0.01_wp) THEN
		emulators%tpot_inv(jk,krow) = -999
		emulators%h2o_inv(jk,krow) = -999
		emulators%h2o_pbl(jk,krow) = -999
		emulators%tpot_pbl(jk,krow) = -999
		emulators%pbl_h(jk,krow) = -999
    	    	emulators%emu_mask(jk,krow) = 0.0_wp
		emulators%cloud_top(jk,krow) = -999
		emulators%cloud_base(jk,krow) = -999
		!emulators%clw_profile(jk,:,krow) = -999
		emulators%hpalim700(jk,krow) = -999
                !emulators%tpot(jk,:,krow) = -999
                !emulators%q_total(jk,:,krow) = -999
                !emulators%ice(jk,:,krow) = -999
                emulators%emu_w(jk,krow) = -999
                emulators%emu_precip(jk,krow) = -999
                CYCLE
	    END IF
	    !Eliminate points where there is now low level clouds There is cloud if there is more than 0.01 g/kg of cloud water below 700hPa
	    IF(COUNT(pxlm1(jk,zfull700(jk):klev)*1000 .GE. 0.01) < 1.0_wp) THEN
		emulators%tpot_inv(jk,krow) = -999
		emulators%h2o_inv(jk,krow) = -999
		emulators%h2o_pbl(jk,krow) = -999
		emulators%tpot_pbl(jk,krow) = -999
		emulators%pbl_h(jk,krow) = -999
    	    	emulators%emu_mask(jk,krow) = 0.0_wp
		emulators%cloud_top(jk,krow) = -999
		emulators%cloud_base(jk,krow) = -999
		!emulators%clw_profile(jk,:,krow) = -999
		emulators%hpalim700(jk,krow) = -999
                !emulators%tpot(jk,:,krow) = -999
                !emulators%q_total(jk,:,krow) = -999
                !emulators%ice(jk,:,krow) = -999
                emulators%emu_w(jk,krow) = -999
                emulators%emu_precip(jk,krow) = -999
                CYCLE
	    END IF
 	    !3. Eliminate all points where there is fog, check if there is cloud water on lowest layer and exluce this points
	    IF(pxlm1(jk,klev)*1000 .GE. 0.01) THEN
		emulators%tpot_inv(jk,krow) = -999
		emulators%h2o_inv(jk,krow) = -999
		emulators%h2o_pbl(jk,krow) = -999
		emulators%tpot_pbl(jk,krow) =- 999
		emulators%pbl_h(jk,krow) = -999
    	    	emulators%emu_mask(jk,krow) = 0.0_wp
		emulators%cloud_top(jk,krow) = -999
		emulators%cloud_base(jk,krow) = -999
		!emulators%clw_profile(jk,:,krow) = -999
		emulators%hpalim700(jk,krow) = -999
                !emulators%tpot(jk,:,krow) = -999
                !emulators%q_total(jk,:,krow) = -999
                !emulators%ice(jk,:,krow) = -999
                emulators%emu_w(jk,krow) = -999
                emulators%emu_precip(jk,krow) = -999
                CYCLE
	    END IF
	    !4. step, locate lowest cloud
	    !find_cloud2(lev,indhi,x,cloud_top,cloud_base)
 	    !CALL find_cloud(pxlm1(jk,zfull700(jk):klev)*1000,klev-zfull700(jk),cloud_top,cloud_base)

            if (emul_cloud == 2) THEN 	
     	    	CALL find_cloud2(klev,zfull700(jk),pxlm1(jk,:)*1000,cloud_top,cloud_base)
	    END IF
            if (emul_cloud == 1) THEN 	
     	    	CALL find_cloud(pxlm1(jk,zfull700(jk):klev)*1000,klev-zfull700(jk),cloud_top,cloud_base)
	    	cloud_top = cloud_top+zfull700(jk)
	    	cloud_base = cloud_base+zfull700(jk)
	    END IF
	    !write(*,*) 'cloud base top',emul_cloud,cloud_top,cloud_base
 	    !5. step calculate LWP inside lowest layer layer, this is used to identify that most of the cloud water in column is in low level cloud
	    CALL calculate_lwp(pxlm1(jk,:)*1000,klev,cloud_top,pdpg(jk,:),zlwp700(jk))
	    CALL calculate_lwp(pxlm1(jk,:)*1000,klev,1,pdpg(jk,:),zlwp) 
	    CALL calculate_lwp(pxim1(jk,:)*1000,klev,cloud_top,pdpg(jk,:),ziwp700)
	    CALL calculate_lwp(pxim1(jk,:)*1000,klev,1,pdpg(jk,:),ziwp) 
	    !check that lowest cloud is below 700hPa
	    !write(*,*) 'CLOUDTOP',cloud_top, '700hpa',zfull700(jk), 'cloud_base',cloud_base
	    !write(*,*) 'LWP700',zlwp700(jk), 'TOTAL LWP+IWP:',zlwp+ziwp,'IWP: ',ziwp700
		!write(*,*) 'LWP700',zlwp700(jk), 'TOTAL LWP+IWP:',(zlwp+ziwp)*0.5,'IWP: ',ziwp700, cloud_top,cloud_base,klev
	        IF((zlwp700(jk) > ((zlwp+ziwp)*0.5)) .AND. ( 0.1*zlwp700(jk)  > ziwp700)) THEN
		        !CALCULATE Emulator inputs
			inv_ind_hi = cloud_top-2
			ztotalq = pxlm1(jk,:)+pqm1(jk,:)
			inv_ind_lo = MIN(cloud_base+2,klev)
			!write(*,*) 'inversion indx, hi: ',inv_ind_hi,' lo: ',inv_ind_lo
			emulators%tpot_inv(jk,krow) = MAXVAL(tpot(jk,inv_ind_hi:inv_ind_lo,krow))-MINVAL(tpot(jk,inv_ind_hi:inv_ind_lo,krow))
			emulators%h2o_inv(jk,krow) = 1000.0*(MAXVAL(ztotalq(inv_ind_hi:inv_ind_lo))-MINVAL(ztotalq(inv_ind_hi:inv_ind_lo)))
		 	emulators%h2o_pbl(jk,krow) = MINVAL(ztotalq(inv_ind_hi:inv_ind_lo))
			emulators%tpot_pbl(jk,krow) = MINVAL(tpot(jk,inv_ind_hi:inv_ind_lo,krow))
		        emulators%pbl_h(jk,krow) = 0.01*(aps(jk,krow)-phalf(jk,cloud_top))
    	    	        emulators%emu_mask(jk,krow) = 1.0_wp
			emulators%cloud_top(jk,krow) = cloud_top
			emulators%cloud_base(jk,krow) = cloud_base
			emulators%emu_mask_3d(jk,cloud_top:cloud_base,krow) = 1.0_wp
			emulators%hpalim700(jk,krow) = REAL(zfull700(jk))
	   		emulators%emu_lwp(jk,krow) = zlwp700(jk)
      ! JPK, 28.10.2019, added +1 to the last argument to prevend division by zero. KN 10.2.2020, change to accorgin to PR
      ! Not sure if cloud calculation is correct.
			!emulators%emu_cdnc(jk,krow) = 0.000001*mean(cdnc(jk,cloud_top-1:cloud_base),MAX(cloud_base-(cloud_top-1),1))
			emulators%emu_cdnc(jk,krow) = 0.000001*mean(pcdnc(jk,cloud_top:cloud_base),cloud_base-cloud_top+1)
                        !>> jpk
                        ! Emulate updraft velocity (JPK, 2019.4.23)
                        ! Emulate precipitation production (JPK, 2019.10.4)   
                        ! Set emulator input
                        emuInputVec(1)=emulators%h2o_inv(jk,krow)
                        emuInputVec(2)=emulators%tpot_inv(jk,krow)
                        emuInputVec(3)=emulators%emu_lwp(jk,krow)
                        emuInputVec(4)=emulators%tpot_pbl(jk,krow)
                        emuInputVec(5)=emulators%pbl_h(jk,krow)
			emuInputVec(6)=emulators%emu_cdnc(jk,krow)
			emuInputVec(7)=cos_mu(jk,krow) !solar zenith angle

                        ! Depending on the solar zenith angle, call day or night emulator
                        IF (cos_mu(jk,krow)<=0.0) THEN
                           ! Emulate updraft
                           meanlt = mean(t1_w_n,500)
                           stdlt  = std(t1_w_n,meanlt,500)
                           CALL standardize_emulator_inputs(emuInputVec(1:6),emuInputVecStandard(1:6),500,6,ex1_w_n,jk,krow) 

			   IF(ANY(emuInputVecStandard(1:6)==-999)) THEN
			   	emulators%emu_mask_3d(jk,cloud_top:cloud_base,krow) = 0._wp
    	    	       	   	emulators%emu_mask(jk,krow) = 0.0_wp
			   ELSE
                           	emulators%emu_w(jk,krow) = unstandardize_s(gp_emu_w_n%predict(emuInputVecStandard(1:6),0),meanlt,stdlt)
			   END IF		


                           ! Emulate precipitation production
                           meanlt = mean(t1_precip_n,500)
                           stdlt  = std(t1_precip_n,meanlt,500)
                           CALL standardize_emulator_inputs(emuInputVec(1:6),emuInputVecStandard(1:6),500,6,ex1_precip_n,jk,krow) 

			   IF(ANY(emuInputVecStandard(1:6)==-999)) THEN
			   	emulators%emu_mask_3d(jk,cloud_top:cloud_base,krow) = 0._wp
    	    	       	   	emulators%emu_mask(jk,krow) = 0.0_wp
			   ELSE
                         	emulators%emu_precip(jk,krow) = unstandardize_s(gp_emu_precip_n%predict(emuInputVecStandard(1:6),0),meanlt,stdlt)
			   	IF(emulators%emu_precip(jk,krow)==-999) THEN
    	    	       	   		emulators%emu_mask(jk,krow) = 0.0_wp
				END IF
 			   END IF

                        ELSE
                           ! Emulate updraft
                           meanlt = mean(t1_w_d,497)
                           stdlt  = std(t1_w_d,meanlt,497)
                           CALL standardize_emulator_inputs(emuInputVec(:),emuInputVecStandard(:),497,7,ex1_w_d,jk,krow)

			   IF(ANY(emuInputVecStandard(:)==-999)) THEN
			   	emulators%emu_mask_3d(jk,cloud_top:cloud_base,krow) = 0._wp
    	    	       	  	emulators%emu_mask(jk,krow) = 0.0_wp
			   ELSE
                           	emulators%emu_w(jk,krow) = unstandardize_s(gp_emu_w_d%predict(emuInputVecStandard(:),0),meanlt,stdlt)

			   END IF


                           ! Emulate precipitation production
                           meanlt = mean(t1_precip_d,500)
                           stdlt  = std(t1_precip_d,meanlt,500)
                           CALL standardize_emulator_inputs(emuInputVec(:),emuInputVecStandard(:),500,7,ex1_precip_d,jk,krow)
			   IF(ANY(emuInputVecStandard(:)==-999)) THEN
			   	emulators%emu_mask_3d(jk,cloud_top:cloud_base,krow) = 0._wp
    	    	       	   	emulators%emu_mask(jk,krow) = 0.0_wp
			   ELSE
                         	emulators%emu_precip(jk,krow) = unstandardize_s(gp_emu_precip_d%predict(emuInputVecStandard(:),0),meanlt,stdlt)
			   	IF(emulators%emu_precip(jk,krow)==-999) THEN
    	    	       	   		emulators%emu_mask(jk,krow) = 0.0_wp
				END IF
			   END IF

                        END IF


		END IF
    102 END DO 

      !DEALLOCATE (gp_emu_w_d)
      !DEALLOCATE (gp_emu_w_n)
      !DEALLOCATE(t1_w_n)
      !DEALLOCATE(t1_w_d)
      !DEALLOCATE(ex1_w_n)
      !DEALLOCATE(ex1_w_d)
      !DEALLOCATE(a_n)
      !DEALLOCATE(a_d)
!      DEALLOCATE(lt)
      DEALLOCATE(emuInputVec)

   END SUBROUTINE emulator

   SUBROUTINE standardize_emulator_inputs(emu_input_vec,emu_input_vec_out,train_n,train_dim,train_out,jk,krow)
     ! Standardise emulator inputs.
     INTEGER, INTENT(IN)                                  :: train_n, train_dim,jk,krow
     REAL(wp), DIMENSION(train_dim), INTENT(IN)           :: emu_input_vec
     REAL(wp), DIMENSION(train_dim), INTENT(out)           :: emu_input_vec_out
     REAL(wp)                                             :: stdx, meanx,min_,max_
     REAL(wp), DIMENSION(train_n,train_dim), INTENT(in)   :: train_out
     INTEGER                                              :: jj
     LOGICAL			                          :: test

     test =.TRUE.

     DO jj = 1,train_dim
        meanx = mean(train_out(:,jj),train_n)
        stdx  = std(train_out(:,jj),meanx,train_n)
	max_ = MAXVAL(train_out(:,jj))
	min_ = MINVAL(train_out(:,jj))
	IF ((emu_input_vec(jj) >= min_) .AND. (emu_input_vec(jj) <= max_)) THEN
        	emu_input_vec_out(jj) = standardize(emu_input_vec(jj),meanx,stdx)
	ELSE
		!write(*,*) jj,'bad value',emu_input_vec(jj),min_,max_
    		emulators%emu_diag2(jk,jj,krow) = emu_input_vec(jj)
        	emu_input_vec_out(jj) = -999
    		emulators%emu_diag(jk,jj,krow) = 1.0_wp	
		
	END IF
     END DO
     
    !IF(test) THEN
    !	 write(*,*), 'ok vector'
    !ELSE
    !	write(*,*) 'bad vector'
    !END IF

  END SUBROUTINE  standardize_emulator_inputs

   !>> jpk 
   ! Some helper functions 
   ! Copied from Kalle's emulator implementation

  FUNCTION standardize(x,meanx,stdx) RESULT(res)
    REAL(wp) :: x
    REAL(wp) :: res
    REAL(wp) :: meanx 
    REAL(wp) :: stdx

    res = (x - meanx)/(stdx)
  END FUNCTION standardize

   FUNCTION unstandardize_s(x,meanx,stdx) RESULT(res)
     REAL(dp) x
     REAL(dp) :: res
     REAL(dp) :: meanx 
     REAL(dp) :: stdx

     res = (x * stdx) + meanx
   END FUNCTION unstandardize_s
   
   FUNCTION mean(x,dmn) RESULT(res)
     INTEGER dmn
     REAL(dp) x(dmn)
     REAL(dp) :: res

     res = SUM(x)/dmn
   END FUNCTION mean

   FUNCTION std(x,meanx,dmn) RESULT(res)
     INTEGER :: dmn
     REAL(dp) :: x(dmn)
     REAL(dp) :: meanx
     REAL(dp) :: res

     res = SQRT(SUM((x - meanx)**2)/dmn)
   END FUNCTION std
   !<<

END MODULE mo_emulator
