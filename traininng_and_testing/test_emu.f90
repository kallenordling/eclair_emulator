!! Eclair emulator module tests

PROGRAM emulator
    use m_util, only : dp
    USE m_gp, ONLY : BaseGP
    IMPLICIT NONE
    ! Emulator training data
    REAL(dp) :: ex2(500,7), t2(500)
    INTEGER :: obs2(500)
    ! Emulator
    CLASS(BaseGP),ALLOCATABLE :: gp
    !
    ! Tests
    ! ====
    ! a) Compare previously trained emulators with training data
    !CALL test_emulator('./emulators/DATA_precip_day_v2','emulators/out.gp.day.precip_v2','res_precip_day.dat','')
    !CALL test_emulator('./emulators/DATA_precip_night_v2','emulators/out.gp.night.precip_v2','res_precip_night.dat','')
    !CALL test_emulator('./emulators/DATA_w_night','emulators/out.gp.night.w_v2','res_w_night.dat','')
    !CALL test_emulator('./emulators/DATA_w_day','emulators/out.gp.day.w_v2','res_w_day.dat','')
    !
    ! b) Train emulator first (and save emulator) and then test
    !CALL test_emulator('DATA_day_w',' ','test_day_w.dat','gp_w_day.out')
    ! ...
    !
    ! c) Leave-one-out tests
    CALL leave_one_out('./data_day.precip/DATA_new','lou_day_precip.dat')
    CALL leave_one_out('./data_night.precip/DATA_new','lou_night_precip.dat')
    CALL leave_one_out('./data_day.w/DATA','lou_day_w.dat')
    CALL leave_one_out('./data_night.w/DATA','lou_night_w.dat')
    ! ...

CONTAINS

    SUBROUTINE test_emulator(dataname,emuname,outname,outemu)
        ! Compare existing or new emulator against training data
        use m_gp_dense
        ! Inputs
        CHARACTER(*), INTENT(IN) :: dataname,emuname ! Training data and emulator files
        CHARACTER(*), INTENT(IN) :: outname ! Output file name
        CHARACTER(*), INTENT(IN) :: outemu  ! Output emulator name
        ! Local variables
        REAL(dp), ALLOCATABLE :: emuInputVec(:)
        INTEGER :: i, j, n, m, obs, iout
        REAL(dp) :: meanx, stdx, meant, stdt
        REAL(dp) :: pred, tmp(500)
        !
        ! Read training data: ex2(500,7), t2(500) and obs2(500)
        CALL read_training(dataname,n,m)
        !
        ! Allocate data
        ALLOCATE(emuInputVec(m))
        !
        ! Standardize training data
        meant = mean(t2(1:n),n)
        stdt  = std(t2(1:n),meant,n)
        DO j=1,m ! cos_mu included
            meanx = mean(ex2(1:n,j),n)
            stdx  = std(ex2(1:n,j),meanx,n)
            ex2(1:n,j)=(ex2(1:n,j) - meanx)/stdx
        ENDDO
        !
        IF (LEN_TRIM(emuname)>1) THEN
            ! Read the emulator
            IF (allocated(gp)) DEALLOCATE(gp)
            allocate(gp, source = DenseGP(emuname))
        ELSE
            ! Train emulator
            tmp(1:n)=(t2(1:n)-meant)/stdt
            CALL train_emulator(n, m, ex2(1:n,1:m), obs2(1:n), tmp(1:n))
        ENDIF
        !
        ! Output
        iout=-1
        IF (LEN_TRIM(outname)>1) THEN
            iout=11
            OPEN(UNIT=iout,FILE=outname,ACTION='WRITE')
        ENDIF
        !
        ! Predict
        DO i=1,n
            emuInputVec(:)=ex2(i,1:m)
            obs=obs2(i)
            ! Emulate, unstandardize and invert logistic
            pred = unstandardize_s(gp%predict(emuInputVec,obs),meant,stdt)
            IF (iout>0) THEN
                WRITE(iout,*) i,',', t2(i),',', pred,',', t2(i)-pred
            ELSE
                WRITE(*,*) i,',', t2(i),',', pred,',', t2(i)-pred
            ENDIF
        ENDDO
        IF (iout>0) CLOSE(iout)
        !
        ! Save emulator only if it is trained here (emuname='')
        IF (LEN_TRIM(outemu)>1 .AND. LEN_TRIM(emuname)==0) call gp%write_out(outemu)
    END SUBROUTINE test_emulator


    SUBROUTINE leave_one_out(dataname,outname)
        ! Leave-one-out tests
        IMPLICIT NONE
        ! Inputs
        CHARACTER(*), INTENT(IN) :: dataname ! Training data file
        CHARACTER(*), INTENT(IN) :: outname ! Output file name
        ! Local variables
        REAL(dp), ALLOCATABLE :: ex1(:,:),t1(:),emuInputVec(:)
        INTEGER, ALLOCATABLE :: obs1(:)
        REAL(dp) :: emu, meanx, stdx, meant, stdt
        INTEGER :: i, j, n, m, obs, iout
        !
        ! Read training data: ex2(500,7), t2(500) and obs2(500)
        CALL read_training(dataname,n,m)
        !
        ! Allocate data
        ALLOCATE(ex1(n-1,m),obs1(n-1),t1(n-1),emuInputVec(m)) 
        !
        ! Output
        iout=-1
        IF (LEN_TRIM(outname)>1) THEN
            iout=11
            OPEN(UNIT=iout,FILE=outname,ACTION='WRITE')
        ENDIF
        !
        ! Leave-one-out tests
        DO i=1,n
            ! Training data (without row i)
            IF (i>1) THEN
                ex1(1:i-1,1:m)=ex2(1:i-1,1:m)
                t1(1:i-1)=t2(1:i-1)
                obs1(1:i-1)=obs2(1:i-1)
            ENDIF
            IF (i<n) THEN
                ex1(i:n-1,1:m)=ex2(i+1:n,1:m)
                t1(i:n-1)=t2(i+1:n)
                obs1(i:n-1)=obs2(i+1:n)
            ENDIF
            ! Emulator input (row i)
            emuInputVec(:)=ex2(i,1:m)
            obs=obs2(i)
            !
            ! Standardize training data and emulator inputs
            meant = mean(t1,n-1)
            stdt  = std(t1,meant,n-1)
            t1(:)=(t1(:) - meant)/(stdt)
            DO j=1,m ! cos_mu included
                meanx = mean(ex1(:,j),n-1)
                stdx  = std(ex1(:,j),meanx,n-1)
                ex1(:,j)=(ex1(:,j) - meanx)/(stdx)
                !
                ! Note: no min/max limits for emulator inputs here (used in mo_emulator.f90)!
                emuInputVec(j) = (emuInputVec(j) - meanx)/(stdx)
            ENDDO
            !
            ! Emulator training
            CALL train_emulator(n-1, m, ex1, obs1, t1)
            !
            ! Emulate, unstandardize and invert logistic
            emu = unstandardize_s(gp%predict(emuInputVec,obs),meant,stdt)
            !
            ! Write output (emulator prediction, the LES output and the difference)
            IF (iout>0) THEN
                WRITE(UNIT=iout,FMT=*) emu, t2(i), t2(i)-emu
            ELSE
                WRITE(UNIT=*,FMT=*) emu, t2(i), t2(i)-emu
            ENDIF
        ENDDO
        ! All done
        IF (iout>0) CLOSE(iout)
    END SUBROUTINE leave_one_out


    SUBROUTINE train_emulator(n, input_dimension, x, obs_type, t)
        ! Emulator training
        use m_util
        use m_gp
        use m_gp_dense
        use m_gp_optim
        use m_cov_all
        use m_noise_all
        implicit none
        ! Inputs
        integer, INTENT(IN) :: n, input_dimension ! Number of training points and dimension of the input
        INTEGER, INTENT(IN) :: obs_type(n) ! Observation type
        real(dp), INTENT(IN) :: x(n,input_dimension), t(n) ! Training data
        ! Emulator parameters
        INTEGER :: nnu, ntheta
        character(len=max_name_len) :: covariance_function = 'LINSQEXP'
        character(len=max_name_len) :: noise_model_name    = 'VAL'
        class(cov_fn), allocatable :: cf
        class(noise_model), allocatable :: nm
        REAL(dp), ALLOCATABLE :: nu(:), theta(:), lbounds(:), ubounds(:)
        LOGICAL :: optimize = .true.
        INTEGER :: optimize_max_iter = 10000
        real(dp) :: optimize_ftol = 1.0d-7

        IF (.NOT.ALLOCATED(cf)) call string_to_cov_fn(covariance_function, cf)
        IF (.NOT.ALLOCATED(nm)) call string_to_noise_model(noise_model_name, nm)

        nnu = nm%nparams_required(input_dimension)
        ntheta = cf%ntheta_required(input_dimension)
  
        ALLOCATE(nu(nnu),theta(ntheta),lbounds(nnu+ntheta),ubounds(nnu+ntheta))

        nu = 0.001 ! Here nnu=1 and ntheta=input_dimension+1
        IF (ntheta==8) THEN
            theta = (/ 0.9010,0.9650,0.6729,3.5576,4.7418,1.2722,4.0612,0.5 /)
        ELSEIF (ntheta==7) THEN
            theta = (/ 0.9010,0.9650,0.6729,3.5576,4.7418,1.2722,4.0612  /)
        ELSE
            theta(:) = 1.
        ENDIF
        lbounds(1) = 0.001
        lbounds(2:) = 0.01
        ubounds(:) = 100.0
        ubounds(1) = 0.001

        ! Update emulator
        IF (allocated(gp)) DEALLOCATE(gp)
        allocate(gp, source=DenseGP(nu, theta, x, obs_type, t, cf, nm))

        if (optimize) then
            call log_lik_optim(nnu + ntheta, gp, lbounds, ubounds, optimize_max_iter, optimize_ftol)
        else
            call gp%update_matrices
        end if
    END SUBROUTINE train_emulator


    SUBROUTINE read_training(fname,n,m)
        ! Function for reading training data files
        CHARACTER(*) :: fname ! Data file name
        INTEGER, INTENT(OUT) :: n, m ! Dimensions
        ! Emulator inputs and outputs
        REAL(dp) :: tmp(9)
        INTEGER :: i

        ! Input file
        open(UNIT=11,file=fname,ACTION="READ")
        !
        ! Number of columns
        tmp(:)=-1000.
        read (11,*,IOSTAT=n) tmp(1:9)
        IF (n/=0) THEN
            STOP 'IO error while reading input file (line 1)!'
        ELSEIF (abs(tmp(8))<1e-10) THEN
            ! 9 columns (6+1 variables, column with zeros, and LES output), so daytime data
            m=7
        ELSEIF (abs(tmp(7))<1e-10) THEN
            ! 8 columns, so nighttime data
            m=6
        ELSE
            ! Bad data
            STOP 'Error reading input file (line 1)!'
        ENDIF
        REWIND(11)
        !
        ! Read all
        n=500 ! Default/maximum size
        DO i=1,500
            read (11,*,IOSTAT=n) ex2(i,1:m), obs2(i), t2(i)
            IF (n .NE. 0) EXIT ! EOF
        ENDDO
        close(11)
        n=i-1
        WRITE(*,*) 'Data loaded, dimensions',n,m
        if (n .le. 10) STOP 'Not enough data!'
    END SUBROUTINE read_training


  FUNCTION standardize(x,meanx,stdx) RESULT(res)
    REAL(dp) :: x
    REAL(dp) :: res
    REAL(dp) :: meanx 
    REAL(dp) :: stdx

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

END PROGRAM emulator
