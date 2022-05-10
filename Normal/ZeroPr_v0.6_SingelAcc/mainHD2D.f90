!=================================================================
PROGRAM HD2D
  !=================================================================
  ! MHD2D code
  !
  ! Numerically integrates the incompressible HD equations
  ! in 2 dimensions using the streamfunction formulation
  ! with an external force.
  ! A pseudo-spectral method is used to compute spatial
  ! derivatives, while variable order Runge-Kutta method
  ! is used to evolve the system in time domain.
  ! To compile, you need the FFTW library installed on
  ! your system. You should link with the FFTP subroutines
  ! and use the FFTPLANS and MPIVARS modules (see the file
  ! 'fftp_mod.f90').
  !
  ! NOTATION: index 'i' is 'x'
  !           index 'j' is 'y'
  !
  ! 2014 Vassilios Dallas
  !      Laboratoire de Physique Statistique
  !      Departement de Physique
  !      Ecole Normale Superieure
  !      e-mail: vdallas@lps.ens.fr
  !
  ! NOTE: 1. Computation of kenergy and penergy balance
  !=================================================================

  USE fprecision
  USE commtypes
  USE mpivars
  USE fft
  USE grid
  USE ali
  USE var
  USE kes
  USE boxsize
  use io
  USE random

  IMPLICIT NONE

  !
  ! Integration parameters
  !     ord  : order of the Runge-Kutta method used

  INTEGER, PARAMETER :: ord = 3
  INTEGER :: ini
  INTEGER :: step
  INTEGER :: tstep
  INTEGER :: cstep
  INTEGER :: sstep
  !
  ! streamfunction, vector potential, z component
  ! of the fields and external force matrixes

  COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:) :: ps
  !
  ! Temporal data storing matrixes

  COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:) :: C1,C2,C3
  REAL(KIND=GP),    ALLOCATABLE, DIMENSION (:,:) :: R1, R2 ! Phil: R2 will be used for temperature calculations
  !
  ! Some auxiliary matrixes

  REAL(KIND=GP) :: dt,CFL
  !!!      REAL(KIND=GP) :: prm1,prm2
  REAL(KIND=GP) :: dump,tmp
  REAL(KIND=GP) :: time,timef
  COMPLEX(KIND=GP) :: rannum
  !REAL(KIND=GP) :: cort
  REAL(KIND=GP) :: nu,kappa, Ra
  REAL(KIND=GP) :: phase1,phase2
  REAL(KIND=GP) :: rmp,rmq
  REAL(KIND=GP) :: Gamma
  REAL(KIND=GP) :: RaC, dRa
  REAL(KIND=GP) :: OEEF, OE, EE
  REAL(KIND=GP) :: enke, enko ! Current even energy, and current odd energy

  INTEGER :: mult,stat
  INTEGER :: t,o
  INTEGER :: i,j,ir,jr
  INTEGER :: ic,id,iu
  INTEGER :: jc,jd,ju,jt
  INTEGER :: timet,timec,times
  INTEGER :: seed,seedf, InC
  !INTEGER :: kfup,kfdn
  INTEGER :: znl
  INTEGER :: size

  CHARACTER        :: c,d,u,th
  CHARACTER(len=3) :: node,ext
  CHARACTER(len=4) :: ext4

  !
  ! Initializes the MPI library

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

  ic = 48+int(myrank/100)
  id = 48+int(myrank/10)-int(myrank/100)*10
  iu = 48+int(myrank)-int(myrank/10)*10
  c = char(ic)
  d = char(id)
  u = char(iu)
  node = c // d // u

  CALL range(1,nx/2+1,nprocs,myrank,ista,iend)
  CALL range(1,ny,nprocs,myrank,jsta,jend)

  !
  ! Initializes the FFT library
  ! Use FFTW_ESTIMATE in short runs and FFTW_MEASURE
  ! in long runs

  CALL fftp2d_create_plan(planrc,(/nx,ny/),FFTW_REAL_TO_COMPLEX,0, &
  FFTW_MEASURE)
  CALL fftp2d_create_plan(plancr,(/nx,ny/),FFTW_COMPLEX_TO_REAL,1, &
  FFTW_MEASURE)
  CALL fftp2d_create_plan(plancc,(/nx,ny/),FFTW_COSINE,2,          &
  FFTW_MEASURE)

  !
  ! Allocates memory for distributed blocks

  ALLOCATE( R1(nx,jsta:jend) )
  ALLOCATE( R2(nx,jsta:jend) )
  ALLOCATE( C1(ny,ista:iend) )
  ALLOCATE( C2(ny,ista:iend) )
  ALLOCATE( C3(ny,ista:iend) )
  ALLOCATE( ps(ny,ista:iend) )
  ALLOCATE( kx(nx), ky(ny), kn2(ny,ista:iend) )
  ALLOCATE( kk2(ny,ista:iend) )


  !
  ! Reads from the external file 'status.txt'
  ! the status of a previous run (if any)
  !     stat: last output of a previous run
  !     mult: time step multiplier

  IF (myrank.eq.0) THEN
    OPEN(1,file='status.prm',status='unknown')
    READ(1,*) stat
    READ(1,*) time
    READ(1,*) znl
    READ(1,*) mult
    CLOSE(1)
  ENDIF
  CALL MPI_BCAST(stat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(time,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(znl ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(mult,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


  Lx = 1.0_GP
  Ly = 1.0_GP
  Dkx = 1.0_GP
  Dky = 1.0_GP
  Dkk = 0.0_GP
  !
  ! Reads from the external file 'parameter.txt' the
  ! parameters that will be used during the integration
  !     dt   : time step size
  !     step : total number of time step to compute
  !     tstep: number of steps between I/O writing
  !     sstep: number of steps between power spectrum I/O
  !     cstep: number of steps between information output
  !     f0   : amplitude of the external kinetic force, is actually theta0 now
  !     e0   : amplitude of the external electric kinetic force
  !     a0   : amplitude of the initial vector potential
  !     kup  : minimum wave number in the external force
  !     kdn  : maximum wave number in the external force
  !     kmup : minimum wave number in the external electric force
  !     kmdn : maximum wave number in the external electric force
  !     nu   : kinematic viscosity
  !     mu   : magnetic diffusivity
  !     seed : seed for random function
  !     prm1 : free parameter 1
  !     prm2 : free parameter 2
  !     ldir : local directory for I/O

  IF (myrank.eq.0) THEN
    OPEN(1,file='input.prm',status='unknown') ! Below code gets inputs from the input.prm file
    READ(1,'(a100)') ldir              ! 25
    READ(1,'(a100)') cdir              ! 26
    READ(1,'(a100)') sdir              ! 27
    READ(1,*) dRa                       ! 17
    READ(1,*) CFL                      ! 1
    READ(1,*) step                     ! 2
    READ(1,*) tstep                    ! 3
    READ(1,*) sstep                    ! 4
    READ(1,*) cstep                    ! 5
    READ(1,*) Lx                       ! 6
    READ(1,*) Ly                       ! 7        ! 8
    READ(1,*) seed                     ! 22
    READ(1,*) InC
    READ(1,*) OEEF                     ! Odd/Even energy fraction. If 0, then we do not consider this
    CLOSE(1)
    CFL = CFL/dble(mult)
    step = step*mult
    tstep = tstep*mult
    sstep = sstep*mult
    cstep = cstep*mult
  ENDIF
  CALL MPI_BCAST(  CFL,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( step,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(tstep,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(sstep,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(cstep,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(   Lx,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(   Ly,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(  Dkk,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(   dRa,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(   OEEF,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( seed,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( InC,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( ldir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( cdir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( sdir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  Dkx = 1.0_GP/Lx
  Dky = 1.0_GP/Ly
  IF (Dkk.lt.1e-5) Dkk = min(Dkx,Dky) ! lt is less than. Dkk usually set to 0 originally. so this will be true, ie it will be sqrt(2) since Lx = sqrt(2)


  size = ny*nx/2
  Gamma = Lx*2.0_GP/Ly
  RaC = pi**4*(4.0_GP+Gamma**2)**3/(4.0_GP*Gamma**4);
  Ra = RaC + dRa
  CALL MPI_BCAST(size,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(   Ra,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  !
  ! Some numerical constants

  ic = 48
  id = 48
  iu = 48
  jt = 48
  jc = 48
  jd = 48
  ju = 48

  !
  ! Some constants for the FFT
  !     kmax: maximum truncation for dealiasing. Dealiasing gets rid of big wavenumbers.

  kmax = 1.0_GP/9.0_GP
  nmax = int(max(nx*Dkx,ny*Dky)/Dkk) ! Currently nmax is something used in modules which we dont use
  IF ((Ly/Lx)**2.ge.(2*ny/nx)**2) THEN
    k2max = (pi**2)*(((nx/3.0d0)**2)*((Ly/Lx)**2-(2*ny/nx)**2)+(2.0d0*ny/3.0d0)**2)
  ELSE
    k2max = (pi**2)*(((2.0d0*ny/3.0d0)**2)*(1-(Ly*nx/(Lx*2*ny))**2)+(Ly*nx/(3.0d0*Lx))**2)
  END IF
  specmax = NINT(sqrt(k2max)/(pi*Ly))+1 ! How big K is wrt to the normalsing length.

  ! Builds arrays with the wavenumbers and the
  ! square wavenumbers. At the end, kx, ky, and kz
  ! have wavenumbers with dimensions, kk2 has the
  ! squared wavenumbers with dimensions, and kn2 has
  ! the dimensionless and normalized squared
  ! wavenumbers used for dealiasing.

  DO i = 1,nx/2
    kx(i) = real(i-1,kind=GP)
    kx(i+nx/2) = real(i-nx/2-1,kind=GP) ! (0, 1,..,nx/2-1,-nx/2,-nx/2+1,..,-1)
  END DO
  DO j = 1,ny
    ky(j) = real(j,kind=GP) ! (1,2,..,ny)
  END DO
  rmp = 1.0_GP/real(nx,kind=GP)**2
  rmq = 1.0_GP/(2.0d0*real(ny,kind=GP))**2
  DO i = ista,iend
    DO j = 1,ny
      kn2(j,i) = rmp*kx(i)**2+rmq*ky(j)**2 ! This together with defintion of kmax above, we get that m < 2ny/3 and n < nx/3. Which is what we want
    END DO
  END DO
  kx = kx*pi*Ly/Lx ! non-dim wavenumbers
  ky = ky*pi
  DO i = ista,iend
    DO j = 1,ny
      kk2(j,i) = kx(i)**2+ky(j)**2
    END DO
  END DO



  !
  ! Sets the initial conditions.
  !
  IF (stat.eq.0) THEN ! Stat is 0, taken from status file
    ! InCs that is good for bursts

  IF (InC.eq.1) THEN
    DO i = ista,iend
      DO j = 1,ny
        IF (i.eq.2.and.j.eq.1) THEN
          ps(j,i) = 1.0*1e-10 ! Small things in 1,1
        ELSE
          ps(j,i) = 0.0
        ENDIF
      END DO
    END DO

  ELSE IF (InC.eq.2) then ! random

  DO i = ista,iend
    DO j = 1,ny
      IF (kn2(j,i).le.kmax) THEN
        rannum = 2*pi*randu(seed)*im
        ps(j,i) = real(nx,kind=GP)*real(ny,kind=GP)*exp(rannum)/kk2(j,i) ! absolulte value is 1 in all of the modes
      ELSE
        ps(j,i) = 0.0
      ENDIF
    END DO
  END DO

ELSE IF (InC.eq.3) then ! from .txt file
  ALLOCATE( Rpsi(size), Ipsi(size))

  IF (myrank.eq.0) THEN
    OPEN(1,file='PsiE.txt',status='unknown') ! Getting Psi ICs
    do i = 1, size, 1
      read(1,*) Rpsi(i)
      read(1,*) Ipsi(i)
    end do
    CLOSE(1)
  END IF
  DO i = ista,iend ! Setting to 0 and then adding ICs
    DO j = 1,ny
      ps(j,i) = 0.0
      IF (i.lt.(nx/2+1).and.j.lt.(ny+1).and.kn2(j,i).le.kmax) THEN ! Checking that mode is even and that we are within the bounds of PsiE and ThetaE
        ps(j,i) = real(Rpsi(i + (ny/2)*(j-1)),kind=GP)*real(nx,kind=GP)*real(ny,kind=GP) ! The number here is just N/4
        IF (i.ne.1) THEN ! If we have an imaginary part also
          ps(j,i) = ps(j,i) + im*real(Ipsi(i + (ny/2)*(j-1)),kind=GP)*real(nx,kind=GP)*real(ny,kind=GP)
        END IF
      END IF
    END DO
  END DO

  DEALLOCATE( Rpsi, Ipsi )

  end if

  ! here we consider the even and odd energy stuff
  IF (OEEF.neq.0.0) THEN
    OE = min(OEEF,1.0)
    EE = min(1.0/OEEF,1.0)
    CALL MPI_BCAST(   OEF,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(   EEF,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)

    ! Calculating current energies
    CALL oeenergy(ps, enke, enko, 1)
    CALL MPI_BCAST(enke,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(enke,1,GC_REAL,0,MPI_COMM_WORLD,ierr)

    ! Fixing even energy
    EEF=sqrt(EEF/enke)
    DO i = ista,iend
      DO j = 1,ny
      IF (MOD(i+j,2).eq.1) THEN! Sum is odd, so we have an even mode
        IF (kn2(j,i).le.kmax) THEN ! Don't actually think we need this here since alising only comes into play when non-linear terms are used
          ps(j,i) = EEF*ps(j,i)
        ELSE
          ps(j,i) = 0.0 ! Gets rid of the wavenumbers that are too large
        END IF
      END IF
      END DO
    END DO

    ! Fixing odd energy
    OEF=OEF/sqrt(enko)
    DO i = ista,iend
      DO j = 1,ny
      IF (MOD(i+j,2).eq.0) THEN! Sum is even, so we have an odd mode
        IF (kn2(j,i).le.kmax) THEN ! Don't actually think we need this here since alising only comes into play when non-linear terms are used
          ps(j,i) = OEF*ps(j,i)
        ELSE
          ps(j,i) = 0.0 ! Gets rid of the wavenumbers that are too large
        END IF
      END IF
      END DO
    END DO




  END IF


  !END IF


  ELSE  !! stat

    ic = 48+int(float(stat)/100)
    id = 48+int(float(stat)/10)-int(float(stat)/100)*10
    iu = 48+int(stat)-int(float(stat)/10)*10
    c = char(ic)
    d = char(id)
    u = char(iu)

    DO j = jsta,jend
      DO i = 1,nx
        R1(i,j) = 0.0
      END DO
    END DO

    OPEN(1,file=trim(ldir) // '/hd2Dps.' // node // '.' &
    // c // d // u //'.dat',form='unformatted')
    READ(1) R1
    CLOSE(1)
    CALL fftp2d_real_to_complex(planrc,R1,ps,MPI_COMM_WORLD)

  END IF  !! stat

  !
  ! Time integration scheme starts here
  ! Uses Runge-Kutta of order 'ord'

  !! Output times
  if (stat.eq.0) then ! Phil: If stat is 0, we will get output straig away basically
    timet = tstep  !! for fields
    times = sstep  !! for spectra
    timec = cstep  !! for checks
    timef = 0.0
    !if (cort.gt.0.0) seedf = seed
    ini = 0
  else
    timet = 0  !! for fields
    times = 0  !! for spectra
    timec = 0  !! for checks
    !if (cort.gt.0.) timef = mod(time,cort)
    !if (cort.gt.0.) seedf = seed
    ini = int((stat-1)*tstep)

    dump = dble(ini)/dble(sstep)+1
    jt = 48+int(dump/1000)
    jc = 48+int(dump/100)-int(dump/1000)*10
    jd = 48+int(dump/10)-int(dump/100)*10
    ju = 48+int(dump)-int(dump/10)*10
  end if

  !#################### MAIN LOOP ######################
  RK : DO t = ini,step

    !
    ! Checks the necessary CFL condition for numerical stability
    CALL CFL_condition(ps,Ly, Lx,CFL,ord,dt)

    !
    ! Every 'cstep' steps, generates external files
    ! to check consistency and convergency. See the
    ! mhdcheck subroutine for details.
    IF (timec.eq.cstep) THEN
      CALL hdcheck(ps,time,znl)
      timec = 0
    ENDIF

    !
    ! Every 'sstep' steps, generates external files
    ! with the power spectrum
    IF (times.eq.sstep) THEN
      times = 0
      ju = ju+1
      IF (ju.eq.58) THEN
        ju = 48
        jd = jd+1
      ENDIF
      IF (jd.eq.58) THEN
        jd = 48
        jc = jc+1
      ENDIF
      IF (jc.eq.58) THEN
        jc = 48
        jt = jt+1
      ENDIF
      th= char(jt)
      c = char(jc)
      d = char(jd)
      u = char(ju)
      ext4 = th // c // d // u
      CALL spectrum(ps,ext4) ! Phil: Writes kspectrum file
      IF (znl.eq.1) THEN
        CALL zonalmean(ps,ext4,1) ! Phil: Writes zonalmean file
        !!!            CALL laplak2(ps,C1)     ! make W
        !!!            CALL vectrans(ps,ps,C1,'euu',ext4)
        !!!            CALL vectrans(C1,ps,C1,'vrt',ext4)
      ENDIF
      IF (myrank.eq.0) THEN
        OPEN(1,file=trim(sdir) // '/spec_times.txt',position='append')
        WRITE(1,*) ext4,time
        !   13      FORMAT( A4,    F12.6)
        CLOSE(1)
      ENDIF

    ENDIF

    !
    ! Every 'tstep' steps, stores the results of the integration
    IF (timet.eq.tstep) THEN

      timet = 0
      iu = iu+1
      IF (iu.eq.58) THEN
        iu = 48
        id = id+1
      ENDIF
      IF (id.eq.58) THEN
        id = 48
        ic = ic+1
      ENDIF
      c = char(ic)
      d = char(id)
      u = char(iu)
      ext = c // d // u
      !            CALL spectrum2D(ps,aa,ext)

      rmp = 1.0_GP/(2.0_GP*real(nx,kind=GP)*real(ny,kind=GP))

      DO i = ista,iend
        DO j = 1,ny
          C1(j,i) = ps(j,i)*rmp
        END DO
      END DO
      CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
      OPEN(1,file=trim(ldir) // '/hd2Dps.' // node // '.' &
      // c // d // u // '.dat',form='unformatted')
      WRITE(1) R1
      CLOSE(1)

      CALL laplak2(C1,C3)
      CALL fftp2d_complex_to_real(plancr,C3,R1,MPI_COMM_WORLD)
      OPEN(1,file=trim(ldir) // '/hd2Dww.' // node // '.' &
      // c // d // u // '.dat',form='unformatted')
      WRITE(1) R1
      CLOSE(1)

      IF (myrank.eq.0) THEN
        OPEN(1,file=trim(ldir)//'/field_times.txt',position='append')
        WRITE(1,*) c//d//u,time,seedf
        CLOSE(1)
      ENDIF

    ENDIF

    ! Copies the streamfunction and theta into the auxiliary matrix C1 and C2

    DO i = ista,iend
      DO j = 1,ny
        C1(j,i) = ps(j,i)
      END DO
    END DO

    ! Runge-Kutta intermediate step(s)

    DO o = ord,2,-1
      CALL laplak2(C1,C3)     ! make W
      CALL poisson(C1,C3,C3)  ! u grad w. Poisson bracket

      rmp = 1.0_GP/real(o,kind=GP)
      DO i = ista,iend
        DO j = 1,ny
          IF (kn2(j,i).le.kmax) THEN
            C1(j,i) = (C1(j,i) + dt*rmp*( C3(j,i)/kk2(j,i)+Ra*(kx(i)**2)*C1(j,i)/(kk2(j,i)**2))) &
            /(1.0 + kk2(j,i)*dt*rmp)

          ELSE
            C1(j,i) = 0.0
          ENDIF
        END DO
      END DO

    END DO

    !
    ! Runge-Kutta final step

    o = 1
    CALL laplak2(C1,C3)     ! make W
    CALL poisson(C1,C3,C3)  ! u grad w. Poisson bracket

    rmp = 1.0_GP/real(o,kind=GP)
    DO i = ista,iend
      DO j = 1,ny
        IF (kn2(j,i).le.kmax) THEN
          ps(j,i) = (ps(j,i) + dt*rmp*( C3(j,i)/kk2(j,i)+Ra*(kx(i)**2)*C1(j,i)/(kk2(j,i)**2))) &
          /(1.0 + kk2(j,i)*dt*rmp)

        ELSE
          ps(j,i) = 0.0
        ENDIF

      END DO
    END DO

    timet = timet+1
    times = times+1
    timec = timec+1
    timef = timef+dt
    time  = time+dt

  END DO RK
  !##############  END OF MAIN LOOP ###################

  !
  ! End of Runge-Kutta

  CALL MPI_FINALIZE(ierr)
  CALL fftp2d_destroy_plan(plancr)
  CALL fftp2d_destroy_plan(planrc)
  DEALLOCATE( R1,R2 )
  DEALLOCATE( ps )
  DEALLOCATE( C1,C2,C3 )
  DEALLOCATE( kx,ky )
  DEALLOCATE( kk2 )

  DEALLOCATE( kn2 )

END PROGRAM HD2D
