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
  INTEGER :: size
  !
  ! streamfunction, vector potential, z component
  ! of the fields and external force matrixes

  COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:) :: ps
  COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:) :: theta
  !
  ! Temporal data storing matrixes

  COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:) :: C1,C3,C2, C4
  REAL(KIND=GP),    ALLOCATABLE, DIMENSION (:,:) :: R1, R2 ! Phil: R2 will be used for temperature calculations
  !
  ! Some auxiliary matrixes

  REAL(KIND=GP) :: enerk
  REAL(KIND=GP) :: dt,CFL
  !!!      REAL(KIND=GP) :: prm1,prm2
  REAL(KIND=GP) :: dump,tmp
  REAL(KIND=GP) :: theta0, u0, perts
  REAL(KIND=GP) :: delta,convlim,deltaold ! Used for convergance
  REAL(KIND=GP) :: time
  REAL(KIND=GP) :: rannum
  !REAL(KIND=GP) :: cort
  REAL(KIND=GP) :: nu,kappa, Ra, Pr
  REAL(KIND=GP) :: rmp,rmq
  REAL(KIND=GP) :: Rp, Ip
  REAL(KIND=GP) :: OEEF, OE, EE
  REAL(KIND=GP) :: enke, enko ! Current even energy, and current odd energy

  INTEGER :: mult,stat, InC, pert, norm, convs, even
  INTEGER :: ordvf,ordvh
  INTEGER :: t,o
  INTEGER :: i,j,ir,jr
  INTEGER :: ic,id,iu,ix
  INTEGER :: jc,jd,ju,jt
  INTEGER :: timet,timec,times
  INTEGER :: seed
  !INTEGER :: kfup,kfdn

  CHARACTER        :: c,d,u,th,x
  CHARACTER(len=3) :: node
  CHARACTER(len=4) :: ext4, ext

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
  ALLOCATE( C4(ny,ista:iend) )
  ALLOCATE( theta(ny,ista:iend) )
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
    READ(1,*) mult
    READ(1,*) convs
    CLOSE(1)
  ENDIF
  CALL MPI_BCAST(convs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(stat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(time,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
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
  !     nu   : kinematic viscosity
  !     mu   : magnetic diffusivity
  !     prm1 : free parameter 1
  !     prm2 : free parameter 2
  !     ldir : local directory for I/O

  IF (myrank.eq.0) THEN
    OPEN(1,file='input.prm',status='unknown') ! Below code gets inputs from the input.prm file
    READ(1,'(a100)') ldir              ! 25
    READ(1,'(a100)') cdir              ! 26
    READ(1,'(a100)') sdir              ! 27
    READ(1,*) Ra                       ! 17
    READ(1,*) Pr
    READ(1,*) CFL
    READ(1,*) step                     ! 2
    READ(1,*) tstep                    ! 3
    READ(1,*) sstep                    ! 4
    READ(1,*) cstep                    ! 5
    READ(1,*) Lx                       ! 6
    READ(1,*) Ly                       ! 7
    READ(1,*) theta0                       ! 9
    READ(1,*) u0
    READ(1,*) InC                     ! 22
    READ(1,*) seed
    READ(1,*) pert
    READ(1,*) norm
    READ(1,*) perts
    READ(1,*) convlim
    READ(1,*) even
    READ(1,*) OEEF                     ! Odd/Even energy fraction. If 0, then we do not consider this
    CLOSE(1)
    ! These are some things which we had in input before, but they dont seem to change, so i have just hard coded here instead
    ordvf = 1
    ordvh = 0
    CFL = CFL/dble(mult)
    step = step*mult
    tstep = tstep*mult
    sstep = sstep*mult
    cstep = cstep*mult
  ENDIF
  delta = 1.0d0
  deltaold = 1.0d0
  CALL MPI_BCAST(  CFL,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( step,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(tstep,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(sstep,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(cstep,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(   Lx,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(   Ly,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(  Dkk,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(   OEEF,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(theta0, 1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(delta, 1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(deltaold, 1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(convlim, 1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(perts, 1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(u0, 1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(ordvf,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(ordvh,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(   Ra,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(   Pr,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( ldir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( cdir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( sdir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(InC,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(pert,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(even,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(norm,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( seed,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)
  Dkx = 1.0_GP/Lx
  Dky = 1.0_GP/Ly
  IF (Dkk.lt.1e-5) Dkk = min(Dkx,Dky) ! lt is less than. Dkk usually set to 0 originally. so this will be true, ie it will be sqrt(2) since Lx = sqrt(2)


  ! Calculating nu and kappa
  nu = sqrt((pi*Ly)**3*Pr/Ra)
  kappa = sqrt((pi*Ly)**3/(Ra*Pr))
  size = ny*nx/4
  CALL MPI_BCAST(   nu,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(kappa,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(size,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)

  !
  ! Some numerical constants

  ic = 48
  id = 48
  iu = 48
  ix = 48
  jt = 48
  jc = 48
  jd = 48
  ju = 48

  !
  ! Some constants for the FFT
  !     kmax: maximum truncation for dealiasing. Dealiasing gets rid of big wavenumbers.

  kmax = 1.0_GP/9.0_GP
  kxmax = int(nx*sqrt(1/9-1/(4*ny^2))) ! int rounds down, which is what we want
  kymax = int(2*ny/3)
  nmax = int(max(nx*Dkx,ny*Dky)/Dkk)
  specmax = int(sqrt(2.0d0)*max(nx,2*ny)/(3*min(Lx,Ly))) + 1 ! largest value NINT(sqrt(kk2(j,i))) can take

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
      kn2(j,i) = rmp*kx(i)**2+rmq*ky(j)**2
    END DO
  END DO
  kx = kx*Dkx
  ky = ky*Dky
  DO i = ista,iend
    DO j = 1,ny
      kk2(j,i) = kx(i)**2+ky(j)**2
    END DO
  END DO



  !
  ! Sets the initial conditions.
  !
  IF (stat.eq.0) THEN

    IF (InC.eq.0) THEN
      DO i = ista,iend
        DO j = 1,ny
          IF (kn2(j,i).le.kmax) THEN
            rannum = 2*pi*randu(seed)
            ps(j,i) = real(nx,kind=GP)*real(ny,kind=GP)*exp(rannum*im)/kk2(j,i) ! absolulte value is 1 in all of the modes
          ELSE
            ps(j,i) = 0.0d0
          ENDIF
        END DO
      END DO

    ELSE IF (InC.eq.1) THEN

    ELSE IF (InC.eq.2) THEN

    ELSE IF (InC.eq.3) THEN

    ELSE IF (InC.eq.4) THEN

    ELSE IF (InC.eq.5) THEN
      ALLOCATE( Rpsi(size), Ipsi(size), Rtheta(size), Itheta(size) ) ! Allocation memory to get the ICs, length of the IC vector in r+ij form

      IF (myrank.eq.0) THEN
        OPEN(1,file='PsiE.txt',status='unknown') ! Getting Psi ICs
        do i = 1, size, 1
          read(1,*) Rpsi(i)
          read(1,*) Ipsi(i)
        end do
        CLOSE(1)
        OPEN(1,file='ThetaE.txt',status='unknown') ! Getting Theta ICs
        do i = 1, size, 1
          read(1,*) Rtheta(i)
          read(1,*) Itheta(i)
        end do
        CLOSE(1)
      END IF

      DO i = ista,iend ! Setting to 0 and then adding ICs
        DO j = 1,ny
          ps(j,i) = 0.0d0
          theta(j,i) = 0.0d0
          IF (MOD(i+j,2).eq.1.and.i.lt.(nx/2+1).and.j.lt.(ny+1)) THEN ! Checking that mode is even and that we are within the bounds of PsiE and ThetaE
            ps(j,i) = real(Rpsi((j-1)*(nx/4)+1+(i-1-mod(i-1,2))/2),kind=GP)*real(nx,kind=GP)*real(ny,kind=GP) ! The number here is just N/4
            theta(j,i) = real(Rtheta((j-1)*(nx/4)+1+(i-1-mod(i-1,2))/2),kind=GP)*real(nx,kind=GP)*real(ny,kind=GP)
            IF (i.ne.1) THEN ! If we have an imaginary part also
              ps(j,i) = ps(j,i) + im*real(Ipsi((j-1)*(nx/4)+1+(i-1-mod(i-1,2))/2),kind=GP)*real(nx,kind=GP)*real(ny,kind=GP)
              theta(j,i) = theta(j,i) + im*real(Itheta((j-1)*(nx/4)+1+(i-1-mod(i-1,2))/2),kind=GP)*real(nx,kind=GP)*real(ny,kind=GP)
            END IF
          END IF


        END DO
      END DO

      DEALLOCATE( Rpsi, Ipsi, Rtheta, Itheta )


    ELSE IF (InC.eq.6) Then
      IF (convs.eq.1) THEN ! Read in ICs from fortran
            ix = 48+int(float(convs)/1000)
      ic = 48+int(float(convs)/100)-int(float(convs)/1000)*10
      id = 48+int(float(convs)/10)-int(float(convs)/100)*10
      iu = 48+int(convs)-int(float(convs)/10)*10
      x = char(ix)
      c = char(ic)
      d = char(id)
      u = char(iu)

      DO j = jsta,jend
        DO i = 1,nx
          R1(i,j) = 0.0d0
        END DO
      END DO

      OPEN(1,file=trim(ldir) // '/hd2Dps.' // node // '.' &
      // x // c // d // u //'.dat',form='unformatted')
      READ(1) R1
      CLOSE(1)
      CALL fftp2d_real_to_complex(planrc,R1,ps,MPI_COMM_WORLD)

      DO j = jsta,jend
        DO i = 1,nx
          R1(i,j) = 0.0d0
        END DO
      END DO
      OPEN(1,file=trim(ldir) // '/hd2Dtheta.' // node // '.' &
      // x // c // d // u //'.dat',form='unformatted')
      READ(1) R1
      CLOSE(1)

      !DO j = jsta,jend
      !  DO i = 1,nx
      !     R1(i,j) = R1(i,j) - 1.0d0-(dble(j)-0.5d0)/dble(ny) ! Calculating theta from temperature
      !  END DO
      !END DO

      CALL fftp2d_real_to_complex(planrc,R1,theta,MPI_COMM_WORLD)
      ElSE IF (convs.eq.2) THEN ! load in from MATLAB
        ! These number have been set for 172
        ALLOCATE( Rpsi(size), Ipsi(size), Rtheta(size), Itheta(size) ) ! Allocation memory to get the ICs, length of the IC vector in r+ij form

        IF (myrank.eq.0) THEN
          OPEN(1,file='PsiE.txt',status='unknown') ! Getting Psi ICs
          do i = 1, size, 1
            read(1,*) Rpsi(i)
            read(1,*) Ipsi(i)
          end do
          CLOSE(1)
          OPEN(1,file='ThetaE.txt',status='unknown') ! Getting Theta ICs
          do i = 1, size, 1
            read(1,*) Rtheta(i)
            read(1,*) Itheta(i)
          end do
          CLOSE(1)
        END IF

        DO i = ista,iend ! Setting to 0 and then adding ICs
          DO j = 1,ny
            ps(j,i) = 0.0d0
            theta(j,i) = 0.0d0
            IF (MOD(i+j,2).eq.1.and.i.lt.(nx/2+1).and.j.lt.(ny+1)) THEN ! Checking that mode is even and that we are within the bounds of PsiE and ThetaE
              ps(j,i) = real(Rpsi((j-1)*(nx/4)+1+(i-1-mod(i-1,2))/2),kind=GP)*real(nx,kind=GP)*real(ny,kind=GP) ! The number here is just N/4
              theta(j,i) = real(Rtheta((j-1)*(nx/4)+1+(i-1-mod(i-1,2))/2),kind=GP)*real(nx,kind=GP)*real(ny,kind=GP)
              IF (i.ne.1) THEN ! If we have an imaginary part also
                ps(j,i) = ps(j,i) + im*real(Ipsi((j-1)*(nx/4)+1+(i-1-mod(i-1,2))/2),kind=GP)*real(nx,kind=GP)*real(ny,kind=GP)
                theta(j,i) = theta(j,i) + im*real(Itheta((j-1)*(nx/4)+1+(i-1-mod(i-1,2))/2),kind=GP)*real(nx,kind=GP)*real(ny,kind=GP)
              END IF
            END IF


          END DO
        END DO

        DEALLOCATE( Rpsi, Ipsi, Rtheta, Itheta )

     ElSE IF (convs.eq.3) THEN ! Do the low-dim approx
        DO i = ista,iend
        DO j = 1,ny
        IF ((j.eq.1).and.(i.eq.2)) THEN
            ps(j,i) = 1.0d0*2*sqrt((Ra-8*pi**4)/(8*pi**4))*kappa*im*real(nx,kind=GP)*real(ny,kind=GP)
            theta(j,i) = -1.0d0*2*sqrt((Ra-8*pi**4)/(8*pi**4))*4*real(nx,kind=GP)*real(ny,kind=GP)*pi**3/Ra
        ELSE IF ((j.eq.2).and.(i.eq.1)) THEN
          theta(j,i) = -1.0d0*((2*sqrt((Ra-8*pi**4)/(8*pi**4)))**2)*2*real(nx,kind=GP)*real(ny,kind=GP)*pi**3/Ra

        ELSE IF (MOD(i+j,2).eq.1.and.kn2(j,i).le.kmax) THEN ! Adding some even modes on top
          rannum = 2*pi*randu(seed)*kx(i)
          ps(j,i) = real(nx,kind=GP)*real(ny,kind=GP)*1.0d0*1e-10*(cos(rannum)+im*sin(rannum))
          theta(j,i) = real(nx,kind=GP)*real(ny,kind=GP)*1.0d0*1e-10*(cos(rannum)-im*sin(rannum))
        ELSE
          ps(j,i) = 0.0d0
          theta(j,i) = 0.0d0
        END IF
    END DO
  END DO


      ELSE ! (within ICs 6, random ICs)
      DO i = ista,iend ! Setting random ICs in even modes
        DO j = 1,ny
          ps(j,i) = 0.0d0
          theta(j,i) = 0.0d0
          IF (MOD(i+j,2).eq.1.and.kn2(j,i).le.kmax) Then
            rannum = randu(seed)
            ps(j,i) = (sin(rannum)+im*cos(rannum))*real(nx,kind=GP)*real(ny,kind=GP)
            theta(j,i) = (sin(rannum)-im*cos(rannum))*real(nx,kind=GP)*real(ny,kind=GP)
          END IF
        END DO
      END DO
              ! Normalising
        CALL energy(theta,enerk,0)
        CALL MPI_BCAST(enerk,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
        tmp=theta0/sqrt(enerk)
        DO i = ista,iend
          DO j = 1,ny
            IF (kn2(j,i).le.kmax) THEN ! Don't actually think we need this here since alising only comes into play when non-linear terms are used
              theta(j,i) = tmp*theta(j,i)
            ELSE
              theta(j,i) = 0.0d0 ! Gets rid of the wavenumbers that are too large
            END IF
          END DO
        END DO


        CALL energy(ps,enerk,1)
        CALL MPI_BCAST(enerk,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
        tmp=u0/sqrt(enerk)
        DO i = ista,iend
          DO j = 1,ny
            IF (kn2(j,i).le.kmax) THEN ! Don't actually think we need this here since alising only comes into play when non-linear terms are used
              ps(j,i) = tmp*ps(j,i)
            ELSE
              ps(j,i) = 0.0d0 ! Gets rid of the wavenumbers that are too large
            END IF
          END DO
        END DO
      END IF
      timec = cstep
      DO WHILE (abs(delta).gt.convlim)        !
        ! Checks the necessary CFL condition for numerical stability
        CALL CFL_condition(ps,nu,kappa,CFL,ord,ordvf,Lx,Ly,dt)

        !
        ! Every 'cstep' steps, generates external files
        ! to check consistency and convergency. See the
        ! mhdcheck subroutine for details.
        IF (timec.eq.cstep) THEN
          IF (myrank.eq.0) THEN ! calculate delta
           delta = abs(ps(1,2))*(1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP))) - deltaold
           deltaold = abs(ps(1,2))*(1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)))
          END IF
          CALL MPI_BCAST(delta, 1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(deltaold, 1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
          CALL hdchecktwo(ps,time,delta)
          timec = 0
        ENDIF

        ! Runge-Kutta step 1
        ! Copies the streamfunction and theta into the auxiliary matrix C1 and C2

        DO i = ista,iend
          DO j = 1,ny
            C1(j,i) = ps(j,i)
            C2(j,i) = theta(j,i)
          END DO
        END DO

        ! Runge-Kutta step 2

        DO o = ord,2,-1
          CALL poisson(C1,C2,C4)  ! make u grad theta
          CALL laplak2(C1,C3)     ! make W
          CALL poisson(C1,C3,C3)  ! u grad w. Poisson bracket

          rmp = 1.0_GP/real(o,kind=GP)
          DO i = ista,iend
            DO j = 1,ny
              IF (kn2(j,i).le.kmax.and.MOD(i+j,2).eq.1) THEN ! We only want even modes
                C1(j,i) = (ps(j,i) + dt*rmp*( C3(j,i)/kk2(j,i)-im*kx(i)*C2(j,i)/kk2(j,i))) &
                /(1.0d0 + (nu*kk2(j,i))*dt*rmp)

                C2(j,i) = (theta(j,i) + dt*rmp*(-C4(j,i)+im*kx(i)*C1(j,i)/pi)) &
                /(1.0d0 + (kappa*kk2(j,i))*dt*rmp)
              ELSE
                C1(j,i) = 0.0d0
                C2(j,i) = 0.0d0
              ENDIF

            END DO
          END DO

        END DO

        !
        ! Runge-Kutta step 3

        o = 1

        CALL poisson(C1,C2,C4)  ! make u grad theta
        CALL laplak2(C1,C3)     ! make W
        CALL poisson(C1,C3,C3)  ! u grad w. Poisson bracket.

        rmp = 1.0_GP/real(o,kind=GP)
        DO i = ista,iend
          DO j = 1,ny
            IF (kn2(j,i).le.kmax.and.MOD(i+j,2).eq.1) THEN
              ps(j,i) = (ps(j,i)+ dt*rmp*( C3(j,i)/kk2(j,i)-im*kx(i)*C2(j,i)/kk2(j,i))) &
              /(1.0d0 + (nu*kk2(j,i))*dt*rmp)

              theta(j,i) = (theta(j,i)+ dt*rmp*(-C4(j,i)+im*kx(i)*C1(j,i)/pi)) &
              /(1.0d0 + (kappa*kk2(j,i))*dt*rmp)
            ELSE
              ps(j,i) = 0.0d0
              theta(j,i) = 0.0d0
            ENDIF

          END DO
        END DO
        CALL MPI_BCAST(delta,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
        timec = timec+1
        time  = time+dt
      end do

    if (InC.eq.6) then
        CALL hdchecktwo(ps,time,delta)
    end if


    ELSE IF (InC.eq.7) THEN ! steady state with odd and even modes from matlab

      ! mmax*(nmax+1)
      ALLOCATE( Rpsi(size*2), Ipsi(size*2), Rtheta(size*2), Itheta(size*2) ) ! Allocation memory to get the ICs, length of the IC vector in r+ij form

      IF (myrank.eq.0) THEN
        OPEN(1,file='PsiE.txt',status='unknown') ! Getting Psi ICs
        do i = 1, size*2, 1
          read(1,*) Rpsi(i)
          read(1,*) Ipsi(i)
        end do
        CLOSE(1)
        OPEN(1,file='ThetaE.txt',status='unknown') ! Getting Theta ICs
        do i = 1, size*2, 1
          read(1,*) Rtheta(i)
          read(1,*) Itheta(i)
        end do
        CLOSE(1)
      END IF

      DO i = ista,iend ! Setting to 0 and then adding ICs
        DO j = 1,ny
          ps(j,i) = 0.0d0
          theta(j,i) = 0.0d0
          IF (i.lt.(nx/2+1).and.j.lt.(ny+1)) THEN ! Checking that mode is even and that we are within the bounds of PsiE and ThetaE
            ps(j,i) = real(Rpsi(i + (nx/2)*(j-1)),kind=GP)*real(nx,kind=GP)*real(ny,kind=GP) ! The number here is just N/4
            theta(j,i) = real(Rtheta(i + (nx/2)*(j-1)),kind=GP)*real(nx,kind=GP)*real(ny,kind=GP)
            IF (i.ne.1) THEN ! If we have an imaginary part also
              ps(j,i) = ps(j,i) + im*real(Ipsi(i + (nx/2)*(j-1)),kind=GP)*real(nx,kind=GP)*real(ny,kind=GP)
              theta(j,i) = theta(j,i) + im*real(Itheta(i + (nx/2)*(j-1)),kind=GP)*real(nx,kind=GP)*real(ny,kind=GP)
            END IF
          END IF
        END DO
      END DO

      DEALLOCATE( Rpsi, Ipsi, Rtheta, Itheta )

    END IF ! end if for ics, there is another one for stat


      IF (norm.eq.1) THEN ! Dont want to normalise when ICs is from MATLAB
        ! Normalising
        CALL energy(theta,enerk,0)
        CALL MPI_BCAST(enerk,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
        tmp=theta0/sqrt(enerk)
        DO i = ista,iend
          DO j = 1,ny
            IF (kn2(j,i).le.kmax) THEN ! Don't actually think we need this here since alising only comes into play when non-linear terms are used
              theta(j,i) = tmp*theta(j,i)
            ELSE
              theta(j,i) = 0.0d0 ! Gets rid of the wavenumbers that are too large
            END IF
          END DO
        END DO


        CALL energy(ps,enerk,1)
        CALL MPI_BCAST(enerk,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
        tmp=u0/sqrt(enerk)
        DO i = ista,iend
          DO j = 1,ny
            IF (kn2(j,i).le.kmax) THEN ! Don't actually think we need this here since alising only comes into play when non-linear terms are used
              ps(j,i) = tmp*ps(j,i)
            ELSE
              ps(j,i) = 0.0d0 ! Gets rid of the wavenumbers that are too large
            END IF
          END DO
        END DO
      END IF


      ! here we consider the even and odd energy stuff
      IF (OEEF.neq.0.0d0) THEN
        OE = min(OEEF,1.0d0)
        EE = min(1.0d0/OEEF,1.0d0)
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
              ps(j,i) = 0.0d0 ! Gets rid of the wavenumbers that are too large
            END IF
          END IF
          END DO
        END DO

        ! Fixing odd energy
        OEF=sqrt(OEF/enko)
        DO i = ista,iend
          DO j = 1,ny
          IF (MOD(i+j,2).eq.0) THEN! Sum is even, so we have an odd mode
            IF (kn2(j,i).le.kmax) THEN ! Don't actually think we need this here since alising only comes into play when non-linear terms are used
              ps(j,i) = OEF*ps(j,i)
            ELSE
              ps(j,i) = 0.0d0 ! Gets rid of the wavenumbers that are too large
            END IF
          END IF
          END DO
        END DO




      END IF



    ELSE  !! stat

      ix = 48+int(float(stat)/1000)
      ic = 48+int(float(stat)/100)-int(float(stat)/1000)*10
      id = 48+int(float(stat)/10)-int(float(stat)/100)*10
      iu = 48+int(stat)-int(float(stat)/10)*10
      x = char(ix)
      c = char(ic)
      d = char(id)
      u = char(iu)

      DO j = jsta,jend
        DO i = 1,nx
          R1(i,j) = 0.0d0
        END DO
      END DO

      OPEN(1,file=trim(ldir) // '/hd2Dps.' // node // '.' &
      // x // c // d // u //'.dat',form='unformatted')
      READ(1) R1
      CLOSE(1)
      CALL fftp2d_real_to_complex(planrc,R1,ps,MPI_COMM_WORLD)

      DO j = jsta,jend
        DO i = 1,nx
          R1(i,j) = 0.0d0
        END DO
      END DO
      OPEN(1,file=trim(ldir) // '/hd2Dtheta.' // node // '.' &
      // x // c // d // u //'.dat',form='unformatted')
      READ(1) R1
      CLOSE(1)

      !DO j = jsta,jend
      !  DO i = 1,nx
      !     R1(i,j) = R1(i,j) - 1.0d0-(dble(j)-0.5d0)/dble(ny) ! Calculating theta from temperature
      !  END DO
      !END DO

      CALL fftp2d_real_to_complex(planrc,R1,theta,MPI_COMM_WORLD)



    END IF  !! stat
    IF (pert.eq.1) THEN
    DO i = ista,iend
      DO j = 1,ny
        IF (MOD(i+j,2).eq.0.and.kn2(j,i).le.kmax) THEN ! Adding a small odd pertubation
          rannum = 2*pi*randu(seed)*kx(i)
          ps(j,i) = real(nx,kind=GP)*real(ny,kind=GP)*perts*(cos(rannum)+im*sin(rannum))
          theta(j,i) = real(nx,kind=GP)*real(ny,kind=GP)*perts*(cos(rannum)-im*sin(rannum))
        END IF
        IF (MOD(i+j,2).eq.1.and.kn2(j,i).le.kmax.and.even.eq.1) THEN ! adding even pert
          rannum = 2*pi*randu(seed)*kx(i)
          ps(j,i) = real(nx,kind=GP)*real(ny,kind=GP)*perts*(cos(rannum)+im*sin(rannum))
          theta(j,i) = real(nx,kind=GP)*real(ny,kind=GP)*perts*(cos(rannum)-im*sin(rannum))
        END IF
      END DO
    END DO

  ELSE IF (pert.eq.2) THEN ! Adding eigenfunction from matlab
    ALLOCATE( Rpsi(size), Ipsi(size), Rtheta(size), Itheta(size) ) ! Allocation memory to get the ICs, length of the IC vector in r+ij form

    IF (myrank.eq.0) THEN
      OPEN(1,file='PsiV.txt',status='unknown') ! Getting Psi ICs
      do i = 1, size, 1
        read(1,*) Rpsi(i)
        read(1,*) Ipsi(i)
      end do
      CLOSE(1)
      OPEN(1,file='ThetaV.txt',status='unknown') ! Getting Theta ICs
      do i = 1, size, 1
        read(1,*) Rtheta(i)
        read(1,*) Itheta(i)
      end do
      CLOSE(1)
    END IF

    DO i = ista,iend ! adding the part
      DO j = 1,ny
        IF (MOD(i+j,2).eq.0.and.i.lt.(nx/2+1).and.j.lt.(ny+1)) THEN ! Checking that mode is odd and that we are within the bounds of PsiV and ThetaV
          ps(j,i) = real(Rpsi((j-1)*(nx/4)+1+(i-1-mod(i-1,2))/2),kind=GP)*real(nx,kind=GP)*real(ny,kind=GP)*perts ! The number here is just N/4
          theta(j,i) = real(Rtheta((j-1)*(nx/4)+1+(i-1-mod(i-1,2))/2),kind=GP)*real(nx,kind=GP)*real(ny,kind=GP)*perts
          IF (i.ne.1) THEN ! If we have an imaginary part also
            ps(j,i) = ps(j,i) + im*real(Ipsi((j-1)*(nx/4)+1+(i-1-mod(i-1,2))/2),kind=GP)*real(nx,kind=GP)*real(ny,kind=GP)*perts
            theta(j,i) = theta(j,i) + im*real(Itheta((j-1)*(nx/4)+1+(i-1-mod(i-1,2))/2),kind=GP)*real(nx,kind=GP)*real(ny,kind=GP)*perts
          END IF
        END IF


      END DO
    END DO

    DEALLOCATE( Rpsi, Ipsi, Rtheta, Itheta )


   END IF
    !
    ! Time integration scheme starts here
    ! Uses Runge-Kutta of order 'ord'

    !! Output times
    if (stat.eq.0) then ! Phil: If stat is 0, we will get output straig away basically
      time = 0.0d0
      timet = tstep  !! for fields
      times = sstep  !! for spectra
      timec = cstep  !! for checks
      !if (cort.gt.0.0d0) seedf = seed
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
      CALL CFL_condition(ps,nu,kappa,CFL,ord,ordvf,Lx,Ly,dt)

      !
      ! Every 'cstep' steps, generates external files
      ! to check consistency and convergency. See the
      ! mhdcheck subroutine for details.
      IF (timec.eq.cstep) THEN
        CALL hdcheck(ps,theta,time,ordvf,ordvh)
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
        !CALL spectrum2D(ps,theta,ext4)
        CALL fieldsfs(ps,theta,ext4)
        CALL spectrum(ps,theta,ext4) ! Phil: Writes kspectrum file
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
        IF (ic.eq.58) THEN
          ic = 48
          ix = ix+1
        ENDIF
        x = char(ix)
        c = char(ic)
        d = char(id)
        u = char(iu)

        ext = x // c // d // u

        rmp = 1.0_GP/(2.0_GP*real(nx,kind=GP)*real(ny,kind=GP))

        DO i = ista,iend
          DO j = 1,ny
            C1(j,i) = ps(j,i)*rmp
          END DO
        END DO
        CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
        OPEN(1,file=trim(ldir) // '/hd2Dps.' // node // '.' &
        // ext // '.dat',form='unformatted')
        WRITE(1) R1
        CLOSE(1)

        CALL laplak2(C1,C3)
        CALL fftp2d_complex_to_real(plancr,C3,R1,MPI_COMM_WORLD)
        OPEN(1,file=trim(ldir) // '/hd2Dww.' // node // '.' &
        // ext // '.dat',form='unformatted')
        WRITE(1) R1
        CLOSE(1)

        ! Printing theta
        DO i = ista,iend
          DO j = 1,ny
            C1(j,i) = theta(j,i)*rmp
          END DO
        END DO
        CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)

        OPEN(1,file=trim(ldir) // '/hd2Dtheta.' // node // '.' &
        // ext // '.dat',form='unformatted')
        WRITE(1) R1
        CLOSE(1)

        IF (myrank.eq.0) THEN
          OPEN(1,file=trim(ldir)//'/field_times.txt',position='append')
          WRITE(1,*) x//c//d//u,time
          CLOSE(1)
        ENDIF

      ENDIF

      !
      ! Runge-Kutta step 1
      ! Copies the streamfunction and theta into the auxiliary matrix C1 and C2

      DO i = ista,iend
        DO j = 1,ny
          C1(j,i) = ps(j,i)
          C2(j,i) = theta(j,i)
        END DO
      END DO

      ! Runge-Kutta step 2

      DO o = ord,2,-1
        CALL poisson(C1,C2,C4)  ! make u grad theta
        CALL laplak2(C1,C3)     ! make W
        CALL poisson(C1,C3,C3)  ! u grad w. Poisson bracket

        rmp = 1.0_GP/real(o,kind=GP)
        DO i = ista,iend
          DO j = 1,ny
            IF (kn2(j,i).le.kmax) THEN
              C1(j,i) = (ps(j,i) + dt*rmp*( C3(j,i)/kk2(j,i)-im*kx(i)*C2(j,i)/kk2(j,i))) &
              /(1.0d0 + (nu*kk2(j,i))*dt*rmp)

              C2(j,i) = (theta(j,i) + dt*rmp*(-C4(j,i)+im*kx(i)*C1(j,i)/pi)) &
              /(1.0d0 + (kappa*kk2(j,i))*dt*rmp)
            ELSE
              C1(j,i) = 0.0d0
              C2(j,i) = 0.0d0
            ENDIF

          END DO
        END DO

      END DO

      !
      ! Runge-Kutta step 3

      o = 1

      CALL poisson(C1,C2,C4)  ! make u grad theta
      CALL laplak2(C1,C3)     ! make W
      CALL poisson(C1,C3,C3)  ! u grad w. Poisson bracket.

      rmp = 1.0_GP/real(o,kind=GP)
      DO i = ista,iend
        DO j = 1,ny
          IF (kn2(j,i).le.kmax) THEN
            ps(j,i) = (ps(j,i)+ dt*rmp*( C3(j,i)/kk2(j,i)-im*kx(i)*C2(j,i)/kk2(j,i))) &
            /(1.0d0 + (nu*kk2(j,i))*dt*rmp)

            theta(j,i) = (theta(j,i)+ dt*rmp*(-C4(j,i)+im*kx(i)*C1(j,i)/pi)) &
            /(1.0d0 + (kappa*kk2(j,i))*dt*rmp)
          ELSE
            ps(j,i) = 0.0d0
            theta(j,i) = 0.0d0
          ENDIF

        END DO
      END DO

      timet = timet+1
      times = times+1
      timec = timec+1
      time  = time+dt

    END DO RK
    !##############  END OF MAIN LOOP ###################

    !
    ! End of Runge-Kutta

    CALL MPI_FINALIZE(ierr)
    CALL fftp2d_destroy_plan(plancr)
    CALL fftp2d_destroy_plan(planrc)
    DEALLOCATE( R1,R2 )
    DEALLOCATE( ps,theta )
    DEALLOCATE( C1,C3,C2,C4 )
    DEALLOCATE( kx,ky )
    DEALLOCATE( kk2 )

    DEALLOCATE( kn2 )

  END PROGRAM HD2D
