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
  INTEGER :: rstep
  INTEGER :: size
  !
  ! streamfunction, vector potential, z component
  ! of the fields and external force matrixes

  COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:) :: ps, ph
  COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:) :: theta2, thetav
  !
  ! Temporal data storing matrixes

  COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:) :: C1,C3,C2, C4, C5, C6, C7, C8, C9, C10
  REAL(KIND=GP),    ALLOCATABLE, DIMENSION (:,:) :: R1, R2 ! Phil: R2 will be used for temperature calculations
  !
  ! Some auxiliary matrixes

  REAL(KIND=GP) :: enerk
  REAL(KIND=GP) :: dt
  REAL(KIND=GP) :: phase1,phase2
  !!!      REAL(KIND=GP) :: prm1,prm2
  REAL(KIND=GP) :: dump,tmp
  REAL(KIND=GP) :: theta20, u0, ph0, thetav0
  REAL(KIND=GP) :: time
  REAL(KIND=GP) :: rannum
  !REAL(KIND=GP) :: cort
  REAL(KIND=GP) :: nu,kappa, Ra, Pr, Wid
  REAL(KIND=GP) :: rmp,rmq
  REAL(KIND=GP) :: enke, enko ! Current even energy, and current odd energy

  INTEGER :: mult,stat, InC, norm,normtwo
  INTEGER :: ordvf,ordvh
  INTEGER :: t,o
  INTEGER :: i,j,ir,jr
  INTEGER :: ic,id,iu,ix
  INTEGER :: jc,jd,ju,jt
  INTEGER :: timet,timec,times,timer
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
  ALLOCATE( C5(ny,ista:iend) )
  ALLOCATE( C6(ny,ista:iend) )
  ALLOCATE( C7(ny,ista:iend) )
  ALLOCATE( C8(ny,ista:iend) )
  ALLOCATE( C9(ny,ista:iend) )
  ALLOCATE( C10(ny,ista:iend) )
  ALLOCATE( theta2(ny,ista:iend) )
  ALLOCATE( ps(ny,ista:iend) )
  ALLOCATE( ph(ny,ista:iend) )
  ALLOCATE( thetav(ny,ista:iend) )
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
    CLOSE(1)
  ENDIF
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
    READ(1,*) Wid
    READ(1,*) dt
    READ(1,*) step                     ! 2
    READ(1,*) tstep                    ! 3
    READ(1,*) sstep                    ! 4
    READ(1,*) cstep                    ! 5
    READ(1,*) rstep                    ! 5
    READ(1,*) Lx                       ! 6
    READ(1,*) Ly                       ! 7
    READ(1,*) theta20                       ! 9
    READ(1,*) u0
    READ(1,*) thetav0                       ! 9
    READ(1,*) ph0
    READ(1,*) InC                     ! 22
    READ(1,*) seed
    READ(1,*) norm
    READ(1,*) normtwo
    CLOSE(1)
    ! These are some things which we had in input before, but they dont seem to change, so i have just hard coded here instead
    ordvf = 1
    ordvh = 0
    step = step*mult
    tstep = tstep*mult
    sstep = sstep*mult
    cstep = cstep*mult
    Wid = Wid*pi
  ENDIF
  CALL MPI_BCAST(  dt,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( step,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(tstep,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(sstep,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(cstep,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(rstep,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(   Lx,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(   Ly,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(  Dkk,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(thetav0, 1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(u0, 1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(theta20, 1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(ph0, 1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(ordvf,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(ordvh,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(   Ra,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(   Pr,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(   Wid,  1,GC_REAL,      0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( ldir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( cdir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( sdir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(InC,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(norm,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(normtwo,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( seed,  1,MPI_INTEGER,  0,MPI_COMM_WORLD,ierr)
  Dkx = 1.0_GP/Lx
  Dky = 1.0_GP/Ly
  IF (Dkk.lt.1e-5) Dkk = min(Dkx,Dky) ! lt is less than. Dkk usually set to 0 originally. so this will be true, ie it will be sqrt(2) since Lx = sqrt(2)


  ! Calculating nu and kappa
  nu = sqrt((pi*Ly)**3*Pr/Ra)
  kappa = sqrt((pi*Ly)**3/(Ra*Pr))
  size = ny*nx/2 ! Think this is the right size, will have to check
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
  kxmax = int(nx*sqrt(1.0d0/9.0d0-1.0d0/(4.0d0*ny**2))) ! int rounds down, which is what we want
  kymax = int(2.0d0*ny/3.0d0)
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

    IF (InC.eq.1) THEN
      DO i = ista,iend
        DO j = 1,ny
          IF (kn2(j,i).le.kmax) THEN
            rannum = 2*pi*randu(seed)
            ps(j,i) = real(nx,kind=GP)*real(ny,kind=GP)*exp(rannum*im)/kk2(j,i) ! absolulte value is 1 in all of the modes
            rannum = 2*pi*randu(seed)
            theta2(j,i) = real(nx,kind=GP)*real(ny,kind=GP)*exp(rannum*im)
            rannum = 2*pi*randu(seed)
            ph(j,i) = real(nx,kind=GP)*real(ny,kind=GP)*exp(rannum*im)/kk2(j,i) ! absolulte value is 1 in all of the modes
            rannum = 2*pi*randu(seed)
            thetav(j,i) = real(nx,kind=GP)*real(ny,kind=GP)*exp(rannum*im)
          ELSE
            ps(j,i) = 0.0d0
            theta2(j,i) = 0.0d0
            ph(j,i) = 0.0d0
            thetav(j,i) = 0.0d0
          ENDIF
        END DO
      END DO


    ELSE IF (InC.eq.2) THEN
      ALLOCATE( Rps(size), Ips(size), Rtheta2(size), Itheta2(size), Rph(size), Iph(size), Rthetav(size), Ithetav(size) )
      IF (myrank.eq.0) THEN
        OPEN(1,file='ps.txt',status='unknown') ! Getting ps ICs
        do i = 1, size, 1
          read(1,*) Rps(i)
          read(1,*) Ips(i)
        end do
        CLOSE(1)
        OPEN(1,file='theta2.txt',status='unknown') ! Getting theta2 ICs
        do i = 1, size, 1
          read(1,*) Rtheta2(i)
          read(1,*) Itheta2(i)
        end do
        CLOSE(1)
        OPEN(1,file='ph.txt',status='unknown') ! Getting theta2 ICs
        do i = 1, size, 1
          read(1,*) Rph(i)
          read(1,*) Iph(i)
        end do
        CLOSE(1)
        OPEN(1,file='thetav.txt',status='unknown') ! Getting theta2 ICs
        do i = 1, size, 1
          read(1,*) Rthetav(i)
          read(1,*) Ithetav(i)
        end do
        CLOSE(1)
      END IF
      tmp=real(nx,kind=GP)*real(ny,kind=GP)
      DO i = ista,iend ! Setting to 0 and then adding ICs
        DO j = 1,ny
          ps(j,i) = 0.0d0
          theta2(j,i) = 0.0d0
          ph(j,i) = 0.0d0
          thetav(j,i) = 0.0d0
          IF (i.lt.(nx/2+1).and.j.lt.(ny+1)) THEN ! Checking that mode is even and that we are within the bounds of psE and ThetaE
            ps(j,i) = real(Rps((j-1)*(nx/2)+i),kind=GP)*tmp
            theta2(j,i) = real(Rtheta2((j-1)*(nx/2)+i),kind=GP)*tmp
            ph(j,i) = real(Rph((j-1)*(nx/2)+i),kind=GP)*tmp
            thetav(j,i) = real(Rthetav((j-1)*(nx/2)+i),kind=GP)*tmp
            IF (i.ne.1) THEN ! If we have an imaginary part also
              ps(j,i) = ps(j,i) + im*real(Ips((j-1)*(nx/2)+i),kind=GP)*tmp
              theta2(j,i) = theta2(j,i) + im*real(Itheta2((j-1)*(nx/2)+i),kind=GP)*tmp
              ph(j,i) = ps(j,i) + im*real(Iph((j-1)*(nx/2)+i),kind=GP)*tmp
              thetav(j,i) = theta2(j,i) + im*real(Ithetav((j-1)*(nx/2)+i),kind=GP)*tmp
            END IF
          END IF


        END DO
      END DO

      DEALLOCATE( Rps, Ips, Rtheta2, Itheta2, Rph, Iph, Rthetav, Ithetav )

    ELSE IF (InC.eq.3) THEN ! random for 3d fields, matlab for 2d
      ALLOCATE( Rps(size), Ips(size), Rtheta2(size), Itheta2(size) )
      IF (myrank.eq.0) THEN
        OPEN(1,file='ps.txt',status='unknown') ! Getting ps ICs
        do i = 1, size, 1
          read(1,*) Rps(i)
          read(1,*) Ips(i)
        end do
        CLOSE(1)
        OPEN(1,file='theta2.txt',status='unknown') ! Getting theta2 ICs
        do i = 1, size, 1
          read(1,*) Rtheta2(i)
          read(1,*) Itheta2(i)
        end do
        CLOSE(1)
      END IF
      tmp=real(nx,kind=GP)*real(ny,kind=GP)
      DO i = ista,iend ! Setting to 0 and then adding ICs
        DO j = 1,ny
          ps(j,i) = 0.0d0
          theta2(j,i) = 0.0d0
          IF (i.lt.(nx/2+1).and.j.lt.(ny+1)) THEN ! Checking that mode is even and that we are within the bounds of psE and ThetaE
            ps(j,i) = real(Rps((j-1)*(nx/2)+i),kind=GP)*tmp
            theta2(j,i) = real(Rtheta2((j-1)*(nx/2)+i),kind=GP)*tmp
            IF (i.ne.1) THEN ! If we have an imaginary part also
              ps(j,i) = ps(j,i) + im*real(Ips((j-1)*(nx/2)+i),kind=GP)*tmp
              theta2(j,i) = theta2(j,i) + im*real(Itheta2((j-1)*(nx/2)+i),kind=GP)*tmp
            END IF
          END IF


        END DO
      END DO
      DEALLOCATE( Rps, Ips, Rtheta2, Itheta2 )
      DO i = ista,iend
        DO j = 1,ny
          IF (kn2(j,i).le.kmax) THEN
            rannum = 2*pi*randu(seed)
            ph(j,i) = real(nx,kind=GP)*real(ny,kind=GP)*exp(rannum*im)/kk2(j,i) ! absolulte value is 1 in all of the modes
            rannum = 2*pi*randu(seed)
            thetav(j,i) = real(nx,kind=GP)*real(ny,kind=GP)*exp(rannum*im)
          ELSE
            ph(j,i) = 0.0d0
            thetav(j,i) = 0.0d0
          ENDIF
        END DO
      END DO

    END IF ! end if for ics, there is another one for stat


    IF (norm.eq.1) THEN ! Dont want to normalise when ICs is from MATLAB
      ! theta
      CALL energy(theta2,enerk,0)
      CALL MPI_BCAST(enerk,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      tmp=theta20/sqrt(enerk)
      DO i = ista,iend
        DO j = 1,ny
          IF (kn2(j,i).le.kmax) THEN ! Don't actually think we need this here since alising only comes into play when non-linear terms are used
            theta2(j,i) = tmp*theta2(j,i)
          ELSE
            theta2(j,i) = 0.0d0 ! Gets rid of the wavenumbers that are too large
          END IF
        END DO
      END DO
      ! psi
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

    IF (normtwo.eq.1) THEN ! Dont want to normalise when ICs is from MATLAB
      CALL energy(thetav,enerk,0)
      CALL MPI_BCAST(enerk,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      tmp=thetav0/sqrt(enerk)
      DO i = ista,iend
        DO j = 1,ny
          IF (kn2(j,i).le.kmax) THEN ! Don't actually think we need this here since alising only comes into play when non-linear terms are used
            thetav(j,i) = tmp*thetav(j,i)
          ELSE
            thetav(j,i) = 0.0d0 ! Gets rid of the wavenumbers that are too large
          END IF
        END DO
      END DO
      ! ph
      CALL energy(ph,enerk,1)
      CALL MPI_BCAST(enerk,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      tmp=ph0/sqrt(enerk)
      DO i = ista,iend
        DO j = 1,ny
          IF (kn2(j,i).le.kmax) THEN ! Don't actually think we need this here since alising only comes into play when non-linear terms are used
            ph(j,i) = tmp*ph(j,i)
          ELSE
            ph(j,i) = 0.0d0 ! Gets rid of the wavenumbers that are too large
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
    OPEN(1,file=trim(ldir) // '/hd2Dtheta2.' // node // '.' &
    // x // c // d // u //'.dat',form='unformatted')
    READ(1) R1
    CLOSE(1)
    CALL fftp2d_real_to_complex(planrc,R1,theta2,MPI_COMM_WORLD)

    DO j = jsta,jend
      DO i = 1,nx
        R1(i,j) = 0.0d0
      END DO
    END DO
    OPEN(1,file=trim(ldir) // '/hd2Dph.' // node // '.' &
    // x // c // d // u //'.dat',form='unformatted')
    READ(1) R1
    CLOSE(1)
    CALL fftp2d_real_to_complex(planrc,R1,ph,MPI_COMM_WORLD)

    DO j = jsta,jend
      DO i = 1,nx
        R1(i,j) = 0.0d0
      END DO
    END DO
    OPEN(1,file=trim(ldir) // '/hd2Dthetav.' // node // '.' &
    // x // c // d // u //'.dat',form='unformatted')
    READ(1) R1
    CLOSE(1)
    CALL fftp2d_real_to_complex(planrc,R1,thetav,MPI_COMM_WORLD)



  END IF  !! stat
  !
  ! Time integration scheme starts here
  ! Uses Runge-Kutta of order 'ord'

  !! Output times
  if (stat.eq.0) then ! Phil: If stat is 0, we will get output straig away basically
    time = 0.0d0
    timet = tstep  !! for fields
    times = sstep  !! for spectra
    timec = cstep  !! for checks
    timer = rstep  !! for random field
    !if (cort.gt.0.0d0) seedf = seed
    ini = 0
  else
    timet = 0  !! for fields
    times = 0  !! for spectra
    timec = 0  !! for checks
    timer = 0  !! for random field
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
    !
    ! Every 'cstep' steps, generates external files
    ! to check consistency and convergency. See the
    ! mhdcheck subroutine for details.
    IF (timec.eq.cstep) THEN
      CALL hdcheck(ps,theta2,ph,thetav,time,ordvf,ordvh)
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
      !CALL fieldsfs(ps,theta,ext4)
      CALL spectrum(ps,theta2,ext4,'spec2d')
      CALL spectrum(ph,thetav,ext4,'spec3d')
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

      ! Printing theta
      DO i = ista,iend
        DO j = 1,ny
          C1(j,i) = theta2(j,i)*rmp
        END DO
      END DO
      CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)

      OPEN(1,file=trim(ldir) // '/hd2Dtheta2.' // node // '.' &
      // ext // '.dat',form='unformatted')
      WRITE(1) R1
      CLOSE(1)


      ! Printing ph
      DO i = ista,iend
        DO j = 1,ny
          C1(j,i) = ph(j,i)*rmp
        END DO
      END DO
      CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)

      OPEN(1,file=trim(ldir) // '/hd2Dph.' // node // '.' &
      // ext // '.dat',form='unformatted')
      WRITE(1) R1
      CLOSE(1)

      ! Printing thetav
      DO i = ista,iend
        DO j = 1,ny
          C1(j,i) = thetav(j,i)*rmp
        END DO
      END DO
      CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)

      OPEN(1,file=trim(ldir) // '/hd2Dthetav.' // node // '.' &
      // ext // '.dat',form='unformatted')
      WRITE(1) R1
      CLOSE(1)

      IF (myrank.eq.0) THEN
        OPEN(1,file=trim(ldir)//'/field_times.txt',position='append')
        WRITE(1,*) x//c//d//u,time
        CLOSE(1)
      ENDIF

    ENDIF



    IF (timer.eq.rstep) THEN ! An instance where we wanna do the random field
      timer = 0
      ! first evolve perp fields with random 2d field
      !
      ! Runge-Kutta step 1
      ! Copies the streamfunction and theta2 into the auxiliary matrix C1 and C2

      DO i = ista,iend
        DO j = 1,ny
          C1(j,i) = ps(j,i)
          C2(j,i) = theta2(j,i)
          C5(j,i) = ph(j,i)
          C6(j,i) = thetav(j,i)
        END DO
      END DO

      ! Randomising fields, this could have been put in subroutine maybe
      DO i = ista,iend
        DO j = 1,ny
          IF (kn2(j,i).le.kmax.and.i.gt.1) THEN ! make sure we don't have a zonal mode
            phase1 = pi*randu(seed)
            phase2 = pi*randu(seed)
            C1(j,i) = C1(j,i)*exp(im*phase1)
            C2(j,i) = C2(j,i)*exp(im*phase2)
          ELSE IF (kn2(j,i).le.kmax.and.i.eq.1) THEN ! We have a zonal modes
            phase1 = randu(seed)
            phase2 = randu(seed)
            C1(j,i) = C1(j,i)*sign(1.0d0,phase1)
            C2(j,i) = C2(j,i)*sign(1.0d0,phase2)
          ENDIF
        END DO
      END DO


      ! Runge-Kutta step 2

      DO o = ord,2,-1
        CALL laplak2(C1,C3)     ! make W
        CALL poisson(C5,C3,C7)  ! {phi, nabla^2 psi}
        CALL laplak2(C5,C8)     ! nabla^2 Phi
        CALL poisson(C1,C8,C8) ! {psi, nabla^2 phi}
        CALL poisson(C1,C6,C9) ! {psi, thetav}
        CALL poisson(C5,C2,C10) ! {phi, theta2}

        rmp = 1.0_GP/real(o,kind=GP)
        DO i = ista,iend
          DO j = 1,ny
            IF (kn2(j,i).le.kmax) THEN
              C5(j,i) = (ph(j,i) + dt*rmp*( C7(j,i)/kk2(j,i)+C8(j,i)/kk2(j,i)-im*kx(i)*C6(j,i)/kk2(j,i))) &
              /(1.0d0 + (nu*(kk2(j,i)+(2.0d0*pi/Wid)**2))*dt*rmp)

              C6(j,i) = (thetav(j,i) + dt*rmp*( -C9(j,i)-C10(j,i)+im*kx(i)*C5(j,i)/pi)) &
              /(1.0d0 + (kappa*(kk2(j,i)+(2.0d0*pi/Wid)**2))*dt*rmp)
            ELSE
              C5(j,i) = 0.0d0
              C6(j,i) = 0.0d0
            ENDIF
          END DO
        END DO
      END DO

      !
      ! Runge-Kutta step 3

      o = 1

      CALL laplak2(C1,C3)     ! make W
      CALL poisson(C5,C3,C7)  ! {phi, nabla^2 psi}
      CALL laplak2(C5,C8)     ! nabla^2 Phi
      CALL poisson(C1,C8,C8) ! {psi, nabla^2 phi}
      CALL poisson(C1,C6,C9) ! {psi, thetav}
      CALL poisson(C5,C2,C10) ! {phi, theta2}

      rmp = 1.0_GP/real(o,kind=GP)
      DO i = ista,iend
        DO j = 1,ny
          IF (kn2(j,i).le.kmax) THEN
            ph(j,i) = (ph(j,i) + dt*rmp*( C7(j,i)/kk2(j,i)+C8(j,i)/kk2(j,i)-im*kx(i)*C6(j,i)/kk2(j,i))) &
            /(1.0d0 + (nu*(kk2(j,i)+(2.0d0*pi/Wid)**2))*dt*rmp)

            thetav(j,i) = (thetav(j,i) + dt*rmp*( -C9(j,i)-C10(j,i)+im*kx(i)*C5(j,i)/pi)) &
            /(1.0d0 + (kappa*(kk2(j,i)+(2.0d0*pi/Wid)**2))*dt*rmp)
          ELSE
            ph(j,i) = 0.0d0
            thetav(j,i) = 0.0d0
          ENDIF

        END DO
      END DO

      ! Now evolve the 2D field normally

      !
      ! Runge-Kutta step 1
      ! Copies the streamfunction and theta2 into the auxiliary matrix C1 and C2

      DO i = ista,iend
        DO j = 1,ny
          C1(j,i) = ps(j,i)
          C2(j,i) = theta2(j,i)
        END DO
      END DO

      ! Runge-Kutta step 2

      DO o = ord,2,-1
        CALL poisson(C1,C2,C4)  ! make u grad theta2
        CALL laplak2(C1,C3)     ! make W
        CALL poisson(C1,C3,C3)  ! u grad w. Poisson bracket

        rmp = 1.0_GP/real(o,kind=GP)
        DO i = ista,iend
          DO j = 1,ny
            IF (kn2(j,i).le.kmax) THEN
              C1(j,i) = (ps(j,i) + dt*rmp*( C3(j,i)/kk2(j,i)-im*kx(i)*C2(j,i)/kk2(j,i))) &
              /(1.0d0 + (nu*kk2(j,i))*dt*rmp)
              C2(j,i) = (theta2(j,i) + dt*rmp*(-C4(j,i)+im*kx(i)*C1(j,i)/pi)) &
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

      CALL poisson(C1,C2,C4)  ! make u grad theta2
      CALL laplak2(C1,C3)     ! make W
      CALL poisson(C1,C3,C3)  ! u grad w. Poisson bracket

      rmp = 1.0_GP/real(o,kind=GP)
      DO i = ista,iend
        DO j = 1,ny
          IF (kn2(j,i).le.kmax) THEN
            ps(j,i) = (ps(j,i)+ dt*rmp*( C3(j,i)/kk2(j,i)-im*kx(i)*C2(j,i)/kk2(j,i))) &
            /(1.0d0 + (nu*kk2(j,i))*dt*rmp)

            theta2(j,i) = (theta2(j,i)+ dt*rmp*(-C4(j,i)+im*kx(i)*C1(j,i)/pi)) &
            /(1.0d0 + (kappa*kk2(j,i))*dt*rmp)
          ELSE
            ps(j,i) = 0.0d0
            theta2(j,i) = 0.0d0
          ENDIF

        END DO
      END DO


    ELSE

      !
      ! Runge-Kutta step 1
      ! Copies the streamfunction and theta2 into the auxiliary matrix C1 and C2

      DO i = ista,iend
        DO j = 1,ny
          C1(j,i) = ps(j,i)
          C2(j,i) = theta2(j,i)
          C5(j,i) = ph(j,i)
          C6(j,i) = thetav(j,i)
        END DO
      END DO

      ! Runge-Kutta step 2

      DO o = ord,2,-1
        CALL poisson(C1,C2,C4)  ! make u grad theta2
        CALL laplak2(C1,C3)     ! make W
        CALL poisson(C5,C3,C7)  ! {phi, nabla^2 psi}
        CALL poisson(C1,C3,C3)  ! u grad w. Poisson bracket
        CALL laplak2(C5,C8)     ! nabla^2 Phi
        CALL poisson(C1,C8,C8) ! {psi, nabla^2 phi}
        CALL poisson(C1,C6,C9) ! {psi, thetav}
        CALL poisson(C5,C2,C10) ! {phi, theta2}

        rmp = 1.0_GP/real(o,kind=GP)
        DO i = ista,iend
          DO j = 1,ny
            IF (kn2(j,i).le.kmax) THEN
              C1(j,i) = (ps(j,i) + dt*rmp*( C3(j,i)/kk2(j,i)-im*kx(i)*C2(j,i)/kk2(j,i))) &
              /(1.0d0 + (nu*kk2(j,i))*dt*rmp)
              C2(j,i) = (theta2(j,i) + dt*rmp*(-C4(j,i)+im*kx(i)*C1(j,i)/pi)) &
              /(1.0d0 + (kappa*kk2(j,i))*dt*rmp)

              C5(j,i) = (ph(j,i) + dt*rmp*( C7(j,i)/kk2(j,i)+C8(j,i)/kk2(j,i)-im*kx(i)*C6(j,i)/kk2(j,i))) &
              /(1.0d0 + (nu*(kk2(j,i)+(2.0d0*pi/Wid)**2))*dt*rmp)

              C6(j,i) = (thetav(j,i) + dt*rmp*( -C9(j,i)-C10(j,i)+im*kx(i)*C5(j,i)/pi)) &
              /(1.0d0 + (kappa*(kk2(j,i)+(2.0d0*pi/Wid)**2))*dt*rmp)
            ELSE
              C1(j,i) = 0.0d0
              C2(j,i) = 0.0d0
              C5(j,i) = 0.0d0
              C6(j,i) = 0.0d0
            ENDIF
          END DO
        END DO
      END DO

      !
      ! Runge-Kutta step 3

      o = 1

      CALL poisson(C1,C2,C4)  ! make u grad theta2
      CALL laplak2(C1,C3)     ! make W
      CALL poisson(C5,C3,C7)  ! {phi, nabla^2 psi}
      CALL poisson(C1,C3,C3)  ! u grad w. Poisson bracket
      CALL laplak2(C5,C8)     ! nabla^2 Phi
      CALL poisson(C1,C8,C8) ! {psi, nabla^2 phi}
      CALL poisson(C1,C6,C9) ! {psi, thetav}
      CALL poisson(C5,C2,C10) ! {phi, theta2}

      rmp = 1.0_GP/real(o,kind=GP)
      DO i = ista,iend
        DO j = 1,ny
          IF (kn2(j,i).le.kmax) THEN
            ps(j,i) = (ps(j,i)+ dt*rmp*( C3(j,i)/kk2(j,i)-im*kx(i)*C2(j,i)/kk2(j,i))) &
            /(1.0d0 + (nu*kk2(j,i))*dt*rmp)

            theta2(j,i) = (theta2(j,i)+ dt*rmp*(-C4(j,i)+im*kx(i)*C1(j,i)/pi)) &
            /(1.0d0 + (kappa*kk2(j,i))*dt*rmp)

            ph(j,i) = (ph(j,i) + dt*rmp*( C7(j,i)/kk2(j,i)+C8(j,i)/kk2(j,i)-im*kx(i)*C6(j,i)/kk2(j,i))) &
            /(1.0d0 + (nu*(kk2(j,i)+(2.0d0*pi/Wid)**2))*dt*rmp)

            thetav(j,i) = (thetav(j,i) + dt*rmp*( -C9(j,i)-C10(j,i)+im*kx(i)*C5(j,i)/pi)) &
            /(1.0d0 + (kappa*(kk2(j,i)+(2.0d0*pi/Wid)**2))*dt*rmp)
          ELSE
            ps(j,i) = 0.0d0
            theta2(j,i) = 0.0d0
            ph(j,i) = 0.0d0
            thetav(j,i) = 0.0d0
          ENDIF

        END DO
      END DO
    END IF ! ends if we want to do random or not

    timet = timet+1
    times = times+1
    timec = timec+1
    timer = timer+1
    time  = time+dt

  END DO RK
  !##############  END OF MAIN LOOP ###################

  !
  ! End of Runge-Kutta

  CALL MPI_FINALIZE(ierr)
  CALL fftp2d_destroy_plan(plancr)
  CALL fftp2d_destroy_plan(planrc)
  DEALLOCATE( R1,R2 )
  DEALLOCATE( ps,theta2 )
  DEALLOCATE( C1,C3,C2,C4,C5,C6,C7,C8,C9,C10 )
  DEALLOCATE( kx,ky )
  DEALLOCATE( kk2 )

  DEALLOCATE( kn2 )

END PROGRAM HD2D
