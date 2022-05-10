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
! 2004 Pablo D. Mininni.
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
!
! NOTE: 1. Full 2D spectrum
!       2. Check ftype = 3 before running a job
!=================================================================

      USE mpivars
      USE fft
      USE ali
      USE var
      USE kes
      USE grid
      USE random
      use io

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

      DOUBLE COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: theta
!
! Temporal data storing matrixes

      DOUBLE COMPLEX,   ALLOCATABLE, DIMENSION (:,:) :: C1
      DOUBLE COMPLEX,   ALLOCATABLE, DIMENSION (:,:) :: C2
      DOUBLE COMPLEX,   ALLOCATABLE, DIMENSION (:,:) :: C4
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: R1
!
! Some auxiliary matrixes

      DOUBLE PRECISION :: enerk
      DOUBLE PRECISION :: dt,CFL
      DOUBLE PRECISION :: kr
      !DOUBLE PRECISION :: kup,kdn,kr
      DOUBLE PRECISION :: prm1,prm2
      DOUBLE PRECISION :: dump,tmp,tmp1
      DOUBLE PRECISION :: theta0
      DOUBLE PRECISION :: time,timef
      DOUBLE PRECISION :: cort
      DOUBLE PRECISION :: hnu,kappa, Ra, Rhnu
      DOUBLE PRECISION :: phase1,phase2

      INTEGER :: mult,stat
      INTEGER :: ordvf,ordvh
      INTEGER :: t,o
      INTEGER :: i,j,ir,jr
      INTEGER :: ic,id,iu
      INTEGER :: jc,jd,ju,jt
      INTEGER :: timet,timec,times
      INTEGER :: seed,seedf
      INTEGER :: kup,kdn

      CHARACTER     :: c,d,u,th
      CHARACTER*3   :: node,ext
      CHARACTER*4   :: ext4

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

!
! Allocates memory for distributed blocks

      CALL range(1,n/2+1,nprocs,myrank,ista,iend)
      CALL range(1,n,nprocs,myrank,jsta,jend)

      ALLOCATE( R1(n,jsta:jend) )
      ALLOCATE( C1(n,ista:iend) )
      ALLOCATE( C2(n,ista:iend) )
      ALLOCATE( C4(n,ista:iend) )
      ALLOCATE( theta(n,ista:iend) )
      ALLOCATE( ka(n), ka2(n,ista:iend) )

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
      CALL MPI_BCAST(time,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mult,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!
! Reads from the external file 'parameter.txt' the
! parameters that will be used during the integration
!     dt   : time step size
!     step : total number of time steps to compute
!     tstep: number of steps between I/O writing
!     sstep: number of steps between power spectrum I/O
!     cstep: number of steps between information output
!     f0   : amplitude of the external kinetic force
!     u0   : amplitude of the initial streamfunction
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
         OPEN(1,file='input.prm',status='unknown')
         READ(1,'(a100)') ldir              ! 27
         READ(1,'(a100)') cdir              ! 28
         READ(1,'(a100)') sdir              ! 29
         READ(1,'(a100)') tdir              ! 30
         READ(1,*) Ra                       ! 15
         READ(1,*) Rhnu                      ! 17
         READ(1,*) step                     ! 2
         READ(1,*) tstep                    ! 3
         READ(1,*) sstep                    ! 4
         READ(1,*) cstep                    ! 5
         READ(1,*) ordvf                    ! 25
         READ(1,*) ordvh                    ! 26
         READ(1,*) theta0                       !
         READ(1,*) CFL                      ! 1
         READ(1,*) kup                      ! 11
         READ(1,*) kdn                      ! 12
         READ(1,*) cort                     ! 21
         READ(1,*) seed                     ! 22
         READ(1,*) prm1                     ! 23
         READ(1,*) prm2                     ! 24
         CLOSE(1)
         CFL = CFL/dble(mult)
         step = step*mult
         tstep = tstep*mult
         sstep = sstep*mult
         cstep = cstep*mult
      ENDIF
      CALL MPI_BCAST(  CFL,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST( step,  1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tstep,  1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sstep,  1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cstep,  1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(   theta0,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(  kup,  1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(  kdn,  1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ordvf,  1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ordvh,  1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(   Ra,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(  Rhnu,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST( cort,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST( seed,  1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST( prm1,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST( prm2,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST( ldir,100,MPI_CHARACTER,       0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST( cdir,100,MPI_CHARACTER,       0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST( sdir,100,MPI_CHARACTER,       0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST( tdir,100,MPI_CHARACTER,       0,MPI_COMM_WORLD,ierr)



        ! Calculating nu and kappa and d
  kappa = (2*pi)**(4*ordvf-1)/Ra
  hnu = Rhnu/(2*pi)**(2*(ordvf+ordvh)) ! have set nu = 1
  CALL MPI_BCAST(   kappa,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(   hnu,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
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
!     kmax: maximum truncation for dealiasing
!     tiny: minimum truncation for dealiasing

      kmax = (dble(n)/3.d0)**2
      tiny =  0.000001d0
      enerk = 1.0d0
!!!      enst =  1.0d0

!
! Builds the wave number and the square wave
! number matrixes

      DO i = 1,n/2
         ka(i) = dble(i-1)
         ka(i+n/2) = dble(i-n/2-1)
      END DO
      DO i = ista,iend
         DO j = 1,n
            ka2(j,i) = ka(i)**2+ka(j)**2
         END DO
      END DO

!
! Initializes the FFT library
! Use FFTW_ESTIMATE in short runs and FFTW_MEASURE
! in long runs

      CALL fftp2d_create_plan(planrc,n,FFTW_REAL_TO_COMPLEX, &
                             FFTW_MEASURE)
      CALL fftp2d_create_plan(plancr,n,FFTW_COMPLEX_TO_REAL, &
                             FFTW_MEASURE)

!
! Sets the initial conditions.




!
!STREAM FUNCTION R1 & VECTOR POTENTIAL R2

      IF (stat.eq.0) THEN

      ini = 0
      timet = tstep
      timec = cstep
      times = sstep
      timef = 0.0d0
      if (cort.gt.0.0d0) seedf = seed

      DO j = jsta,jend
         DO i = 1,n
            R1(i,j) = 0.0d0
         END DO
      END DO

      do ir = kup,kdn
         phase1 = pi*randu(seed)
         phase2 = pi*randu(seed)

      DO j = jsta,jend
         DO i = 1,n
            R1(i,j) = R1(i,j) + sin(2*pi*ir*(dble(i)-1)/dble(n)+phase1) &
                              * sin(2*pi*ir*(dble(j)-1)/dble(n)+phase2)
         END DO
      END DO

      end do


      ! Normalise theta and normalise
      CALL fftp2d_real_to_complex(planrc,R1,theta,MPI_COMM_WORLD)
      CALL energy(theta,enerk,0)
      call MPI_BCAST(enerk,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      tmp1=theta0/sqrt(enerk)
      DO i = ista,iend
         DO j = 1,n
            IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
               theta(j,i) = tmp1*theta(j,i)
            ELSE
               theta(j,i) = 0.0d0
            ENDIF
         END DO
      END DO

      ELSE  !! stat

         !print*,'READING...',stat
         ini = int((stat-1)*tstep)
         dump = dble(ini)/dble(sstep)+1
         times = 0
         timet = 0
         timec = 0
         if (cort.gt.0.0d0) timef = mod(time,cort)
         if (cort.gt.0.0d0) seedf = seed

         jt = 48+int(dump/1000)
         jc = 48+int(dump/100)-int(dump/1000)*10
         jd = 48+int(dump/10)-int(dump/100)*10
         ju = 48+int(dump)-int(dump/10)*10

         ic = 48+int(float(stat)/100)
         id = 48+int(float(stat)/10)-int(float(stat)/100)*10
         iu = 48+int(stat)-int(float(stat)/10)*10
         c = char(ic)
         d = char(id)
         u = char(iu)

         OPEN(1,file=trim(ldir) // '/hd2Dtheta.' // node // '.' &
                           // c // d // u //'.dat',form='unformatted')
         READ(1) R1
         CALL fftp2d_real_to_complex(planrc,R1,theta,MPI_COMM_WORLD)

      ENDIF  !! stat

!
! Time integration scheme starts here
! Uses Runge-Kutta of order 'ord'
!#################### MAIN LOOP ######################
 RK : DO t = ini,step
      CALL makeps(theta,hnu,ordvf,ordvh,C1) ! make ps inn C1
      CALL CFL_condition(C1,kappa,CFL,ord,ordvf,dt)
!      print*,t

!
! Every 'cstep' steps, generates external files
! to check consistency and convergency. See the
! mhdcheck subroutine for details.
          IF (timec.eq.cstep) THEN
            timec = 0
              CALL hdcheck(C1,theta,time,ordvf)
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
            ! At this stage, I think C1 is still ps
            CALL spectrum(theta,0,ext4,'PEspec')
            CALL vectrans(theta,C1,theta,'puu',ext4) ! Potential Flux
            CALL derivk2(C1,C2,1)       ! v, want to keep C1 as ps for fields output
            CALL vectrans2(theta,C2,0,'fenktrans',ext4) ! Temperature Gradient Forcing
            CALL vectrans2(theta,theta,ordvf,'denptrans',ext4) ! Potential diss
            ! Dont do normal forcing spec for now
            !CALL vectrans(C1,ps,C1,'vrt',ext4)
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
            DO i = ista,iend
               DO j = 1,n
                  C1(j,i) = C1(j,i)/dble(n)**2
               END DO
            END DO
            CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
            OPEN(1,file=trim(ldir) // '/hd2Dps.' // node // '.' &
                 // c // d // u // '.dat',form='unformatted')
            WRITE(1) R1
            CLOSE(1)

            DO i = ista,iend
               DO j = 1,n
                  C1(j,i) = theta(j,i)/dble(n)**2
               END DO
            END DO
            CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
            OPEN(1,file=trim(ldir) // '/hd2Dtheta.' // node // '.' &
                 // c // d // u // '.dat',form='unformatted')
            WRITE(1) R1
            CLOSE(1)

            IF (myrank.eq.0) THEN
            OPEN(1,file=trim(ldir)//'/field_times.txt',position='append')
            WRITE(1,*) c//d//u,time,seedf
            CLOSE(1)
            ENDIF

         ENDIF

!
! Runge-Kutta step 1
! Copies the streamfunction into the auxiliary matrix C1

         DO i = ista,iend
            DO j = 1,n
               C2(j,i) = theta(j,i)
            END DO
         END DO

!
! Runge-Kutta step 2

         DO o = ord,2,-1
         CALL makeps(C2,hnu,ordvf,ordvh,C1) ! make ps
         CALL poisson(C1,C2,C4)  ! make u grad theta

         DO i = ista,iend
            DO j = 1,n
            IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
              C2(j,i) = (theta(j,i) + dt*(-C4(j,i)+im*ka(i)*C1(j,i)/(2*pi))/dble(o)) &
              /(1.0d0 + (kappa*ka2(j,i)**ordvf)*dt/dble(o))
            ELSE
               C2(j,i) = 0.0d0
            ENDIF

            END DO
         END DO

         END DO

!
! Runge-Kutta step 3
! Copies the result from the auxiliary matrixes into ps, az

         o = 1
         CALL makeps(C2,hnu,ordvf,ordvh,C1) ! make ps
         CALL poisson(C1,C2,C4)  ! make u grad theta

         DO i = ista,iend
            DO j = 1,n
            IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
              theta(j,i) = (theta(j,i)+ dt*(-C4(j,i)+im*ka(i)*C1(j,i)/(2*pi))/dble(o)) &
              /(1.0d0 + (kappa*ka2(j,i)**ordvf)*dt/dble(o))
            ELSE
               theta(j,i) = 0.0d0
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
      DEALLOCATE( R1 )
      DEALLOCATE( theta )
      DEALLOCATE( C1,C2,C4 )
      DEALLOCATE( ka,ka2 )

      END PROGRAM HD2D
