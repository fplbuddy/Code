!=================================================================
! MODULES for 2D codes
!
! 2014 Vassilios Dallas
!      Laboratoire de Physique Statistique
!      Departement de Physique
!      Ecole Normale Superieure
!      e-mail: vdallas@lps.ens.fr
!=================================================================

!=================================================================

  MODULE grid
!
! n: number of points in the spatial grid
! choose: 128   256  512  1024  2048  4096  8192  16384
      INTEGER :: nx = 64
      INTEGER :: ny = 64
      SAVE

  END MODULE grid
!=================================================================

  MODULE boxsize
!
! Lx, Ly, Lz, and Dk are the lengths of the sides
! of the box, and the width of Fourier shells
      USE fprecision
      REAL(KIND=GP)    :: Lx,Ly
      REAL(KIND=GP)    :: Dkx,Dky,Dkk
      SAVE

  END MODULE boxsize
!=================================================================

  MODULE fft
!
      USE fftplans
      TYPE (FFTPLAN) :: planrc, plancr, plancc
      SAVE

  END MODULE fft
!=================================================================

  MODULE ali
      USE fprecision
      REAL(KIND=GP) :: kmax
      SAVE

  END MODULE ali
!=================================================================

  MODULE var
      USE fprecision
      REAL(KIND=GP)    :: pi = 3.14159265358979323846_GP
      COMPLEX(KIND=GP) :: im = (0.0_GP,1.0_GP)
      SAVE

  END MODULE var
!=================================================================

  MODULE kes
      USE fprecision
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:)           :: kx,ky
      REAL(KIND=GP), TARGET, ALLOCATABLE, DIMENSION (:,:) :: kn2
      REAL(KIND=GP), POINTER, DIMENSION (:,:)             :: kk2
      INTEGER :: nmax, specmax, kxmax, kymax!,nmaxperp
      SAVE

  END MODULE kes
!=================================================================

  MODULE io
      CHARACTER(len=100) :: ldir,cdir,sdir,tdir
      SAVE

  END MODULE io
!=================================================================

  MODULE random
      USE fprecision
      CONTAINS
!
! Uniform distributed random numbers between -1 and
! 1. The seed idum must be between 0 and the value
! of mask

       REAL(KIND=GP) FUNCTION randu(idum)

       INTEGER, PARAMETER :: iq=127773,ir=2836,mask=123459876
       INTEGER, PARAMETER :: ia=16807,im=2147483647
       INTEGER            :: k,idum
       REAL(KIND=GP), PARAMETER :: am=1./im

       idum = ieor(idum,mask)
       k = idum/iq
       idum = ia*(idum-k*iq)-ir*k
       IF (idum.lt.0) idum = idum+im
       randu = am*idum
       randu = (randu-.5)*2
       idum = ieor(idum,mask)
       END FUNCTION randu

!
! Normally distributed random numbers with zero mean
! and unit variance. The seed idum must be between 0
! and the value of mask in randu.

       REAL(KIND=GP) FUNCTION randn(idum)

       REAL(KIND=GP)      :: v1,v2,ran1
       REAL(KIND=GP)      :: fac,rsq
       REAL(KIND=GP),SAVE :: gset
       INTEGER       :: idum
       INTEGER, SAVE :: iset

       IF ((iset.ne.0).or.(iset.ne.1)) iset=0
       IF (idum.lt.0) iset=0
       IF (iset.eq.0) THEN
          rsq = 2.
          DO WHILE ((rsq.ge.1.).or.(rsq.eq.0.))
             v1 = randu(idum)
             v2 = randu(idum)
             rsq = v1**2+v2**2
          END DO
          fac = sqrt(-2.*log(rsq)/rsq)
          gset = v1*fac
          randn = v2*fac
          iset = 1
       ELSE
          randn = gset
          iset = 0
       ENDIF
       END FUNCTION randn

  END MODULE random
!=================================================================
