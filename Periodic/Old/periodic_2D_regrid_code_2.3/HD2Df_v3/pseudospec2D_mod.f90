!=================================================================
! MODULES for 2D codes
!
! 2003 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar 
!=================================================================

!=================================================================

  MODULE grid
!
! n: number of points in the spatial grid
! choose: 128   256  512  1024  2048  4096  8192  16384
      INTEGER :: n = 256
      SAVE

  END MODULE grid
!=================================================================

  MODULE fft
!
      USE fftplans
      TYPE (FFTPLAN) :: planrc, plancr
      SAVE

  END MODULE fft
!=================================================================

  MODULE ali
      DOUBLE PRECISION :: kmax
      DOUBLE PRECISION :: tiny
      SAVE

  END MODULE ali
!=================================================================

  MODULE var  
      DOUBLE PRECISION    :: pi = 3.1415926535897932384d0
      DOUBLE COMPLEX :: im = (0.,1.)
      SAVE

  END MODULE var
!=================================================================

  MODULE kes
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)   :: ka
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: ka2
      SAVE

  END MODULE kes
!=================================================================

  MODULE io
      CHARACTER(len=100) :: ldir,cdir,sdir,tdir
      SAVE

  END MODULE io
!=================================================================

  MODULE random
      CONTAINS
       REAL FUNCTION randu(idum)
!
! Uniform distributed random numbers between -1 and 
! 1. The seed idum must be between 0 and the value 
! of mask

       INTEGER, PARAMETER :: iq=127773,ir=2836,mask=123459876
       INTEGER, PARAMETER :: ia=16807,im=2147483647
       INTEGER            :: k,idum
       REAL, PARAMETER    :: am=1./im

       idum = ieor(idum,mask)
       k = idum/iq
       idum = ia*(idum-k*iq)-ir*k
       IF (idum.lt.0) idum = idum+im
       randu = am*idum
       randu = (randu-.5)*2
       idum = ieor(idum,mask)

       END FUNCTION randu

  END MODULE random
!=================================================================
