!=================================================================
! MODULES for 2D and 3D codes
!
! 2009 Duane Rosenberg and Pablo D. Mininni.
!      National Center for Atmospheric Research.
!=================================================================

!=================================================================

  MODULE fprecision
!
! 
      INTEGER, PARAMETER :: GP = KIND(0.0D0)
      !INTEGER, PARAMETER :: GP = KIND(0.0)
      SAVE

  END MODULE fprecision
!=================================================================

!=================================================================

  MODULE commtypes
      INCLUDE 'mpif.h'
!
! 
      INTEGER, SAVE :: GC_REAL    = MPI_DOUBLE_PRECISION
      INTEGER, SAVE :: GC_COMPLEX = MPI_DOUBLE_COMPLEX
      !INTEGER, SAVE :: GC_REAL    = MPI_REAL
      !INTEGER, SAVE :: GC_COMPLEX = MPI_COMPLEX

  END MODULE commtypes
!=================================================================
