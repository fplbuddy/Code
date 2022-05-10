!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Subroutines for computing spatial derivatives and nonlinear
! terms in incompressible HD, MHD and Hall-MHD equations in 2D
! using a pseudo-spectral method. You should use the FFTPLANS
! and MPIVARS modules (see the file 'fftp2D_mod.f90') in each
! program that call any of the subroutines in this file.
!
! NOTATION: index 'i' is 'x'
!           index 'j' is 'y'
!
! 2014 Vassilios Dallas and Kannabiran Seshasayanan
!      Laboratoire de Physique Statistique
!      Departement de Physique
!      Ecole Normale Superieure
!      e-mail: vdallas@lps.ens.fr
!=================================================================

!*****************************************************************
SUBROUTINE derivk2(a,b,dir)
  !-----------------------------------------------------------------
  !
  ! Two-dimensional derivative of the matrix 'a'
  !
  ! Parameters
  !     a  : input matrix
  !     b  : at the output contains the derivative da/dk_dir
  !     dir: =1 derivative in the x-direction
  !          =2 derivative in the y-direction
  !
  USE fprecision
  USE mpivars
  USE grid
  USE kes
  USE var
  IMPLICIT NONE

  COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(ny,ista:iend) :: a
  COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(ny,ista:iend) :: b
  INTEGER, INTENT(IN) :: dir
  INTEGER :: i,j

  !
  ! Derivative in the x-direction
  !
  IF (dir.eq.1) THEN
    DO i = ista,iend
      DO j = 1,ny
        b(j,i) = im*kx(i)*a(j,i)
      END DO
    END DO
    !
    ! Derivative in the y-direction
    !
  ELSE
    DO i = ista,iend
      DO j = 1,ny
        b(j,i) = ky(j)*a(j,i)
      END DO
    END DO
  ENDIF

  RETURN
END SUBROUTINE derivk2

!*****************************************************************
SUBROUTINE laplak2(a,b)
  !-----------------------------------------------------------------
  !
  ! Two-dimensional Laplacian of the matrix 'a'
  !
  ! Parameters
  !     a: input matrix
  !     b: at the output contains the Laplacian d2a/dka2
  !
  USE fprecision
  USE kes
  USE grid
  USE mpivars
  IMPLICIT NONE

  COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(ny,ista:iend) :: a
  COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(ny,ista:iend) :: b
  INTEGER :: i,j

  DO i = ista,iend
    DO j = 1,ny
      b(j,i) = -kk2(j,i)*a(j,i)
    END DO
  END DO

  RETURN
END SUBROUTINE laplak2

!*****************************************************************
SUBROUTINE poisson(a,b,c)
  !-----------------------------------------------------------------
  !
  ! Poisson bracket of the scalar fields A and B
  ! in real space.
  !
  ! Parameters
  !     a: input matrix
  !     b: input matrix
  !     c: Poisson bracket {A,B} [output]
  !
  USE fprecision
  USE commtypes
  USE mpivars
  USE grid
  USE fft
  IMPLICIT NONE

  COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(ny,ista:iend) :: a,b
  COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(ny,ista:iend) :: c
  COMPLEX(KIND=GP), DIMENSION(ny,ista:iend) :: c1,c2,c3
  REAL(KIND=GP), DIMENSION(nx,jsta:jend)    :: r1,r2,r3
  REAL(KIND=GP) :: tmp
  INTEGER :: i,j

  !
  ! Computes dA/dx.dB/dy
  !
  CALL derivk2(a,c1,1)
  CALL derivk2(b,c2,2)
  DO i = ista,iend
    c3(1,i) = 0.0_GP
    DO j = 2,ny
      c3(j,i) = c2(j-1,i) ! moving stuff up one
    END DO
  END DO
  CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
  CALL fftp2d_complex_to_real(plancc,c3,r2,MPI_COMM_WORLD) ! Guessing this is a different fourier transform
  DO j = jsta,jend
    DO i = 1,nx
      r3(i,j) = r1(i,j)*r2(i,j)
    END DO
  END DO

  !
  ! Computes dA/dy.dB/dx
  !
  CALL derivk2(a,c1,2)
  CALL derivk2(b,c2,1)
  DO i = ista,iend
    c3(1,i) = 0.0_GP
    DO j = 2,ny
      c3(j,i) = c1(j-1,i)
    END DO
  END DO
  CALL fftp2d_complex_to_real(plancc,c3,r1,MPI_COMM_WORLD)
  CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
  tmp = 1.0_GP/(2.0_GP*real(nx,kind=GP)*real(ny,kind=GP))**2
  DO j = jsta,jend
    DO i = 1,nx
      r3(i,j) = (r3(i,j)-r1(i,j)*r2(i,j))*tmp
    END DO
  END DO
  CALL fftp2d_real_to_complex(planrc,r3,c,MPI_COMM_WORLD)

  RETURN
END SUBROUTINE poisson

!*****************************************************************
SUBROUTINE energy(a,b,kin)
  !-----------------------------------------------------------------
  !
  ! Computes the mean kinetic or magnetic energy in 2D,
  ! and the mean square current density or vorticity.
  ! The output is valid only in the first node.
  !
  ! Parameters
  !     a  : input matrix with the scalar field
  !     b  : at the output contains the energy
  !     kin: =2 computes the square of the scalar field
  !          =1 computes the energy
  !          =0 computes the current or vorticity
  !
  !
  USE fprecision
  USE commtypes
  USE mpivars
  USE grid
  USE kes
  IMPLICIT NONE

  COMPLEX(KIND=GP), INTENT(IN), DIMENSION(ny,ista:iend) :: a
  DOUBLE PRECISION, INTENT(OUT) :: b
  DOUBLE PRECISION              :: bloc
  REAL(KIND=GP)                 :: tmp, tmp1
  INTEGER, INTENT(IN) :: kin
  INTEGER             :: i,j

  bloc = 0.0_GP
  tmp = 1.0_GP/(2.0_GP*real(nx,kind=GP)*real(ny,kind=GP))**2

  !
  ! Computes the square of the scalar field
  !
  IF (kin.eq.0) THEN
    IF (ista.eq.1) THEN
      DO j = 1,ny
        bloc = bloc+abs(a(j,1))**2*tmp
      END DO
      DO i = 2,iend
        DO j = 1,ny
          bloc = bloc+2*abs(a(j,i))**2*tmp
        END DO
      END DO
    ELSE
      DO i = ista,iend
        DO j = 1,ny
          bloc = bloc+2*abs(a(j,i))**2*tmp
        END DO
      END DO
    ENDIF
    !
    ! Computes the energy
    !
  ELSE IF (kin.eq.1) THEN
    IF (ista.eq.1) THEN
      DO j = 1,ny
        bloc = bloc+kk2(j,1)*abs(a(j,1))**2*tmp
      END DO
      DO i = 2,iend
        DO j = 1,ny
          bloc = bloc+2*kk2(j,i)*abs(a(j,i))**2*tmp
        END DO
      END DO
    ELSE
      DO i = ista,iend
        DO j = 1,ny
          bloc = bloc+2*kk2(j,i)*abs(a(j,i))**2*tmp
        END DO
      END DO
    ENDIF
    !
    ! Computes the current or vorticity
    !
  ELSE
    IF (ista.eq.1) THEN
      DO j = 1,ny
        bloc = bloc+kk2(j,1)**kin*abs(a(j,1))**2*tmp ! Phil: An extra factor of 2 below.
      END DO
      DO i = 2,iend
        DO j = 1,ny
          bloc = bloc+2*kk2(j,i)**kin*abs(a(j,i))**2*tmp
        END DO
      END DO
    ELSE
      DO i = ista,iend
        DO j = 1,ny
          bloc = bloc+2*kk2(j,i)**kin*abs(a(j,i))**2*tmp
        END DO
      END DO
    ENDIF
  ENDIF

  ! Computes the reduction between nodes
  !
  CALL MPI_REDUCE(bloc,b,1,GC_REAL,MPI_SUM,0, &
  MPI_COMM_WORLD,ierr) ! Phil: sums up all the blocs in each processor

  RETURN
END SUBROUTINE energy

!*****************************************************************
SUBROUTINE oeenergy(a,b,c,kin)
  !-----------------------------------------------------------------
  !
  ! Computes the mean kinetic or magnetic energy in 2D,
  ! and the mean square current density or vorticity.
  ! The output is valid only in the first node.
  !
  ! Parameters
  !     a  : input matrix with the scalar field
  !     b,c  : at the output contains the e/o energy
  !     kin: =1 computes the energy
  !          =0 computes the current or vorticity
  !
  !
  USE fprecision
  USE commtypes
  USE mpivars
  USE grid
  USE kes
  IMPLICIT NONE

  COMPLEX(KIND=GP), INTENT(IN), DIMENSION(ny,ista:iend) :: a
  DOUBLE PRECISION, INTENT(OUT) :: b,c
  DOUBLE PRECISION              :: bloce,  bloco
  REAL(KIND=GP)                 :: tmp, tmp1
  INTEGER, INTENT(IN) :: kin
  INTEGER             :: i,j

  bloce = 0.0_GP
  bloco = 0.0_GP
  tmp = 1.0_GP/(2.0_GP*real(nx,kind=GP)*real(ny,kind=GP))**2

  !
  ! Computes the square of the scalar field
  !
  IF (kin.eq.0) THEN
    IF (ista.eq.1) THEN
      DO j = 2,ny,2
        bloce = bloce + abs(a(j,1))**2*tmp
      END DO
      DO j = 1,ny,2
        bloco = bloco + abs(a(j,1))**2*tmp
      END DO
      DO i = 2,iend
        IF (MOD(i,2).eq.1) THEN ! i is odd, so kx is even
          DO j = 2,ny,2
            bloce = bloce + 2*abs(a(j,i))**2*tmp
          END DO
          DO j = 1,ny,2
            bloco = bloco + 2*abs(a(j,i))**2*tmp
          END DO
        ELSE
          DO j = 1,ny,2
            bloce = bloce + 2*abs(a(j,i))**2*tmp
          END DO
          DO j = 2,ny,2
            bloco = bloco + 2*abs(a(j,i))**2*tmp
          END DO
        END IF
      END DO
    ELSE
      DO i = ista,iend
        IF (MOD(i,2).eq.1) THEN ! i is odd, so kx is even
          DO j = 2,ny,2
            bloce = bloce + 2*abs(a(j,i))**2*tmp
          END DO
          DO j = 1,ny,2
            bloco = bloco + 2*abs(a(j,i))**2*tmp
          END DO
        ELSE
          DO j = 1,ny,2
            bloce = bloce + 2*abs(a(j,i))**2*tmp
          END DO
          DO j = 2,ny,2
            bloco = bloco + 2*abs(a(j,i))**2*tmp
          END DO
        END IF
      END DO
    ENDIF
    !
    ! Computes the energy
    !
  ELSE IF (kin.eq.1) THEN
    IF (ista.eq.1) THEN
      DO j = 2,ny,2
        bloce = bloce + kk2(j,1)*abs(a(j,1))**2*tmp
      END DO
      DO j = 1,ny,2
        bloco = bloco + kk2(j,1)*abs(a(j,1))**2*tmp
      END DO
      DO i = 2,iend
        IF (MOD(i,2).eq.1) THEN ! i is odd, so kx is even
          DO j = 2,ny,2
            bloce = bloce + 2*kk2(j,i)*abs(a(j,i))**2*tmp
          END DO
          DO j = 1,ny,2
            bloco = bloco + 2*kk2(j,i)*abs(a(j,i))**2*tmp
          END DO
        ELSE
          DO j = 1,ny,2
            bloce = bloce + 2*kk2(j,i)*abs(a(j,i))**2*tmp
          END DO
          DO j = 2,ny,2
            bloco = bloco + 2*kk2(j,i)*abs(a(j,i))**2*tmp
          END DO
        END IF
      END DO
    ELSE
      DO i = ista,iend
        IF (MOD(i,2).eq.1) THEN ! i is odd, so kx is even
          DO j = 2,ny,2
            bloce = bloce + 2*kk2(j,i)*abs(a(j,i))**2*tmp
          END DO
          DO j = 1,ny,2
            bloco = bloco + 2*kk2(j,i)*abs(a(j,i))**2*tmp
          END DO
        ELSE
          DO j = 1,ny,2
            bloce = bloce + 2*kk2(j,i)*abs(a(j,i))**2*tmp
          END DO
          DO j = 2,ny,2
            bloco = bloco + 2*kk2(j,i)*abs(a(j,i))**2*tmp
          END DO
        END IF
      END DO
    ENDIF
  ENDIF

  ! Computes the reduction between nodes
  !
  CALL MPI_REDUCE(bloce,b,1,GC_REAL,MPI_SUM,0, &
  MPI_COMM_WORLD,ierr)
  CALL MPI_REDUCE(bloco,c,1,GC_REAL,MPI_SUM,0, &
  MPI_COMM_WORLD,ierr)

  RETURN
END SUBROUTINE oeenergy

!*****************************************************************
SUBROUTINE inerprod(a,b,kin,rslt)
  !-----------------------------------------------------------------
  ! Parameters
  !     a  : first  input matrix
  !     b  : second input matrix
  !     kin: = multiplies by the laplacian to this power

  USE fprecision
  USE commtypes
  USE mpivars
  USE grid
  USE kes

  IMPLICIT NONE

  COMPLEX(KIND=GP), INTENT(IN), DIMENSION(ny,ista:iend) :: a,b
  DOUBLE PRECISION, INTENT(OUT) :: rslt
  REAL(KIND=GP)                 :: tmp,tmq
  INTEGER, INTENT(IN) :: kin
  INTEGER :: i,j

  tmp = 0.0_GP
  tmq = 1.0_GP/(2.0_GP*real(nx,kind=GP)*real(ny,kind=GP))**2

  IF (kin.eq.0) THEN   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (ista.eq.1) THEN
      DO j = 1,ny
        tmp = tmp+  dble(b(j,1)*conjg(a(j,1)))*tmq
      END DO
      DO i = 2,iend
        DO j = 1,ny
          tmp = tmp+2*dble(b(j,i)*conjg(a(j,i)))*tmq
        END DO
      END DO
    ELSE
      DO i = ista,iend
        DO j = 1,ny
          tmp = tmp+2*dble(b(j,i)*conjg(a(j,i)))*tmq
        END DO
      END DO
    ENDIF
  ELSE IF (kin.eq.1) THEN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (ista.eq.1) THEN
      DO j = 1,ny
        tmp = tmp+  kk2(j,1)*dble(b(j,1)*conjg(a(j,1)))*tmq
      END DO
      DO i = 2,iend
        DO j = 1,ny
          tmp = tmp+2*kk2(j,i)*dble(b(j,i)*conjg(a(j,i)))*tmq
        END DO
      END DO
    ELSE
      DO i = ista,iend
        DO j = 1,ny
          tmp = tmp+2*kk2(j,i)*dble(b(j,i)*conjg(a(j,i)))*tmq
        END DO
      END DO
    ENDIF
  ELSE      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (ista.eq.1) THEN
      DO j = 1,ny
        tmp = tmp+  (kk2(j,1)**kin)*dble(b(j,1)*conjg(a(j,1)))*tmq
      END DO
      DO i = 2,iend
        DO j = 1,ny
          tmp = tmp+2*(kk2(j,i)**kin)*dble(b(j,i)*conjg(a(j,i)))*tmq
        END DO
      END DO
    ELSE
      DO i = ista,iend
        DO j = 1,ny
          tmp = tmp+2*(kk2(j,i)**kin)*dble(b(j,i)*conjg(a(j,i)))*tmq
        END DO
      END DO
    ENDIF
  ENDIF
  CALL MPI_REDUCE(tmp,rslt,1,GC_REAL,MPI_SUM,0, &
  MPI_COMM_WORLD,ierr)

  RETURN
END SUBROUTINE inerprod

!*****************************************************************
SUBROUTINE hdcheck(a,b,d,e,t,ordvf,ordvh)
  !-----------------------------------------------------------------
  !
  ! Consistency check for the conservation of energy in HD 2D
  !
  ! Parameters
  !     a  : streamfunction
  !     t  : time
  !     b  : temp pertubation (theta)
  !
  USE fprecision
  USE commtypes
  USE kes
  USE grid
  USE mpivars
  USE io
  USE fft

  IMPLICIT NONE

  COMPLEX(KIND=GP), INTENT(IN), DIMENSION(ny,ista:iend) :: a,b,d,e
  COMPLEX(KIND=GP), DIMENSION(ny,ista:iend) :: c1,c2
  DOUBLE PRECISION  :: enk,denk,fenk
  DOUBLE PRECISION  :: enp,denp
  DOUBLE PRECISION  :: enk2,denk2,fenk2
  DOUBLE PRECISION  :: enp2,denp2
  DOUBLE PRECISION :: enzm,enfl
  REAL(KIND=GP) :: t
  REAL(KIND=GP) :: nrm,nrm2,tmp1,tmp2
  INTEGER :: ordvf,ordvh
  INTEGER :: i,j

  !
  ! Computes the mean energy
  !
  CALL energy(a,enk,1)
  CALL energy(a,denk,1+ordvf)
  CALL derivk2(a,c1,1)       ! v
  CALL inerprod(c1,b,0,fenk)  ! integrates v*theta

  CALL energy(d,enk2,1)
  CALL energy(d,denk2,1+ordvf)
  CALL derivk2(d,c1,1)
  CALL inerprod(c1,e,0,fenk2)

  !
  ! Computes the mean potential energy
  !

  CALL energy(b,enp,0)
  CALL energy(b,denp,ordvf)

  CALL energy(e,enp2,0)
  CALL energy(e,denp2,ordvf)

  !
  ! For Analysis of 0,1
  !
  CALL laplak2(a,c1)     ! nabla^2 psi
  CALL poisson(d,c1,c1)  ! {phi, nabla^2 psi}
  CALL laplak2(d,c2)     ! nabla^2 phi
  CALL poisson(a,c2,c2)  ! {psi, nabla^2 phi}

  !
  ! Creates external files to store the results
  !
  IF (myrank.eq.0) THEN
    nrm = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP))
    OPEN(1,file=trim(cdir)//'/kenergy.txt',position='append')
    WRITE(1,20) t,enk,denk,fenk,REAL(a(1,1))*nrm
    20    FORMAT( E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3 )
    CLOSE(1)
    OPEN(1,file=trim(cdir)//'/penergy.txt',position='append')
    WRITE(1,20) t,enp,denp
    CLOSE(1)

    OPEN(1,file=trim(cdir)//'/kenergy2.txt',position='append')
    WRITE(1,20) t,enk2,denk2,fenk2
    CLOSE(1)
    OPEN(1,file=trim(cdir)//'/penergy2.txt',position='append')
    WRITE(1,20) t,enp2,denp2
    CLOSE(1)

    OPEN(1,file=trim(cdir)//'/Z0An.txt',position='append')
    WRITE(1,20) t,REAL(d(1,1))*nrm, REAL(c1(1,1))*nrm, REAL(c2(1,1))*nrm
    CLOSE(1)
  ENDIF

  RETURN
END SUBROUTINE hdcheck

!*****************************************************************
SUBROUTINE hdchecktwo(a,t,delta)
  !-----------------------------------------------------------------
  !
  ! Consistency check for the conservation of energy in HD 2D
  !
  ! Parameters
  !     a  : streamfunction
  !     t  : time
  !     b  : temp pertubation (theta)
  !
  USE fprecision
  USE commtypes
  USE kes
  USE grid
  USE mpivars
  USE io
  USE fft

  IMPLICIT NONE

  COMPLEX(KIND=GP), INTENT(IN), DIMENSION(ny,ista:iend) :: a
  REAL(KIND=GP) :: t
  REAL(KIND=GP) :: nrm
  REAL(KIND=GP) :: delta


  IF (myrank.eq.0) THEN
    nrm = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP))
    OPEN(1,file=trim(cdir)//'/conv.txt',position='append')
    WRITE(1,20) t, delta,REAL(a(1,2))*nrm,AIMAG(a(1,2))*nrm,REAL(a(1,4))*nrm,AIMAG(a(1,4))*nrm
    20    FORMAT( E33.24E3,E33.24E3,E33.24E3,E33.24E3,E33.24E3,E33.24E3 )
    CLOSE(1)
  ENDIF

  RETURN
END SUBROUTINE hdchecktwo

!*****************************************************************
SUBROUTINE zonalmean(a,ext,kin)
  !-----------------------------------------------------------------
  !
  ! Computes the zonal mean velocity
  !
  ! Parameters
  !     a  : thing we want to take zonal mean of.
  !     kin: 1 Computes the zonal mean u-velocity
  !          2 Computes the zonal temperature pertubation
  !          else Computes the zonal mean stream function

  USE fprecision
  USE commtypes
  USE kes
  USE grid
  USE mpivars
  USE io
  USE fft

  IMPLICIT NONE

  COMPLEX(KIND=GP), DIMENSION(ny,ista:iend) :: a
  COMPLEX(KIND=GP), DIMENSION(ny,ista:iend) :: c1,c2
  REAL(KIND=GP),    DIMENSION(nx,jsta:jend) :: r1
  REAL(KIND=GP),    DIMENSION(ny) :: bloc,b
  DOUBLE PRECISION :: nrm,nrm2
  INTEGER :: i,j
  CHARACTER*4 :: ext
  INTEGER :: kin

  nrm  = 1.0_GP/(2.0_GP*real(nx,kind=GP)*real(ny,kind=GP))
  nrm2 = 1.0_GP/(real(nx,kind=GP))


  IF (kin.eq.1) Then
    DO j = 1,ny
      bloc(j) = 0.0d0
    END DO

    DO i = ista,iend
      DO j = 1,ny
        c1(j,i) = a(j,i)*nrm
      END DO
    END DO
    CALL derivk2(c1,c2,2)  !NB: psi_y = -u
    DO i = ista,iend
      C1(1,i) = 0.0_GP
      DO j = 2,ny
        C1(j,i) = C2(j-1,i)
      END DO
    END DO
    CALL fftp2d_complex_to_real(plancc,c1,r1,MPI_COMM_WORLD)

    DO j = jsta,jend
      DO i = 1,nx
        bloc(j) = bloc(j) - r1(i,j)*nrm2  !NB: u = -psi_y
      END DO
    END DO

    !
    ! Computes the reduction between nodes
    ! and exports the result to a file
    !
    CALL MPI_REDUCE(bloc,b,ny,GC_REAL,MPI_SUM,0, &
    MPI_COMM_WORLD,ierr)
    IF (myrank.eq.0) THEN
      OPEN(1,file=trim(sdir) // '/zonalmean.' // ext // '.txt')
      WRITE(1,20) b
      20    FORMAT( E23.15 )
      CLOSE(1)
    ENDIF

  ELSE IF (kin.eq.2) Then
    DO j = 1,ny
      bloc(j) = 0.0d0
    END DO

    DO i = ista,iend
      DO j = 1,ny
        c1(j,i) = a(j,i)*nrm
      END DO
    END DO

    CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)

    DO j = jsta,jend
      DO i = 1,nx
        bloc(j) = bloc(j) + r1(i,j)*nrm2
      END DO
    END DO

    !
    ! Computes the reduction between nodes
    ! and exports the result to a file
    !
    CALL MPI_REDUCE(bloc,b,ny,GC_REAL,MPI_SUM,0, &
    MPI_COMM_WORLD,ierr)
    IF (myrank.eq.0) THEN
      OPEN(1,file=trim(sdir) // '/Zonaltheta.' // ext // '.txt')
      WRITE(1,20) b
      CLOSE(1)
    ENDIF

  ELSE


    ! Phil: Normalising
    DO i = ista,iend
      DO j = 1,ny
        c1(j,i) = a(j,i)*nrm
      END DO
    END DO



    CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)

    DO j = jsta,jend
      DO i = 1,nx
        bloc(j) = bloc(j) + r1(i,j)*nrm2
      END DO
    END DO

    !
    ! Computes the reduction between nodes
    ! and exports the result to a file
    !
    CALL MPI_REDUCE(bloc,b,ny,GC_REAL,MPI_SUM,0, &
    MPI_COMM_WORLD,ierr)
    IF (myrank.eq.0) THEN
      OPEN(1,file=trim(sdir) // '/psmean.' // ext // '.txt')
      WRITE(1,20) b
      CLOSE(1)
    ENDIF

  ENDIF

  RETURN
END SUBROUTINE zonalmean

!*****************************************************************
SUBROUTINE spectrum(a,b,ext,name)
  !-----------------------------------------------------------------
  !
  ! Computes the energy power spectrum in 2D.
  ! The output is written to a file by the first node.
  !
  ! Parameters
  !     a  : streamfunction
  !     b  : theta
  !     ext: the extension used when writting the file
  !
  USE fprecision
  USE commtypes
  USE mpivars
  USE grid
  USE kes
  USE io
  USE boxsize

  IMPLICIT NONE

  COMPLEX(KIND=GP), INTENT(IN), DIMENSION(ny,ista:iend) :: a,b
  DOUBLE PRECISION, DIMENSION(specmax) :: Ek,Ektot,Ep, Eptot
  REAL(KIND=GP) :: tmp,two
  INTEGER       :: kin
  INTEGER       :: kmn
  INTEGER       :: i,j
  CHARACTER(len=4), INTENT(IN) :: ext
  CHARACTER(len=6), INTENT(IN) :: name

  !
  ! Sets Ek to zero
  !
  DO i = 1,specmax
    Ek(i) = 0.0d0
    Ep(i) = 0.0d0
  END DO
  !
  ! Computes the kinetic energy spectrum
  !
  tmp = 1.0_GP/(2.0_GP*real(nx,kind=GP)*real(ny,kind=GP))**2
  DO i = ista,iend
    two=2.0d0
    IF (i.eq.1) two=1.0d0
    DO j = 1,ny
      kmn = NINT(sqrt(kk2(j,i))) ! nint round to the nearest integare, we go up in steps of 1, nint rounds upwards. Maybe should change this so that kk2 is just integers? At the moment is is normalised wrt d and L.
      IF (kmn.eq.0) kmn = 1 ! Making sure kmn is not 0
      IF (kmn.le.specmax) THEN
        Ek(kmn) = Ek(kmn)+two*kk2(j,i)*abs(a(j,i))**2*tmp
        Ep(kmn) = Ep(kmn)+two*abs(b(j,i))**2*tmp
      ENDIF
    END DO
  END DO
  !
  ! Computes the reduction between nodes
  ! and exports the result to a file
  !
  CALL MPI_REDUCE(Ek,Ektot,specmax,GC_REAL,MPI_SUM,0, &
  MPI_COMM_WORLD,ierr)
  CALL MPI_REDUCE(Ep,Eptot,specmax,GC_REAL,MPI_SUM,0, &
  MPI_COMM_WORLD,ierr)
  IF (myrank.eq.0) THEN
    OPEN(1,file=trim(sdir) // '/Ek' // name // '.'  // ext // '.txt')
    DO i = 1,specmax
      WRITE(1,20) 1.0D0*i,Ektot(i) ! Devided Ektot by Dkk before? It was for formatting reasons i guess. Now we just mulitply by 1.0D0 so that we go up in integer steps. Equallly, could get rid of 1.0D0 and format the below as an integer
    ENDDO
    20    FORMAT( E23.15, E23.15 )
    CLOSE(1)
    OPEN(1,file=trim(sdir) // '/Ep' // name // '.' // ext // '.txt')
    DO i = 1,specmax
      WRITE(1,20) 1.0d0*i,Eptot(i)
    ENDDO
    CLOSE(1)
  ENDIF

  RETURN
END SUBROUTINE spectrum

! This should work not, although Ek and Ep are alot bigger than they need to be due to allising. Does not really matter though
! could change the defintion for nmax so that it is as big as it needs to be.

!*****************************************************************
SUBROUTINE spectrum2D(a,b,ext)
  !-----------------------------------------------------------------
  !
  ! Computes the energy power spectrum in 2D.
  ! The output is written to a file by the first node.
  !
  ! Parameters
  !     a  : streamfunction
  !     b  : vector potential
  !     ext: the extension used when writting the file
  !     kin: =1 computes the kinetic spectrum
  !          =0 computes the magnetic spectrum
  !
  USE fprecision
  USE commtypes
  USE mpivars
  USE grid
  USE kes
  use io
  USE boxsize

  IMPLICIT NONE

  DOUBLE PRECISION, DIMENSION((kxmax + 1),kymax) :: Ek,Ektot,Ep,Eptot !The + 1 is there since we have 0 in kx
  COMPLEX(KIND=GP), INTENT(IN), DIMENSION(ny,ista:iend) :: a, b
  DOUBLE PRECISION        :: tmp,two
  INTEGER     :: kin
  INTEGER     :: kmnx,kmny
  INTEGER     :: i,j
  CHARACTER*4 :: ext

  !
  ! Sets Ek to zero
  !
  DO i = 1,kxmax + 1
    DO j = 1,kymax
      Ek(i,j) = 0.0d0
      Ep(i,j) = 0.0d0
    END DO
  END DO
  tmp = 1.0_GP/(2.0_GP*real(nx,kind=GP)*real(ny,kind=GP))**2
  DO i = ista,iend
    two=2.0d0
    IF (i.eq.1) two=1.0d0
    DO j = 1,ny
      kmnx = int(abs(kx(i))/Dkx+1.5d0) !int rounds dows, so 0 -> 1, 1 -> 2 here.
      kmny = int(abs(ky(j))/Dky+.5d0) !Dky NOT sure
      IF ((kmnx.gt.0).and.(kmnx.le.(nx/3+1))) THEN
        IF ((kmny.gt.0).and.(kmny.le.(2*ny/3+1))) THEN
          Ek(kmnx,kmny) = Ek(kmnx,kmny)+two*kk2(j,i)*abs(a(j,i))**2*tmp
          Ep(kmnx,kmny) = Ep(kmnx,kmny)+two*abs(b(j,i))**2*tmp
        ENDIF
      ENDIF
    END DO
  END DO
  CALL MPI_REDUCE(Ek,Ektot,(kxmax + 1)*kymax,GC_REAL,MPI_SUM, &
  0,MPI_COMM_WORLD,ierr)
  CALL MPI_REDUCE(Ep,Eptot,(kxmax + 1)*kymax,GC_REAL,MPI_SUM, &
  0,MPI_COMM_WORLD,ierr)
  IF (myrank.eq.0) THEN
    OPEN(1,file=trim(sdir) // '/spectrum2D_UU.' // ext // '.out',form='unformatted')
    WRITE(1) Ektot
    CLOSE(1)
    OPEN(1,file=trim(sdir) // '/spectrum2D_PP.' // ext // '.out',form='unformatted')
    WRITE(1) Eptot
    CLOSE(1)
  ENDIF

  RETURN
END SUBROUTINE spectrum2D

!*****************************************************************
SUBROUTINE vectrans(a,b,c,ext1,ext2)
  !-----------------------------------------------------------------
  !
  ! Computes the square vector potential transfer in
  ! Fourier space in 2D MHD. The output is written
  ! to a file by the first node.
  !
  ! Parameters
  !     a  : streamfunction
  !     b  : vector potential
  !     ext: the extension used when writting the file
  !
  USE fprecision
  USE commtypes
  USE mpivars
  USE grid
  USE kes
  use io
  USE boxsize

  IMPLICIT NONE

  COMPLEX(KIND=GP), INTENT(IN), DIMENSION(ny,ista:iend) :: a,b,c
  COMPLEX(KIND=GP), DIMENSION(ny,ista:iend) :: d
  DOUBLE PRECISION, DIMENSION(nmax/2+1) :: Ek,Ektot
  REAL(KIND=GP)    :: tmp
  INTEGER          :: kmn,kin
  INTEGER          :: i,j
  CHARACTER(len=3) :: ext1
  CHARACTER(len=4) :: ext2

  !
  ! Sets Ek to zero
  !
  DO i = 1,nmax/2+1
    Ek(i) = 0.0_GP
  END DO
  !
  ! Computes the square vector potential flux
  !
  tmp = 1.0_GP/(2.0_GP*real(nx,kind=GP)*real(ny,kind=GP))**2
  CALL poisson(b,c,d)
  IF (ista.eq.1) THEN
    DO j = 1,ny
      kmn = int(sqrt(kk2(j,1))/Dkk+.5)
      IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
        Ek(kmn) = Ek(kmn)+dble(a(j,1)*conjg(d(j,1)))*tmp
      ENDIF
    END DO
    DO i = 2,iend
      DO j = 1,ny
        kmn = int(sqrt(kk2(j,i))/Dkk+.5)
        IF (kmn.le.nmax/2+1) THEN
          Ek(kmn) = Ek(kmn)+2*dble(a(j,i)*conjg(d(j,i)))*tmp
        ENDIF
      END DO
    END DO
  ELSE
    DO i = ista,iend
      DO j = 1,ny
        kmn = int(sqrt(kk2(j,i))/Dkk+.5)
        IF (kmn.le.nmax/2+1) THEN
          Ek(kmn) = Ek(kmn)+2*dble(a(j,i)*conjg(d(j,i)))*tmp
        ENDIF
      END DO
    END DO
  ENDIF
  !
  ! Computes the reduction between nodes
  ! and exports the result to a file
  !
  CALL MPI_REDUCE(Ek,Ektot,nmax/2+1,GC_REAL,MPI_SUM,0, &
  MPI_COMM_WORLD,ierr)
  IF (myrank.eq.0) THEN
    OPEN(1,file=trim(tdir) // '/transfer_' // ext1 // '.' // ext2 // '.txt')
    DO i=1,nmax/2
      WRITE(1,20) Dkk*i,Ektot(i)/Dkk
    ENDDO
    20    FORMAT( E23.15, E23.15 )
    CLOSE(1)
  ENDIF

  RETURN
END SUBROUTINE vectrans

!*****************************************************************
SUBROUTINE CFL_condition(a,nu,kappa,cfl,ord,vord,Lx,Ly,dt)
  !-----------------------------------------------------------------
  !
  ! Checks the CFL condition
  !
  ! Parameters
  !     cfl : cfl factor
  !      a : stream fun
  !      b : vector pot
  !      b0 : external

  USE fprecision
  USE commtypes
  USE mpivars
  USE kes
  USE grid
  USE ali
  USE fft
  USE var

  IMPLICIT NONE

  COMPLEX(KIND=GP), DIMENSION(ny,ista:iend) :: a
  COMPLEX(KIND=GP), DIMENSION(ny,ista:iend) :: c1,c2,c3
  REAL(KIND=GP), DIMENSION(nx,jsta:jend)    :: r1,r2!,r3
  INTEGER :: i,j
  INTEGER :: ord
  INTEGER :: vord
  DOUBLE PRECISION :: nu,kappa
  DOUBLE PRECISION :: dt,cfl
  DOUBLE PRECISION :: umax,dx
  DOUBLE PRECISION :: tmp,tmq
  DOUBLE PRECISION :: Lx,Ly

  tmq = 1.0_GP/(2.0_GP*real(nx,kind=GP)*real(ny,kind=GP))

  CALL derivk2(a,c1,1)
  CALL derivk2(a,c2,2)
  DO i = ista,iend
    c3(1,i) = 0.0_GP
    DO j = 2,nx
      c3(j,i) = c2(j-1,i)
    END DO
  END DO
  CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
  CALL fftp2d_complex_to_real(plancc,c3,r2,MPI_COMM_WORLD)
  !u^2
  tmp = 0.0d0
  DO j = jsta,jend
    DO i = 1,ny
      !r3(i,j) = r1(i,j)*r1(i,j)+r2(i,j)*r2(i,j)
      umax = r1(i,j)*r1(i,j)+r2(i,j)*r2(i,j)
      if (tmp.lt.umax) tmp = umax
    END DO
  END DO
  !tmp=maxval(r3)
  call MPI_REDUCE(tmp,umax,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(umax,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
  umax = sqrt(umax)*tmq ! Phil: So umax is the greatest value of sqrt(u^2) * tmq in any square

  ! See Lele JCP v103 p16 (1992) p.32 for the CFL factors
  !!!      kcut = max(real(nx,kind=GP)/lx, 2.0_GP*real(ny,kind=GP)/ly)/3.0_GP
  !dx = min(Ly*pi*3.0_GP/(2.0_GP*real(ny,kind=GP)), Lx*pi*6.0_GP/(real(nx,kind=GP)))
  dx = 3.0_GP/real(nx,kind=GP) ! back to what we had before...
  !kcut = real(nx,kind=GP)/3.0_GP
  ! If it doesnt work multiply by a cfl factor < 1
  if (ord.eq.2) then
    dt = cfl*min((sqrt(3.0d0)*dx/pi)/umax, &
    (2.0d0*dx**(2*vord)/pi**2)/nu, &
    (2.0d0*dx**(2*vord)/pi**2)/kappa)
  else if (ord.eq.3) then
    dt = cfl*min((sqrt(3.0d0)*dx/pi)/umax, &
    (2.5d0*dx**(2*vord)/pi**2)/nu, &
    (2.5d0*dx**(2*vord)/pi**2)/kappa)
  end if

  RETURN
END SUBROUTINE CFL_condition
