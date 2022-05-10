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
! 2003 Pablo D. Mininni.
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
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
      USE ali
      USE kes
      USE var
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b
      INTEGER :: dir
      INTEGER :: i,j

!
! Derivative in the x-direction
!
      IF (dir.eq.1) THEN
         DO i = ista,iend
            DO j = 1,n
               IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
                  b(j,i) = im*ka(i)*a(j,i)
               ELSE
                  b(j,i) = 0.0d0
               ENDIF
            END DO
         END DO
!
! Derivative in the y-direction
!
      ELSE
         DO i = ista,iend
            DO j = 1,n
               IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
                  b(j,i) = im*ka(j)*a(j,i)
               ELSE
                  b(j,i) = 0.0d0
               ENDIF
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
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b
      INTEGER :: i,j

      DO i = ista,iend
         DO j = 1,n
            b(j,i) = -ka2(j,i)*a(j,i)
         END DO
      END DO

      RETURN
      END SUBROUTINE laplak2

!*****************************************************************
      SUBROUTINE add_bb0(b0,a,b)
!-----------------------------------------------------------------
!
! Parameters
!     a: input matrix
!     b: at the output contains the Laplacian d2a/dka2
!
      USE kes
      USE grid
      USE var
      USE ali
      USE mpivars
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b
      INTEGER :: i,j
      DOUBLE PRECISION :: b0

      DO i = ista,iend
         DO j = 1,n
            IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
            b(j,i) = b(j,i)-im*ka(i)*a(j,i)*b0
            ELSE
            b(j,i) = 0.0d0
            ENDIF
         END DO
      END DO

      RETURN
      END SUBROUTINE add_bb0


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
      USE mpivars
      USE kes
      USE ali
      USE grid
      USE fft
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b,c
      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: c1,c2
      DOUBLE PRECISION, DIMENSION(n,jsta:jend)    :: r1,r2,r3
      INTEGER :: i,j

!
! Computes dA/dx.dB/dy
!
      CALL derivk2(a,c1,1)
      CALL derivk2(b,c2,2)
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      DO j = jsta,jend
         DO i = 1,n
            r3(i,j) = r1(i,j)*r2(i,j)
         END DO
      END DO

!
! Computes dA/dy.dB/dx
!
      CALL derivk2(a,c1,2)
      CALL derivk2(b,c2,1)
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      DO j = jsta,jend
         DO i = 1,n
            r3(i,j) = (r3(i,j)-r1(i,j)*r2(i,j))/dble(n)**4
         END DO
      END DO
      CALL fftp2d_real_to_complex(planrc,r3,c,MPI_COMM_WORLD)
      DO i = ista,iend
         DO j = 1,n
            IF ((ka2(j,i).ge.kmax).and.(ka2(j,i).le.tiny)) THEN
               c(j,i) = 0.0d0
            ENDIF
         END DO
      END DO
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
      USE kes
      USE grid
      USE mpivars
      USE ali
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a
      DOUBLE PRECISION  :: b
      DOUBLE PRECISION  :: bloc
      DOUBLE PRECISION  :: tmp
      INTEGER :: kin
      INTEGER :: i,j

      bloc = 0.0d0
      tmp = 1.0d0/dble(n)**4 ! Phil: So we are returning 2<a^2> (See parseval theorem notes)

!
! Computes the square of the scalar field
!
      IF (kin.eq.0) THEN
         IF (ista.eq.1) THEN
            DO j = 1,n
               bloc = bloc+abs(a(j,1))**2*tmp ! Phil: this is running throuhg the modes (0,-n), (0,0), (0,n). No extra c.c to consider.
            END DO
            DO i = 2,iend
               DO j = 1,n
                  bloc = bloc+2*abs(a(j,i))**2*tmp ! Phil: this is running through the modes (kx,-n), (kx,0), (kx,n), kx \neq 0. Have energy contribution from (-kx,-n), (-kx,0) (-kx,n). Hence the 2.
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  bloc = bloc+2*abs(a(j,i))**2*tmp
               END DO
            END DO
         ENDIF
!
! Computes the energy
!
      ELSE IF (kin.eq.1) THEN
         IF (ista.eq.1) THEN
            DO j = 1,n
               bloc = bloc+ka2(j,1)*abs(a(j,1))**2*tmp
            END DO
            DO i = 2,iend
               DO j = 1,n
                  bloc = bloc+2*ka2(j,i)*abs(a(j,i))**2*tmp
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  bloc = bloc+2*ka2(j,i)*abs(a(j,i))**2*tmp
               END DO
            END DO
         ENDIF
!
! Computes the current or vorticity
!
      ELSE
         IF (ista.eq.1) THEN
            DO j = 1,n
               i=1
               IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
               bloc = bloc+ka2(j,1)**kin*abs(a(j,1))**2*tmp
               ENDIF
            END DO
            DO i = 2,iend
               DO j = 1,n
                  IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
                  bloc = bloc+2*ka2(j,i)**kin*abs(a(j,i))**2*tmp
                  ENDIF
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
                  bloc = bloc+2*ka2(j,i)**kin*abs(a(j,i))**2*tmp
                  ENDIF
               END DO
            END DO
         ENDIF
      ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(bloc,b,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE energy

!###############################
!*****************************************************************
      SUBROUTINE inerprod(a,b,kin,rslt)
!-----------------------------------------------------------------
! Parameters
!     a  : first  input matrix
!     b  : second input matrix
!     kin: = multiplies by the laplacian to this power

      USE kes
      USE grid
      USE mpivars
      USE ali

      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b
      DOUBLE PRECISION  :: tmp,tmq
      DOUBLE PRECISION  :: rslt
      INTEGER :: kin
      INTEGER :: i,j

      tmp = 0.0d0
      tmq = 1./dble(n)**4

      IF (kin.eq.0) THEN   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (ista.eq.1) THEN
         DO j = 1,n
            tmp = tmp+  dble(b(j,1)*conjg(a(j,1)))*tmq
         END DO
         DO i = 2,iend
            DO j = 1,n
            tmp = tmp+2*dble(b(j,i)*conjg(a(j,i)))*tmq
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
            tmp = tmp+2*dble(b(j,i)*conjg(a(j,i)))*tmq
            END DO
         END DO
      ENDIF
      ELSE IF (kin.eq.1) THEN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (ista.eq.1) THEN
         DO j = 1,n
            tmp = tmp+  (ka2(j,1))*dble(b(j,1)*conjg(a(j,1)))*tmq
         END DO
         DO i = 2,iend
            DO j = 1,n
            tmp = tmp+2*(ka2(j,i))*dble(b(j,i)*conjg(a(j,i)))*tmq
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
            tmp = tmp+2*(ka2(j,i))*dble(b(j,i)*conjg(a(j,i)))*tmq
            END DO
         END DO
      ENDIF
      ELSE      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (ista.eq.1) THEN
         DO j = 1,n
            i=1
            IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
            tmp = tmp+  (ka2(j,1)**kin)*dble(b(j,1)*conjg(a(j,1)))*tmq
            ENDIF
         END DO
         DO i = 2,iend
            DO j = 1,n
            IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
            tmp = tmp+2*(ka2(j,i)**kin)*dble(b(j,i)*conjg(a(j,i)))*tmq
            ENDIF
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
            IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
            tmp = tmp+2*(ka2(j,i)**kin)*dble(b(j,i)*conjg(a(j,i)))*tmq
            ENDIF
            END DO
         END DO
      ENDIF
      ENDIF
      CALL MPI_REDUCE(tmp,rslt,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE inerprod
!################################


!*****************************************************************
      SUBROUTINE hdcheck(a,b,c,t,ordvf,ordvh)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of energy in HD 2D
!
! Parameters
!     a  : streamfunction
!     b  : theta
!     c  : external force
!     t  : time
!
      USE kes
      USE grid
      USE mpivars
      use io

      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b,c1,c
      DOUBLE PRECISION  :: enk,denk,fenk,henk
      DOUBLE PRECISION  :: enp,denp,fenp
      DOUBLE PRECISION  :: t
      DOUBLE PRECISION  :: tmq,tmp
      INTEGER :: i,j,ordvf,ordvh

!
! Computes the mean energy
!
      CALL energy(a,enk,1)
      CALL energy(a,denk,1+ordvf)
      CALL derivk2(a,c1,1)       ! v
      CALL inerprod(c1,b,0,fenk)  ! integrates v*theta
      CALL energy(a,henk,1-ordvh)

!
! Computes the mean potential energy
!

      CALL energy(b,enp,0)       ! integrates tt^2
      CALL energy(b,denp,ordvf)
      CALL inerprod(b,c,0,fenp)  ! integrates tt*f


!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(cdir)//'/kenergy.txt',position='append')
         WRITE(1,20) t,enk,denk,fenk,henk
   20    FORMAT( E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3 )
         CLOSE(1)
         OPEN(1,file=trim(cdir)//'/penergy.txt',position='append')
         WRITE(1,20) t,enp,denp,fenp
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE hdcheck

!*****************************************************************
      SUBROUTINE mhdcheck(a,b,c,d,t,eng,kup,ordvf,ordvh)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of energy in MHD 2D
!
! Parameters
!     a  : streamfunction
!     b  : vector potential
!     c  : external kinetic force
!     d  : external magnetic force
!     t  : time
!
      USE kes
      USE grid
      USE mpivars
      use io

      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b,c,d
      DOUBLE PRECISION  :: eng, kup
      DOUBLE PRECISION  :: enk, denk, fenk, henk
      DOUBLE PRECISION  :: enm, denm, fenm, henm
      DOUBLE PRECISION  :: ens, dens, fens, hens
      DOUBLE PRECISION  :: crs, dcrs, fcrs ,ecrs, hcrs
      DOUBLE PRECISION  :: asq, dasq, fasq, hasq
      DOUBLE PRECISION  :: t,tmp1,tmp2,tmp3
      DOUBLE PRECISION  :: euf,ebf,two
      INTEGER :: i,j,ordvf,ordvh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ENERGY BALANCE
      CALL energy(a, enk, 1)
      CALL energy(a,denk, 1+ordvf)
      !CALL energy(a,henk, 1-ordvh)
!!$      CALL energy(b, enm, 1)
!!$      CALL energy(b,denm, 1+ordvf)
!!$      !CALL energy(b,henm, 1-ordvh)
      CALL inerprod(c,a,1,fenk)
!!$      CALL inerprod(d,b,1,fenm)
!!$      eng = enk+enm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ENSTROPHY BALANCE
      ens=denk
      CALL energy(a,dens,2+ordvf)
      CALL inerprod(c,a,2,fens)
      !CALL energy(a,hens,2-ordvh)
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!     CROSS HELICITY BALANCE
!!$      CALL inerprod(b,a, 1, crs)
!!$      CALL inerprod(b,a, 1+ordvf,dcrs)
!!$      !CALL inerprod(b,a, 1-ordvh,hcrs)
!!$      CALL inerprod(a,d, 1,fcrs)
!!$      CALL inerprod(b,c, 1,ecrs)
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!     SQUARE VECTOR POTENTIAL
!!$      CALL energy(b, asq,0)
!!$      CALL energy(b,dasq, ordvf)
!!$      CALL inerprod(b,d,0,fasq)
!!$      !CALL energy(b,hasq,-ordvh)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!     CALCULATES ENERGIES AT FORCING SCALE
!!!      tmp1=0.0d0
!!!      tmp2=0.0d0
!!!      tmp3=1.0d0/dble(n)**4
!!!         DO i = ista,iend
!!!            DO j = 1,n
!!!            IF ((ka2(j,i).le.2.01*kup*kup).and.(ka2(j,i).ge.kup*kup )) THEN
!!!!           print*,"DBG",myrank,kup,2.01*kup*kup,ka2(j,i),kup*kup,ka2(j,i)*abs(a(j,i))**2*tmp3
!!!            two=2.0d0
!!!            IF (i.eq.1) two=1.0d0
!!!            tmp1 = tmp1+  two*ka2(j,i)*abs(a(j,i))**2*tmp3
!!!            tmp2 = tmp2+  two*ka2(j,i)*abs(b(j,i))**2*tmp3
!!!            ENDIF
!!!            END DO
!!!         END DO
!!!         CALL MPI_REDUCE(tmp1,euf,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
!!!                      MPI_COMM_WORLD,ierr)
!!!         CALL MPI_REDUCE(tmp2,ebf,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
!!!                      MPI_COMM_WORLD,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(cdir)//'/kenergy.txt',position='append')
         WRITE(1,10) t,enk,denk,fenk!,henk
   10    FORMAT( E23.14E3,E23.14E3,E23.14E3,E23.14E3)
!   10    FORMAT( E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3)
         CLOSE(1)
!!$         OPEN(1,file=trim(cdir)//'/menergy.txt',position='append')
!!$         WRITE(1,10) t,enm,denm,fenm!,henm
!!$         CLOSE(1)
         OPEN(1,file=trim(cdir)//'/enstrophy.txt',position='append')
         WRITE(1,10) t,ens,dens,fens!,hens
         CLOSE(1)
!!$         OPEN(1,file=trim(cdir)//'/sqr_vecpot.txt',position='append')
!!$         WRITE(1,10) t,asq,dasq,fasq!,hasq
!!$         CLOSE(1)
!!$         OPEN(1,file=trim(cdir)//'/cross.txt',position='append')
!!$         WRITE(1,11) t,crs,dcrs,fcrs,ecrs!,hcrs
!!$   11    FORMAT( E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3)
!!$!   11    FORMAT( E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3)
!!$         CLOSE(1)
!!!         OPEN(1,file=trim(cdir)//'/engf.txt',position='append')
!!!         WRITE(1,12) t,euf,ebf
!!!   12    FORMAT( E23.14E3,E23.14E3,E23.14E3)
!!!         CLOSE(1)

      ENDIF

      RETURN
      END SUBROUTINE mhdcheck

!*****************************************************************
      SUBROUTINE spectrum(a,kin,ext,ext2)
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
      USE kes
      USE grid
      USE mpivars
      use io

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1)        :: Ek,Ektot
      DOUBLE COMPLEX, DIMENSION(n,ista:iend)    :: a!,b
      DOUBLE PRECISION        :: tmp,two
      INTEGER     :: kin
      INTEGER     :: kmn
      INTEGER     :: i,j
      CHARACTER*4 :: ext
      CHARACTER*6 :: ext2


      tmp = 1.0d0/dble(n)**4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Computes the VELOCITY energy spectrum!!!!!!!!!!!!!!
      DO i = 1,n/2+1
         Ek(i) = 0.0d0
      END DO
      DO i = ista,iend
         two=2.0d0
         IF (i.eq.1) two=1.0d0
         DO j = 1,n
            kmn = int(sqrt(ka2(j,i))+.5d0)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
               Ek(kmn) = Ek(kmn)+two*ka2(j,i)**kin*abs(a(j,i))**2*tmp
            ENDIF
         END DO
      END DO
      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
            OPEN(1,file=trim(sdir) // '/' // ext2 // ext // '.txt')
         WRITE(1,20) Ektot
         CLOSE(1)
      ENDIF
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!!!!!!!!! Computes the MAGNETIC energy spectrum!!!!!!!!!!!!!!
!!$      DO i = 1,n/2+1
!!$         Ek(i) = 0.0d0
!!$      END DO
!!$      DO i = ista,iend
!!$         two=2.0d0
!!$         IF (i.eq.1) two=1.0d0
!!$         DO j = 1,n
!!$            kmn = int(sqrt(ka2(j,i))+.5d0)
!!$            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!!$               Ek(kmn) = Ek(kmn)+two*ka2(j,i)*abs(b(j,i))**2*tmp
!!$            ENDIF
!!$         END DO
!!$      END DO
!!$      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
!!$                      MPI_COMM_WORLD,ierr)
!!$      IF (myrank.eq.0) THEN
!!$            OPEN(1,file=trim(sdir) // '/spectrum_bb.' // ext // '.txt')
!!$         WRITE(1,20) Ektot
!!$         CLOSE(1)
!!$      ENDIF
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!!!!!!!!! Computes the ELS1 energy spectrum!!!!!!!!!!!!!!
!!$      DO i = 1,n/2+1
!!$         Ek(i) = 0.0d0
!!$      END DO
!!$      DO i = ista,iend
!!$         two=2.0d0
!!$         IF (i.eq.1) two=1.0d0
!!$         DO j = 1,n
!!$            kmn = int(sqrt(ka2(j,i))+.5d0)
!!$            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!!$               Ek(kmn)= Ek(kmn)+two*ka2(j,i)*(abs(b(j,i)+a(j,i)))**2*tmp
!!$            ENDIF
!!$         END DO
!!$      END DO
!!$      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
!!$                      MPI_COMM_WORLD,ierr)
!!$      IF (myrank.eq.0) THEN
!!$            OPEN(1,file=trim(sdir) // '/spectrum_zp.' // ext // '.txt')
!!$         WRITE(1,20) Ektot
!!$         CLOSE(1)
!!$      ENDIF
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!!!!!!!!! Computes the ELS2     energy spectrum!!!!!!!!!!!!!!
!!$      DO i = 1,n/2+1
!!$         Ek(i) = 0.0d0
!!$      END DO
!!$      DO i = ista,iend
!!$         two=2.0d0
!!$         IF (i.eq.1) two=1.0d0
!!$         DO j = 1,n
!!$            kmn = int(sqrt(ka2(j,i))+.5d0)
!!$            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!!$               Ek(kmn)= Ek(kmn)+two*ka2(j,i)*(abs(b(j,i)-a(j,i)))**2*tmp
!!$            ENDIF
!!$         END DO
!!$      END DO
!!$      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
!!$                      MPI_COMM_WORLD,ierr)
!!$      IF (myrank.eq.0) THEN
!!$            OPEN(1,file=trim(sdir) // '/spectrum_zm.' // ext // '.txt')
!!$         WRITE(1,20) Ektot
!!$         CLOSE(1)
!!$      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Computes the ELS1 PERP energy spectrum!!!!!!!!!!!!!!
!      DO i = ista,iend
!        IF (i.le.5) THEN
!          DO j = 1,n/2+1
!           Ek(j) = 0.0d0
!          END DO
!          two=2.0d0
!          IF (i.eq.1) two=1.0d0
!          DO j = 1,n
!             kmn = int(abs(ka(j))+.5d0)
!             IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!               Ek(kmn)= Ek(kmn)+two*ka2(j,i)*(abs(b(j,i)+a(j,i)))**2*tmp
!             ENDIF
!          END DO
!          OPEN(1,file='spectrum_m' // char(i+48) // '.'// ext // '.txt')
!          WRITE(1,20) Ek
!          CLOSE(1)
!        ENDIF
!      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Computes the ELS2 PERP energy spectrum!!!!!!!!!!!!!
!      DO i = ista,iend
!        IF (i.le.5) THEN
!          DO j = 1,n/2+1
!           Ek(j) = 0.0d0
!          END DO
!          two=2.0d0
!          IF (i.eq.1) two=1.0d0
!          DO j = 1,n
!             kmn = int(abs(ka(j))+.5d0)
!             IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!               Ek(kmn)= Ek(kmn)+two*ka2(j,i)*(abs(b(j,i)-a(j,i)))**2*tmp
!             ENDIF
!          END DO
!          OPEN(1,file='spectrum_p' // char(i+48) // '.'// ext // '.txt')
!          WRITE(1,20) Ek
!          CLOSE(1)
!        ENDIF
!      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   20    FORMAT( E23.14E3 )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      RETURN
      END SUBROUTINE spectrum

!*****************************************************************
      SUBROUTINE spectrum2D(a,ext)
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
      USE kes
      USE grid
      USE mpivars
      use io

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n,n)        :: Ek,Ektot
      DOUBLE COMPLEX, DIMENSION(n,ista:iend)  :: a
      DOUBLE PRECISION        :: tmp,two
      INTEGER     :: kin
      INTEGER     :: kx,ky
      INTEGER     :: i,j
      CHARACTER*3 :: ext


      tmp = 1.0d0/dble(n)**4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Computes the VELOCITY energy spectrum!!!!!!!!!!!!!!
      DO i = 1,n
      DO j = 1,n
         Ek(j,i) = 0.0d0
      END DO
      ENDDO
      DO i = ista,iend
         DO j = 1,n
            Ek(j,i) = Ek(j,i)+ka2(j,i)*abs(a(j,i))**2*tmp
         END DO
      END DO
      CALL MPI_REDUCE(Ek,Ektot,n*n,MPI_DOUBLE_PRECISION,&
                                       MPI_SUM,0,MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(sdir) // '/spectrum2D_UU.' // ext // '.dat',form='unformatted')
         WRITE(1) Ektot
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
      USE kes
      USE grid
      USE mpivars
      use io

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1)        :: Ek,Ektot
      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b,c,d
      DOUBLE PRECISION        :: tmp
      INTEGER     :: kmn,kin
      INTEGER     :: i,j
      CHARACTER*3 :: ext1
      CHARACTER*4 :: ext2

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.
      END DO
!
! Computes the square vector potential flux
!
      tmp = 1.0d0/dble(n)**4
      CALL poisson(b,c,d)
      IF (ista.eq.1) THEN
         DO j = 1,n
            kmn = int(sqrt(ka2(j,1))+.5)
            Ek(kmn) = Ek(kmn)+dble(a(j,1)*conjg(d(j,1)))*tmp
         END DO
         DO i = 2,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.5)
               IF (kmn.le.n/2+1) THEN
                  Ek(kmn) = Ek(kmn)+2*dble(a(j,i)*conjg(d(j,i)))*tmp
               ENDIF
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.5)
               IF (kmn.le.n/2+1) THEN
                  Ek(kmn) = Ek(kmn)+2*dble(a(j,i)*conjg(d(j,i)))*tmp
               ENDIF
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(sdir) // '/vectrans_' // ext1 // '.' // ext2 // '.txt')
         WRITE(1,20) Ektot
   20    FORMAT( E23.14E3 )
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE vectrans

      !*****************************************************************
      SUBROUTINE vectrans2(a,b,kin,ext1,ext2)
        !-----------------------------------------------------------------
        !
        ! Computes the square vector potential transfer in
        ! Fourier space in 2D MHD. The output is written
        ! to a file by the first node.
        !
        ! Parameters
        !     a  : field1
        !     b  : field2
        !     kin: laplacian to this power
        !     ext: the extension used when writting the file
        !
        USE kes
        USE grid
        USE mpivars
        use io

        IMPLICIT NONE

        DOUBLE PRECISION, DIMENSION(n/2+1)        :: Ek,Ektot
        DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b
        DOUBLE PRECISION        :: tmp
        INTEGER     :: kmn,kin
        INTEGER     :: i,j
        CHARACTER*9 :: ext1
        CHARACTER*4 :: ext2

        !
        ! Sets Ek to zero
        !
        DO i = 1,n/2+1
          Ek(i) = 0.
        END DO
        !
        ! Computes the square vector potential flux
        !
        tmp = 1.0d0/dble(n)**4
        IF (ista.eq.1) THEN
          DO j = 1,n
            kmn = int(sqrt(ka2(j,1))+.5)
            Ek(kmn) = Ek(kmn)+ka2(j,1)**kin*dble(a(j,1)*conjg(b(j,1)))*tmp
          END DO
          DO i = 2,iend
            DO j = 1,n
              kmn = int(sqrt(ka2(j,i))+.5)
              IF (kmn.le.n/2+1) THEN
                Ek(kmn) = Ek(kmn)+2*ka2(j,i)**kin*dble(a(j,i)*conjg(b(j,i)))*tmp
              ENDIF
            END DO
          END DO
        ELSE
          DO i = ista,iend
            DO j = 1,n
              kmn = int(sqrt(ka2(j,i))+.5)
              IF (kmn.le.n/2+1) THEN
                Ek(kmn) = Ek(kmn)+2*ka2(j,i)**kin*dble(a(j,i)*conjg(b(j,i)))*tmp
              ENDIF
            END DO
          END DO
        ENDIF
        !
        ! Computes the reduction between nodes
        ! and exports the result to a file
        !
        CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
        MPI_COMM_WORLD,ierr)
        IF (myrank.eq.0) THEN
          OPEN(1,file=trim(sdir) // '/vectrans_' // ext1 // '.' // ext2 // '.txt')
          WRITE(1,20) Ektot
          20    FORMAT( E23.14E3 )
          CLOSE(1)
        ENDIF

        RETURN
      END SUBROUTINE vectrans2

!*****************************************************************
      SUBROUTINE CFL_condition(c1,nu,kappa,cfl,ord,vord,dt)
!-----------------------------------------------------------------
!
! Checks the CFL condition
!
! Parameters
!     cfl : cfl factor
!      c1 : stream fun
!      c2 : vector pot
!      b0 : external

      USE mpivars
      USE kes
      USE ali
      USE grid
      USE fft
      USE var

      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend)   :: c1,c3,c4
      DOUBLE PRECISION, DIMENSION(n,jsta:jend) :: r1,r2!,r3
      INTEGER :: i,j
      INTEGER :: ord
      INTEGER :: vord
      DOUBLE PRECISION :: nu,kappa
      DOUBLE PRECISION :: dt,cfl
      DOUBLE PRECISION :: umax,dx
      DOUBLE PRECISION :: tmp,nrm

      nrm=(dble(n))**2

      CALL derivk2(c1,c3,1)
      CALL derivk2(c1,c4,2)
      CALL fftp2d_complex_to_real(plancr,c3,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c4,r2,MPI_COMM_WORLD)

      !u^2
      tmp = 0.0d0
      DO j = jsta,jend
         DO i = 1,n
            !r3(i,j) = r1(i,j)*r1(i,j)+r2(i,j)*r2(i,j)
            umax = r1(i,j)*r1(i,j)+r2(i,j)*r2(i,j)
            if (tmp.lt.umax) tmp = umax
         END DO
      END DO
      !tmp=maxval(r3)
      call MPI_REDUCE(tmp,umax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(umax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      umax = sqrt(umax)/nrm

      ! See Lele JCP v103 p16 (1992) p.32 for the CFL factors
      dx = 2.0d0*pi/dble(n)
      ! If it doesnt work multiply by a cfl factor < 1
      if (ord.eq.2) then
         dt = cfl*min((sqrt(3.0d0)/pi)*dx/umax, &
                       (2.0d0/pi**2)*dx**(2*vord)/nu, &
                       (2.0d0/pi**2)*dx**(2*vord)/kappa)
      else if (ord.eq.3) then
         dt = cfl*min((sqrt(3.0d0)/pi)*dx/umax, &
                      (2.5d0/pi**2)*dx**(2*vord)/nu, &
                      (2.5d0/pi**2)*dx**(2*vord)/kappa)
      end if

      RETURN
      END SUBROUTINE CFL_condition
