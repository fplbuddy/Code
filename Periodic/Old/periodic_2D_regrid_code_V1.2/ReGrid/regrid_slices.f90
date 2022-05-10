!=================================================================
      PROGRAM rd_slices
!=================================================================
!=================================================================
      IMPLICIT NONE
!
!     fields 

      REAL(8), ALLOCATABLE, DIMENSION (:,:)    :: Vold ! a slice of the field
      REAL(8), ALLOCATABLE, DIMENSION (:,:)    :: Vnew ! the total NxNxN field

!
! Some auxiliary parameters

      REAL :: rfile,rnode

      INTEGER :: Nnodes,Mnodes,inode
      INTEGER :: n,ns
      INTEGER :: m,ms
      INTEGER :: is
      INTEGER :: iname
      INTEGER :: ifile
      INTEGER :: i ,k
      INTEGER :: i2,k2
      INTEGER :: ic,id,iu
      INTEGER :: jc,jd,ju

      CHARACTER     :: c,d,u
      CHARACTER*3   :: node,ext
      CHARACTER*100 :: ldir,fdir,odir
      
      OPEN(1,file='input.prm',status='unknown')
      READ(1,*) ifile  ! filenumber
      READ(1,*) n      ! old grid size
      READ(1,*) Nnodes ! old number of nodes (ie number of slices)
      READ(1,*) m      ! new grid size (m>n)
      READ(1,*) Mnodes ! new number of nodes (ie number of slices)
      READ(1,'(a100)') fdir
      READ(1,'(a100)') odir
      CLOSE(1)
      
      ns     = n/Nnodes  ! Old Slice thickness
      ms     = m/Mnodes  ! New Slice thickness (ms=<ns)

      ALLOCATE( Vold(n,ns) )
      ALLOCATE( Vnew(m,ms) )

      print*,'OLD FILES N=',n,' Nodes=',Nnodes,' slice=',ns
      print*,'NEW FILES N=',m,' Nodes=',Mnodes,' slice=',ms
      print*,'ONE OLD SLICE GIVES ',Mnodes/Nnodes,' SLICES'
      print*,'ONE GRID POINT GIVES',m/n,' POINTS'
      !pause

!!@      !print*, 10/0
!!@      print*, 10./0  ! Inf
!!@      print*, 10/0.  ! Inf
!!@      print*, 10./0. ! Inf
!!@      !print*, 0/0
!!@      print*, 0./0   ! NaN
!!@      print*, 0/0.   ! NaN
!!@      print*, 0./0.  ! NaN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DO iname=1,1!2 ! LOOP OVER HEADERS
 
      IF (iname.eq.1) ldir  = 'hd2Dps.'
      IF (iname.eq.2) ldir  = 'hd2Daa.'

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DO inode=0,Nnodes-1  ! LOOP OVER ALL OLD SLICES %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!        Make the file name of slice # inode
         rnode = float(inode)

         ic = 48+int(rnode/100)
         id = 48+int(rnode/10 )-int(rnode/100)*10
         iu = 48+int(rnode)-int(rnode/10)*10
         c = char(ic)
         d = char(id)
         u = char(iu)
         node = c // d // u

         rfile = float(ifile)
         
         jc = 48+int(rfile/100)
         jd = 48+int(rfile/10)-int(rfile/100)*10
         ju = 48+int(rfile)-int(rfile/10)*10
         c = char(jc)
         d = char(jd)
         u = char(ju)
         ext = c // d // u

         print*," "
         print*,"READ SLICE ",trim(fdir)//trim(ldir) & 
                                           //node//'.'//ext//'.dat'
         OPEN(1,file=trim(fdir)//trim(ldir)//node//'.'//ext//'.dat', &
                form='unformatted')
         READ(1) Vold
         CLOSE(1)

         DO is = 1,Mnodes/Nnodes
         !print*,"DBG making slice #",is,is+inode*ns/ms,ns/ms,ns,ms
         !pause

         DO k=1,ms
            k2=(k+1+(is-1)*ms)*n/m
            if (k2.eq.0) k2 = 1
            !print*,"k=",k," k2=",k2
            DO i=1,m
               i2=(i+1)*n/m
               if (i2.eq.0) i2 = 1
               !print*,"i=",i," i2=",i2
               !pause
               Vnew(i,k)=Vold(i2,k2) 
            END DO
         END DO 

         ! INDEX FOR NEW SLICES
         rnode = float(inode*Mnodes/Nnodes+is-1)
         ic = 48+int(rnode/100)
         id = 48+int(rnode/10 )-int(rnode/100)*10
         iu = 48+int(rnode)-int(rnode/10)*10
         c = char(ic)
         d = char(id)
         u = char(iu)
         node = c // d // u

         print*,'WRITE SLICE ', trim(odir) // trim(ldir) &
                                             // node // '.001.dat'
!!!         if (m.gt.1024) then 
!!!            OPEN(1,file= trim(odir)//trim(ldir) // node // '.001.dat', &
!!!                   form='unformatted',convert='BIG_ENDIAN') !JUQUEEN
!!!         else
            OPEN(1,file= trim(odir)//trim(ldir) // node // '.001.dat', &
                   form='unformatted') !Intel machines
!!!         end if
         WRITE(1) Vnew
         CLOSE(1)

         END DO !is

      END DO !inode 

      END DO !iname

      print*, 'DONE'

      DEALLOCATE( Vold )
      DEALLOCATE( Vnew )

      END PROGRAM rd_slices
