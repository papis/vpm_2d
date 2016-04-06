
Program test_pm
use vpm_lib
use test_mod
use pmgrid
use MPI
Implicit None
double precision :: Vref,NI_in,DT_in,RMETM,OMET,OG,FACDEF,T,XMIN,XMAX,UINF(3)
double precision,allocatable ::velsavex(:,:,:)
double precision,allocatable ::XPDUM(:,:),QPDUM(:,:)
integer          :: Noutput, NDimoutput,NVR_turb,NVR_all
integer :: my_rank,np,ierr,i,neq,j

call MPI_INIT(ierr)
call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
read (*,*) DT_in
neq=1
if (my_rank.eq.0) then 
    open(1,file='particles.dat',form='formatted')
    read(1,*) NVR_ext
    allocate(XPR(3,NVR_ext),QPR(neq+1,NVR_ext))
    QPR=0;XPR=0
    write(*,*) 'NVR=',NVR_ext
    do i=1,NVR_ext
       read(1,*) XPR(1,i),XPR(2,i),XPR(3,i)
    enddo
close(1)
endif

mrem=1
UINF=0
 !-Iwhattodo
 call vpm(XPR,QPR,UPR,GPR,NVR_ext,neq,0,RHS_pm_in,velx,vely,velz,0,NI_in,NVR_ext)

 if (my_rank.eq.0) st=MPI_WTIME()
 call remesh_particles_2d(0)
 if (my_rank.eq.0)then
     et= MPI_WTIME()
     write(*,*) 'remeshing',int((et-st)/60),'m',mod(et-st,60.d0),'s'
 endif
T=0
do i=1,6000
!get velocities and deformations
 T = DT_in 
if (my_rank.eq.0) allocate(UPR(3,NVR_ext),GPR(3,NVR_ext))
call vpm(XPR,QPR,UPR,GPR,NVR_ext,neq,2,RHS_pm_in,velx,vely,velz,i,NI_in,NVR_ext)
if (my_rank.eq.0) then 
     st=MPI_WTIME()
     do j= 1,NVR_ext
         XPR(1:3,j) = XPR(1:3,j)  + (UPR(1:3,j)+UINF(1:3)) * DT_in
     enddo
     !!$omp enddo
     et= MPI_WTIME()
     write(*,*) 'Convection',int((et-st)/60),'m',mod(et-st,60.d0),'s'
     deallocate (UPR,GPR)
      
endif
 call vpm(XPR,QPR,UPR,GPR,NVR_ext,neq,0,RHS_pm_in,velx,vely,velz,i,NI_in,NVR_ext)
 if (my_rank.eq.0) st=MPI_WTIME()
 if (mod(i,1).eq.0) call remesh_particles_2d(1)
 if (my_rank.eq.0)then
     et= MPI_WTIME()
     write(*,*) 'remeshing',int((et-st)/60),'m',mod(et-st,60.d0),'s'
 endif
!call vpm(XPR,QPR,UPR,GPR,NVR_ext,neq,5,RHS_pm_in,velx,vely,velz,i,NI_in,NVR_ext)
!if (my_rank.eq.0) then 
!   write(*,*) maxval(GPR(:,:))
!    do j= 1,NVR_ext
!        QPR(1:3,j) = QPR(1:3,j)  -GPR(1:3,j) * DT_in
!    enddo
!endif
!get velocities and deformation
enddo

!if (my_rank.eq.0) then 
!allocate (QPDUM(7,NVR_ext),XPDUM(3,NVr_ext))
!XPDUM =XPR
!QPDUM(1:3,:)=QPR(1:3,:)
!QPDUM(4:6,:)=QPR(1:3,:)
!QPDUM(7,:)  =QPR(4,:/)
!endif
! call vpm(XPDUM,QPDUM,UPR,GPR,NVR_ext,6,1,RHS_pm_in,velx,vely,velz,i,NI_in,NVR_ext)
!  if (mod(i,1).eq.0) call remesh_particles_3d(1)
 call MPI_FINALIZE(ierr)
End Program test_pm

 Subroutine writepar(NTIME,XPR,NVR)
    Implicit None
     integer,intent(in) :: NTIME,NVR
     double precision,intent(in):: XPR(3,NVR)
     integer ::i
    character*80 :: filout1
    write(filout1,'(i5.5,a)') NTIME,'vr.dat'
    open(10,file=filout1)
    WRITE(10,*)'VARIABLES = "X" "Y" "Z" '
    do  i=1, NVR
        write(10,'(3(F20.10,1x))') XPR(1,i), XPR(2,i),XPR(3,i)
    enddo
    call system('~/bin/preplot '&
        //filout1//' >/dev/null')
    call system('rm '//filout1)
    close(10)

 End Subroutine writepar
