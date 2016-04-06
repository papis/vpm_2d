Subroutine definevort(RHS_pm,Xbound,Dpm,NN,NN_bl)
  Implicit None
 
  double precision, intent(in)   :: Xbound(6),Dpm(3)
  integer,intent(in)             :: NN(3),NN_bl(6)
  double precision, intent(inout):: RHS_pm(2,NN(1),NN(2),NN(3))
  double precision               :: X(4), XC(4), Y(4), YC(4), ksi1,ksi2,th1,th2,xi,yi,w1,w2,Vol,ANG,dens1,dens2
  integer                        :: i ,j ,k, NXpm1, NYpm1, NZpm1, NCELLSpm, npar,ndumc
  double precision               :: xc1,xc2,yc1,yc2,rmax,r,th,dth,dr,ths,rs,PI
  logical          :: I_EXIST
    xc1=0;xc2=0;
    yc1=-1.5;yc2=1.5
    PI = 4.*datan(1.d0)
    do j = NN_bl(2),NN_bl(5)
        do i = NN_bl(1),NN_bl(4)
            xi=Xbound(1) + (i-1)*Dpm(1)
            yi=Xbound(2) + (j-1)*Dpm(2)
            ksi1=sqrt((xi-xc1)**2+(yi-yc1)**2)
            ksi2=sqrt((xi-xc2)**2+(yi-yc2)**2)
            th1=atan2((yi-yc1),xi-xc1)
            th2=atan2((yi-yc2),xi-xc2)
            if (th1.lt.0.d0) th1=th1+2.d0*PI
            if (th2.lt.0.d0) th2=th2+2.d0*PI
            w1   = (2.d0 - ksi1**2) * exp (0.5d0*(1.d0-ksi1**2))
            w2   =-(2.d0 - ksi2**2 )* exp (0.5d0*(1.d0-ksi2**2))
            RHS_pm(1,i,j,1)=-(w1+w2)
         enddo
    enddo

 End Subroutine definevort

