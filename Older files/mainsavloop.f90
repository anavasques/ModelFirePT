!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                      PROGRAM TO MODEL "SAVANNAH"                     !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

! FIRE IS MODELED AS AN INSTANTANEOUS EVENT THAT LEAVES JUST A FRACTION
! OF THE ADULT TREES AND BRINGS GRASSES AND SEEDLINGS TO A THRESHOLD EPS 
! (UNLESS THEY ARE ALREADY LOWER THAN THE THRESHOLD!!)

! PROGRAM TO MODEL LOOP OVER PARAMETERS c3 and c2

program mainsav
use parasav

implicit none

double precision,dimension(3) ::  b,db,f,df,fout,bav, b0
real,dimension(:), allocatable::firev, c2vec!, firec 
real              ::  t, dummy, ran3
integer           ::  i, i3, i2, ii, dimfv, iday, idummy, maxTd, NN, iii
character (len=55)::  file12



!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! CODE STARTS

!open (21,file=file1,status='old')
!open (22,file='firedays')
!open (23,file=file4)

maxTd=maxT*365
dimfv=maxTd
NN=maxTd*1/h
ALLOCATE(firev(dimfv))
!ALLOCATE(firec(dimfv))
!call firedays(firev,firec,dimfv)
call firedays(firev,dimfv)
!ALLOCATE(rainv(NN))
!call raindays(rainv,NN)
!ALLOCATE(c2vec(NN))

!idummy=7896

!write(22,'(1f20.3)') firev
!write(23,'(1f20.3)') firec


hy=h/365              ! integration interval in years

!           stochastic amplitude            
!do i=1,NN
!      iday=i                                 ! COINCIDE SE h=1d
!      dummy=ran3(idummy)
!      c2vec(i)=cc2*(rainv(i))
!end do

!write(*,*) '<c2>=',sum(c2vec)/NN, '<c3>=',sum(c2vec)/NN/cc2*cc3

! loop over initial conditions
do iii=1,1
 if (iii.eq.1) then
  file12='bloopc3c2LowTree.dat'
  b0(1)=b01lt
  b0(2)=b02lt
 else
  file12='bloopc3c2LowGrass.dat'
  b0(1)=b01lg
  b0(2)=b02lg
 end if
 write(*,*) file12
 open (21,file=file12)
! loop over values of parameters
 do i3=0,100,5
 do i2=0,100,5

!     variables initialization
      b(1)=b0(1)
      b(2)=b0(2)
      b(3)=b03

      bav=0*b
      c3=i3*0.1
      c2=i2*0.1
      
      
!     temporal loop
      do i=1,(NN/2)
!            c3=i3*0.1*(1-firec(i)+1/2*firec(i))
            f=b
            call derivs(i,f,df)
            call rk4(i,f,df,fout)
            b=fout
            b(1)=b(1)*(1-firev(i)+fractree*firev(i))
!            b(2)=b(2)*(1-firev(i)+0.6*firev(i))
!            b(3)=b(3)*(1-firev(i)+0.6*firev(i))
            if (b(2)>eps) b(2)=b(2)*(1-firev(i))+eps*firev(i)
            if (b(3)>eps) b(3)=b(3)*(1-firev(i))+eps*firev(i)
 
!            if (b(2)>eps) then
!              b(2)=b(2)*(1-firev(i))+eps*firev(i)
!            else
!              b(2)=b(2)*(1-firev(i)+eps*firev(i))
!            end if
!              b(3)=b(3)*(1-firev(i))+0.9*firev(i)
!            else
!              b(3)=b(3)*(1-firev(i)+0.9*firev(i))
!            end if
            
            

                                  
!           if (mod(i,100).eq.0 .or. firev(iday)>0) then
!                  write (23,'(1i15,4f15.4)') iday,b, c3
!           end if
!            cav=cav+c3/1000
      end do 
      ii=0
      do i=(NN/2)+1,NN
!            c3=i3*0.1*(1-firec(i)+1/2*firec(i))            
            ii=ii+1
                        
            f=b
            call derivs(i,f,df)
            call rk4(i,f,df,fout)
            b=fout
            b(1)=b(1)*(1-firev(i)+fractree*firev(i))
!            b(2)=b(2)*(1-firev(i)+0.6*firev(i))
!            b(3)=b(3)*(1-firev(i)+0.6*firev(i))
            if (b(2)>eps) b(2)=b(2)*(1-firev(i))+eps*firev(i)
            if (b(3)>eps) b(3)=b(3)*(1-firev(i))+eps*firev(i)
 
!            if (b(2)>eps) then
!              b(2)=b(2)*(1-firev(i))+eps*firev(i)
!            else
!              b(2)=b(2)*(1-firev(i)+eps*firev(i))
!            end if
!              b(3)=b(3)*(1-firev(i))+0.9*firev(i)
!            else
!              b(3)=b(3)*(1-firev(i)+0.9*firev(i))
!            end if
                        
!            if (mod(i,100).eq.0 .or. firev(iday)>0) then
!                  write (23,'(1i15,4f15.4)') iday,b, c3
!            end if
            bav=bav+b/1000
      end do 
      bav=bav/ii*1000
      write (21,'(2f6.2,1i15,4f10.4)') c2, c3, i,bav
      write (* ,'(2f6.2,1i15,4f10.4)') c2, c3, i,bav
end do
end do
end do

end program mainsav


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

        subroutine rk4(tt,y,dydx,yout)
        use parasav
        
        implicit none
        
        real                          :: tt,h6
        double precision,dimension(3) :: dydx,y,yout,k1,k2,k3,k4,y1,y2,y3
        h6=hy/6.

        k1=dydx
        y1=y+0.5*k1*hy
        call derivs(tt,y1,k2)
        y2=y+0.5*k2*hy
        call derivs(tt,y2,k3)
        y3=y+k3*hy
        call derivs(tt,y3,k4)

        yout=y+h6*(k1+2*k2+2*k3+k4)

        return
        end


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

        subroutine derivs(tt,v,dv)
!       v vegetation cover, dv derivatives
        use parasav
        
        implicit none
        real               :: tt
        double precision,dimension(3)  :: v,dv
      
        dv(1)=g1*v(3)-m1*v(1)
        dv(2)=c2*v(2)*(1-v(1)-v(2))-m2*v(2)
        dv(3)=c3*v(1)*(1-v(1)-v(2)-v(3))-m3*v(3)-g1*v(3)-c2*v(2)*v(3)

        return
        end
        
        
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

        subroutine derivs2(tt,v,dv)
!       v vegetation cover, dv derivatives
        use parasav
        
        implicit none
        real               :: tt
        double precision,dimension(3)  :: v,dv
        
        if (m2*hy>2) then
            m2=2/hy
        end if
        
        if (m3*hy>2) then
            m3=2/hy
        end if
      
        dv(1)=g1*v(3)-m1*v(1)
        dv(2)=c2*v(2)*(1-v(1)-v(2))-m2*v(2)
        dv(3)=c3*v(1)*(1-v(1)-v(2)-v(3))-m3*v(3)-g1*v(3)-c2*v(2)*v(3)

        return
        end
        
        
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!      subroutine firedays(firevec,firecol,dimdum)
      subroutine firedays(firevec,dimdum)
      use parasav
      
      implicit none
      
      real                    :: dum, ran3
      integer                 :: i, ji, jji, idum,dimdum, numday
      real,dimension(dimdum)  :: firevec!, firecol
      
      firevec=0.
      idum=23459
      ji=1
!      jji=1

      do 
!          do
            dum=ran3(idum)
            numday=nint(-log(dum)*365*Tvect(ifire))
!            if (numday>ndf) exit
!          end do
            ji=ji+numday
            if (ji>maxT*365) exit
            firevec(ji:ji)=1
!            if (ji+ndf-1>maxT*365) then
!              firecol(ji:maxT*365)=1
!              exit
!            else  
!              firecol(ji:ji+ndf-1)=1
!            end if
      end do
      

      end subroutine
      
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

      subroutine raindays(rainvec,dimdum)
      use parasav
      
      implicit none
      
      real                    :: numday, dum, ran3
      integer                 :: i, ji, idum,dimdum
      real,dimension(dimdum)  :: rainvec
      
      rainvec=0.
      idum=23459
      ji=1

      do 
         do	
            dum=ran3(idum)
            numday=nint(-log(dum)*Tmoist)       ! time btw two rain events
            if (numday>0) exit
         end do  
            rainvec(ji+numday-1:ji+numday-1+ndf-1)=1   
            ji=ji+numday-1+ndf
            if (ji>maxT*365) exit
      end do
      

      end subroutine

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

      FUNCTION RAN3(IDUM)
      SAVE
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.E-9)
      DIMENSION MA(55)
      DATA IFF /0/
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        MJ=MSEED-IABS(IDUM)
        MJ=MOD(MJ,MBIG)
        MA(55)=MJ
        MK=1
        DO 11 I=1,54
          II=MOD(21*I,55)
          MA(II)=MK
          MK=MJ-MK
          IF(MK.LT.MZ)MK=MK+MBIG
          MJ=MA(II)
11      CONTINUE
        DO 13 K=1,4
          DO 12 I=1,55
            MA(I)=MA(I)-MA(1+MOD(I+30,55))
            IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
12        CONTINUE
13      CONTINUE
        INEXT=0
        INEXTP=31
        IDUM=1
      ENDIF
      INEXT=INEXT+1
      IF(INEXT.EQ.56)INEXT=1
      INEXTP=INEXTP+1
      IF(INEXTP.EQ.56)INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF(MJ.LT.MZ)MJ=MJ+MBIG
      MA(INEXT)=MJ
      RAN3=MJ*FAC
      RETURN
      END
