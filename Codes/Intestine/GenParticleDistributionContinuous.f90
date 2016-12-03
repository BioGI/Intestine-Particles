
!===================================================================================================
!This subroutine creates the particle distribution list 
!===================================================================================================
IMPLICIT NONE
REAL*8      :: Dose, nu_m, M_d
REAL*8			:: dmin,dmax,dcen,dg,dd,dgmin,dgmax
REAL*8			:: sigma,miu,pi
REAL*8			:: vd,vdmin,vdmax,vtot,deltd,ntot
REAL*8			:: testnp,testv
REAL*8			:: vf,V_c,V_P,C_tot_Over_Cs, C_s, nu_P, C_tot, Vp_tot,Ctot_Over_Cs_test
REAL*8      :: d_end,d_up,d1,d2,eps,A,B1,B2,dB1,dB2,V_integ,V_par_i
INTEGER			:: Counter,nptot, ngrp, i,j
real,    allocatable 	:: vfrac(:),nfrac(:),np(:)
integer, allocatable 	:: intnp(:)

Dose  = 500.0             ! \mu g
nu_m  = 268.0e6  	         ! (\mu m)^3 / (\mu mol)  		
M_d   = 206.285            ! (\mu g)   / (\mu mol)
Vp_tot= (Dose/M_d) * nu_m  ! (\mu m) ^3

ngrp	= 20
dmin	= 5.0
dmax	= 195.0
miu   = 100.0
sigma	= 25.0

allocate(vfrac(ngrp),nfrac(ngrp),np(ngrp))

pi      = 4.0 * atan(1.0)
deltd   = (dmax-dmin)/(ngrp-1) 		
vtot    = 0.0
ntot    = 0.0
vfrac(:)= 0.0
nfrac(:)= 0.0

DO i= 1,ngrp
   dg= dmin + (i-1)*deltd
   !----- Integrating inside each bin (descritized to 100 sections)
   DO j= 1, 100
      dd = (dg-(deltd/2.0)) + ((j-1)+0.5)*(deltd/100)
      vf = (1.0/sqrt(2.*pi*sigma**2)) * exp(-(dd-miu)**2/(2.*sigma**2))
      vfrac(i) = vfrac(i) + (deltd/100) * vf
   END DO
   vtot= vtot + vfrac(i)
END DO

!----- Computing the number of particles in each bin -----------------------------------------------
open(49,file='Plot_file.dat')
do i= 1,ngrp
   dg       = dmin + (i-1)*deltd        	
   vfrac(i) = vfrac(i) * (Vp_tot/vtot)
   nfrac(i) = vfrac(i) / ( (4.0*pi/3.0)* (dg/2.0)**3)
   ntot     = ntot + nfrac(i)
   write(49,*) dg,  (vfrac(i)/deltd), ((int(nfrac(i)+0.5)*(4.0/3.0*pi*(dg/2.0)**3))/deltd), int(nfrac(i)+0.5) 
end do

V_P= 0.0

do i= 1,ngrp
   dg       = dmin + (i-1)*deltd
   nfrac(i) = vfrac(i) / ( (4.0*pi/3.0)* (dg/2.0)**3)
   V_P = V_P +  (int(nfrac(i)+0.5)) * (4.0/3.0*pi*(dg/2.0)**3) 
end do

!----- Creating the Par_Dist file with particle distribution data ----------------------------------
open(50,file='Par_Dist.dat')
do i= 1,ngrp
   dg  	= dmin + (i-1)*deltd        
   write(50,*) dg, 1./sqrt(2.*pi*sigma**2)*exp(-(dg-miu)**2/(2.*sigma**2)), vfrac(i), nfrac(i), int(nfrac(i)+0.5)  
end do

!----- Testing total number of particles and total volume ------------------------------------------
testnp=0.0
testv=0.0
do i= 1,ngrp
   TestNp= TestNp + int(nfrac(i)+0.5) 
   TestV = TestV  + vfrac(i)
end do

write(*,*) '----------------------------------------------------------------------------------------'
write(*,*) 'Desired V_P:  ', Vp_tot
write(*,*) 'Real V_P:     ', V_P
Write(*,*) 'Desired Dose: ', Vp_tot*M_d/nu_m!Dose
Write(*,*) 'Real Dose:    ', V_P*M_d/nu_m

!-----OutPut file ----------------------------------------------------------------------------------
open(51,file='Plot_file_Vexpected.dat')
write(51,2) Vp_tot, Vp_tot
write(51,2) TestNp,TestNp 
write(51,2) dmin, dmin
write(51,2) miu, miu
write(51,2) dmax,dmax
write(51,2) sigma,sigma
write(51,*) ngrp,ngrp

DO i= 1, 2200
   dg= dmin
   dd = (dg-(deltd/2.0)) + ((i-1)+0.5) * ((dmax-dmin)/2000)
   vf = (1.0/sqrt(2.*pi*sigma**2)) * exp(-(dd-miu)**2/(2.*sigma**2))
   write(51,2) dd,  vf
END DO
2 format (F25.8,F25.8)
!-------------

d_end= 170
d_up = d_end

DO i=1,100000
   eps = 0.0
   V_integ= 0.0
   V_par_i= 1.0
   Counter =0
   d2 = d_up
   DO WHILE (V_par_i .GE. V_integ)
      Counter =Counter+1
      eps = eps + 0.0000001  
      d1 =  d2 - eps
      IF (d1 .LE. 2.5) stop
      V_integ = V_integ + (1.0/sqrt(2.*pi*sigma**2)) * exp(-(((d2+d1)/2.0)-miu)**2/(2.*sigma**2)) * (d2-d1) *Vp_tot
!     A  =  1.0 / sqrt(2.0*PI*(sigma**2.0))
!     B1 = -( (d1-miu)**2.0 / (2.0*(sigma**2.0)) )
!     B2 = -( (d2-miu)**2.0 / (2.0*(sigma**2.0)) )
!     dB1= -( 2.0*(d1-miu)  / (2.0*(sigma**2.0)) )
!     dB2= -( 2.0*(d2-miu)  / (2.0*(sigma**2.0)) )
!     V_integ = A*(exp(B2)/dB2 - exp(B1)/dB1 ) * Vp_tot
      V_par_i = (4.0/3.0)* PI* ((d1+d2)/4.0)**3.0
!     write(*,*) 'd1,d2i,V_integ,V_par', d1,d2, V_integ,V_par_i
      d2 = d1
   END DO
   Vf=(V_par_i/Vp_tot)/(d_up-d1) 
   d_up= d1
   write(*,*) (d1+d_up)/2.0, V_integ,V_par_i,Vf
END DO

stop

end
