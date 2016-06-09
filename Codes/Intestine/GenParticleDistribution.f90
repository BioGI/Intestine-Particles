
!===================================================================================================
!This subroutine creates the particle distribution list 
!===================================================================================================
IMPLICIT NONE

REAL*8			:: dmin,dmax,dcen,dg,dd,dgmin,dgmax
REAL*8			:: sigma,miu,pi
REAL*8			:: vd,vdmin,vdmax,vtot,deltd,ntot
REAL*8			:: testnp,testv
REAL*8			:: vf,V_c,V_P,C_tot_Over_Cs, C_s, nu_P, C_tot, Vp_tot,Ctot_Over_Cs_test
INTEGER			:: nptot, ngrp, i,j
real,    allocatable 	:: vfrac(:),nfrac(:),np(:)
integer, allocatable 	:: intnp(:)

V_c		= 9.424778e12 		! (\mu m) ^3
C_tot_Over_Cs	= 0.2000000000000	! dimensionless
C_s		= 3.3e-19              	! mole / (\mu m ^3)
nu_P		= 268.0e12 		! (\mu m ^3) / mole  		
C_tot		= C_s * C_tot_Over_Cs   ! mole / (\mu m ^3)
Vp_tot  	= C_tot * nu_P * V_c    ! (\mu m) ^3

ngrp	= 20
dmin	= 5.0
dmax	= 195.0
miu	= 100.0
sigma	= 25.0

allocate(vfrac(ngrp),nfrac(ngrp),np(ngrp))

pi	= 4.0 * atan(1.0)
deltd 	= (dmax-dmin)/(ngrp-1) 		
vtot	= 0.0
ntot	= 0.0
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

Ctot_Over_Cs_test = V_P/(V_c * nu_P * C_s)  

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
write(*,*) 'Desired V_P    :',Vp_tot
write(*,*) 'Real V_P       :',V_P
write(*,*) 'Desired Ctot/Cs:', C_tot/C_s
write(*,*) 'Real Ctot/Cs   :',Ctot_Over_Cs_test


!-----OutPut file ----------------------------------------------------------------------------------
open(51,file='Plot_file_Vexpected.dat')
write(51,2) C_tot_Over_Cs, C_tot_Over_Cs
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

stop

end
