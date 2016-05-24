
!===================================================================================================
!This subroutine creates the particle distribution list 
!===================================================================================================

IMPLICIT NONE

real   			::   dmin,dmax,dcen,dg,dgmin,dgmax
real   			::   sigma,miu,pi
real   			::   vd,vdmin,vdmax,vtot,deltd,ntot,Vp_tot
real   			::   testnp,testv
integer			::   nptot,ngrp, i,j,k
real,    allocatable 	::  vfrac(:),nfrac(:),np(:)
integer, allocatable 	:: intnp(:)


!----- Read the input parameters -------------------------------------------------------------------
!write(*,*) "please enter the total number of particles"				!nptot=	500
!read(*,*) nptot
!write(*,*) "please enter the total number of bins" 				!ngrp= 	20
!read(*,*) ngrp 
!write(*,*) "please enter the minimum diameter"					!dmin= 	3.0
!read(*,*) dmin
!write(*,*) "please enter the maximum diameter"					!dmax= 	92.1
!read(*,*) dmax 
!write(*,*) "please enter the diameter at the center of the distribution"	!dcen= 52.0
!read(*,*) dcen
!write(*,*) "please enter the standard deviation"   				!sigma= 17.32
!read(*,*) sigma
!write(*,*) "please enter the medium diameter"	                                !miu= 52.0
!read(*,*) miu
!---------------------------------------------------------------------------------------------------

nptot	= 500
Vp_tot	= 833.5e5 ! 1667.0e5 
ngrp	= 20
dmin	= 10.0
dmax	= 200.0
dcen	= 100.0
miu	= 100.0
sigma	= 30.0

allocate(vfrac(ngrp),nfrac(ngrp),np(ngrp))
allocate(intnp(ngrp))

pi	= 4.0 * atan(1.0)
deltd 	= (dmax-dmin)/(ngrp-1) 		! (dmax-dmin)/(ngrp-1)
vtot	= 0.0
ntot	= 0.0
vfrac(:)= 0.0
nfrac(:)= 0.0

do k= 1,ngrp
   dgmin    = dmin- deltd/2.0 + (k-1)*deltd
   dgmax    = dmin+ deltd/2.0 + (k-1)*deltd
   vdmin    = 1.0/sqrt(2.*pi*sigma**2) * exp(-(dgmin-miu)**2/(2.*sigma**2))                    
   vdmax    = 1.0/sqrt(2.*pi*sigma**2) * exp(-(dgmax-miu)**2/(2.*sigma**2))                    
   vfrac(k) = (vdmin+vdmax)*deltd/2.0
   vtot	    = vtot + vfrac(k)
!  nfrac(k) = vfrac(k) / (4.0/3.0*pi* ((dgmin+dgmax)/4.0)**3)
!  ntot	    = ntot + nfrac(k)
end do

!vfrac = vfrac / vtot
!nfrac = nfrac / ntot
!np    = nfrac * nptot

do k= 1,ngrp
   dg    	= dmin + (k-1)*deltd        	
   vfrac(k) 	= vfrac(k) * (Vp_tot/vtot)
   nfrac(k) 	= vfrac(k) / (4.0/3.0*pi* ((dg)/2.0)**3)
   ntot     	= ntot + nfrac(k)
   write(*,*) 'k,vfrac,nfrac',k,vfrac(k),nfrac(k)
end do


open(50,file='Par_Dist.dat')
do i= 1,ngrp
   dg    	= dmin + (i-1)*deltd        	
   write(50,*) dg, 1./sqrt(2.*pi*sigma**2)*exp(-(dg-miu)**2/(2.*sigma**2)),      &
               vfrac(i), nfrac(i), int(nfrac(i)+0.5)  
end do

!----- Testing total number of particles and total volume
testnp=0.0
testv=0.0
do i= 1,ngrp
   TestNp= TestNp + int(nfrac(i)+0.5) 
   
   TestV = TestV  + vfrac(i)
end do
write(*,*) testnp,testv

stop

end
