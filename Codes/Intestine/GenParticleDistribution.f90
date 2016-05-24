
!===================================================================================================
!This subroutine creates the particle distribution list 
!===================================================================================================

IMPLICIT NONE

real   			::   dmin,dmax,dcen,dgmin,dgmax
real   			::   sigma,miu,pi
real   			::   vd,vdmin,vdmax,vtot,deltd,ntot
real   			::   testnp,testv
integer			::   nptot,ngrp, i,j,k
real,    allocatable 	::  vfrac(:),nfrac(:),np(:)
integer, allocatable 	:: intnp(:)


!----- Read the input parameters -------------------------------------------------------------------
!write(*,*) "please enter the total number of particles"				!nptot=	500
!read(*,*) nptot
!
!write(*,*) "please enter the total number of bins" 				!ngrp= 	20
!read(*,*) ngrp 
!
!write(*,*) "please enter the minimum diameter"					!dmin= 	3.0
!read(*,*) dmin
!
!write(*,*) "please enter the maximum diameter"					!dmax= 	92.1
!read(*,*) dmax 
!
!write(*,*) "please enter the diameter at the center of the distribution"	!dcen= 52.0
!read(*,*) dcen
!
!write(*,*) "please enter the standard deviation"   				!sigma= 17.32
!read(*,*) sigma
! 
!write(*,*) "please enter the medium diameter"	                                !miu= 52.0
!read(*,*) miu
!
!---------------------------------------------------------------------------------------------------
allocate(vfrac(ngrp),nfrac(ngrp),np(ngrp))
allocate(intnp(ngrp))

nptot	= 500
ngrp	= 20
dmin	= 5
dmax	= 195
dcen	= 50
miu	= 50
sigma	= 30

pi	= 4.0 * atan(1.0)
deltd 	= (dmax-dmin)/(ngrp-1)
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
   nfrac(k) = vfrac(k) / (4.0/3.0*pi* ((dgmin+dgmax)/4.0)**3)
   ntot	    = ntot + nfrac(k)
end do

vfrac = vfrac / vtot
nfrac = nfrac / ntot
np    = nfrac * nptot

do k= 1, ngrp
   write(*,*) 'k,vfrac,nfrac,np',k,vfrac(k),nfrac(k),np(k)
end do

open(50,file='Par_Dist.dat')
do i= 1,ngrp
   dgmin= (dmin-deltd)/2.0 + (i-1)*deltd
   dgmax= (dmin+deltd)/2.0 + (i-1)*deltd
   write(50,*) (dgmin+dgmax)/2.,1./sqrt(2.*pi*sigma**2)*exp(-(dgmin/2.+dgmax/2.-miu)**2/(2.*sigma**2)),      &
               vfrac(i),np(i),int(np(i)+0.5)
end do

!----- Testing total number of particles and total volume
testnp=0.0
testv=0.0
do i= 1,ngrp
   testnp= testnp+np(i)
   testv = testv+vfrac(i)
end do
write(*,*) testnp,testv

stop

end
