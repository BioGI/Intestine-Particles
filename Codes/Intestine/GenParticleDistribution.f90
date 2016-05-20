       IMPLICIT NONE

       real   ::   dmin,dmax,dcen
       real   ::   dgmin,dgmax
       real   ::   sigma,miu
       real   ::   vd,vdmin,vdmax,vtot,deltd,ntot
       real   ::   pi
       real   ::   testnp,testv

       integer   ::   nptot,ngrp
       integer   ::   i,j,k

       real, allocatable  ::  vfrac(:),nfrac(:),np(:)
       integer, allocatable  :: intnp(:)

       pi=4.*atan(1.0)

       dmin = 3.0
       dmax = 92.1
       dcen = 52.0

       sigma = 17.32
       miu=52.0

       write(*,*) "please enter the total number of particles"
       read(*,*) nptot
       write(*,*) "please enter the total number of bins"
       read(*,*) ngrp 

       allocate(vfrac(ngrp),nfrac(ngrp),np(ngrp))
       allocate(intnp(ngrp))

       vtot=0.0
       ntot=0.0
       vfrac(:)=0.0
       nfrac(:)=0.0

       deltd=(dmax-dmin)/(ngrp-1)
       do k=1,ngrp
         dgmin=dmin-deltd/2.+(k-1)*deltd
         dgmax=dmin+deltd/2.+(k-1)*deltd
         vdmin=1./sqrt(2.*pi*sigma**2)*exp(-(dgmin-miu)**2/(2.*sigma**2))                    
         vdmax=1./sqrt(2.*pi*sigma**2)*exp(-(dgmax-miu)**2/(2.*sigma**2))                    
         vfrac(k)=(vdmin+vdmax)*deltd/2.
         vtot=vtot+vfrac(k)
         nfrac(k)=vfrac(k)/(4./3.*pi*((dgmin+dgmax)/4.)**3)
         ntot=ntot+nfrac(k)
       end do

       vfrac=vfrac/vtot
       nfrac=nfrac/ntot
       np=nfrac*nptot

       open(50,file='Par_Dist.dat')
       do i=1,ngrp
         dgmin=dmin-deltd/2.+(i-1)*deltd
         dgmax=dmin+deltd/2.+(i-1)*deltd
         write(50,*) (dgmin+dgmax)/2.,1./sqrt(2.*pi*sigma**2)*exp(-(dgmin/2.+dgmax/2.-miu)**2/(2.*sigma**2)),      &
                  vfrac(i),np(i),int(np(i)+0.5)
       end do

       testnp=0.0
       testv=0.0
       do i=1,ngrp
         testnp=testnp+np(i)
         testv=testv+vfrac(i)
       end do
       write(*,*) testnp,testv

       stop
       end
