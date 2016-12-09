!===================================================================================================
!This subroutine creates the particle distribution list 
!===================================================================================================
IMPLICIT NONE
INTEGER			:: Counter,i,j
REAL*8      :: Dose, nu_m, M_d, dmin,dmax,sigma,miu,pi
REAL*8      :: d_end,d_up,d1,d2,eps,V_integ,V_par_i,Vf,Vp_tot,Vp_tot_Achieved

Write(*,*) 'Please enter the dosage (\mu g):'
read(*,*) Dose
Write(*,*) 'Please enter the maximum cut-off diameter, D_max (micron):'
read(*,*) dmax
Write(*,*) 'Please enter the average particle diameter, D_ave (micron):'
read(*,*) miu
Write(*,*) 'Please enter the minimum cut-off diameter, D_min (micron):'
read(*,*) dmin
Write(*,*) 'Please enter the standard deviation of the particle distribution, Sigma (micron):'
read(*,*) sigma

pi              = 4.0 * atan(1.0)
nu_m            = 268.0e6  	           ! (\mu m)^3 / (\mu mol)  		
M_d             = 206.285              ! (\mu g)   / (\mu mol)
Vp_tot_Achieved = 0.0
Vp_tot          = (Dose/M_d) * nu_m    ! (\mu m) ^3
d_up            = dmax
d1              = dmax
Counter         = 0
eps=1e-6

open(50,file='Particle_Sizes.dat')

Write(50,*) 'Total dosage                 ', Dose
write(50,*) 'Total Particle volume desired', Vp_tot 
write(50,*) 'Minimum cut-off diameter     ', dmin 
write(50,*) 'Average particle diameter    ', miu
Write(50,*) 'Maximum cut-off diameter     ', dmax

DO WHILE (d1.GE.dmin) 
   Counter = Counter + 1
   V_integ = 0.0
   V_par_i = 1.0
   d2      = d_up
   DO WHILE (V_par_i .GE. V_integ)
!      eps = eps + 1e-11  
      d1 =  d2 - eps
      V_integ = V_integ + (1.0/sqrt(2.*pi*sigma**2)) * exp(-(((d2+d1)/2.0)-miu)**2/(2.*sigma**2)) * (d2-d1) *Vp_tot
      V_par_i = (4.0/3.0)* PI* ((d1+d2)/4.0)**3.0
      d2 = d1
   END DO
   Vf = (V_par_i/Vp_tot)/(d_up-d1) 
   Vp_tot_Achieved= Vp_tot_Achieved + V_Par_i
   write(50,*) Counter, (d1+d_up)/2.0, V_integ, V_par_i, Vf 
   d_up= d1
END DO

write(*,*) 'Number of particle generated:', Counter
write(*,*) 'Total particle volume desired and achieved', Vp_tot, Vp_tot_Achieved
write(*,*) 'Achieved Dosage', (1.0- ((Vp_tot-Vp_tot_Achieved)/Vp_tot)) * Dose

stop

end
