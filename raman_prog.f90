PROGRAM le
IMPLICIT NONE
INTEGER::k,i,tot_struct = 19999,j1 = 2
REAL::eps_data(1000000,10),ac(1000000,10)
REAL::s = 0,sum_avg,sum_up,sum_down,x,sum_fourier_cos,sum_fourier_sin,dt,integration
REAL::h_k,r,w,k_in,w_in,temp,exp_part,wavenum,f,pi,n_in

pi = (22.0/7.0)

OPEN(UNIT=9,FILE='./tables/CO2_acorr_sep.csv')
OPEN(UNIT=16,FILE='Raman_11_model')

! PRINT*,"Give wavelength (in nm) of incident light:  "
! READ*,n_in 

n_in = 532

w_in = (45.613245464/n_in)                  !!! in a.u 

PRINT*,w_in

! PRINT*,"Give the temp(in K): "
! READ*,temp

temp = 300

DO k = 1,tot_struct
  READ(9,*)(ac(k,i),i=1,j1)
ENDDO

dt = ac(2,1) - ac(1,1)

PRINT*, dt

h_k = 3.15775*100000     !!!!!! h_bar/K value (in Sec.K )  (h = planck const, K = boltzmann)

DO k= -(tot_struct/2),(tot_struct/2)

  sum_fourier_cos = 0
  sum_fourier_sin = 0

  IF(k .ne. 0)THEN

       f = float(k)/(float(tot_struct)*dt)

 
       wavenum = (33.333333*f)   !!! wavenumber   in cm-1

       w = 1.5204417*0.0001*f    !!! in a.u.

     DO i= 1,tot_struct
       
       sum_fourier_cos = sum_fourier_cos + ac(i,2)*COS((2*pi*float(i)*float(k))/float(tot_struct))  
       sum_fourier_sin = sum_fourier_sin + ac(i,2)*SIN((2*pi*float(i)*float(k))/float(tot_struct))

     ENDDO
  ENDIF

       integration = SQRT(sum_fourier_cos**2 + sum_fourier_sin**2)  

   exp_part = exp(-(h_k*w)/temp)

   r = ((w_in - w)**4/(w*(1 - exp_part)))*integration

   WRITE(16,*)wavenum,r !-((w_in - w)**4/(w*(1 - exp_part)))
ENDDO 

END PROGRAM le
