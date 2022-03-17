MODULE nrtype
  INTEGER, PARAMETER :: ITT=SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: RT=SELECTED_REAL_KIND(P=15)
  INTEGER, PARAMETER :: RT2=SELECTED_REAL_KIND(P=31)
  INTEGER, PARAMETER :: RTC=SELECTED_REAL_KIND(P=15)
  INTEGER, PARAMETER :: CT=RT
  REAL(rt),PARAMETER :: PI=3.14159265359
  REAL(rt),PARAMETER :: Euler=2.71828182846
  REAL(rt),PARAMETER :: MACH_INFINITY=HUGE(0.0_rt)
END MODULE nrtype

!======================IMPORTANT======================================!
!this module calls the Cosmological model wich we are going to
!use, all the cosmologies are in the C_models.f90 file!!!
!====================================================================!

MODULE Global_cosmology
  USE nrtype
  !USE two_form_cosmology
  !Use one_form_cosmology
  !USE non_abelian_cosmology
  !Use EDS_cosmology
   USE LCDM_cosmology


  !========Global-variables=======!
  REAL(rt) :: G_m0_density, G_r0_density  ! *intial densities for all
  REAL(rt) :: G_k0_density, G_de0_density ! the models

  !========SPLINE_NUMBER GLOBAL====!
  INTEGER :: G_spline_numbe
  REAL(rt), ALLOCATABLE :: Vec_EOS(:), Vec_time(:), Vec_EOSpp(:)
  REAL(rt), ALLOCATABLE :: Vec_EOS_int(:),Vec_EOS_int_pp(:)

 CONTAINS

  SUBROUTINE Model_initialization

  !================spline-global variables================!

   CALL present_data(G_m0_density,G_r0_density,G_de0_density,G_k0_density)
   WRITE(*,*) "Densities today...(Intial data)"
   WRITE(*,*) "Matter:",G_m0_density
   WRITE(*,*) "Radiation:",G_r0_density
   WRITE(*,*) "Dark Energy:",G_de0_density
   WRITE(*,*) "Total:",G_m0_density+&
              g_r0_density+G_de0_density

   IF(Spline_model==1) THEN

   WRITE(*,*)"              "
   WRITE(*,*)"==========SPLINE INITIALIZATIONS================"
    CALL timestamp()

    !========Assign the number of splines=======!
    G_spline_numbe = Number_of_splines
    ALLOCATE(Vec_EOS(Number_of_splines),Vec_time(Number_of_splines))
    ALLOCATE(Vec_EOSpp(Number_of_splines),Vec_EOS_int(Number_of_splines))
    ALLOCATE(Vec_EOS_int_pp(Number_of_splines))

    !========EOS spline=================!
    CALL EOS_Spline_initialization(Vec_time,Vec_EOS,Vec_EOSpp)

    !=======int_EOS_spline_initialization====!
    CALL int_EOS_spline_initialization(Vec_time,Vec_EOS_int,Vec_EOS_int_pp)
    !===============================================!

   CALL timestamp()

   !======allow radiation and curvature.....
   G_k0_density=0.0_rt !no curvature
   G_r0_density=0.0_rt !no radiation,
   
   !IMPORTANT: this density of radiation contribution can be also be removed 
   !from the S.C equations itself
   
   WRITE(*,*)"================================================"
   WRITE(*,*)"              "

 ELSE IF(spline_model==0)THEN

    WRITE(*,*)"============================="
    WRITE(*,*)"Model without splines!!!"
    WRITE(*,*)"with present data..."
    G_k0_density=0.0_rt !no curvature
    G_r0_density=0.0_rt !no radiation
    WRITE(*,*) G_m0_density, G_r0_density,G_k0_density, G_de0_density
    WRITE(*,*)"============================="

   END IF

 END SUBROUTINE Model_initialization

END MODULE Global_cosmology

!==============================================
!Here, will be all the functions that are auxiliar
!but depend on th equation of state used and things like that
!================================================

MODULE Cosmology_func
  USE nrtype
  USE Global_cosmology ! the global model!!!!
  IMPLICIT NONE
  PUBLIC

  !=========Tags===========!
  !Omega_M0----> matter
  !Omega_R0----> radiation
  !Omega_K0----> curvature
  !Omega_DE----> DE
  !Omega_CASS--->Cassimir (vaccum)

CONTAINS

! Exp_func := computes the expansion function, and its derivative, notice
! that like the evaluation of g_func is expensive, turns better to call it the
! minima amount of times

  SUBROUTINE Exp_func(a,E_a,dEda)

    !==========input-values===================!
    REAL(rt) :: a, E_a,dEda
    REAL(rt) :: Omega_M0,Omega_K0,Omega_DE0,Omega_R0

    !==========subroutine-variables===========!
    REAL(rt) :: ga,dgda
    REAL(rt) :: Eos,wp,wpp

    !========Assignaments of the parameters======!
    !====According to the corresponding model:

    Omega_M0=G_m0_density
    Omega_R0=G_r0_density
    Omega_K0=G_k0_density
    Omega_DE0=G_de0_density

    IF(Spline_model==1) THEN

    !==============================================================================
    !Waning: in this definition we neglect any radiation and curvature contributions
    !for cosmologies with radiation and curvature uncomment bellow.
    !===============================================================================

    CALL g_function(a,ga) !dark energy component
    E_a=DSqrt((Omega_M0/(a**3))+(Omega_DE0*ga))

    !uncomment this for cosmologies with radiation and curvature.
    !E_a=DSqrt((Omega_M0/(a**3))+(Omega_K0/(a**2))+(Omega_R0/(a**4))+(Omega_DE0*ga))

    CALL spline_cubic_val(G_spline_numbe,Vec_time,Vec_EOS,Vec_EOSpp,&
                          Dlog(a),Eos,wp,wpp)

    dgda=-3.0_rt*ga*(EOS+1.0_rt)/a ! derivative of the dark-energy component

    dEda=(1.0_rt/(2*E_a))*((-3.0_rt*Omega_M0/(a**4))+(Omega_DE0*dgda))

    !uncomment this for cosmologies with radiation and curvature.    
    !dEda=(1.0_rt/(2*E_a))*((-3.0_rt*Omega_M0/(a**4))+(-2.0_rt*Omega_K0/(a**3)) &
    !+(-4.0*Omega_R0/(a**5))+(Omega_DE0*dgda))

    ELSE IF(Spline_model==0) THEN

       CALL E_model_function(a,E_a,dEda)

    END IF

  END SUBROUTINE Exp_func

  !SUBROUTINE g_function: computes the dark energy density contribution of the model
  !notice that if the EOS is complicated we must call the splines
  !funtion for avoid the multiple integrations.....

  SUBROUTINE g_function(a,ga)

    REAL(rt) :: a,ga
    REAL(rt) :: y,yp,ypp,t

    IF (Spline_model==1) THEN

    t=Dlog(a)
    CALL spline_cubic_val(G_spline_numbe,Vec_time,Vec_EOS_int,&
                          Vec_EOS_int_pp,t,y,yp,ypp)

    ga=Dexp(-3.0_rt*(t+ga))

    END IF

  END SUBROUTINE g_function


END MODULE Cosmology_func

!======================================================
!In this module are a number of different implementations of the
!S.C model, notice that different implementations bring some
!advantages
!======================================================

MODULE Spherical_collapse
  USE nrtype
  USE Cosmology_func
  !use fgsl
  !use, intrinsic :: iso_c_binding
  
  IMPLICIT NONE

  !====Global-Quantities=======================!
  REAL(rt) :: slope_G, lna0_G, Lnac_G,x_G,Target_global
  Real(rt) :: delta_ini ! only for the radious implementation
  !=================================================

  !===parameters like the numerical infinity and things like that
  REAL(rt), PARAMETER :: Eps_infinity=1.0e10_rt ! numerical inifinity
  REAL(rt), PARAMETER :: a_ini=1.0e-5_rt !initial scale factor
  REAL(rt), PARAMETER :: Eps_Zero=1.0_rt/Eps_infinity
  REAL(rt) :: abserr=1.0e-7_rt!DSqrt(EPSILON(abserr)) !acepted errors for the propagation
  REAL(rt) :: relerr=1.0e-7_rt!Dsqrt(EPSILON(relerr)) !acepted errors for th propagation
  !===========================

CONTAINS

  !====================Log_implementations======================!
  !In tis equations the time is set up as the number of e-foldings
  !i.e : t=log(a)!!!
  !=============================================================!

  SUBROUTINE Snd_NL_log(lna,x,dx_out)
    REAL(rt), DIMENSION(2) :: x,dx_out
    REAL(rt) :: E_lna, dE_lnada
    REAL(rt) :: lna,a

    E_lna=0.0_rt
    dE_lnada=0.0_rt
    dx_out=0.0_rt
    a=DExp(lna)

    CALL Exp_func(a,E_lna,dE_lnada)

    dx_out(1)=x(2)
    dx_out(2)=-((2.0_rt+(dE_lnada*a/E_lna))*x(2)) &
              +((4.0_rt/3)*((x(2)**2.0_rt)/(1.0_rt+x(1)))) &
              +((3.0_rt/2)*(G_m0_density/((a**3.0_rt)*(E_lna**2.0_rt)))*x(1)*(x(1)+1.0_rt))

     END SUBROUTINE Snd_NL_log


   SUBROUTINE Snd_L_log(lna,x,dx_out)
     REAL(rt), DIMENSION(2) :: x,dx_out
     REAL(rt) :: E_lna, dE_lnada
     REAL(rt) :: lna,a

     E_lna=0.0_rt
     dE_lnada=0.0_rt
     dx_out=0.0_rt
     a=DExp(lna)

     CALL Exp_func(a,E_lna,dE_lnada)

     dx_out(1)=x(2)
     dx_out(2)=-((2.0_rt+(dE_lnada*a/E_lna))*x(2)) &
               +((3.0_rt/2)*(G_m0_density/((a**3.0_rt)*(E_lna**2.0_rt)))*x(1))

    END SUBROUTINE Snd_L_log

  !===============Inverse_log_implementations==================!
  !the previous subroutiones evolve the overdensity
  !here we will evolve the inverse of the overdensity f=1/d
  !============================================================

  SUBROUTINE F_Nl_log(lna,x,dx_out)
    REAL(rt), DIMENSION(2) :: x,dx_out
    REAL(rt) :: E_lna, dE_lnada
    REAL(rt) :: lna,a
   
    E_lna=0.0_rt
    dE_lnada=0.0_rt
    dx_out=0.0_rt
    a=DExp(lna)

    CALL Exp_func(a,E_lna,dE_lnada)

    dx_out(1)=x(2)
    dx_out(2)=-((2.0_rt+(dE_lnada*a/E_lna))*x(2)) &
              +((2.0_rt-(4.0_rt/3*(1+x(1))))*(x(2)**2.0_rt)/x(1)) &
              -((3.0_rt/2)*(G_m0_density/((a**3)*(E_lna**2)))*(x(1)+1))

  END SUBROUTINE F_Nl_log

  
  !=====================================================================!
  !================Inverse-log/divergence implementations===============!
  !=====================================================================!

  Subroutine Ftheta(lna,x,dx_out)
  
   Real(rt), DIMENSION(2) :: x, dx_out
    REAL(rt) :: E_lna, dE_lnada
    REAL(rt) :: lna,a
      
    E_lna=0.0_rt
    dE_lnada=0.0_rt
    dx_out=0.0_rt
    a=DExp(lna)
    
    CALL Exp_func(a,E_lna,dE_lnada)
   
    dx_out(1)=x(1)*(1.0_rt+x(1))*x(2)
    dx_out(2)=-(2.0_rt+(dE_lnada*a/E_lna))*x(2) - ((x(2)**2.0_rt)/3.0_rt) &  
    -((3/2.0_rt)*(G_m0_density/((a**3.0_rt)*(E_lna**2.0_rt)))*(1.0_rt/x(1)))
   
 END SUBROUTINE Ftheta  
 
  
 Subroutine Ftheta_linear(lna,x,dx_out)
  
   Real(rt), DIMENSION(2) :: x, dx_out
    REAL(rt) :: E_lna, dE_lnada
    REAL(rt) :: lna,a
      
    E_lna=0.0_rt
    dE_lnada=0.0_rt
    dx_out=0.0_rt
    a=DExp(lna)
    
    CALL Exp_func(a,E_lna,dE_lnada)
   
    dx_out(1)=(x(1)**2.0_rt)*x(2)
    dx_out(2)=-(2.0_rt+(dE_lnada*a/E_lna))*x(2) &
    -((3/2.0_rt)*(G_m0_density/((a**3)*(E_lna**2.0_rt)))*(1.0_rt/x(1)))
   
 END SUBROUTINE Ftheta_linear
 

 !===================================================================
 !Radious implementations
 !=================================================================

 subroutine new_radious_evol(lna,z,dz)
   real(rt), DIMENSION(2) :: z, dz 
   real(rt) :: E, dE 
   real(rt) :: lna, a

   E=0
   dE=0
   a=Dexp(lna)
   
   call Exp_func(a,E,dE)
   
   delta=(delta_ini+1.0_rt)/((((a_ini/a)*z(2))+1)**3) - 1.0_rt

   dz(1)=(-a*dE/E)z(1)+((1.0_rt+(-a*dE/E)-(1.0_rt/2)*Omega_M0*delta/((a**3)*(E**2)))*z(2)) &
   -(1.0_rt/2)*Omega_M0*delta/((a**2)*(E**2)*a_ini)
   dz(2)=z(1)

 End subroutine
 
 !==========================================================
 !==========================================================
 !==========================================================
 
 !=========Power-law=========================================!
 subroutine power_law(a_0,n1,n2)
 
 real(rt) :: a_0 ! d~a_0^n
 real(rt) :: n1,n2 ! law exponents
 REAL(rt) :: E_lna, dE_lnada
 
 real(rt) :: discriminant
 
 
    E_lna=0.0_rt
    dE_lnada=0.0_rt
    CALL Exp_func(a_0,E_lna,dE_lnada)
    
discriminant=Dsqrt((2+(a_0*dE_lnada/E_lna))**2 +(4*(3.0_rt/2)*(G_m0_density/((a_0**3)*(E_lna**2)))))
 
 n1=(-(2+(a_0*dE_lnada/E_lna))+discriminant)/2.0_rt
 n2=(-(2+(a_0*dE_lnada/E_lna))-discriminant)/2.0_rt
  
 end subroutine
 
 
 subroutine SC_Observables(n_partitions,z_min,z_max)
 
 integer :: mu
 
 Integer :: n_partitions,j_step
 real(rt) :: z_min,z_max,z
 real(rt) :: t,t_out,t_tr
 
 real(rt) :: fl_collapse,fnl_tr,fl_tr
 real(rt) :: y_vir,vir_collapse
 
  fl_collapse=0.0_rt
  fl_tr=0.0_rt
  fnl_tr=0.0_rt
  y_vir=0.0_rt
  vir_collapse=0.0_rt
  
   
 CALL timestamp()  
  OPEN(newunit=mu, file="SP_Observables.dat", status="replace")
 
do j_step=1,n_Partitions

   z=z_min+j_step*((z_max-z_min)/n_Partitions)
   
   write(*,*)"======",j_step,"========="
   call Collapse_time(z,t_out,t_tr)
   write(*,*)"========================="
   write(*,*) " "
   if(t_out.eq.0.0_rt) exit   
   
   call propagate(z,t_out,fl_collapse,0) !linear
   call propagate(z,t_tr,fnl_tr,1) !non-linear
   call virial_radious(fnl_tr,t_out,t_tr,y_vir)
   call virial_quantities(fnl_tr,t_out,t_tr,y_vir,vir_collapse)
   
   write(mu,*) z,t_out,t_tr,(1.0_rt/fl_collapse),(1.0_rt/fnl_tr)+1.0_rt,y_vir,vir_collapse  

 end do
  
  close(mu)
  
CALL timestamp()  
  
  end subroutine SC_Observables  
  
  
 !subroutine Target_conditions
  
 ! real(rt) :: result, abserr
 ! integer(fgsl_int) :: status
 ! type(fgsl_function) :: pwr
!  type(c_ptr) :: param_unused
 ! real(rt) initial_inverse,final_inverse
  
 ! Target_global=0.0_rt
 ! initial_inverse=100.0_rt
 ! final_inverse=0.0_rt
  !param_unused=0
  
 ! pwr = fgsl_function_init(f, c_null_ptr)
  
 ! status = fgsl_deriv_central (pwr, initial_inverse, 1.E-8_fgsl_double,result, abserr)
 !  write(*,*) "initial f:",initial_inverse
 ! final_inverse=f(initial_inverse,param_unused)
 ! write(*,*) "collapse ini_f: ",final_inverse
 ! final_inverse=initial_inverse-(final_inverse/result)
!   write(*,*) "final invese: ", final_inverse
 !  write(*,*) "collapse final invese: ", f(final_inverse,param_unused) 
  
 ! write(*,*) "the results"
 ! write(*,*) 
  
 ! end subroutine Target_conditions
  
  
  subroutine Collapse_time(f_in,lna_out,lna_tr)
  
  real(rt) :: f_in !initial conditions
  real(rt) :: lna_out !collapse point
  real(rt) :: lna_tr  !turn around point 
  real(rt) :: yr,yr_aux ! inverse of radius
  
  real(rt) :: n1,n2 !power-like exponents
  real(rt) :: t ! running time
  
  INTEGER :: flag
  REAL(rt),ALLOCATABLE :: z(:),dz(:)
  
  ALLOCATE(z(2),dz(2))
  
  call power_law(a_ini,n1,n2)
  
     n1=max(n1,n2)
     z(1)=f_in
     z(2)=-n1*z(1)/(z(1)*(1.0_rt+z(1)))
     t=Dlog(a_ini)
          
     flag=-1 !single step mode (-1)
     yr=MACH_INFINITY
     
     call Ftheta(t,z,dz)
        
     call timestamp()
     
     do while(flag.lt.2)
     
       CALL r8_rkf45(Ftheta,2,z,dz,t,0.0_rt,relerr,abserr,flag)
       if(flag.eq.4) flag=-2
       
       yr_aux=DLog(((1.0_rt/z(1))+1.0_rt)/(Dexp(t)**3.0_rt))
       
       if(yr_aux.lt.yr) then 
       yr=yr_aux
       lna_tr=t
       end if 
       
       if(z(1).lt.1.0e-8_rt) exit
       
     end do
     
     call timestamp()
  
     lna_out=t
     write(*,*) "the flag", flag
     
DEALLOCATE(z,dz)
      
  
  end subroutine collapse_time
  
  
  !============================================================
  !VIRIAL_RADIOUS : computes the radious of the sphere at the 
  !virialization scale factor normilized at the turn around
  !============================================================
  
  subroutine virial_radious(f_tr,lna_collapse,lna_tr,y_vir)
  
  real(rt) :: f_tr ! f non-linear at turn around, its initialized 
  !initial condition with lna_c=lna_collapse
  real(rt) :: lna_collapse,lna_tr ! collase and turn_around logs
  real(rt) :: y_vir ! y_virialization
  
  real(rt) :: delta_tr,nu_tr,nu_v
  
  real(rt) :: m_density,r_density,de_density
  real(rt) :: w,dlnw
  
  delta_tr=(1.0_rt/f_tr)+1.0_rt ! overdensity at tr
  y_vir=0.0_rt
  
  !CALL evol_var(lna_tr,m_density,r_density,de_density,w,dlnw) !call from the diferential equations
  
  m_density=G_m0_density/(Dexp(lna_tr)**3) !for LCDM_cosmology
  de_density=1.0_rt-m_density
  
  
  nu_tr=(2.0_rt/delta_tr)*(de_density/m_density)
  
  !call evol_var(lna_collapse,m_density,r_density,de_density,w,dlnw)
  
  m_density=G_m0_density/(Dexp(lna_collapse)**3) !for LCDM_cosmology
  de_density=1.0_rt-m_density
  
  
  nu_v=(2.0_rt/delta_tr)*(de_density/m_density)*((Dexp(lna_tr)/Dexp(lna_collapse))**3)
  
  y_vir=(1.0_rt-nu_v/2.0_rt)/(2.0_rt+nu_tr-(3/2.0_rt)*nu_v) 
  
    end subroutine virial_radious
    
 !======================================================================  
 !subroutine: virial_quantities......
 !returns the virial overdensity at the collapse_time
 !======================================================================
 
 subroutine virial_quantities(f_tr,lna_collapse,lna_tr,y_vir,vir_collapse)
 
 real(rt) :: f_tr,y_vir 
 real(rt) :: lna_collapse,lna_tr 
 real(rt) :: vir_collapse
 
 real(rt) :: a_collapse,a_tr
 real(rt) :: E_lna,dE_lna
 real(rt) :: delta_tr
 real(rt) :: omega_m_collapse
 
 a_collapse=Dexp(lna_collapse)
 a_tr=Dexp(lna_tr)
 delta_tr=(1.0_rt/f_tr)+1.0_rt
 vir_collapse=0.0_rt
 
 call Exp_func(a_collapse,E_lna,dE_lna)
 
 omega_m_collapse=(G_m0_density/((a_collapse**3)*(E_lna**2)))
 
 vir_collapse=delta_tr*omega_m_collapse*((a_collapse/a_tr)**3)/(y_vir**3) 
 
 end subroutine virial_quantities 
  
  
!=====================================================
!subroutine: propagate
!propagates the implamentations
!======================================================
  
  subroutine propagate(f_ini, t_end,f_end, mode)
   
   real(rt) :: f_ini, t_end !initial condition and collapse_time
   real(rt) :: f_end ! value of the inverse overdensity at t_end
   integer :: mode ! 0 for linear and 1 for non-linear
   
   real(rt) :: n1,n2 !power-like exponents
   real(rt) :: t_0 ! initial time
  
  INTEGER :: flag
  REAL(rt),ALLOCATABLE :: z(:),dz(:)
  
   ALLOCATE(z(2),dz(2))
  
    call power_law(a_ini,n1,n2)
  
     n1=max(n1,n2)
     z(1)=f_ini
     z(2)=-n1*z(1)/(z(1)*(1.0_rt+z(1)))
     t_0=Dlog(a_ini)
     
     flag=1
     
     if (mode.eq.0) then
     
       call Ftheta_linear(t_0,z,dz)
     
         do while((flag.lt.2).or.(flag.eq.4))
          CALL r8_rkf45(Ftheta_linear,2,z,dz,t_0,t_end,relerr,abserr,flag)
         end do
  
      f_end=z(1) 
      
      else if (mode.eq.1) then
      
      call Ftheta(t_0,z,dz)
     
         do while((flag.lt.2).or.(flag.eq.4))
          CALL r8_rkf45(Ftheta,2,z,dz,t_0,t_end,relerr,abserr,flag)
         end do
  
      f_end=z(1) 
      
      else 
      
      write(*,*) "ERROR, MODE ISNT 0 OR 1"
      
      end if
           
    DEALLOCATE(z,dz)
    
  end subroutine propagate
  
  
   !==========================================
   !==========derivative fgsl=================
   !==========================================
  
   !function f(x,parameters):= returns the reciprocal time evolved initial
   !overdensity x:=1/d, from a_ini to t_out.

     FUNCTION f(x,c_params)
   
     !=====input======!
     REAL(rt), value :: x
     type(c_ptr), value :: c_params
     real(rt) :: f
     !================!

     INTEGER :: flag
     REAL(rt),ALLOCATABLE :: z(:),dz(:)
     REAL(rt) :: t0,tf
     
     real(rt) :: n1,n2

     ALLOCATE(z(2),dz(2))

     call power_law(a_ini,n1,n2)
     n1=max(n1,n2)
     z(1)=x
     z(2)=-n1*z(1)/(z(1)*(1.0_rt+z(1)))
     t0=Dlog(a_ini)
     tf=Target_global

        
     flag=1 !single step mode
     
     call Ftheta(t0,z,dz)
     
     do while((flag.lt.2).or.(flag.eq.4))
       CALL r8_rkf45(Ftheta,2,z,dz,t0,tf,relerr,abserr,flag)  
     end do
     
     f=z(1)
              
     DEALLOCATE(z,dz)

   END FUNCTION f

 !subroutine a vs t: computes the time evolution of the model
   
   subroutine scale_vs_time(a0,dt)
     
    integer :: w
    real(rt) :: a0,dt !initial scale factor and resolution
    real(rt) :: E_a, dE_da,t,a

    a=a0 
    t=0.0_rt   
    open(newunit=w, file="scale_vs_time.dat", status="replace")

     do while(a.lt.1.0_rt)
    
      call Exp_func(a,E_a,dE_da)
      write(w, *) a,t,E_a
      a=a+(a*E_a*dt)
      t=t+dt
     end do 

    close(w)

  
   end subroutine scale_vs_time


END MODULE Spherical_collapse



PROGRAM MAIN_KERNEL
  USE nrtype
  USE Global_cosmology
  USE Cosmology_func
  USE Spherical_collapse
  IMPLICIT NONE

   INTEGER :: number,iter,o
   REAL(rt) :: t0,tf,dt
   REAL(rt) :: m,r,de,w,dlnw
   
   Real(rt) :: Eos_ws
   Real(rt) :: Eosp_ws
   Real(rt) :: Eospp_ws
   
   real(rt) :: c_time,tr_time
   real(rt) :: f_botom, f_top

  !===========WARNIN=================!
  !before all, we must first initialize the cosmology

  WRITE(*,*)"Initialize Cosmology..."
  WRITE(*,*)"    "
  CALL Model_initialization
  WRITE(*,*)"    "
  WRITE(*,*)"Done: the cosmology is ready!"
  
  f_botom=10.0_rt
  f_top=48000.0_rt
  c_time=0.0_rt
  tr_time=0.0_rt

WRITE(*,*)"Initial conditions: "
!call Collapse_time(f_top,c_time,tr_time)
!Call SC_Observables(200,f_botom,f_top)
write(*,*) "initial scale factor:",a_ini
write(*,*) "initial condition:",f_top
write(*,*) "turn around:",Dexp(tr_time)
write(*,*) "collapse_time:", Dexp(c_time)
WRITE(*,*)"End initial conditions:  "
WRITE(*,*)"       "
WRITE(*,*)"=====Main Input======"

number=5000
t0=DLog(1.0_rt/(1.0_rt+1.0e8_rt)) ! only plot form radiation era
tf=0.0_rt
dt=(tf-t0)/number

OPEN(newunit=o, file="EOS.dat", status="replace")

DO iter=1,number
t0=t0+dt

!CALL evol_var(t0,m,r,de,w,dlnw) !call from the diferential equations
!CALL spline_cubic_val(G_spline_numbe,Vec_time,Vec_EOS,Vec_EOSpp,&
                         !t0,Eos_ws,Eosp_ws,Eospp_ws) !call from the spline
!WRITE(o,*) 1.0_rt/Dexp(t0),m,r,de,w,dlnw,Eos_ws,abs(Eos_ws-w)
END DO

CLOSE(o)

call scale_vs_time(0.1_rt,0.001_rt)

WRITE(*,*) "Densities today..."
WRITE(*,*) "Matter:",G_m0_density
WRITE(*,*) "Radiation:",G_r0_density
WRITE(*,*) "Dark Energy:",G_de0_density
WRITE(*,*) "Total:",G_m0_density+&
           g_r0_density+G_de0_density

END PROGRAM MAIN_KERNEL
