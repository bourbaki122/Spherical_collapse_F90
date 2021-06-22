
!================cosmology-models=====================!
!for a better presentation of the code, each cosmology
!will have their own module, were their respective 
!equation of satete and the needed variables will be 
!computed, note that all modules must have the same structure:
!
!*) Equation of state
!*) Present data
!*) Splines
!=====================================================!

MODULE non_abelian_cosmology
  USE nrtype
  IMPLICIT NONE
  PUBLIC

  !the variable spline model, is to chek is the model has
  !spline or not, 1=yes, 0=no

  !=======Global_spine_variable=========!
  INTEGER, PARAMETER :: Spline_model=1
  Integer, PARAMETER :: Number_of_splines=200

  !========Eos-spline-varibles==========!
  INTEGER, PARAMETER :: EOS_nsplines=Number_of_splines
  REAL(rt), DIMENSION(EOS_nsplines) :: vec_w, vec_lna, vec_wpp
  INTEGER :: spline_flag0=1,spline_flagf=1 ! flag for the boundaries

  !=======w_integral-spline-variables=====!
  INTEGER, PARAMETER :: intEOS_nsplines=Number_of_splines
  REAL(rt), DIMENSION(intEOS_nsplines) :: vec_IntEos, vec_t!values and times of evaluation
  REAL(rt), DIMENSION(intEOS_nsplines) :: vec_IntEos_p,vec_IntEos_pp ! first and second derivatives!

  !========model-parameters=========!
  !
  REAL(rt) :: alpha=1.0e9_rt
  Real(rt) :: omegai=3.0e-27_rt !MASSLESS case on
  Real(rt) :: xi=3.0e-30_rt
  Real(rt) :: yi=0.0010_rt
  real(rt) :: zi=1.1e5_rt !~field
  real(rt) :: omega_ri=0.99998_rt !initial radiation density
  real(rt) :: z0=1.74e8_rt !initial redshift


CONTAINS

  !EOS := equation of state for the model
  ! althought evol_var returns those quantities
  ! we need this for the integration of the g(a)
  ! process

   Real(rt) FUNCTION EOS(lna)
    REAL(rt) :: lna
    REAL(rt) :: matter,radiation,dark_e,W,wp ! ancila data

    CALL evol_var(lna,matter,radiation,dark_e,w,wp)
    Eos=w
    
  END FUNCTION EOS


  !w_integral: implementation of the integral of the equation
  !of state, notice that we are integrating w(lna), so the equation of state
  !must be put in this way, w_integral=int(w/a da)

  SUBROUTINE w_integral(lna,wa_int)

    REAL(rt) :: lna,wa_int

    !======variables for the integration=====!
    INTEGER :: kf, iflag !flag and number of steps
    REAL(rt) :: t0,tf
    REAL(rt) :: eps=DSqrt(EPSILON(eps)) ! precission
    REAL(rt) :: wa_int_error ! error of the integration

    t0=0.0_rt  !ln(a=1), collapse now!!
    tf=lna ! evaluation of g(a)

    CALL q1da(EOS,t0,tf,eps,wa_int,wa_int_error,kf,iflag)

  END SUBROUTINE w_integral


  !Dynamic_eq:= equations of motion for the cosmological model

  SUBROUTINE Dynamic_eq(lna,z,dz_out)

    !====Variables of the subroutine======!
    REAL(rt) :: lna
    REAL(rt), DIMENSION(5) :: z,dz_out
    Real(rt) :: q

    !conventions:
    !z(1)-->x, z(2)-->y
    !z(3)-->w, z(4)-->z, z(5)-->rho_r
    
    !IMPORTANT: IF WE WANT THE MASSLESS CASE, TURN dz_out(3):=dw=0!
    !           AND SET THE INITIAL VALUE OF W=0
    !The equations are those in the mathematica nb.file

    q=(1 + z(2)**2 - z(3)**2 + (z(1)**2)*(1 - 12*alpha*(z(4)**4)) + z(5))/2.

    dz_out(1)=q*z(1) - ((2*(z(2)**2))/z(4) + 8*alpha*(z(1)**2)*(z(4)**3) +&
                         z(1)*(1 - 12*alpha*(z(4)**4)) )/(1 + 4*alpha*(z(4)**4))
                          
    dz_out(2)=z(2)*(-1.0_rt + q + (2*z(1))/z(4))

    dz_out(3)=z(3)*(q + (z(1)/z(4))) !MASSLESS case off
    !dz_out(3)=0.0_rt !MASSLESS case on

    dz_out(4)=z(1)-z(4)

    dz_out(5)=2*(-1.0_rt + q)*z(5)

  END SUBROUTINE Dynamic_eq


  !presnt data:= computes de values of m0_density, r0_density, and
  !de0_density, in the diff_model the curvature is set to zero

  SUBROUTINE present_data(m0_rho,r0_rho,de0_rho,k0_rho)

    REAL(rt) :: m0_rho,r0_rho,de0_rho,k0_rho
    REAL(rt) :: w,wp !auxiliar variable for the call
    CALL evol_var(0.0_rt,m0_rho,r0_rho,de0_rho,w,wp)
    k0_rho=0.0_rt

  END SUBROUTINE present_data

  !evol_var: computes de evolution of the density parameters and the
  !equantion of state at the desired time lna

  SUBROUTINE evol_var(lna,m_density,r_density,de_density,w,dw_dlna)

    REAL(rt) :: lna
    REAL(rt) :: w,q !w(lna) output value of the eos
    REAL(rt) :: dw_dlna ! derivative of the eos
    REAL(rt) :: m_density,r_density,de_density ! output values of densities

    !====varibles for the call of the diff_equation===!
    REAL(rt), ALLOCATABLE :: z(:),dz(:)
    REAL(rt) :: t0,tf
    INTEGER :: flag
    REAL(rt) :: r_error,abs_error

    !==auxial_variables derivative of W==========!
    real(rt) :: num,dem
    real(rt) :: d_num,d_dem

    !=====initialization

    ALLOCATE(z(5),dz(5))

    r_error=DSqrt(EPSILON(r_error))
    abs_error=DSqrt(EPSILON(abs_error))

    flag=1

    !=====initialization values for one of the proposed in the paper========!
    t0=DLog(1.0_rt/(1.0_rt+z0))
    tf=lna
    z(1)=xi
    z(2)=yi
    z(3)=DSqrt(omegai)*zi
    z(4)=zi
    z(5)=omega_ri

    CALL Dynamic_eq(t0,z,dz)

    CALL r8_rkf45(Dynamic_eq,5,z,dz,t0,tf,r_error,abs_error,flag)

    m_density=1 - z(2)**2 - z(3)**2 - z(1)**2*(1 + 4*alpha*z(4)**4) - z(5)
    r_density=z(5)
    de_density=z(1)**2 + z(2)**2 + z(3)**2 + (4*alpha*(z(1)**2)*(z(4)**4))
    q=(1 + z(2)**2 - z(3)**2 + z(1)**2*(1 - 12*alpha*z(4)**4) + z(5))/2.

    !w=(z(1)**2 - 12*alpha*(z(1)**2)*(z(4)**4) + z(2)**2 -z(3)**2)/de_density
    !w=w/3

    w=(1/3.)*((z(1)**2)*(1-(12*alpha*(z(4)**4)))+(z(2)**2)-(z(3)**2))/&
       ((z(1)**2)*(1+(4*alpha*(z(4)**4)))+(z(2)**2)+(z(3)**2))

    d_num=(2*z(2)*dz(2))-(2*z(3)*dz(3))+(2*z(1)*dz(1)*(1 - 12*alpha*(z(4)**4)))&
           -((z(1)**2)*48*alpha*(z(4)**3)*dz(4))

    d_dem=(2*z(2)*dz(2))+(2*z(3)*dz(3))+(2*z(1)*dz(1)*(1 + 4*alpha*(z(4)**4)))&
           +((z(1)**2)*16*alpha*(z(4)**3)*dz(4))

    num=((z(1)**2)*(1-(12*alpha*(z(4)**4)))+(z(2)**2)-(z(3)**2))

    dem=((z(1)**2)*(1+(4*alpha*(z(4)**4)))+(z(2)**2)+(z(3)**2))

    dw_dlna=(1/3.0_rt)*((d_num*dem)-(d_dem*num))/(dem**2)

    !check the constraints....
    !write(*,*) "check: ", w,z(3)

    DEALLOCATE(z,dz)

  END SUBROUTINE evol_var

  !============splined=================!

  SUBROUTINE EOS_Spline_initialization(Aux_vec_t,Aux_vec,Aux_vec_pp)

    INTEGER :: i_val
    REAL(rt) :: t0,tf,dt
    REAL(rt) :: dx,dy,dw,w,wp ! ancila data

    REAL(rt) :: wp_0,wp_f !derivatives of the eos at the edges

    !==============definition of the input vectors========!
    REAL(rt),DIMENSION(EOS_nsplines) :: Aux_vec_t,Aux_vec,Aux_vec_pp
    Aux_vec_t=0.0_rt
    Aux_vec=0.0_rt
    Aux_vec_pp=0.0_rt
    !====================================================!

    !initialization of vec_w
    WRITE(*,*)"           "
    WRITE(*,*)"EOS_SPLINE INITIALIZATION..."

    t0=DLog(1.0_rt/(1.0_rt+z0))
    tf=0.0_rt
    dt=(tf-t0)/REAL(EOS_nsplines-1,kind=8)

    DO i_val=1,EOS_nsplines
       vec_lna(i_val)=t0+REAL(i_val-1,kind=8)*dt
       CALL evol_var(vec_lna(i_val),dx,dy,dw,w,wp)
       vec_w(i_val)=w
    END DO

    spline_flag0=1
    spline_flagf=1

    !Boundary conditions at the edges of the interval
    !================================================!
    CALL evol_var(t0,dx,dy,dw,w,wp)
    wp_0=wp
    CALL evol_var(tf,dx,dy,dw,w,wp)
    wp_f=wp
    !================================================!

    CALL spline_cubic_set(EOS_nsplines,vec_lna,vec_w,spline_flag0, &
         wp_0,spline_flagf,wp_f,vec_wpp)

    !==========this is where we take the information back=======!
    Aux_vec_t=vec_lna
    Aux_vec=vec_w
    Aux_vec_pp=vec_wpp
    !==========================================================!

    WRITE(*,*)"EOS_SPLINE: check"

  END SUBROUTINE EOS_Spline_initialization


  !int_EOS_spline_initialization:= spline for the integral...
  SUBROUTINE int_EOS_spline_initialization(Aux_vec_t,Aux_vec,Aux_vec_pp)

    INTEGER :: i_val
    REAL(rt) :: t0,tf,dt
    REAL(rt) :: int_w

    INTEGER :: intEOs_flag0,intEOs_flagf
    REAL(rt) :: Eos_f, Eos_0 ! edges derivatives

    !==============definition of the input vectors========!
    REAL(rt),DIMENSION(intEOS_nsplines) :: Aux_vec_t,Aux_vec,Aux_vec_pp
    Aux_vec_t=0.0_rt
    Aux_vec=0.0_rt
    Aux_vec_pp=0.0_rt
    !=====================================================!

    WRITE(*,*)"           "
    WRITE(*,*)"Int_EOS_SPLINE INITIALIZATION..."

    t0=DLog(1.0_rt/(1.0_rt+z0))
    tf=0.0_rt
    dt=(tf-t0)/REAL(intEOS_nsplines-1,kind=8)

    DO i_val=1,intEOS_nsplines
       vec_t(i_val)=t0+REAL(i_val-1,kind=8)*dt
       CALL w_integral(vec_t(i_val),int_w)
       vec_IntEos(i_val)=int_w
    END DO

    intEOs_flag0=1
    intEOs_flagf=1

    Eos_0=Eos(t0)
    Eos_f=Eos(tf)

    CALL spline_cubic_set(intEOS_nsplines,vec_t,vec_IntEos,intEos_flag0, &
         Eos_0,intEos_flagf,Eos_f,vec_IntEos_pp)

    !==========Assignaments==================!
    Aux_vec_t=vec_t
    Aux_vec=vec_IntEos
    Aux_vec_pp=vec_IntEos_pp
    !========================================!

    WRITE(*,*)"Int_EOS_SPLINE: check"
    WRITE(*,*)"             "

  END SUBROUTINE int_EOS_spline_initialization

END MODULE non_abelian_cosmology





MODULE two_form_cosmology
  USE nrtype
  IMPLICIT NONE
  PUBLIC

  !the variable spline model, is to chek is the model has
  !spline or not, 1=yes, 0=no

  !=======Global_spine_variable=========!
  INTEGER, PARAMETER :: Spline_model=1
  Integer, PARAMETER :: Number_of_splines=200

  !========Eos-spline-varibles==========!
  INTEGER, PARAMETER :: EOS_nsplines=Number_of_splines
  REAL(rt), DIMENSION(EOS_nsplines) :: vec_w, vec_lna, vec_wpp
  INTEGER :: spline_flag0=1,spline_flagf=1 ! flag for the boundaries

  !=======w_integral-spline-variables=====!
  INTEGER, PARAMETER :: intEOS_nsplines=Number_of_splines
  REAL(rt), DIMENSION(intEOS_nsplines) :: vec_IntEos, vec_t!values and times of evaluation
  REAL(rt), DIMENSION(intEOS_nsplines) :: vec_IntEos_p,vec_IntEos_pp ! first and second derivatives!

  !========model-parameters=========!
  REAL(rt) :: lambda=2.0_rt
  REAL(rt) :: mu=5.5_rt
  Real(rt) :: x1=1.0e-13_rt
  Real(rt) :: x2=1.0e-14_rt
  real(rt) :: sigma=0
  real(rt) :: omega_b=1.0e-10_rt
  real(rt) :: omega_r=0.999961
  real(rt) :: z0=4.5e7_rt !initial redshift for the lmanda and mu, values

CONTAINS

  !EOS := equation of state for the model
  ! althought evol_var returns those quantities
  ! we need this for the integration of the g(a)
  ! process

   Real(rt) FUNCTION EOS(lna)
    REAL(rt) :: lna
    REAL(rt) :: matter,radiation,dark_e,W,wp ! ancila data

    CALL evol_var(lna,matter,radiation,dark_e,w,wp)

    Eos=w
  END FUNCTION EOS


  !w_integral: implementation of the integral of the equation
  !of state, notice that we are integrating w(lna), so the equation of state
  !must be put in this way, w_integral=int(w/a da)

  SUBROUTINE w_integral(lna,wa_int)

    REAL(rt) :: lna,wa_int

    !======variables for the integration=====!
    INTEGER :: kf, iflag !flag and number of steps
    REAL(rt) :: t0,tf
    REAL(rt) :: eps=DSqrt(EPSILON(eps)) ! precission
    REAL(rt) :: wa_int_error ! error of the integration

    t0=0.0_rt  !ln(a=1), collapse now!!
    tf=lna ! evaluation of g(a)

    CALL q1da(EOS,t0,tf,eps,wa_int,wa_int_error,kf,iflag)

  END SUBROUTINE w_integral


  !Dynamic_eq:= equations of motion for the cosmological model

  SUBROUTINE Dynamic_eq(lna,z,dz_out)

    !====Variables of the subroutine======!
    REAL(rt) :: lna
    REAL(rt), DIMENSION(5) :: z,dz_out

    !lna=0.0_rt !the condition for an autonomous system

    !conventions:
    !z(1)-->x_1, z(2)-->x_2
    !z(3)-->Sigma, z(4)-->rho_b, z(5)-->rho_r

    dz_out(1)=(3.0_rt/2)*z(1)*((z(1)**2)-(z(2)**2)+(z(3)**2)-1.0_rt &
         -((1.0_rt/3)*z(4))+((1.0_rt/3)*z(5))) &
         +(DSqrt(6.0_rt)/2)*((lambda*(z(2)**2))-(mu*z(4)))

    dz_out(2)=(1.0_rt/2)*z(2)*(3.0_rt*((z(1)**2)-(z(2)**2)+(z(3)**2)+1.0_rt) &
         -(DSqrt(6.0_rt)*lambda*z(1))-z(4)+z(5))

    dz_out(3)=((1.0_rt/2)*z(3)*(3.0_rt*((z(1)**2)-(z(2)**2)+(z(3)**2)-1.0_rt) &
         -z(4)+z(5)))-2.0_rt*z(4)

    dz_out(4)=z(4)*((3.0_rt*((z(1)**2)-(z(2)**2)+(z(3)**2)))+(4.0_rt*z(3)) &
         +1.0_rt+(DSqrt(6.0_rt)*mu*z(1))-z(4)+z(5))

    dz_out(5)=z(5)*((3.0_rt*((z(1)**2)-(z(2)**2)+(z(3)**2)))-1.0_rt-z(4)+z(5))

  END SUBROUTINE Dynamic_eq


  !presnt data:= computes de values of m0_density, r0_density, and
  !de0_density, in the diff_model the curvature is set to zero

  SUBROUTINE present_data(m0_rho,r0_rho,de0_rho,k0_rho)

    REAL(rt) :: m0_rho,r0_rho,de0_rho,k0_rho
    REAL(rt) :: w,wp !auciliary variable for the call
    CALL evol_var(0.0_rt,m0_rho,r0_rho,de0_rho,w,wp)
    k0_rho=0.0_rt

  END SUBROUTINE present_data

  !evol_var: computes de evolution of the density parameters and the
  !equantion of state at the desired time lna

  SUBROUTINE evol_var(lna,m_density,r_density,de_density,w,dw_dlna)

    REAL(rt) :: lna
    REAL(rt) :: w !w(lna) output value of the eos
    REAL(rt) :: dw_dlna ! derivative of the eos
    REAL(rt) :: m_density,r_density,de_density ! output values of densities

    !====varibles for the call of the diff_equation===!
    REAL(rt), ALLOCATABLE :: z(:),dz(:)
    REAL(rt) :: t0,tf
    INTEGER :: flag
    REAL(rt) :: r_error,abs_error

    !=====initialization

    ALLOCATE(z(5),dz(5))

    r_error=DSqrt(EPSILON(r_error))
    abs_error=DSqrt(EPSILON(abs_error))

    flag=1

    !=====initialization values for one of the proposed in the paper========!
    t0=DLog(1.0_rt/(1.0_rt+z0))
    
    tf=lna
    z(1)=x1
    z(2)=x2
    z(3)=sigma
    z(4)=omega_b
    z(5)=omega_r

    CALL Dynamic_eq(t0,z,dz)
      
    CALL r8_rkf45(Dynamic_eq,5,z,dz,t0,tf,r_error,abs_error,flag)
    
    
    m_density=1.0_rt-(z(1)**2)-(z(2)**2)-(z(3)**2)-z(4)-z(5)
    r_density=z(5)
    de_density=1.0_rt-m_density-r_density

    w=(3.0_rt*((z(1)**2)-(z(2)**2)+(z(3)**2)))-z(4)
    w=w/(3.0_rt*((z(1)**2)+(z(2)**2)+(z(3)**2)+z(4)))

    dw_dlna=(6.0_rt*((z(1)*dz(1))-(z(2)*dz(2))+(z(3)*dz(3))))-dz(4)
    dw_dlna=dw_dlna*(3.0_rt*((z(1)**2)+(z(2)**2)+(z(3)**2)+z(4)))
    dw_dlna=dw_dlna-((6.0_rt*((z(1)*dz(1))+(z(2)*dz(2))+(z(3)*dz(3)))+3.0_rt*dz(4))*&
         (3.0_rt*((z(1)**2)-(z(2)**2)+(z(3)**2))-z(4)))
    dw_dlna=dw_dlna/((3.0_rt*((z(1)**2)+(z(2)**2)+(z(3)**2)+z(4)))**2)

    DEALLOCATE(z,dz)

  END SUBROUTINE evol_var

  !============splined=================!

  SUBROUTINE EOS_Spline_initialization(Aux_vec_t,Aux_vec,Aux_vec_pp)

    INTEGER :: i_val
    REAL(rt) :: t0,tf,dt
    REAL(rt) :: dx,dy,dw,w,wp ! ancila data

    REAL(rt) :: wp_0,wp_f !derivatives of the eos at the edges

    !==============definition of the input vectors========!
    REAL(rt),DIMENSION(EOS_nsplines) :: Aux_vec_t,Aux_vec,Aux_vec_pp
    Aux_vec_t=0.0_rt
    Aux_vec=0.0_rt
    Aux_vec_pp=0.0_rt
    !====================================================!

    !initialization of vec_w
    WRITE(*,*)"           "
    WRITE(*,*)"EOS_SPLINE INITIALIZATION..."

    t0=DLog(1.0_rt/(1.0_rt+z0))
    tf=0.0_rt
    dt=(tf-t0)/REAL(EOS_nsplines-1,kind=8)

    DO i_val=1,EOS_nsplines
       vec_lna(i_val)=t0+REAL(i_val-1,kind=8)*dt
       CALL evol_var(vec_lna(i_val),dx,dy,dw,w,wp)
       vec_w(i_val)=w
    END DO

    spline_flag0=1
    spline_flagf=1

    !Boundary conditions at the edges of the interval
    !================================================!
    CALL evol_var(t0,dx,dy,dw,w,wp)
    wp_0=wp
    CALL evol_var(tf,dx,dy,dw,w,wp)
    wp_f=wp
    !================================================!

    CALL spline_cubic_set(EOS_nsplines,vec_lna,vec_w,spline_flag0, &
         wp_0,spline_flagf,wp_f,vec_wpp)

    !==========this is where we take the information back=======!
    Aux_vec_t=vec_lna
    Aux_vec=vec_w
    Aux_vec_pp=vec_wpp
    !==========================================================!

    WRITE(*,*)"EOS_SPLINE: check"

  END SUBROUTINE EOS_Spline_initialization


  !int_EOS_spline_initialization:= spline for the integral...
  SUBROUTINE int_EOS_spline_initialization(Aux_vec_t,Aux_vec,Aux_vec_pp)

    INTEGER :: i_val
    REAL(rt) :: t0,tf,dt
    REAL(rt) :: int_w

    INTEGER :: intEOs_flag0,intEOs_flagf
    REAL(rt) :: Eos_f, Eos_0 ! edges derivatives

    !==============definition of the input vectors========!
    REAL(rt),DIMENSION(intEOS_nsplines) :: Aux_vec_t,Aux_vec,Aux_vec_pp
    Aux_vec_t=0.0_rt
    Aux_vec=0.0_rt
    Aux_vec_pp=0.0_rt
    !=====================================================!

    WRITE(*,*)"           "
    WRITE(*,*)"Int_EOS_SPLINE INITIALIZATION..."

    t0=DLog(1.0_rt/(1.0_rt+z0))
    tf=0.0_rt
    dt=(tf-t0)/REAL(intEOS_nsplines-1,kind=8)

    DO i_val=1,intEOS_nsplines
       vec_t(i_val)=t0+REAL(i_val-1,kind=8)*dt
       CALL w_integral(vec_t(i_val),int_w)
       vec_IntEos(i_val)=int_w
    END DO

    intEOs_flag0=1
    intEOs_flagf=1

    Eos_0=Eos(t0)
    Eos_f=Eos(tf)

    CALL spline_cubic_set(intEOS_nsplines,vec_t,vec_IntEos,intEos_flag0, &
         Eos_0,intEos_flagf,Eos_f,vec_IntEos_pp)

    !==========Assignaments==================!
    Aux_vec_t=vec_t
    Aux_vec=vec_IntEos
    Aux_vec_pp=vec_IntEos_pp
    !========================================!

    WRITE(*,*)"Int_EOS_SPLINE: check"
    WRITE(*,*)"             "

  END SUBROUTINE int_EOS_spline_initialization

END MODULE two_form_cosmology





MODULE one_form_cosmology
  USE nrtype
  IMPLICIT NONE
  PUBLIC

  !the variable spline model, is to chek is the model has
  !spline or not, 1=yes, 0=no

    
  !=======Global_spine_variable=========!
  INTEGER, PARAMETER :: Spline_model=1
  Integer, PARAMETER :: Number_of_splines=500  
  
  !========Eos-spline-varibles==========!
  INTEGER, PARAMETER :: EOS_nsplines=Number_of_splines !number of knots
  REAL(rt), DIMENSION(EOS_nsplines) :: vec_w, vec_lna, vec_wpp
  INTEGER :: spline_flag0=1,spline_flagf=1 ! flag for the boundaries

  !=======w_integral-spline-variables=====!
  INTEGER, PARAMETER :: intEOS_nsplines=Number_of_splines ! number of splines
  REAL(rt), DIMENSION(intEOS_nsplines) :: vec_IntEos, vec_t!values and times of evaluation
  REAL(rt), DIMENSION(intEOS_nsplines) :: vec_IntEos_p,vec_IntEos_pp ! first and second derivatives!

  !========model-parameters=========!
  REAL(rt) :: lambda=2.0_rt
  REAL(rt) :: mu=5.0_rt
  real(rt) :: z0=5.5e7_rt

CONTAINS

  !EOS := equation of state for the model
  ! althought evol_var returns those quantities
  ! we need this for the integration of the g(a)
  ! process

  REAL(rt) FUNCTION EOS(lna)
    REAL(rt) :: lna
    REAL(rt) :: dx,dy,dw,w,wp ! ancila data

    CALL evol_var(lna,dx,dy,dw,w,wp)

    EOS=w
  END FUNCTION EOS

  !w_integral: implementation of the integral of the equation
  !of state, notice that we are integrating w(lna), so the equation of state
  !must be put in this way, w_integral=int(w/a da)

  SUBROUTINE w_integral(lna,wa_int)

    REAL(rt) :: lna,wa_int

    !======variables for the integration=====!
    INTEGER :: kf, iflag !flag and number of steps
    REAL(rt) :: t0,tf
    REAL(rt) :: eps=DSqrt(EPSILON(eps)) ! precission
    REAL(rt) :: wa_int_error ! error of the integration

    t0=0.0_rt  !ln(a=1), collapse now!!
    tf=lna ! evaluation of g(a)

    CALL q1da(EOS,t0,tf,eps,wa_int,wa_int_error,kf,iflag)

  END SUBROUTINE w_integral


  !Dynamic_eq:= equations of motion for the cosmological model

  SUBROUTINE Dynamic_eq(lna,z,dz_out)

    !====Variables of the subroutine======!
    REAL(rt) :: lna
    REAL(rt), DIMENSION(5) :: z,dz_out

    lna=0.0_rt !the condition for an autonomous system

    !conventions:
    !z(1)-->X, z(2)-->Y
    !z(3)-->Sigma, z(4)-->Omega_1, z(5)-->Omega_r

    dz_out(1)=(3.0_rt/2)*z(1)*((z(1)**2)-(z(2)**2)+(z(3)**2)-1.0_rt &
         +((1.0_rt/3)*z(4))+((1.0_rt/3)*z(5))) &
         +(DSqrt(6.0_rt)/2)*((lambda*(z(2)**2))-(2.0_rt*mu*z(4)))

    dz_out(2)=(1.0_rt/2)*z(2)*(3.0_rt*((z(1)**2)-(z(2)**2)+(z(3)**2)+1.0_rt) &
         -(DSqrt(6.0_rt)*lambda*z(1))+z(4)+z(5))

    dz_out(3)=((1.0_rt/2)*z(3)*(3.0_rt*((z(1)**2)-(z(2)**2)+(z(3)**2)-1.0_rt) &
         +z(4)+z(5)))+2.0_rt*z(4)

    dz_out(4)=z(4)*((3.0_rt*((z(1)**2)-(z(2)**2)+(z(3)**2)))+(4.0_rt*z(3)) &
         -1.0_rt+(2.0_rt*DSqrt(6.0_rt)*mu*z(1))+z(4)+z(5))

    dz_out(5)=z(5)*((3.0_rt*((z(1)**2)-(z(2)**2)+(z(3)**2)))-1.0_rt+z(4)+z(5))

  END SUBROUTINE Dynamic_eq


  !presnt data:= computes de values of m0_density, r0_density, and
  !de0_density, in the diff_model the curvature is set to zero

  SUBROUTINE present_data(m0_rho,r0_rho,de0_rho,k0_rho)

    REAL(rt) :: m0_rho,r0_rho,de0_rho,k0_rho
    REAL(rt) :: w,wp !auciliary variable for the call
    CALL evol_var(0.0_rt,m0_rho,r0_rho,de0_rho,w,wp)
    k0_rho=0.0_rt

  END SUBROUTINE present_data

  !evol_var: computes de evolution of the density parameters and the
  !equantion of state at the desired time lna

  SUBROUTINE evol_var(lna,m_density,r_density,de_density,w,dw_dlna)

    REAL(rt) :: lna
    REAL(rt) :: w !w(lna) output value of the eos
    REAL(rt) :: dw_dlna ! derivative of the eos
    REAL(rt) :: m_density,r_density,de_density ! output values of densities

    !====varibles for the call of the diff_equation===!
    REAL(rt), ALLOCATABLE :: z(:),dz(:)
    REAL(rt) :: t0,tf
    INTEGER :: flag
    REAL(rt) :: r_error,abs_error

    !=====initialization

    ALLOCATE(z(5),dz(5))

    r_error=DSqrt(EPSILON(r_error))
    abs_error=DSqrt(EPSILON(abs_error))

    flag=1

    !=====initialization values for one of the proposed in the paper========!
    t0=DLog(1.0_rt/(1.0_rt+z0))
    tf=lna
    z(1)=1.0e-13_rt
    z(2)=1.0e-14_rt
    z(3)=0
    z(4)=1.0e-5_rt
    z(5)=0.999961
    CALL Dynamic_eq(t0,z,dz)

    CALL r8_rkf45(Dynamic_eq,5,z,dz,t0,tf,r_error,abs_error,flag)

    de_density=(z(1)**2)+(z(2)**2)+(z(3)**2)+z(4)
    r_density=z(5)
    m_density=1.0_rt-de_density-r_density
    
    w=(3.0_rt*((z(1)**2)-(z(2)**2)+(z(3)**2)))+z(4)
    w=w/(3.0_rt*((z(1)**2)+(z(2)**2)+(z(3)**2)+z(4)))

    dw_dlna=(6.0_rt*((z(1)*dz(1))-(z(2)*dz(2))+(z(3)*dz(3))))+dz(4)
    dw_dlna=dw_dlna*(3.0_rt*((z(1)**2)+(z(2)**2)+(z(3)**2)+z(4)))
    dw_dlna=dw_dlna-((6.0_rt*((z(1)*dz(1))+(z(2)*dz(2))+(z(3)*dz(3)))+3.0_rt*dz(4))*&
         (3.0_rt*((z(1)**2)-(z(2)**2)+(z(3)**2))+z(4)))
    dw_dlna=dw_dlna/((3.0_rt*((z(1)**2)+(z(2)**2)+(z(3)**2)+z(4)))**2)

    DEALLOCATE(z,dz)

  END SUBROUTINE evol_var

  !============splined=================!

  SUBROUTINE EOS_Spline_initialization(Aux_vec_t,Aux_vec,Aux_vec_pp)

    INTEGER :: i_val
    REAL(rt) :: t0,tf,dt
    REAL(rt) :: dx,dy,dw,w,wp ! ancila data

    REAL(rt) :: wp_0,wp_f !derivatives of the eos at the edges

    !==============definition of the input vectors========!
    REAL(rt),DIMENSION(EOS_nsplines) :: Aux_vec_t,Aux_vec,Aux_vec_pp
    Aux_vec_t=0.0_rt
    Aux_vec=0.0_rt
    Aux_vec_pp=0.0_rt
    !====================================================!

    !initialization of vec_w
    WRITE(*,*)"           "
    WRITE(*,*)"EOS_SPLINE INITIALIZATION..."

    t0=DLog(1.0_rt/(1.0_rt+7.9e7_rt))
    tf=0.0_rt
    dt=(tf-t0)/REAL(EOS_nsplines-1,kind=8)

    DO i_val=1,EOS_nsplines
       vec_lna(i_val)=t0+REAL(i_val-1,kind=8)*dt
       CALL evol_var(vec_lna(i_val),dx,dy,dw,w,wp)
       vec_w(i_val)=w
    END DO

    spline_flag0=1
    spline_flagf=1

    !Boundary conditions at the edges of the interval
    !================================================!
    CALL evol_var(t0,dx,dy,dw,w,wp)
    wp_0=wp
    CALL evol_var(tf,dx,dy,dw,w,wp)
    wp_f=wp
    !================================================!

    CALL spline_cubic_set(EOS_nsplines,vec_lna,vec_w,spline_flag0, &
         wp_0,spline_flagf,wp_f,vec_wpp)

    !==========this is where we take the information back=======!
    Aux_vec_t=vec_lna
    Aux_vec=vec_w
    Aux_vec_pp=vec_wpp
    !==========================================================!

    WRITE(*,*)"EOS_SPLINE: check"

  END SUBROUTINE EOS_Spline_initialization


  !int_EOS_spline_initialization:= spline for the integral...
  SUBROUTINE int_EOS_spline_initialization(Aux_vec_t,Aux_vec,Aux_vec_pp)

    INTEGER :: i_val
    REAL(rt) :: t0,tf,dt
    REAL(rt) :: int_w

    INTEGER :: intEOs_flag0,intEOs_flagf
    REAL(rt) :: Eos_f, Eos_0 ! edges derivatives

    !==============definition of the input vectors========!
    REAL(rt),DIMENSION(intEOS_nsplines) :: Aux_vec_t,Aux_vec,Aux_vec_pp
    Aux_vec_t=0.0_rt
    Aux_vec=0.0_rt
    Aux_vec_pp=0.0_rt
    !=====================================================!

    WRITE(*,*)"           "
    WRITE(*,*)"Int_EOS_SPLINE INITIALIZATION..."

    t0=DLog(1.0_rt/(1.0_rt+7.9e7_rt))
    tf=0.0_rt
    dt=(tf-t0)/REAL(intEOS_nsplines-1,kind=8)

    DO i_val=1,intEOS_nsplines
       vec_t(i_val)=t0+REAL(i_val-1,kind=8)*dt
       CALL w_integral(vec_t(i_val),int_w)
       vec_IntEos(i_val)=int_w
    END DO

    intEOs_flag0=1
    intEOs_flagf=1

    Eos_0=Eos(t0)
    Eos_f=Eos(tf)

    CALL spline_cubic_set(intEOS_nsplines,vec_t,vec_IntEos,intEos_flag0, &
         Eos_0,intEos_flagf,Eos_f,vec_IntEos_pp)

    !==========Assignaments==================!
    Aux_vec_t=vec_t
    Aux_vec=vec_IntEos
    Aux_vec_pp=vec_IntEos_pp
    !========================================!

    WRITE(*,*)"Int_EOS_SPLINE: check"
    WRITE(*,*)"             "

  END SUBROUTINE int_EOS_spline_initialization

END MODULE one_form_cosmology


!====================================================
!.....EINSTEIN´S DE SITTER UNIVERSE.....
!====================================================

MODULE EDS_cosmology
  USE nrtype

  INTEGER, PARAMETER :: Spline_model=0

CONTAINS

  SUBROUTINE present_data(m0_rho,r0_rho,de0_rho,k0_rho)

    REAL(rt) :: m0_rho,r0_rho,de0_rho,k0_rho
    k0_rho=0.0_rt
    r0_rho=0.0_rt
    de0_rho=0.0_rt
    m0_rho=1.0_rt

  END SUBROUTINE present_data

  SUBROUTINE E_model_function(a,E_f,de_da)

    !==========initial_densities=======!
    REAL(rt) :: matter
    !=========input and output-parameters=========!
    REAL(rt) :: a,E_f,De_da

    matter=1.0_rt

    E_f=DSQRT(matter/(a**3.0_rt))
    de_da=(1.0_rt/2)*(1.0_rt/E_f)*(-3.0_rt*matter*(1/(a**4.0_rt)))

  END SUBROUTINE E_model_function

END MODULE EDS_cosmology

!====================================================
!....Lambda CDM cosmology.....
!====================================================

MODULE LCDM_cosmology
  USE nrtype

  INTEGER, PARAMETER :: spline_model=0
  real(rt), PARAMETER :: DE_today=0.6842_rt
  real(rt), PARAMETER :: M_today=0.3158_rt

  !0.750_rt  !0.6842_rt !planck 2018 values

CONTAINS

  SUBROUTINE present_data(m0_rho,r0_rho,de0_rho,k0_rho)

    REAL(rt) :: m0_rho,r0_rho,de0_rho,k0_rho
    k0_rho=0.0_rt
    r0_rho=0.0_rt
    de0_rho=DE_today
    m0_rho=M_today

  END SUBROUTINE present_data

  SUBROUTINE E_model_function(a,E_f,de_da)

    !==========initial_densities=======!
    REAL(rt) :: matter,d_energy

    !=========input and output-parameters=========!
    REAL(rt) :: a,E_f,De_da

    matter=M_today
    d_energy=DE_today

    E_f=DSQRT((matter/(a**3.0_rt))+d_energy)
    de_da=(1.0_rt/2)*(1.0_rt/E_f)*((-3.0_rt*matter*(1/(a**4.0_rt))))

  END SUBROUTINE E_model_function

END MODULE LCDM_cosmology

!====================================================
!====================================================
!....CASSIMIR cosmology.....
!====================================================
!====================================================

MODULE Cassimir_cosmology
  USE nrtype

  INTEGER, PARAMETER :: Spline_model=0
  real(rt), PARAMETER :: M_today=0.25_rt!0.3158_rt
  real(rt), PARAMETER :: DE_today=0.75_rt!0.6842_rt
  real(rt), PARAMETER :: cass_today=-0.00035

CONTAINS

  SUBROUTINE present_data(m0_rho,r0_rho,de0_rho,k0_rho)

    REAL(rt) :: m0_rho,r0_rho,de0_rho,k0_rho
    k0_rho=0.0_rt
    r0_rho=0.0_rt
    de0_rho=DE_today
    m0_rho=M_today

  END SUBROUTINE present_data

  real(rt) function Eos(lna)
    real(rt) :: lna,aux
    real(rt) :: d_energy, cassimir
    d_energy=DE_today
    cassimir=cass_today

    aux=(-1.0_rt/3)*((3.0_rt*d_energy*(Dexp(lna)**4))+Cassimir)/ &
        ((d_energy*(Dexp(lna)**4))-cassimir)

    Eos=aux

  end function eos

  SUBROUTINE E_model_function(a,E_f,de_da)

    !==========initial_densities=======!
    REAL(rt) :: matter,d_energy
    REAL(rt) :: Cassimir

    !=========input and output-parameters=========!
    REAL(rt) :: a,E_f,De_da

    matter=M_today
    d_energy=DE_today
    Cassimir=cass_today

    E_f=DSQRT((matter/(a**3.0_rt))+d_energy-(Cassimir/(a**4.0_rt)))
    de_da=(1.0_rt/2)*(1.0_rt/E_f)*((-3.0_rt*matter*(1/(a**4.0_rt)))+&
    (4.0_rt*cassimir/(a**5.0_rt)))

  END SUBROUTINE E_model_function

END MODULE Cassimir_cosmology
