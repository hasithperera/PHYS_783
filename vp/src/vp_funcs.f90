!=============================================================================!
!=============================================================================!
!       Functions for   VP: 1D1V Vlasov-Poisson Nonlinear Solver WAVES:       !
!=============================================================================!
!==================================================================== !
!
!------------------------------------------------------------------------------
!                             Copyright, 2015
!                                Greg Howes
!
!------------------------------------------------------------------------------
!  
!------------------------------------------------------------------------------
module vp_funcs
  use vp_data, only: nx,nv,pi
  implicit none

  public :: get_phi, ddx, ddv
  public :: initial_conditions
  public :: adams_bashforth3
  public :: calc_dfsbdt,calc_dfswldt,calc_dfswnbdt,calc_dfswnwdt,calc_dfscdt
  public :: fparity
  public :: probe_output,periodic_output,energy_output,phase_output

  private :: energy_full_ie,number_cons
  
 contains
!------------------------------------------------------------------------------
!                           Greg Howes, 2015
!------------------------------------------------------------------------------
!Compute initial conditions for distribution function
!------------------------------------------------------------------------------
   subroutine initial_conditions
     use vp_data, only: fs,fs0,fsb,fsw,fswl,fswnb,fswnw,fswnc,fsc
     use vp_data, only: v,nl,x,dx, spec, nspec, kle
     use vp_data, only: nk1,k1,dn,ic,eig_opt,w_ic,phi0
     use vp_data, only: t,CN_smooth,CN_N,runname
     use vp_disp_funcs, only: vp_disp,rtsec
    implicit none
    integer :: i,j, ispec
    complex :: w0                 !Initial guess
    complex :: omega                                !Refined Omega
    complex :: dfc                                  !Complex delta f
    complex :: om1,om2                              !Bracket Values
    integer :: iflag                               !Flag for Root search
    real, parameter :: prec=0.01                  !Root Finding precision
    real, parameter :: tol=1.0E-13                  !Root Search Tolerance   
    complex, parameter :: ii=(0.,1.)                !Imaginary Unit
    character(150) :: writeName
    
		write(*,*) "ahe init -0" 
		write(*,*) "ic",ic
    select case(ic)
    !-----------------------------------------
    case(1) !Electron density perturbation
       !NOTE: KGK SEPT 2016: Updated to allow for arbitrary number of velocity distributions
       !Perturbation imposed on first (reference) species
       !Set wavenumber of initial density perturbation
       k1=real(nk1)/nl
       write(*,*) "ahe init -1" 
       !SET uniform Maxwellian equilibrium initial conditions
       ! NOTE that there is assumed no spatial variation of Maxwellian fs0
       do ispec=1,nspec
					write(*,*) "ahe bg f init bump at +/- u=7.5"
					write(*,*) ispec
          do i=-nv/2,nv/2
             !Set background
             fs0(ispec,:,i)=spec(ispec)%D_s*exp(-(v(i)-spec(ispec)%V_s)**2./2.)/sqrt(2.*pi)
						! + 0.5*exp(-(v(i)-(7.5))**2./2.)/sqrt(2.*pi)+ 0.5*exp(-(v(i)-(-7.5))**2./2.)/sqrt(2.*pi)


          enddo
       enddo
       
       !SET initial conditions for perturbed distribution functions
       do ispec=1,nspec
          if (spec(ispec)%p_S) then
             write(*,*)'spec: ',ispec,' perturb!'
             write(*,*)k1,dn
          endif
          fsb(ispec,:,:)=0. !Ballistic term is zero initially
          fswnb(ispec,:,:)=0. !NL Wave ballistic term is zero initially
          fswnw(ispec,:,:)=0. !NL Wave wave term is zero initially
          fswnc(ispec,:,:)=0. !NL Wave coll term is zero initially
          fsc(ispec,:,:)=0. !Collision term is zero initially
          if (spec(ispec)%p_S) then
             do i=1,nx
                fswl(ispec,i,:)=fs0(ispec,i,:)*dn*sin(k1*x(i))  
             enddo
          else
             fswl(ispec,:,:)=0.
          endif
          !Sum to get full distribution function 
          fsw(ispec,:,:)= fswl(ispec,:,:) + fswnb(ispec,:,:) + fswnw(ispec,:,:)+ fswnc(ispec,:,:)
          fs(ispec,:,:)=fs0(ispec,:,:)+ fsb(ispec,:,:)+ fsw(ispec,:,:)+ fsc(ispec,:,:)!!
       enddo
       
    !-----------------------------------------
    case(2) !Exact Eigenfunction initialization
       !NOTE: GGH 18 FEB 2016: Wave term separation not coded here yet!
       !NOTE: KGK 05 MAR 2016: Wave term separation Added. Testing for rho = 0
       !NOTE: KGK SEPT 2016: Updated to allow for arbitrary number of velocity distributions
       
       !USES DIMENSIONLESS PARAMETERS: kle & spec(1:npsec)
       kle=real(nk1)/nl

       !Choose wave frequency to initialize
       select case(eig_opt)
       case(0) !Specify initial guess
          w0= w_ic
       case(1) !Langmuir Wave
          w0=cmplx(sqrt(1.+3.*kle**2.), -sqrt(pi/8)/(abs(kle)**3.)*exp(-1.*(1./(2.*kle**2.)+1.5)) )
       case(2) !Ion Acoustic Wave
          w0=cmplx(kle/sqrt(spec(2)%mu_s),-kle/sqrt(spec(2)%mu_s)* sqrt(pi/8.)*spec(2)%tau_s**-1.5*exp(-1./(2.*spec(2)%tau_s)) )
       end select

       write(writeName,'(2a,i0)')trim(runname),'.vp_ics.s',ispec
       open(19,file=trim(writeName),status='replace')
       !Refine initial guess w0=============================================
       !Set first two guesses bracketing old value
       om1=w0*(1.-prec)
       om2=w0*(1.+prec)
       !Check Parameters
       write(*,'(a,f6.3,a,2es15.7,a,2es15.7)') &
            'kle= ',kle,'  om= ',w0,'    disp= ',vp_disp(w0)
       write(19,'(a,f6.3,a,2es15.7,a,2es15.7)') &
            'kle= ',kle,'  om= ',w0,'    disp= ',vp_disp(w0)
       do ispec=1,nspec
          write(*,'(a,i0,a,es11.3,a,es11.3,a,es11.3,a,es11.3,a,es11.3)') &
               'Species ',ispec,': T_s/T_R=',spec(ispec)%tau_s,&
               ' m_s/m_R=',spec(ispec)%mu_s,&
               ' q_s/q_R=',spec(ispec)%Q_s,&
               ' n_s/n_R=',spec(ispec)%D_s,&
               ' v_s/v_R=',spec(ispec)%V_s
          write(19,'(a,i0,a,es11.3,a,es11.3,a,es11.3,a,es11.3,a,es11.3)') &
               'Species ',ispec,': T_s/T_R=',spec(ispec)%tau_s,&
               ' m_s/m_R=',spec(ispec)%mu_s,&
               ' q_s/q_R=',spec(ispec)%Q_s,&
               ' n_s/n_R=',spec(ispec)%D_s,&
               ' v_s/v_R=',spec(ispec)%V_s
       enddo

       !Find root of dispersion relation-----------------------------------
       omega=cmplx(0.,0.)
       omega=rtsec(vp_disp,om1,om2,tol,iflag)

       !SET uniform Maxwellian equilibrium initial conditions on fe
       ! NOTE that there is assumed no spatial variation of Maxwellian fs0
       do i=-nv/2,nv/2
          !Set background
          do ispec=1,nspec
						 ! added an initial distribution with a bump at u=3.5
             fs0(ispec,:,i)=spec(ispec)%D_s*exp(-(v(i)-spec(ispec)%V_s)**2./2.)/sqrt(2.*pi) 
! + 0.5*exp(-(v(i)-(3.5))**2./2.)/sqrt(2.*pi)
          enddo
       enddo

       !Compute the initial distribution functions fswl
       do i=-nv/2,nv/2
          do j=1,nx
             do ispec=1,nspec
                dfc=phi0*(spec(ispec)%Q_s/spec(ispec)%tau_s)*fs0(ispec,j,i)*&
                     kle*(v(i)-spec(ispec)%V_s)*sqrt(spec(ispec)%tau_s/spec(ispec)%mu_s)/&
                     (omega-kle*v(i)*sqrt(spec(ispec)%tau_s/spec(ispec)%mu_s)) *exp(-ii*kle*x(j))                
                fswl(ispec,j,i)=real(dfc)
             enddo
          enddo
       enddo

! !remove bump init on sp=1
!       !AHE - added a bump on the tail after the wave purtabation
!				write(*,*) "AHE init bump: v=3.5 vth"
!       ! NOTE that there is assumed no spatial variation of Maxwellian fs0
!       do i=-nv/2,nv/2
!          !Set background
!          do ispec=1,nspec
!						 ! added an initial distribution with a bump at u=3.5
!             fs0(ispec,:,i)=fs0(ispec,:,i) + 0.5*exp(-(v(i)-(7.5))**2./2.)/sqrt(2.*pi)
!          enddo
!       enddo



       !Smooth the initial function using a simple diffusion operator
       if (CN_smooth) &
            !smooths fswl
            call Crank_Nicholson_smoothing(CN_N)

       do ispec=1,nspec
          fsb(ispec,:,:)=0. !Ballistic term is zero initially
          fswnb(ispec,:,:)=0. !NL Wave ballistic term is zero initially
          fswnw(ispec,:,:)=0. !NL Wave wave term is zero initially
          fswnc(ispec,:,:)=0. !NL Wave coll term is zero initially
          fsc(ispec,:,:)=0. !Collision term is zero initially
           !Sum to get full distribution function 
          fsw(ispec,:,:)= fswl(ispec,:,:) + fswnb(ispec,:,:) + fswnw(ispec,:,:)+ fswnc(ispec,:,:)
          fs(ispec,:,:)=fs0(ispec,:,:)+ fsb(ispec,:,:)+ fsw(ispec,:,:)+ fsc(ispec,:,:)
       enddo
       
       
       write(*,'(a,es15.7,a,es15.7,a,2es15.7)')&
            'Initial Guess:   wr= ',real(w0),'   wi= ',aimag(w0),'    disp= ',vp_disp(w0)
       write(*,'(a,es15.7,a,es15.7,a,2es15.7)')&
            'Refined Guess:   wr= ',real(omega),'   wi= ',aimag(omega),'    disp= ',vp_disp(omega)
       write(19,'(a,es15.7,a,es15.7,a,2es15.7)')&
            'Initial Guess:   wr= ',real(w0),'   wi= ',aimag(w0),'    disp= ',vp_disp(w0)
       write(19,'(a,es15.7,a,es15.7,a,2es15.7)')&
            'Refined Guess:   wr= ',real(omega),'   wi= ',aimag(omega),'    disp= ',vp_disp(omega)
       close(19)
       
    case default; write(*,*)'ERR: Unknown sweep parameter'; stop
    end select

    !Initialize time
    t=0.
    
  end subroutine initial_conditions

!------------------------------------------------------------------------------
!                           Kristopher Klein, 2016
!------------------------------------------------------------------------------
! Smooth the initial wave distribution using a Crank-Nicholson diffusion operator
!-=-=-=-==-=-
subroutine Crank_Nicholson_smoothing(ntau)
  use vp_data, only: fsw,v,nx,nv,dt,dv,fs0,fswl,nspec,runname
  implicit none
  !Passed
  integer :: ntau !number of time steps for evaulation
  !Local
  !Vectors for Thomas Algorithm for Tridiagonal Matrix Solution
  real, dimension(-nv/2:nv/2) :: aa,bb,cc,dd
  real :: m !Holding variable for Thomas Algorithm
  real, dimension(:,:,:,:), allocatable :: ff_temp
  integer :: i,j,k,ispec
  !rr = nu dt / (2 dv^2)
  real,dimension(-nv/2:nv/2) :: rr
  real :: nu = 1.5625E-3
  character(150) :: writeName
  !In general, the diffusion coefficent can be velocity dependent

  write(*,'(a,i0,a)') &
       'Smoothing Velocity Space with Crank-Nicholson Diffusion Operator; ',ntau,' steps'

  !set up diffusion matrix
  allocate(ff_temp(1:nx,-nv/2:nv/2,1:ntau,1:nspec)); ff_temp = 0.
  do i = -nv/2,nv/2
     rr(i) = ((nu)/(dv**2.))
     do j = 1,nx
        do ispec=1,nspec
           ff_temp(j,i,1,ispec)=fswl(ispec,j,i)
        enddo
     enddo
  enddo
  
  do ispec = 1, nspec
     do j = 1, nx
        do k = 2, ntau
        !-=-=-=-           
           !define diffusion terms
           i=-nv/2 
           aa(i) = 0; bb(i) = 1. + 2.*rr(i); cc(i) = -rr(i)
           dd(i)=(1.-2.*rr(i))*ff_temp(j,i,k-1,ispec) + &
                rr(i+1)*ff_temp(j,i+1,k-1,ispec)
           do i = -nv/2+1,nv/2-1
              aa(i) = -rr(i); bb(i) = 1. + 2. * rr(i); cc(i) = -rr(i)
              dd(i)=rr(i-1)*ff_temp(j,i-1,k-1,ispec)+&
                   (1.-2.*rr(i))*ff_temp(j,i,k-1,ispec) + &
                   rr(i+1)*ff_temp(j,i+1,k-1,ispec)
           enddo
           i=nv/2 
           dd(i)=rr(i-1)*ff_temp(j,i-1,k-1,ispec)+&
                (1.-2.*rr(i))*ff_temp(j,i,k-1,ispec)                     
           aa(i) = -rr(i); bb(i) = 1. + 2.*rr(i); cc(i) = 0.
           !-=-=-=-=-
           !Thomas Algorithm: Forward Reduction
           do i = -nv/2+1,nv/2
              m = aa(i)/bb(i-1)
              bb(i) = bb(i) - m*(cc(i-1))
              dd(i) = dd(i) - m*(dd(i-1))
           enddo
           !Thomas Algorithm: Backward Reduction
           !ff_temp(j,nv/2,k,ispec) = ff_temp(j,nv/2,k-1,ispec)
           ff_temp(j,nv/2,k,ispec) = dd(nv/2)/bb(nv/2)
           do i = nv/2-1,-nv/2,-1
              ff_temp(j,i,k,ispec) = &
                   (dd(i)-cc(i)*ff_temp(j,i+1,k,ispec))/bb(i)
           enddo
        enddo
        !-=-=-=-
        enddo !end spatial loop
     enddo !end species loop
  
     do ispec = 1,nspec
        write(writeName,'(2a,i0)')trim(runname),'.vp_f_smooth.s',ispec
        open(unit=11,file=trim(writeName),status='replace')
        k = 1
        do i = -nv/2,nv/2
           write(11,'(2i5,2es14.4)')&
                k,i,v(i),ff_temp(j,i,k,ispec)
        enddo
        write(11,*);write(11,*)
        k = ntau
        do i = -nv/2,nv/2
           write(11,'(2i5,2es14.4)')&
                k,i,v(i),ff_temp(j,i,k,ispec)
        enddo
        write(11,*);write(11,*)
        close(11)
     enddo
  !send back smoothed matrix
     do ispec=1,nspec
        do i = -nv/2,nv/2
           do j = 1,nx
              fswl(ispec,j,i)=ff_temp(j,i,ntau,ispec)
           enddo
        enddo
     enddo
  
end subroutine Crank_Nicholson_smoothing

!------------------------------------------------------------------------------
!Greg Howes, 2018
!------------------------------------------------------------------------------
! Compute total conservation of number
subroutine number_cons(ns1,ns01,nsb1,nsw1,nsc1)
  use vp_data, only: nspec,dv,dx,fs,fs0,fsb,fsw,fsc
  implicit none
 !Passed
  real, intent(out), dimension(1:nspec) :: ns1
  real, intent(out), dimension(1:nspec) :: ns01
  real, intent(out), dimension(1:nspec) :: nsb1
  real, intent(out), dimension(1:nspec) :: nsw1
  real, intent(out), dimension(1:nspec) :: nsc1
  !Local
  integer :: i,j,ispec
  real :: dFactor

  !Initialize
  ns1=0.
  ns01=0.
  nsb1=0.
  nsw1=0.
  nsc1=0.
  
  !Integrate over velocity and position to compute total particle number
  dFactor=dv*dx
  do ispec=1,nspec
     ns1(ispec)=sum(fs(ispec,:,:))*dFactor
     ns01(ispec)=sum(fs0(ispec,:,:))*dFactor
     nsb1(ispec)=sum(fsb(ispec,:,:))*dFactor
     nsw1(ispec)=sum(fsw(ispec,:,:))*dFactor
     nsc1(ispec)=sum(fsc(ispec,:,:))*dFactor
  enddo
    
end subroutine number_cons

!------------------------------------------------------------------------------
!Greg Howes, 2016
!------------------------------------------------------------------------------
!Compute energy of equilibrium, ballistic, and wave terms
!(ion/electron)
!------------------------------------------------------------------------------
!K.G. Klein; 2016
!Edited to allow for an arbitrary number of velocity distribution
!functions
subroutine energy_full_ie(fs1,fs01,fsb1,fsw1,fswl1,fswnb1,fswnw1,fswnc1,fsc1,dfs0dv1,decs1,intdecs1,qcs1,v1,dv1,dphidx1,dx1,etot1,ephi1,efs1,efs01,efsb1,efsw1,efswl1,efswnb1,efswnw1,efswnc1,efsc1,edecs1,eintdecs1,eqcs1,ekos1,ekos01,ekos11,ekos21)
    use vp_data, only: spec,nspec,coll
    implicit none
    !Passed
    real, intent(in), dimension(1:nspec,1:nx,-nv/2:nv/2) :: fs1
    real, intent(in), dimension(1:nspec,1:nx,-nv/2:nv/2) :: fs01
    real, intent(in), dimension(1:nspec,1:nx,-nv/2:nv/2) :: fsb1
    real, intent(in), dimension(1:nspec,1:nx,-nv/2:nv/2) :: fsw1
    real, intent(in), dimension(1:nspec,1:nx,-nv/2:nv/2) :: fswl1
    real, intent(in), dimension(1:nspec,1:nx,-nv/2:nv/2) :: fswnb1
    real, intent(in), dimension(1:nspec,1:nx,-nv/2:nv/2) :: fswnw1
    real, intent(in), dimension(1:nspec,1:nx,-nv/2:nv/2) :: fswnc1
    real, intent(in), dimension(1:nspec,1:nx,-nv/2:nv/2) :: fsc1
    real, intent(in), dimension(1:nspec,1:nx,-nv/2:nv/2) :: dfs0dv1
    real, intent(in), dimension(1:nspec,1:nx,-nv/2:nv/2) :: decs1    
    real, intent(in), dimension(1:nspec,1:nx,-nv/2:nv/2) :: intdecs1  
    real, intent(in), dimension(1:nspec,1:nx,-nv/2:nv/2) :: qcs1      
    real, intent(in), dimension(-nv/2:) :: v1
    real, intent(in) :: dv1
    real, intent(in), dimension(1:nx) :: dphidx1
    real, intent(in) :: dx1 
    real, intent(out) :: etot1
    real, intent(out) :: ephi1
    real, intent(out), dimension(1:nspec) :: efs1
    real, intent(out), dimension(1:nspec) :: efs01
    real, intent(out), dimension(1:nspec) :: efsb1
    real, intent(out), dimension(1:nspec) :: efsw1
    real, intent(out), dimension(1:nspec) :: efswl1
    real, intent(out), dimension(1:nspec) :: efswnb1
    real, intent(out), dimension(1:nspec) :: efswnw1
    real, intent(out), dimension(1:nspec) :: efswnc1
    real, intent(out), dimension(1:nspec) :: efsc1
    real, intent(out), dimension(1:nspec) :: edecs1
    real, intent(out), dimension(1:nspec) :: eintdecs1
    real, intent(out), dimension(1:nspec) :: eqcs1
    real, intent(out), dimension(1:nspec) :: ekos1
    real, intent(out), dimension(1:nspec) :: ekos01
    real, intent(out), dimension(1:nspec) :: ekos11
    real, intent(out), dimension(1:nspec) :: ekos21
    !Local
    integer :: i,j,ispec
    real :: dFactor

    !Initialize
    etot1=0.
    ephi1=0.
    efs1=0.
    efs01=0.
    efsb1=0.
    efsw1=0.
    efswl1=0.
    efswnb1=0.
    efswnw1=0.
    efswnc1=0.
    efsc1=0.
    edecs1=0.
    eintdecs1=0.
    eqcs1=0.
    ekos1=0.
    ekos01=0.
    ekos11=0.
    ekos21=0.
    
    !Integrate over velocity and position to obtain energy parts
    !Energy normalized to k_B*T_reference
    do i=1,nx
       do ispec=1,nspec
          dFactor=dv1*dx1*spec(ispec)%tau_S!*sqrt(spec(ispec)%D_s)
          do j=-nv/2,nv/2
             efs1(ispec)=efs1(ispec)+&
                  v1(j)**2.*fs1(ispec,i,j)*dFactor
             efs01(ispec)=efs01(ispec)+&
                  v1(j)**2.*fs01(ispec,i,j)*dFactor
             efsb1(ispec)=efsb1(ispec)+&
                  v1(j)**2.*fsb1(ispec,i,j)*dFactor
             efsw1(ispec)=efsw1(ispec)+&
                  v1(j)**2.*fsw1(ispec,i,j)*dFactor
             efswl1(ispec)=efswl1(ispec)+&
                  v1(j)**2.*fswl1(ispec,i,j)*dFactor
             efswnb1(ispec)=efswnb1(ispec)+&
                  v1(j)**2.*fswnb1(ispec,i,j)*dFactor
             efswnw1(ispec)=efswnw1(ispec)+&
                  v1(j)**2.*fswnw1(ispec,i,j)*dFactor
             efswnc1(ispec)=efswnc1(ispec)+&
                  v1(j)**2.*fswnc1(ispec,i,j)*dFactor
             efsc1(ispec)=efsc1(ispec)+&
                  v1(j)**2.*fsc1(ispec,i,j)*dFactor
             !Compute integrated collisional energy diagnostics
             !NOTE: dFactor factor already included in earlier calculation in vp.f90
             if (coll) then
                edecs1(ispec)= edecs1(ispec)+ &
                     decs1(ispec,i,j)
                eintdecs1(ispec)= eintdecs1(ispec)+ &
                     intdecs1(ispec,i,j)
                eqcs1(ispec)= eqcs1(ispec)+ &
                     qcs1(ispec,i,j)
             endif
         
          !Kruskal-Oberman Energy=================================
             !NEED TO HANDLE v=0 point specially!
             !NOTEL: COLLISION TERMS HAVEN"T BEEN ADDED INTO THIS CALCULATION YET!
          if (j .ne. 0) then

             ekos1(ispec)=ekos1(ispec)+&
                  v1(j)*fs1(ispec,i,j)**2./&
                  (-1.*dfs0dv1(ispec,i,j))*dFactor
             !fs0
             ekos01(ispec)=ekos01(ispec)+&
                  v1(j)*fs01(ispec,i,j)**2./&
                  (-1.*dfs0dv1(ispec,i,j))*dFactor
             !fs1
             ekos11(ispec)=ekos11(ispec)+&
                  v1(j)*2.*fs01(ispec,i,j)*&
                  (fsb1(ispec,i,j)+fsw1(ispec,i,j))/&
                  (-1.*dfs0dv1(ispec,i,j))*dFactor
             !fs1^2
             ekos21(ispec)=ekos21(ispec)+&
                  v1(j)*(fsb1(ispec,i,j)+fsw1(ispec,i,j))**2./&
                  (-1.*dfs0dv1(ispec,i,j))*dFactor
          else !Handle 0/0 for v=0 and df0/dv=0.
             !NOTE: Here I am assuming a Maxwellian f0!!!
             !Total
             ekos1(ispec)=ekos1(ispec)+&
                  fs1(ispec,i,j)**2./fs01(ispec,i,j)*dFactor
             !fs0
             ekos01(ispec)=ekos01(ispec)+&
                  fs01(ispec,i,j)**2./fs01(ispec,i,j)*dFactor
             !fs1
             ekos11(ispec)=ekos11(ispec)+&
                  2.*fs01(ispec,i,j)*&
                  (fsb1(ispec,i,j)+fsw1(ispec,i,j))/&
                  fs01(ispec,i,j)*dFactor
             !fs1^2
             ekos21=ekos21+(fsb1(ispec,i,j)+fsw1(ispec,i,j))**2./&
                  fs01(ispec,i,j)*dFactor
          endif
       enddo
    enddo
    ephi1=ephi1+dphidx1(i)**2.*dx1
 enddo
 etot1=ephi1+sum(efs1(:))
 
  end subroutine energy_full_ie

!------------------------------------------------------------------------------
!                           Greg Howes, 2015
!------------------------------------------------------------------------------
! Compute potential from rho using Green's function approach
!===
! NOTE: 09/21/2016 K.G. Klein: Potential now normalized by q_R/T_R,
!       where the reference species is the first species in the input namelist
  subroutine get_phi(rho1,nl1,x1,dx1,phi1)
    implicit none
    !Passed
    real, intent(in), dimension(1:nx) :: rho1
    real, intent(in) :: nl1
    real, intent(in), dimension(1:nx) :: x1
    real, intent(in) :: dx1
    real, intent(out), dimension(1:nx) :: phi1
    !Local
    real, dimension(1:nx) :: il
    real, dimension(1:nx) :: iu
    integer :: i,j

    !Intermediate integrals to get the potential using a Green's function
    il(:)=0.
    iu(:)=0.
    do i=1,nx
       do j=1,i-1
          il(i)=il(i)+(pi*nl1+x1(j))*rho1(j)*dx1
       enddo
       il(i)=il(i)+0.5*(pi*nl1+x1(i))*rho1(i)*dx1
       iu(i)=iu(i)+0.5*(pi*nl1-x1(i))*rho1(i)*dx1
       do j=i+1,nx
         iu(i)=iu(i)+(pi*nl1-x1(j))*rho1(j)*dx1
      enddo
   enddo
   
   !Compute the potential and electric field
   phi1(:)=-1./(2*pi*nl1)*((pi*nl1-x1(:))*il(:)+(pi*nl1+x1(:))*iu(:))

 end subroutine get_phi
!------------------------------------------------------------------------------
!                           Greg Howes, 2015
!------------------------------------------------------------------------------
! NOTE: Density must be initialized such that SUM_s int q_s f_s dv =0
! This means overall neutrality, otherwise Green's function for phi fails.
! Compute Charge density rho(x)
! ===
! 09/2016- K.G.Klein: Edited to allow for 'nspec' distribution funtions rather than 2
  subroutine get_rho(fs1,dv1,rho1)
    use vp_data, only: spec,nspec
    implicit none
    !Passed
    real, intent(in), dimension(nspec,1:nx,-nv/2:nv/2) :: fs1
    real, intent(in) :: dv1
    real, intent(out), dimension(1:nx) :: rho1
    !Local
    integer :: i, is

    !Integrate over velocity at each position to obtain rho
    rho1=0.
    do is=1,nspec
       do i=1,nx
          rho1(i)=rho1(i)-spec(is)%Q_s*sum(fs1(is,i,:))*dv1
       enddo
    enddo

  end subroutine get_rho
!------------------------------------------------------------------------------
!                           Greg Howes, 2015
!------------------------------------------------------------------------------
!Compute Spatial Derivative
!------------------------------------------------------------------------------
  subroutine ddx(g,dx1,dgdx)
    implicit none
    real, intent(in), dimension(1:nx) ::g             !Function
    real, intent(in) :: dx1                         !Spacing dx
    real, intent(out), dimension(1:nx) :: dgdx      !Derivative
    integer :: k
    
    !Compute dg/dx assuming periodicity
    do k=2,nx-1
       dgdx(k)=(g(k+1)-g(k-1))/(2.*dx1)
    enddo
    !Assume periodicity to get endpoints
    dgdx(1)=(g(2)-g(nx))/(2.*dx1)
    dgdx(nx)=(g(1)-g(nx-1))/(2.*dx1)
    
  end subroutine ddx
!------------------------------------------------------------------------------
!Compute Velocity Derivative
!------------------------------------------------------------------------------
  subroutine ddv(g,dv1,dgdv)
    implicit none
    real, intent(in), dimension(-nv/2:nv/2) ::g
    real, intent(in) :: dv1
    real, intent(out), dimension(-nv/2:nv/2) :: dgdv
    integer :: k

    !Compute dg/dv, using uncentered, 1st-order derivatives at endpoints
   do k=-nv/2+1,nv/2-1
        dgdv(k)=(g(k+1)-g(k-1))/(2.*dv1)
   enddo
   dgdv(-nv/2)=(g(-nv/2+1)-g(-nv/2))/(dv1)
   dgdv(nv/2)=(g(nv/2)-g(nv/2-1))/(dv1)
  end subroutine ddv
!------------------------------------------------------------------------------
!3rd order Adams-Bashforth Method
!------------------------------------------------------------------------------
  subroutine adams_bashforth3(fe3,dt1,fe2,dfedt2,dfedt1,dfedt0)
    implicit none
    real, intent(out) :: fe3
    real, intent(in) :: dt1
    real, intent(in) :: fe2
    real, intent(in) :: dfedt2
    real, intent(in) :: dfedt1
    real, intent(in) :: dfedt0

    fe3=fe2+dt1/12.*(23.*dfedt2-16.*dfedt1+5.*dfedt0)
    
  end subroutine adams_bashforth3
!------------------------------------------------------------------------------
!Split odd/even components of f(v)
!------------------------------------------------------------------------------
  subroutine fparity(g,godd,geven)
    implicit none
    real, intent(in), dimension(-nv/2:nv/2) ::g
    real, intent(out), dimension(-nv/2:nv/2) :: godd
    real, intent(out), dimension(-nv/2:nv/2) :: geven
    integer :: k

    !Compute dg/dv, using uncentered, 1st-order derivatives at endpoints
   do k=-nv/2,nv/2
      godd(k)=(g(k)-g(-k))/2.
      geven(k)=(g(k)+g(-k))/2.
   enddo
 end subroutine fparity
!------------------------------------------------------------------------------
!Calculate Species Ballistic Term Change, dfsb/dt
!------------------------------------------------------------------------------
  real function calc_dfsbdt(is,v,dfsbdx1,dfswdx1,dfscdx1)
    use vp_data, only: spec
    implicit none
    integer :: is
    real :: v
    real :: dfsbdx1
    real :: dfswdx1
    real :: dfscdx1

    calc_dfsbdt=spec(is)%sqtm*(-v*dfsbdx1-v*dfswdx1-v*dfscdx1)

  end function calc_dfsbdt
  
  
!------------------------------------------------------------------------------
!Calculate Linear Wave Term Change, dfswl/dt
!------------------------------------------------------------------------------
  real function calc_dfswldt(is,dphidx1,dfs0dv1)
    use vp_data, only: spec
    implicit none
    integer :: is
    real :: dphidx1
    real :: dfs0dv1

    calc_dfswldt=dphidx1*dfs0dv1*spec(is)%Qinvtm

  end function calc_dfswldt

!------------------------------------------------------------------------------
!Calculate Nonlinear Wave Ballistic Term Change, dfswnb/dt
!------------------------------------------------------------------------------
  real function calc_dfswnbdt(is,dphidx1,dfsbdv1)
    use vp_data, only: spec
    implicit none
    integer :: is
    real :: dphidx1
    real :: dfsbdv1

    calc_dfswnbdt=dphidx1*dfsbdv1*spec(is)%Qinvtm

  end function calc_dfswnbdt
  
!------------------------------------------------------------------------------
!Calculate Nonlinear Wave Wave Term Change, dfswnw/dt
!------------------------------------------------------------------------------
  real function calc_dfswnwdt(is,dphidx1,dfswdv1)
    use vp_data, only: spec
    implicit none
    integer :: is
    real :: dphidx1
    real :: dfswdv1

    calc_dfswnwdt=dphidx1*dfswdv1*spec(is)%Qinvtm

  end function calc_dfswnwdt
!------------------------------------------------------------------------------
!Calculate Nonlinear Wave Coll Term Change, dfswnc/dt
!------------------------------------------------------------------------------
  real function calc_dfswncdt(is,dphidx1,dfscdv1)
    use vp_data, only: spec
    implicit none
    integer :: is
    real :: dphidx1
    real :: dfscdv1

    calc_dfswncdt=dphidx1*dfscdv1*spec(is)%Qinvtm

  end function calc_dfswncdt
!------------------------------------------------------------------------------
!Calculate Species Collision Term Change, dfsc/dt
!------------------------------------------------------------------------------
  real function calc_dfscdt(is,j1,fsb1,fsw1,fsc1)
    !NEEDS TO BE UPDATED!
    use vp_data, only: nuspec
    implicit none
    integer :: is
    integer :: j1
    real :: fsb1
    real :: fsw1
    real :: fsc1

    calc_dfscdt=-1.*nuspec(is,j1)*(fsb1+fsw1+fsc1)

  end function calc_dfscdt
!------------------------------------------------------------------------------
!                           Kristopher G Klein, 2016
! Routine for opening output files. Moved from vp.f90 to this subroutine to 
! ease readability of code.
!------------------------------------------------------------------------------
 subroutine open_units
   use vp_data, only: runname,full_out,full_unit
   use vp_data, only: q_unit,e_unit,p_unit, ql_out,num_unit
   use vp_data, only: nx_out,f_unit,phi_unit
   use vp_data, only: moving_probes,nProbe,nspec
   use vp_data, only: probe_unit,probe_f_unit,probe_phi_unit
   use vp_data, only: xout, jE_unit
   use vp_data, only: phase_out,phase_unit, current_unit, current_out
   use vp_data, only: get_unused_unit
   implicit none

   !Local
   integer :: im, ip, ispec
   character(150) :: writeName
   
  !Output files have option of full output of f
  ! as well as an arbitrary number of slices at prescribed x values  
  if (full_out) then
     do ispec = 1, nspec
        write(writeName,'(2a,i0)')trim(runname),'.vp_f_full.s',ispec
        write(*,'(2a)')' => ',trim(writeName)
        call get_unused_unit(full_unit(ispec))
        open(full_unit(ispec),file=trim(writeName),status='replace')
     enddo
  endif

   do im=1,nx_out
      do ispec = 1,nspec
         write(writeName,'(2a,i0,a,i0)')trim(runname),'.vp_f.',xout(im),'.s',ispec
         write(*,'(2a)')' =>',trim(writeName)
         call get_unused_unit(f_unit(ispec,im))
         open(f_unit(ispec,im),file=trim(writeName),status='replace')
      enddo

      write(writeName,'(2a,i0)')trim(runname),'.vp_phi.',xout(im)
      write(*,'(2a)')' =>',trim(writeName)
      call get_unused_unit(phi_unit(im))
      open(phi_unit(im),file=trim(writeName),status='replace')
   enddo
   
   if (moving_probes) then
      do ip = 1, nProbe
         do ispec = 1, nspec
            write(writeName,'(2a,i0,a,i0)')trim(runname),'.vp_f_probe.',ip,'.s',ispec
            write(*,'(2a)')' =>',trim(writeName)
            call get_unused_unit(probe_f_unit(ispec,ip))
            open(probe_f_unit(ispec,ip),file=trim(writeName),status='replace')
         enddo
         
         write(writeName,'(2a,i0)')trim(runname),'.vp_phi_probe.',ip
         write(*,'(2a)')' =>',trim(writeName)
         call get_unused_unit(probe_phi_unit(ip))
         open(probe_phi_unit(ip),file=trim(writeName),status='replace')
      enddo
      write(writeName,'(2a)')trim(runname),'.probe_position'
      write(*,'(2a)')' =>',trim(writeName)
      call get_unused_unit(probe_unit)
      open(probe_unit,file=trim(writeName),status='replace')
   endif

   if (phase_out) then
      do ispec = 1,nspec
         write(writeName,'(2a,i0)')trim(runname),'.phase_space.',ispec
         write(*,'(2a)')' =>',trim(writeName)
         call get_unused_unit(phase_unit(ispec))
         open(phase_unit(ispec),file=trim(writeName),status='replace')
      enddo
   endif

   if (current_out) then
         write(writeName,'(2a)')trim(runname),'.current'
         write(*,'(2a)')' =>',trim(writeName)
         call get_unused_unit(current_unit)
         open(current_unit,file=trim(writeName),status='replace')

         write(writeName,'(2a)')trim(runname),'.jE'
         write(*,'(2a)')' =>',trim(writeName)
         call get_unused_unit(jE_unit)
         open(jE_unit,file=trim(writeName),status='replace')
   endif

   write(writeName,'(2a)')trim(runname),'.vp_phi'
   write(*,'(2a)')' =>',trim(writeName)
   call get_unused_unit(p_unit)
   open(p_unit,file=trim(writeName),status='replace')

   write(writeName,'(2a)')trim(runname),'.vp_energy'
   write(*,'(2a)')' =>',trim(writeName)
   call get_unused_unit(e_unit)
   open(e_unit,file=trim(writeName),status='replace')

   write(writeName,'(2a)')trim(runname),'.vp_number'
   write(*,'(2a)')' =>',trim(writeName)
   call get_unused_unit(num_unit)
   open(num_unit,file=trim(writeName),status='replace')

   if (ql_out) then
      do ispec=1,nspec
         write(writeName,'(2a,i0)')trim(runname),'.vp_ql.s',ispec
         write(*,'(2a)')' =>',trim(writeName)
         call get_unused_unit(q_unit(ispec))
         open(q_unit(ispec),file=trim(writeName),status='replace')
      enddo
   endif

 end subroutine open_units

!------------------------------------------------------------------------------
!                           Kristopher G Klein, 2016
! Routine for periodic distribution and field output. 
!Moved from vp.f90 to this subroutine to ease readability of code.
!------------------------------------------------------------------------------

 subroutine periodic_output
   use vp_data, only: t, x, v, nv, nx, dv, dx
   use vp_data, only: full_out,full_unit, ql_out
   use vp_data, only: q_unit,p_unit, spec
   use vp_data, only: nx_out,f_unit,phi_unit, jE_unit
   use vp_data, only: probe_unit,probe_f_unit,probe_phi_unit
   use vp_data, only: nspec, xout, current_unit, current_out
   use vp_data, only: fs,fs0,fsb,fsbodd,fsbeven,fsw,fswodd,fsweven
   use vp_data, only: fsc,fscodd,fsceven
   use vp_data, only: fswl,fswnb,fswnw,fswnc,phi,rho,dphidx, js, jE_dx
   use vp_data, only: qlfs,qlfs0,qlfsb,qlfsw,qlfswl,qlfswnb,qlfswnw,qlfswnc,qlfsc
   implicit none

   integer :: ip, ispec
   integer :: im, j, i
   character(50) :: current_fmt, jE_fmt

   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   !DISTRIBUTION FUNCTION OUTPUT
   !write(*,'(a)')'VDF DIAGNOSTIC'
   do im = 1,nx_out
      j = xout(im)
      do ispec = 1, nspec
         call fparity(fsb(ispec,j,:),fsbodd(:),fsbeven(:))
         call fparity(fsw(ispec,j,:),fswodd(:),fsweven(:))
         call fparity(fsc(ispec,j,:),fscodd(:),fsceven(:))
         do i=-nv/2,nv/2
            write(f_unit(ispec,im),'(3f12.4,15es15.7)') &
                 t,x(j),v(i),&
                 fs(ispec,j,i),fs0(ispec,j,i),&
                 fsb(ispec,j,i),fsw(ispec,j,i),fsc(ispec,j,i),&
                 fsbodd(i),fsbeven(i),fswodd(i),fsweven(i),fscodd(i),fsceven(i),&
                 fswl(ispec,j,i),fswnb(ispec,j,i),fswnw(ispec,j,i),fswnc(ispec,j,i)
         enddo
      enddo
   enddo
   
   if (full_out) then
      !write(*,'(a)')'FULL VDF DIAGNOSTIC'
      do j = 1,nx
         do ispec = 1,nspec
            call fparity(fsb(ispec,j,:),fsbodd(:),fsbeven(:))
            call fparity(fsw(ispec,j,:),fswodd(:),fsweven(:))
            call fparity(fsc(ispec,j,:),fscodd(:),fsceven(:))
            do i=-nv/2,nv/2
               write(full_unit(ispec),'(3f12.4,15es15.7)') &
                    t,x(j),v(i),&
                    fs(ispec,j,i),fs0(ispec,j,i),&
                    fsb(ispec,j,i),fsw(ispec,j,i),fsc(ispec,j,i),&
                    fsbodd(i),fsbeven(i),fswodd(i),fsweven(i),fscodd(i),fsceven(i),&
                    fswl(ispec,j,i),fswnb(ispec,j,i),fswnw(ispec,j,i),fswnc(ispec,j,i)
            enddo
            write(full_unit(ispec),*)
         enddo
      enddo
   endif

   do ispec=1,nspec
      do im = 1,nx_out
         write(f_unit(ispec,im),*)      
         !write(f_unit(ispec,im),*)      
      enddo
      if (full_out) then
         write(full_unit(ispec),*)
      endif
   enddo

   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   !FIELDS OUTPUT
   !write(*,'(a)')'FIELDS DIAGNOSTIC'
   do im = 1,nx_out
      j = xout(im)
      write(phi_unit(im),'(2f12.4,3es15.7)') &
           t,x(j),phi(j),rho(j),dphidx(j)
   enddo
   do j=1,nx
      write(p_unit,'(2f12.4,3es15.7)') &
           t,x(j),phi(j),rho(j),dphidx(j)
   enddo
   write(p_unit,*)

   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   !CURRENT OUTPUT
   !write(*,'(a)')'CURRENT DIAGNOSTIC'
   if (current_out) then
      write(current_fmt,'(a,i0,a)') '(3es14.4,',nspec,'es14.4)'
      write(jE_fmt,'(a,i0,a)') '(es14.4,',nspec,'es14.4)'
      jE_dx=0.
      do j=1,nx
         js=0.
         do ispec=1,nspec
            do i=-nv/2,nv/2
               js(ispec)=js(ispec)+dv*v(i)*spec(ispec)%Q_s*fs(ispec,j,i)
            enddo
            jE_dx(ispec)=jE_dx(ispec)-dphidx(j)*js(ispec)*dx
         enddo
         write(current_unit,current_fmt)&
              t,x(j),-dphidx(j),js(1:nspec)         
      enddo
      write(jE_unit,jE_fmt)t,jE_dx(1:nspec)
      write(current_unit,*);write(current_unit,*)
   endif
   

   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   !QUASILINEAR DIAGNOSTIC
   !write(*,'(a)')'QL DIAGNOSTIC'
   if (ql_out) then
      do ispec=1,nspec      
         do j=-nv/2,nv/2
            qlfs(j)=sum(fs(ispec,:,j))/real(nx)
            qlfs0(j)=sum(fs0(ispec,:,j))/real(nx)
            qlfsb(j)=sum(fsb(ispec,:,j))/real(nx)
            qlfsw(j)=sum(fsw(ispec,:,j))/real(nx)
            qlfswl(j)=sum(fswl(ispec,:,j))/real(nx)
            qlfswnb(j)=sum(fswnb(ispec,:,j))/real(nx)
            qlfswnw(j)=sum(fswnw(ispec,:,j))/real(nx)
            qlfswnc(j)=sum(fswnc(ispec,:,j))/real(nx)
            qlfsc(j)=sum(fsc(ispec,:,j))/real(nx)
         enddo
         do j=-nv/2,nv/2
            write(q_unit(ispec),'(2f12.4,8es15.7)')&
                 t,v(j),qlfs(j),qlfs0(j),qlfsb(j),qlfsw(j),&
                 qlfswl(j),qlfswnb(j),qlfswnw(j),qlfswnc(j),qlfsc(j)
         enddo
         write(q_unit(ispec),*)
      enddo
   endif
   !write(*,'(a)')'DIAGNOSTICS COMPLETE'
 end subroutine periodic_output

!------------------------------------------------------------------------------
!                           Kristopher G Klein, 2016
!Routine for outputing system energies
 !----------------------------------------------------------------
 subroutine probe_output
   use vp_data, only: dt,ntout,x,dx,nspec,nx,t,x,v,phi,rho,dphidx
   use vp_data, only: fs,fs0,fsb,fsw,fsbodd,fsbeven,fswodd,fsweven,fswl,fswnb,fswnw
   use vp_data, only: fswnc
   use vp_data, only: fsc,fscodd,fsceven
   use vp_data, only: nProbe,probe_unit,probe_f_unit,probe_phi_unit
   use vp_data, only: xProbe, vProbe, nProbe
   use vp_data, only: fs_inter,fs0_inter,fsb_inter,fsw_inter,fsc_inter
   use vp_data, only: fswl_inter,fswnb_inter,fswnw_inter,fswnc_inter
   implicit none
   integer :: ip,ispec,j,i,il,inter,im
   real :: xtmp
   character(50) :: probe_fmt
   real :: inter_step
   real :: phi_inter, rho_inter, dphidx_inter

   write(probe_fmt,'(a,i0,a)') '(',nProbe+1,'es14.4)'

   !Update moving probe positions
   do ip = 1, nProbe
      xtmp= xProbe(ip) + dt*real(ntout)*vProbe(ip)               
      if (xtmp.le.(x(1)-dx)) then                  
                  xtmp = xtmp + (2.*x(nx))
               elseif (xtmp.gt.x(nx)) then
                  xtmp = xtmp - (2.*x(nx))
               endif
               xProbe(ip) = xtmp
               !!!linear interpolation
               !find nearest neighbors
               do il = 1,nx-1
                  inter = 0
                  if (xProbe(ip).eq.x(il)) then
                     !No need to interpolate                     
                     j = il
                     write(probe_phi_unit(ip),'(i12,2f12.4,3es15.7)') &
                          ip,t,x(j),phi(j),rho(j),dphidx(j)
                     do ispec=1,nspec
                        call fparity(fsb(ispec,j,:),fsbodd(:),fsbeven(:))
                        call fparity(fsw(ispec,j,:),fswodd(:),fsweven(:))
                        call fparity(fsc(ispec,j,:),fscodd(:),fsceven(:))
                        do i=-nv/2,nv/2
                           write(probe_f_unit(ispec,ip),'(3f12.4,15es15.7)') &
                                t,x(j),v(i),&
                                fs(ispec,j,i),fs0(ispec,j,i),&
                                fsb(ispec,j,i),fsw(ispec,j,i),fsc(ispec,j,i),&
                                fsbodd(i),fsbeven(i),fswodd(i),fsweven(i), &
                                fscodd(i),fsceven(i),&
                                fswl(ispec,j,i),fswnb(ispec,j,i),fswnw(ispec,j,i),fswnc(ispec,j,i)
                        enddo
                        write(probe_f_unit(ispec,ip),*)
                     enddo
                  elseif(xProbe(ip).eq.x(il+1)) then
                     !No need to interpolate
                     j = il+1
                     write(probe_phi_unit(ip),'(i12,2f12.4,3es15.7)') &
                          ip,t,x(j),phi(j),rho(j),dphidx(j)
                     do ispec=1,nspec
                        call fparity(fsb(ispec,j,:),fsbodd(:),fsbeven(:))
                        call fparity(fsw(ispec,j,:),fswodd(:),fsweven(:))
                        call fparity(fsc(ispec,j,:),fscodd(:),fsceven(:))
                        do i=-nv/2,nv/2
                           write(probe_f_unit(ispec,ip),'(3f12.4,15es15.7)') &
                                t,x(j),v(i),&
                                fs(ispec,j,i),fs0(ispec,j,i),&
                                fsb(ispec,j,i),fsw(ispec,j,i),fsc(ispec,j,i),&
                                fsbodd(i),fsbeven(i),fswodd(i),fsweven(i),&
                                fscodd(i),fsceven(i),&
                                fswl(ispec,j,i),fswnb(ispec,j,i),fswnw(ispec,j,i),fswnc(ispec,j,i)
                        enddo
                        write(probe_f_unit(ispec,ip),*)
                     enddo
                  elseif ((xProbe(ip).gt.x(il)).and.(xProbe(ip).lt.x(il+1))) then
                     !save index for interpolation
                     inter = il
                     exit 
                  elseif (xProbe(ip).lt.x(1)) then
                     inter = -1
                  endif
               enddo
               !if interpolation necessary
               !add in periodic interpolation...

               if (inter.gt.0) then
                  !interpolate fields
                  inter_step = (xProbe(ip) - x(inter))/(x(inter+1)-x(inter))
                  phi_inter = phi(inter) + &
                       (phi(inter+1)-phi(inter))*inter_step
                  rho_inter = rho(inter) + &
                       (rho(inter+1)-rho(inter))*inter_step
                  dphidx_inter = dphidx(inter) + &
                       (dphidx(inter+1)-dphidx(inter))*inter_step

                  write(probe_phi_unit(ip),'(i12,2f12.4,3es15.7)') &
                       ip,t,xProbe(ip),phi_inter,rho_inter,dphidx_inter

                  !interpolate distributions
                  do i = -nv/2,nv/2
                     do ispec=1,nspec
                        fs_inter(ispec,i)=fs(ispec,inter,i)   + (fs(ispec,inter+1,i) -fs(ispec,inter,i)) *inter_step                     
                        fs0_inter(ispec,i)=fs0(ispec,inter,i) + (fs0(ispec,inter+1,i)-fs0(ispec,inter,i))*inter_step                      
                        fsb_inter(ispec,i)=fsb(ispec,inter,i) + (fsb(ispec,inter+1,i)-fsb(ispec,inter,i))*inter_step                     
                        fsw_inter(ispec,i)=fsw(ispec,inter,i) + (fsw(ispec,inter+1,i)-fsw(ispec,inter,i))*inter_step        
                        fswl_inter(ispec,i)= fswl(ispec,inter,i) +  (fswl(ispec,inter+1,i)- fswl(ispec,inter,i))*inter_step 
                        fswnb_inter(ispec,i)=fswnb(ispec,inter,i) + (fswnb(ispec,inter+1,i)-fswnb(ispec,inter,i))*inter_step
                        fswnw_inter(ispec,i)=fswnw(ispec,inter,i) + (fswnw(ispec,inter+1,i)-fswnw(ispec,inter,i))*inter_step
                        fswnc_inter(ispec,i)=fswnc(ispec,inter,i) + (fswnc(ispec,inter+1,i)-fswnc(ispec,inter,i))*inter_step
                        fsc_inter(ispec,i)=fsc(ispec,inter,i) + (fsc(ispec,inter+1,i)-fsc(ispec,inter,i))*inter_step                     
                     enddo
                  enddo

                  do ispec=1,nspec
                     call fparity(fsb_inter(ispec,:),fsbodd(:),fsbeven(:))
                     call fparity(fsw_inter(ispec,:),fswodd(:),fsweven(:))
                     call fparity(fsc_inter(ispec,:),fscodd(:),fsceven(:))
                     !HERE?
                     do i = -nv/2,nv/2
                        write(probe_f_unit(ispec,ip),'(3f12.4,15es15.7)') &
                             t,xProbe(iP),v(i),&
                             fs_inter(ispec,i),fs0_inter(ispec,i),&
                             fsb_inter(ispec,i),fsw_inter(ispec,i),fsc_inter(ispec,i),&
                             fsbodd(i),fsbeven(i),fswodd(i),fsweven(i),&
                             fscodd(i),fsceven(i),&
                             fswl_inter(ispec,i),fswnb_inter(ispec,i),fswnw_inter(ispec,i),fswnc_inter(ispec,i)
                     enddo
                     write(probe_f_unit(ispec,ip),*)
                  enddo




               !Handle periodic Boundaries
               elseif (inter==-1) then
                  !interpolate fields
                  inter_step = (xProbe(ip) - x(1)+dx)/(dx)
                  phi_inter = phi(nx) + (phi(1)-phi(nx))*inter_step                       
                  rho_inter = rho(nx) + (rho(1)-rho(nx))*inter_step                       
                  dphidx_inter = dphidx(nx) + (dphidx(1)-dphidx(nx))*inter_step
                  
                  write(probe_phi_unit(ip),'(i12,2f12.4,3es15.7)') &
                       ip,t,xProbe(ip),phi_inter,rho_inter,dphidx_inter

                  !interpolate distributions
                  do i = -nv/2,nv/2
                     do ispec=1,nspec
                        fs_inter(ispec,i)=fs(ispec,nx,i)   + (fs(ispec,1,i) -fs(ispec,nx,i)) *inter_step                     
                        fs0_inter(ispec,i)=fs0(ispec,nx,i) + (fs0(ispec,1,i)-fs0(ispec,nx,i))*inter_step                      
                        fsb_inter(ispec,i)=fsb(ispec,nx,i) + (fsb(ispec,1,i)-fsb(ispec,nx,i))*inter_step                     
                        fsw_inter(ispec,i)=fsw(ispec,nx,i) + (fsw(ispec,1,i)-fsw(ispec,nx,i))*inter_step        
                        fswl_inter(ispec,i)= fswl(ispec,nx,i) +  (fswl(ispec,1,i)- fswl(ispec,nx,i))*inter_step 
                        fswnb_inter(ispec,i)=fswnb(ispec,nx,i) + (fswnb(ispec,1,i)-fswnb(ispec,nx,i))*inter_step
                        fswnw_inter(ispec,i)=fswnw(ispec,nx,i) + (fswnw(ispec,1,i)-fswnw(ispec,nx,i))*inter_step
                        fswnc_inter(ispec,i)=fswnc(ispec,nx,i) + (fswnc(ispec,1,i)-fswnc(ispec,nx,i))*inter_step
                        fsc_inter(ispec,i)=fsc(ispec,nx,i) + (fsc(ispec,1,i)-fsc(ispec,nx,i))*inter_step                     
                     enddo
                  enddo

                  do ispec=1,nspec
                     call fparity(fsb_inter(ispec,:),fsbodd(:),fsbeven(:))
                     call fparity(fsw_inter(ispec,:),fswodd(:),fsweven(:))
                     call fparity(fsc_inter(ispec,:),fscodd(:),fsceven(:))
                     do i = -nv/2,nv/2
                        write(probe_f_unit(ispec,ip),'(3f12.4,15es15.7)') &
                             t,xProbe(iP),v(i),&
                             fs_inter(ispec,i),fs0_inter(ispec,i),&
                             fsb_inter(ispec,i),fsw_inter(ispec,i),fsc_inter(ispec,i),&
                             fsbodd(i),fsbeven(i),fswodd(i),fsweven(i),&
                             fscodd(i),fsceven(i),&
                             fswl_inter(ispec,i),fswnb_inter(ispec,i),fswnw_inter(ispec,i),fswnc_inter(ispec,i)
                     enddo
                     write(probe_f_unit(ispec,ip),*)
                  enddo
               endif

               

            enddo
            !write probe positions.
            write(probe_unit,trim(probe_fmt)) t,xProbe

 end subroutine probe_output

!------------------------------------------------------------------------------
!                          Gregory G. Howes, 2018
!Routine for outputing Intergated system particle number by species
!------------------------------------------------------------------------------
subroutine number_output(jt)
  use vp_data, only: num_unit,nspec,ns,ns0,nsb,nsw,nsc,t
  implicit none
  !PASSED
  integer :: jt
  !LOCAL
  character(50) :: num_fmt 
  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !NUMBER OUTPUT AND DIAGNOSTICS

  call number_cons(ns,ns0,nsb,nsw,nsc)

  write(num_fmt,'(a,i0,a)')'(i6,1f12.4,',5*nspec,'es15.7)'
  write(num_unit,num_fmt) &
       jt,t,ns(1:nspec),ns0(1:nspec),nsb(1:nspec),nsw(1:nspec),nsc(1:nspec)

  
  
end subroutine number_output
!------------------------------------------------------------------------------
!                           Kristopher G Klein, 2016
!Routine for outputing system energies
!----------------------------------------------------------------
subroutine energy_output(jt)
  use vp_data, only:e_unit,t,fs,fs0,fsb,fsw,fswl,fswnb,fswnw,fswnc,fsc,dfs0dv,decs,intdecs,qcs,v,dv,dphidx,dx,etot,ephi,efs,efs0,efsb,efsw,efswl,efswnb,efswnw,efswnc,efsc,edecs,eintdecs,eqcs,ekos,ekos0,ekos1,ekos2,nspec
  implicit none
  !PASSED
  integer :: jt
  !LOCAL
  character(50) :: e_fmt 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !ENERGY OUTPUT AND DIAGNOSTICS
  
  call energy_full_ie(fs,fs0,fsb,fsw,fswl,fswnb,fswnw,fswnc,fsc,dfs0dv,decs,intdecs,qcs,v,dv,dphidx,dx,etot,ephi,efs,efs0,efsb,efsw,efswl,efswnb,efswnw,efswnc,efsc,edecs,eintdecs,eqcs,ekos,ekos0,ekos1,ekos2)

  write(e_fmt,'(a,i0,a)')'(i6,1f12.4,2es15.7,',16*nspec,'es15.7)'
  write(e_unit,e_fmt) &
       jt,t,etot,ephi,efs(1:nspec),&
       efs0(1:nspec),efsb(1:nspec),efsw(1:nspec),efsc(1:nspec),&
       efswl(1:nspec),efswnb(1:nspec),efswnw(1:nspec),efswnc(1:nspec),&
       edecs(1:nspec),eintdecs(1:nspec),eqcs(1:nspec), &
       ekos(1:nspec),ekos0(1:nspec),ekos1(1:nspec),ekos2(1:nspec) 
 end subroutine energy_output

!------------------------------------------------------------------------------
!                           Kristopher G Klein, 2016
! Routine for outputting phase space plots output files.
!------------------------------------------------------------------------------
 subroutine phase_output
   use vp_data, only: nx,nv,phase_unit,t,x,v,fs,fs0,fsb,fsw,nspec,fsc,coll,decs,intdecs,qcs
   implicit none
   integer :: i,j,ispec
   do ispec = 1,nspec
      do j = 1,nx
         do i = -nv/2,nv/2
            if (coll) then
               write(phase_unit(ispec),'(11es14.4)')&
                    t,x(j),v(i),fs(ispec,j,i),fs0(ispec,j,i),fsb(ispec,j,i),fsw(ispec,j,i),fsc(ispec,j,i),decs(ispec,j,i),intdecs(ispec,j,i),qcs(ispec,j,i)
               else
                  write(phase_unit(ispec),'(8es14.4)')&
                       t,x(j),v(i),fs(ispec,j,i),fs0(ispec,j,i),fsb(ispec,j,i),fsw(ispec,j,i),fsc(ispec,j,i)
               endif
         enddo
         write(phase_unit(ispec),*)
      enddo
      write(phase_unit(ispec),*)
   enddo
 end subroutine phase_output

!------------------------------------------------------------------------------
!                           Kristopher G Klein, 2016
! Routine for closing output files. Moved from vp.f90 to this subroutine to 
! ease readability of code.
!------------------------------------------------------------------------------
 subroutine close_units
   use vp_data, only: full_out,full_unit, jE_unit
   use vp_data, only: q_unit,e_unit,p_unit, ql_out,num_unit
   use vp_data, only: nx_out,f_unit,phi_unit
   use vp_data, only: moving_probes,nProbe,nspec
   use vp_data, only: probe_unit,probe_f_unit,probe_phi_unit
   use vp_data, only: phase_out, phase_unit, current_unit, current_out
   implicit none   

   integer :: im, ip, ispec
  
   close(e_unit); close(p_unit); close(num_unit)

   do im=1,nx_out      
      do ispec=1,nspec
         close(f_unit(ispec,im))
      enddo
      close(phi_unit(im))
   enddo

   do ispec=1,nspec
      if (ql_out) &
           close(q_unit(ispec))
      if (full_out) &
           close(full_unit(ispec))
      if (phase_out) &
           close(phase_unit(ispec))
   enddo


   if (moving_probes) then
      do ip=1,nProbe
         do ispec=1,nspec
            close(probe_f_unit(ispec,ip))
         enddo
         close(probe_phi_unit(ip))
      enddo
      close(probe_unit)
   endif

   if (current_out) then
      close(current_unit)
      close(jE_unit)
   endif
 end subroutine close_units
!------------------------------------------------------------------------------
end module vp_funcs
 
