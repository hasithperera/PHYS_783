!=============================================================================!
!=============================================================================!
!       Data      for   VP: 1D1V Vlasov-Poisson Nonlinear Solver WAVES:       !
!=============================================================================!
!=============================================================================!
!
!------------------------------------------------------------------------------
!                             Copyright, 2015
!                                Greg Howes
!
!------------------------------------------------------------------------------
!  
!------------------------------------------------------------------------------
module vp_data
  implicit none
  !Parameters
  !KGK---These are default values. 
  !   ---the subroutine 'read_in_params' allows run-time parameter input
  !   --- via appending a *.in file after the executable
  integer :: nspec=2 !number of plasma species
  integer :: nx, nv            !nx and nv values
  real :: vmax                    !Max of velocity grid
  real :: nl                       !Size of domain: L=2*pi*nl*lambda_e
  logical :: nonlin=.true.             !Include NL terms
  logical :: coll=.false.              !Include Collisions
  logical :: vdepnu=.false.            !Velocity Dependent Collision Operator
  real :: mu=100.                      !mi/me
  real :: tau=1.0                       !Ti/Te
  real :: cfl=0.0625                    !CFL fraction
  integer :: nt                  !Number of time steps
  integer :: ntout                 !Number of steps per save
  integer :: ntout_phase                 !Number of steps per save
  !Initial Condition Parameters
  integer :: nk1                   !Integral wavenumber w.r.t. domain
  !kle = real(nk1)/nl
  real :: dn=0.02                       !Initial density perturbation
  integer :: ic=1                      !ICs: 1=Density, 2=eigenfunction
  integer :: eig_opt=0                 !ICs: 1=Langmuir, 2=Ion Acoustic, 0=Specify
! OMEGA: Ion Acoustic kle=1.0, tau=0.1, mu=16
!  complex, parameter :: w_ic=cmplx(0.244,-0.0291)!IC: Eigenfunction freq guess
! OMEGA: Ion Acoustic kle=0.5, tau=0.1, mu=16
!  complex, parameter :: w_ic=cmplx(1.39E-01,-1.60E-02)!IC: Eigenfunction freq guess
! OMEGA: Ion Acoustic kle=0.5, tau=0.1, mu=100
!  complex, parameter :: w_ic=cmplx(5.55E-02,-3.56E-03)!IC: Eigenfunction freq guess
!  OMEGA: LANGMUIR kle=0.75, tau=1, mu=100
!  complex, parameter :: w_ic=cmplx(-1.74,-0.465)!IC: Eigenfunction freq guess
!*  OMEGA: LANGMUIR kle=0.5, tau=1, mu=100
!  complex, parameter :: w_ic=cmplx(1.41685,-0.151758)!IC: Eigenfunction freq guess
!  OMEGA: LANGMUIR kle=0.3, tau=1, mu=100
!  complex, parameter :: w_ic=cmplx(1.16,-0.0122)!IC: Eigenfunction freq guess
!  OMEGA: LANGMUIR kle=0.2, tau=1, mu=100
  complex :: w_ic=cmplx(1.06,-4.99E-05)!IC: Eigenfunction freq guess
!  complex, parameter :: w_ic=cmplx(0.06,-0.00353)!IC: Eigenfunction freq guess
  real :: phi0=0.1                   !Initial Wave Potential Amplitude
  


  !Variables- Distribution Functions
  real,  dimension(:,:,:), allocatable :: fs        !Full Dist Func
  real,  dimension(:,:,:), allocatable :: fs0       !Initial equilibrium
  real,  dimension(:,:,:), allocatable :: fsb       !Delta fs ballistic
  real,  dimension(:,:,:), allocatable :: fsw       !Delta fs wave, total
  real,  dimension(:,:,:), allocatable :: fswl      !Delta fs wave, linear
  real,  dimension(:,:,:), allocatable :: fswnb     !Delta fs wave, nl ballistic
  real,  dimension(:,:,:), allocatable :: fswnw     !Delta fs wave, nl wave
  real,  dimension(:,:,:), allocatable :: fswnc     !Delta fs wave, nl coll
  real,  dimension(:,:,:), allocatable :: fsc       !Delta fs collisional
  real,  dimension(:,:,:,:), allocatable :: dfsbdt   !dfs/dt (for 3 timesteps)
  real,  dimension(:,:,:,:), allocatable :: dfswdt   !dfs/dt (for 3 timesteps)
  real,  dimension(:,:,:,:), allocatable :: dfswldt   !dfs/dt (for 3 timesteps)
  real,  dimension(:,:,:,:), allocatable :: dfswnbdt   !dfs/dt (for 3 timesteps)
  real,  dimension(:,:,:,:), allocatable :: dfswnwdt   !dfs/dt (for 3 timesteps)
  real,  dimension(:,:,:,:), allocatable :: dfswncdt   !dfs/dt (for 3 timesteps)
  real,  dimension(:,:,:,:), allocatable :: dfscdt   !dfs/dt (for 3 timesteps)
  real,  dimension(:,:,:), allocatable :: fsbnew     !New ballistic delta fs
  real,  dimension(:,:,:), allocatable :: fswnew     !New wave delta fs
  real,  dimension(:,:,:), allocatable :: fswlnew     !New wave delta fs
  real,  dimension(:,:,:), allocatable :: fswnbnew     !New wave delta fs
  real,  dimension(:,:,:), allocatable :: fswnwnew     !New wave delta fs
  real,  dimension(:,:,:), allocatable :: fswncnew     !New wave delta fs
  real,  dimension(:,:,:), allocatable :: fscnew     !New collisional delta fs
  real,  dimension(:,:,:), allocatable :: dfsbdx      !d/dx(Ball Dist Func)
  real,  dimension(:,:,:), allocatable :: dfsbdv      !d/dv(Ball Dist Func)
  real,  dimension(:,:,:), allocatable :: dfswdx      !d/dx(Wave Dist Func)
  real,  dimension(:,:,:), allocatable :: dfswdv      !d/dv(Wave Dist Func)
  real,  dimension(:,:,:), allocatable :: dfs0dv      !d/dv(Equilib Dist Func)
  real,  dimension(:,:,:), allocatable :: dfscdx      !d/dx(Coll Dist Func)
  real,  dimension(:,:,:), allocatable :: dfscdv      !d/dv(Coll Dist Func)
  real,  dimension(:),allocatable :: dfsbdt_tmp,dfswdt_tmp,dfswldt_tmp,dfswnbdt_tmp,dfswnwdt_tmp,dfscdt_tmp,dfswncdt_tmp
  real,  dimension(:,:), allocatable :: nuspec   !Velocity-dependent collision operator for each species

  !Variables-potential and charge density
  real,  dimension(:), allocatable :: phi                   !Potential
  real,  dimension(:), allocatable :: dphidx                !Potential Gradient
  real,  dimension(:), allocatable :: d2phidx2              !Potential 2nd Deriv
  real,  dimension(:), allocatable :: rho                   !Charge Density
  real :: dXrho                   !Sum_x Charge Density KGK
  !Grid variables
  real,  dimension(:), allocatable  :: x                    !Position grid
  real,  dimension(:), allocatable :: v               !Velocity grid (species dependent normalization)
  real :: t=0.                                    !Time
  real :: dt                                      !Time Step Size
  real :: pi
  real :: dx                                       !Spatial interval
  real :: dv                                       !Velocity-space interval
  !Numerical Implementation Parameters
  real :: sqtm                                     !Sqrt(tau/mu)
  real :: sqinvtm                                  !1./Sqrt(tau*mu)
  ! Initial conditions variables
  real :: k1                                       !Wavevector
  !Diagnostics
  real  ::  etot                                   !Total Energy
  real :: ephi                                     !Total Field Energy
  real,  dimension(:), allocatable :: efs        !Total Dist Fn Energy
  real,  dimension(:), allocatable :: efs0       !Equil Dist Fn Energy
  real,  dimension(:), allocatable :: efsb       !Ballistic Dist Fn Energy
  real,  dimension(:), allocatable :: efsw       !Wave Dist Fn Energy
  real,  dimension(:), allocatable :: efswl      !Wave Dist Fn Energy
  real,  dimension(:), allocatable :: efswnb     !Wave Dist Fn Energy
  real,  dimension(:), allocatable :: efswnw     !Wave Dist Fn Energy
  real,  dimension(:), allocatable :: efswnc     !Wave Dist Fn Energy
  real,  dimension(:), allocatable :: efsc       !Collisional Dist Fn Energy

  real,  dimension(:,:,:), allocatable :: decs     !Delta E_coll_s (x,v)
  real,  dimension(:,:,:), allocatable :: intdecs  !Int dt Delta E_coll_s (x,v)
  real,  dimension(:,:,:), allocatable :: qcs      !Delta E_coll_s (x,v)/Delta t
  real,  dimension(:), allocatable :: edecs        ! E_coll_s energy
  real,  dimension(:), allocatable :: eintdecs     !Int  E_coll_s energy
  real,  dimension(:), allocatable :: eqcs         !d E_coll_s (x,v)/dt energy rate

  real,  dimension(:), allocatable :: ekos        !Kruskal-Oberman Energy,e
  real,  dimension(:), allocatable :: ekos0       !K-O Energy,e, Equilib
  real,  dimension(:), allocatable :: ekos1       !K-O Energy,e, Linear
  real,  dimension(:), allocatable :: ekos2       !K-O Energy,e, Quadratic
  real,  dimension(:), allocatable :: fsbeven,fsbodd !Delta fs ballistic Odd/even
  real,  dimension(:), allocatable :: fsweven,fswodd  !Delta fs wave Odd/even
  real,  dimension(:), allocatable :: fsceven,fscodd !Delta fs collisional Odd/even

  real, dimension(:), allocatable :: ns        !Total species particle number
  real, dimension(:), allocatable :: ns0       !Total Equilib species particle number
  real, dimension(:), allocatable :: nsb       !Total delta number: ballistic
  real, dimension(:), allocatable :: nsw       !Total delta number: wave
  real, dimension(:), allocatable :: nsc       !Total delta number: collisional
  
  

  real,  dimension(:), allocatable :: qlfs              !QL (x-avg) fs
  real,  dimension(:), allocatable :: qlfs0             !QL (x-avg) fs0
  real,  dimension(:), allocatable :: qlfsb             !QL (x-avg) fsb
  real,  dimension(:), allocatable :: qlfsw             !QL (x-avg) fsw
  real,  dimension(:), allocatable :: qlfswl            !QL (x-avg) fsw
  real,  dimension(:), allocatable :: qlfswnb           !QL (x-avg) fsw
  real,  dimension(:), allocatable :: qlfswnw           !QL (x-avg) fsw
  real,  dimension(:), allocatable :: qlfswnc           !QL (x-avg) fsw
  real,  dimension(:), allocatable :: qlfsc             !QL (x-avg) fsc

  real,  dimension(:), allocatable :: js !Current from species s [q_s*v*f_s]
  real,  dimension(:), allocatable :: jE_dx !Volume integrated js.E

  !Vlasov-Poission Parameters
  real :: kle                                             !k lambda_dRef

  !SPECIES PARAMETERS
  public :: specie
     type :: specie
        real :: tau_s     !T_s/T_R
        real :: mu_s      !m_s/m_R
        real :: Q_s       !q_s/q_R
        real :: D_s       !n_s/n_R
        real :: V_s       !V_s/v_tR
        real :: nu_s      !nu_s/nu_R
        logical :: P_s    !Turns on (T) or off (F) initial perturbations
        real :: L2        !D_s*Q_s^2/tau_s
        real :: sqtm      !sqrt(tau_s/mu_s)
        real :: Qinvtm   !Q_s/sqrt(tau_s*mu_s)
     end type specie
  type (specie), dimension (:), allocatable, target :: spec 


  logical :: CN_smooth=.true.
  integer :: CN_N=15

  !interpolated probe distribution
  real,  dimension(:,:), allocatable :: fs_inter !Full Dist Func
  real,  dimension(:,:), allocatable :: fs0_inter   !Initial equilibrium
  real,  dimension(:,:), allocatable :: fsb_inter   !Delta fs ballistic
  real,  dimension(:,:), allocatable :: fsw_inter   !Delta fs wave
  real,  dimension(:,:), allocatable :: fswl_inter  !Delta fs wave, linear
  real,  dimension(:,:), allocatable :: fswnb_inter !Delta fs wave, nl ballistic
  real,  dimension(:,:), allocatable :: fswnw_inter !Delta fs wave, nl wave
  real,  dimension(:,:), allocatable :: fswnc_inter !Delta fs wave, nl wave
  real,  dimension(:,:), allocatable :: fsc_inter   !Delta fs collisional

  !I/O values for namelist
  integer :: unit
  integer, parameter :: stdout_unit=6
  integer, save :: input_unit_no, error_unit_no=stdout_unit
  character(50) :: runname   !string for parameter input file
  integer :: nx_out=3 !number of slices in x to output
  logical :: full_out=.false. !Flag for outputting entire (x,v,t) distribution function
  logical :: ql_out=.false. !Flag for outputting ql distribution function
  integer, dimension(:), allocatable :: xout     !indices associated with nx_out slices
  integer, dimension(:), allocatable :: phi_unit !output units for phi slices
  integer, dimension(:,:), allocatable :: f_unit   !output units for f slices
  integer, dimension(:), allocatable :: full_unit  !output units for full distribution
  integer, dimension(:), allocatable :: q_unit  !output units for quasilinear distribution
  integer, dimension(:), allocatable :: phase_unit  !output units for phase space plots
  logical :: phase_out=.false.
  logical :: current_out=.false. !output or supress periodic output of the current
  integer :: e_unit, p_unit, current_unit, jE_unit, num_unit
  logical :: moving_probes=.false. !Flag for outputting distributions and fields moving along
  integer :: probe_unit
  !trajectories in x(t)
  integer :: nP_x=1 !number of initial positions for probes
  integer :: nP_v=1 !number of velocities for probes
  integer :: nProbe=1 !total number of probes
  real :: Px_min, Px_max !min and max of initial probe positions, in units of 1/lamdba_e  
  real :: Pv_min, Pv_max !min and max of initial probe positions, in units of 1/lamdba_e  
  real, dimension(:), allocatable :: xProbe !current position of each probe
  real, dimension(:), allocatable :: vProbe !velocity of each probe
  integer, dimension(:), allocatable :: probe_phi_unit   !output units for moving probe fields
  integer, dimension(:,:), allocatable :: probe_f_unit     !output units for moving probe vdf
  
  !PARAMETERS
  public :: nx,nv,nl,nonlin,vmax,cfl
  public :: coll,nuspec,spec,kle
  public :: nt,ntout,ntout_phase,phase_out, ql_out
  public :: nspec,CN_smooth,CN_N
  !VARIABLES
  public :: fs,fs0,fsb,fsw,fsc,dfsbdt,dfswdt,dfscdt,fsbnew,fswnew,fscnew 
  public :: fswl,fswnb,fswnw,dfswldt,dfswnbdt,dfswnwdt,fswlnew,fswnbnew,fswnwnew
  public :: fswnc,dfswncdt,fswncnew
  public :: dfsbdx,dfsbdv,dfswdx,dfswdv,dfs0dv,dfscdx,dfscdv
  public :: phi,dphidx,d2phidx2,rho,dXrho,js, jE_dx
  public :: x,dx,v,dv
  public :: t,dt
  public :: pi,sqtm,sqinvtm
  !INITIAL CONDITIONS
  public :: nk1,k1,dn,ic,eig_opt,w_ic,phi0
  !ENERGY DIAGNOSTICS
  public :: etot,ephi,efs,efs0,efsb,efsw,efsc
  public :: efswl,efswnb,efswnw,efswnc
  public :: decs,intdecs,qcs
  public :: ekos,ekos0,ekos1,ekos2
  public :: fsbeven,fsbodd,fsweven,fswodd,fsceven,fscodd 
  public :: qlfs,qlfs0,qlfsb,qlfsw,qlfsc
  public :: qlfswl,qlfswnb,qlfswnw,qlfswnc
  !I/O VARIABLES
  public :: runname, get_unused_unit
  public :: nx_out, full_out, xout 
  public :: phi_unit, f_unit, probe_f_unit, probe_phi_unit
  public :: current_unit, current_out, jE_unit
  public :: full_unit, p_unit,e_unit,q_unit,probe_unit,phase_unit,num_unit
  public :: moving_probes, nP_x, nP_v, Px_min, Px_max, Pv_min, Pv_max
  public :: fs_inter,fs0_inter,fsb_inter,fsw_inter,fswl_inter,fswnb_inter,fswnw_inter
  public :: fsc_inter,fswnc_inter
contains
!------------------------------------------------------------------------------
!                           Greg Howes, 2015
!------------------------------------------------------------------------------
! Initialize Variables
!------------------------------------------------------------------------------
  subroutine initialize_variables
    implicit none
    integer :: i,ispec
    integer :: ix,iv

    real :: dv_probe, dx_probe

    call read_in_params
    call allocate_arrays

    !Calculate pi
    pi=4.0*atan(1.0)

    !SET velocity grid
    do i=-nv/2,nv/2
       v(i)=vmax *real(i)/real(nv/2)
    enddo
    dv=vmax/real(nv/2)
    !SET spatial grid
    dv=vmax/real(nv/2)
    do i=1,nx
       x(i)=(real(i)/real(nx)*2.-1.)*pi*nl
    enddo
    dx=2.*pi*nl/real(nx)
    !Set Krook collision operator coefficient (possibly velocity dependent)
    if (coll) then
       do ispec=1,nspec
          if (vdepnu) then
             write(*,'(a)') "WARNING: Velocity dependent collision operator not yet implemented"
          else
             nuspec(ispec,:)=spec(ispec)%nu_s
          endif
       enddo
    endif

    if (moving_probes) then
       !insure that all probes are within simulation box
       if (Px_min.le.x(1)) Px_min = x(1)
       if (Px_max.gt.x(nx)) Px_max = x(nx)
       !set initial value of probe position
       dv_probe = (Pv_max-Pv_min)/(real(nP_v-1))
       dx_probe = (Px_max-Px_min)/(real(nP_x-1))
       i = 1
       do iv=1,nP_V; do ix=1,nP_x          
          xProbe(i) = Px_min + real(ix-1)*dx_probe
          vProbe(i) = Pv_min + real(iv-1)*dv_probe
          i = i + 1
       enddo; enddo
    endif

 end subroutine initialize_variables
!------------------------------------------------------------------------------
!K.G. Klein- 2016
  subroutine allocate_arrays
    implicit none

 !Dist Fns
    allocate(fs(nspec,1:nx,-nv/2:nv/2));fs=0.
    allocate(fs0(nspec,1:nx,-nv/2:nv/2));fs0=0.
    allocate(fsb(nspec,1:nx,-nv/2:nv/2));fsb=0.
    allocate(fsw(nspec,1:nx,-nv/2:nv/2));fsw=0.
    allocate(fswl(nspec,1:nx,-nv/2:nv/2));fswl=0.
    allocate(fswnb(nspec,1:nx,-nv/2:nv/2));fswnb=0.
    allocate(fswnw(nspec,1:nx,-nv/2:nv/2));fswnw=0.
    allocate(fswnc(nspec,1:nx,-nv/2:nv/2));fswnc=0.
    allocate(fsc(nspec,1:nx,-nv/2:nv/2));fsc=0.
    allocate(dfsbdt(nspec,1:nx,-nv/2:nv/2,3));dfsbdt=0.
    allocate(dfswdt(nspec,1:nx,-nv/2:nv/2,3));dfswdt=0.
    allocate(dfswldt(nspec,1:nx,-nv/2:nv/2,3));dfswldt=0.
    allocate(dfswnbdt(nspec,1:nx,-nv/2:nv/2,3));dfswnbdt=0.
    allocate(dfswnwdt(nspec,1:nx,-nv/2:nv/2,3));dfswnwdt=0.
    allocate(dfswncdt(nspec,1:nx,-nv/2:nv/2,3));dfswncdt=0.
    allocate(dfscdt(nspec,1:nx,-nv/2:nv/2,3));dfscdt=0.
    allocate(fsbnew(nspec,1:nx,-nv/2:nv/2));fsbnew=0.
    allocate(fswnew(nspec,1:nx,-nv/2:nv/2));fswnew=0.
    allocate(fswlnew(nspec,1:nx,-nv/2:nv/2));fswlnew=0.
    allocate(fswnbnew(nspec,1:nx,-nv/2:nv/2));fswnbnew=0.
    allocate(fswnwnew(nspec,1:nx,-nv/2:nv/2));fswnwnew=0.
    allocate(fswncnew(nspec,1:nx,-nv/2:nv/2));fswncnew=0.
    allocate(fscnew(nspec,1:nx,-nv/2:nv/2));fscnew=0.
    allocate(dfsbdx(nspec,1:nx,-nv/2:nv/2));dfsbdx=0.
    allocate(dfsbdv(nspec,1:nx,-nv/2:nv/2));dfsbdv=0.
    allocate(dfswdx(nspec,1:nx,-nv/2:nv/2));dfswdx=0.
    allocate(dfswdv(nspec,1:nx,-nv/2:nv/2));dfswdv=0.
    allocate(dfs0dv(nspec,1:nx,-nv/2:nv/2));dfs0dv=0.
    allocate(dfscdx(nspec,1:nx,-nv/2:nv/2));dfscdx=0.
    allocate(dfscdv(nspec,1:nx,-nv/2:nv/2));dfscdv=0.

    allocate(dfsbdt_tmp(nspec))
    allocate(dfswdt_tmp(nspec))
    allocate(dfswldt_tmp(nspec))
    allocate(dfswnbdt_tmp(nspec))
    allocate(dfswnwdt_tmp(nspec))
    allocate(dfswncdt_tmp(nspec))
    allocate(dfscdt_tmp(nspec))

    !NOTE: THESE VARIABLES NEED TO BE INITIALIZED BECAUSE THEY ARE PASSED
!    if (coll) then
       allocate(nuspec(nspec,-nv/2:nv/2)); nuspec=0.
       allocate(decs(nspec,1:nx,-nv/2:nv/2)); decs=0.
       allocate(intdecs(nspec,1:nx,-nv/2:nv/2)); intdecs=0.
       allocate(qcs(nspec,1:nx,-nv/2:nv/2)); qcs=0.
       allocate(edecs(nspec)); edecs=0.
       allocate(eintdecs(nspec)); eintdecs=0.
       allocate(eqcs(nspec)); eqcs=0.
!    endif
       
    allocate(efs(nspec));  efs=0.
    allocate(efs0(nspec)); efs0=0.
    allocate(efsb(nspec)); efsb=0.
    allocate(efsw(nspec)); efsw=0.
    allocate(efswl(nspec)); efswl=0.
    allocate(efswnb(nspec)); efswnb=0.
    allocate(efswnw(nspec)); efswnw=0.
    allocate(efswnc(nspec)); efswnc=0.
    allocate(efsc(nspec)); efsc=0.

    allocate(ekos(nspec));   ekos=0.
    allocate(ekos0(nspec));  ekos0=0.
    allocate(ekos1(nspec));  ekos1=0.
    allocate(ekos2(nspec));  ekos2=0.

    allocate(ns(nspec));  ns=0.
    allocate(ns0(nspec));  ns0=0.
    allocate(nsb(nspec));  nsb=0.
    allocate(nsw(nspec));  nsw=0.
    allocate(nsc(nspec));  nsc=0.


    allocate(js(nspec)); js=0.
    allocate(jE_dx(nspec)); jE_dx=0.

    !Fields
    allocate(phi(1:nx));phi=0.
    allocate(dphidx(1:nx));dphidx=0.
    allocate(d2phidx2(1:nx));d2phidx2=0.
    allocate(rho(1:nx));rho=0.
    !Grid
    allocate(x(1:nx)); x=0. !inverse reference Debye length
    allocate(v(-nv/2:nv/2));v=0.
    !Probes
    if (moving_probes) then
       allocate(xProbe(1:nProbe)); xProbe=0.
       allocate(vProbe(1:nProbe)); vProbe=0.
    endif
    !Diagnostics
    if (moving_probes) then
       allocate(fs_inter(nspec,-nv/2:nv/2));fs_inter=0.
       allocate(fs0_inter(nspec,-nv/2:nv/2));fs0_inter=0.
       allocate(fsb_inter(nspec,-nv/2:nv/2));fsb_inter=0.
       allocate(fsw_inter(nspec,-nv/2:nv/2));fsw_inter=0.
       allocate(fswl_inter(nspec,-nv/2:nv/2));fswl_inter=0.
       allocate(fswnb_inter(nspec,-nv/2:nv/2));fswnb_inter=0.
       allocate(fswnw_inter(nspec,-nv/2:nv/2));fswnw_inter=0.
       allocate(fswnc_inter(nspec,-nv/2:nv/2));fswnc_inter=0.
       allocate(fsc_inter(nspec,-nv/2:nv/2));fsc_inter=0.
    endif

    allocate(fsbeven(-nv/2:nv/2));fsbeven=0.
    allocate(fsbodd(-nv/2:nv/2));fsbodd=0.
    allocate(fsweven(-nv/2:nv/2));fsweven=0.
    allocate(fswodd(-nv/2:nv/2));fswodd=0.
    allocate(fsceven(-nv/2:nv/2));fsceven=0.
    allocate(fscodd(-nv/2:nv/2));fscodd=0.

    if (ql_out) then
       allocate(qlfs(-nv/2:nv/2));qlfs=0.
       allocate(qlfs0(-nv/2:nv/2));qlfs0=0.
       allocate(qlfsb(-nv/2:nv/2));qlfsb=0.
       allocate(qlfsw(-nv/2:nv/2));qlfsw=0.
       allocate(qlfswl(-nv/2:nv/2));qlfswl=0.
       allocate(qlfswnb(-nv/2:nv/2));qlfswnb=0.
       allocate(qlfswnw(-nv/2:nv/2));qlfswnw=0
       allocate(qlfswnc(-nv/2:nv/2));qlfswnc=0
       allocate(qlfsc(-nv/2:nv/2));qlfsc=0.
    endif

  end subroutine allocate_arrays

!------------------------------------------------------------------------------
!K.G. Klein- 2016
  subroutine deallocate_arrays
    implicit none

    !Dist Fns
    deallocate(fs,fs0,fsb,fsw,fswl,fswnb,fswnw,fswnc,fsc)
    deallocate(dfsbdt,dfswdt,dfswldt,dfswnbdt,dfswnwdt,dfswncdt,dfscdt)
    deallocate(fsbnew,fswnew,fswlnew,fswnbnew,fswnwnew,fswncnew,fscnew)
    deallocate(dfsbdx,dfsbdv,dfswdx,dfswdv,dfs0dv,dfscdx,dfscdv)
    deallocate(dfsbdt_tmp,dfswdt_tmp,dfswldt_tmp,dfswnbdt_tmp,dfswnwdt_tmp,dfswncdt_tmp,dfscdt_tmp)

    !Fields
    deallocate(phi,dphidx,d2phidx2,rho)

    !Energy
    deallocate(efs,efs0,efsb,efsw,efswl,efswnb,efswnw,efswnc,efsc)
    deallocate(decs,intdecs,qcs,edecs,eintdecs,eqcs)
    deallocate(ekos,ekos0,ekos1,ekos2)
    deallocate(js,jE_dx)

    !Grid
    deallocate(x,v)
    if (moving_probes) then
       deallocate(xProbe,vProbe)
       deallocate(fs_inter,fs0_inter,fsb_inter,fsw_inter,fswl_inter,fswnb_inter,fswnw_inter,fswnc_inter,fsc_inter)
    endif
    !Diagnostics
    deallocate(fsbeven,fsbodd,fsweven,fswodd,fsceven,fscodd)
    if (ql_out) &
         deallocate(qlfs,qlfs0,qlfsb,qlfsw,qlfswl,qlfswnb,qlfswnw,qlfswnc,qlfsc)
    deallocate(q_unit,phi_unit,f_unit)
    if (full_out) deallocate(full_unit)
    if (phase_out) deallocate(phase_unit)
    deallocate(spec)
  end subroutine deallocate_arrays

!------------------------------------------------------------------------------
!K.G. Klein- 2016
  subroutine read_in_params
    !Read in system parameters
    !input file is argument after executable:
    !$ ./vp.e test.in
    implicit none
    integer :: l
    integer :: ispec !species index
    integer :: im !slide index
    integer :: x_in !input indices for x
    real :: om, gam !read in values for complex frequency guess
    real :: tauS,muS,Qs,Ds,Vs,nuS
    logical :: perturbS

    nameList /params/ &
         nspec,nx,nv,vmax,nl,nonlin,coll,vdepnu,&
         cfl,nt,ntout,ntout_phase,phase_out,CN_smooth,CN_N,&
         nk1,dn,ic,eig_opt,om,gam,phi0,&
         full_out,ql_out,nx_out,moving_probes,current_out

     nameList /probes/  nP_x,nP_v,Px_min,Px_max,Pv_min,Pv_max

     nameList /species/ tauS,muS,Qs,Ds,Vs,perturbS,nuS
        
     nameList /slice/ x_in
     
    call get_unused_unit (input_unit_no)
    call get_runname(runname)
    runname=trim(runname)//".in"
    unit=input_unit_no
    open (unit=unit,file=runname,status='old',action='read')
    read (unit=unit,nml=params)
    
    allocate(phi_unit(1:nx_out));phi_unit=0
    allocate(f_unit(1:nspec,1:nx_out));  f_unit=0
    if (full_out) then
       allocate(full_unit(1:nspec));  full_unit=0
    endif
    allocate(q_unit(1:nspec));  q_unit=0
    allocate(xout(1:nx_out));    xout=1
    if (phase_out) then
       allocate(phase_unit(1:nspec));  phase_unit=0
    endif

    w_ic=cmplx(om,gam)

    if (moving_probes) then
!       nameList /probes/ &
!            nP_x,nP_v,Px_min,Px_max,Pv_min,Pv_max
       read (unit=unit,nml=probes)
    endif

    nProbe = nP_x*nP_v
    allocate(probe_phi_unit(1:nProbe));probe_phi_unit=0
    allocate(probe_f_unit(1:nspec,1:nProbe));  probe_f_unit=0

    
    allocate(spec(nspec))
    do ispec =1,nspec
       call get_indexed_namelist_unit (unit, "species", ispec)
!       nameList /species/ &
!            tauS,muS,Qs,Ds,Vs,perturbS,nuS
       read (unit=unit,nml=species)
       spec(ispec)%tau_s=tauS
       spec(ispec)%mu_s=muS
       spec(ispec)%Q_s=Qs
       spec(ispec)%D_s=Ds
       spec(ispec)%V_s=Vs
       spec(ispec)%P_s=perturbS
       spec(ispec)%nu_s=nuS
       
       spec(ispec)%L2=Ds*Qs**2./tauS
       spec(ispec)%sqtm=sqrt(tauS/muS)
       spec(ispec)%Qinvtm=Qs/sqrt(tauS*muS)

       write(*,'(a,i4)')    'Species:  ', ispec
       write(*,'(a,es11.3)')'T_s/T_R:  ', spec(ispec)%tau_s
       write(*,'(a,es11.3)')'m_s/m_R:  ', spec(ispec)%mu_s
       write(*,'(a,es11.3)')'q_s/q_R:  ', spec(ispec)%q_s
       write(*,'(a,es11.3)')'n_s/n_R:  ', spec(ispec)%D_s
       write(*,'(a,es11.3)')'V_s/v_ts: ', spec(ispec)%V_s
       write(*,'(a,es11.3)')'nu_s/nu_R: ', spec(ispec)%nu_s

       !GGH Add close(unit) to try to avoid crash in gfortran
       close(unit)
       
    enddo
    

    if (nx_out.ge.1) then
       do im = 1,nx_out
          call get_indexed_namelist_unit (unit, "slice", im)
!          nameList /slice/ &
!               x_in
          read (unit=unit,nml=slice)
          xout(im)=x_in
          !GGH Add close(unit) to try to avoid crash in gfortran
          close(unit)
        enddo
   endif

!    close(unit)

    !remove *.in from runname for data output
    l = len_trim (runname)
    runname = runname(1:l-3)

  end subroutine read_in_params

!---------------------------------------------------------------
! Get runname for output files from input argument
  subroutine get_runname(runname)
    implicit none
    integer       :: l
    character(50) :: arg
    character(50), intent(out) :: runname

    !Get the first argument of the program execution command
    call getarg(1,arg)

    !Check if this is the input file and trim .in extension to get runname
    l = len_trim (arg)
    if (l > 3 .and. arg(l-2:l) == ".in") then
       runname = arg(1:l-3)
    end if
  end subroutine get_runname
!------------------------------------------------------------------------------

!-=-=-=-=-=-
!The following routines:
!    get_indexed_namelist_unit
!    input_unit_exist
!    get_unused_unit
!    input_unit
!were all adopted from the Astrophysical Gyrokinetic Code (AGK)
!as a means of allowing arbitrary namelist group name input.
!A bit of hassle, but worth the effort.
!-=-=-=-=-=-
  subroutine get_indexed_namelist_unit (unit, nml, index_in)
    implicit none
    integer, intent (out) :: unit
    character (*), intent (in) :: nml
    integer, intent (in) :: index_in
    character(500) :: line
    integer :: iunit, iostat, in_file
    integer :: ind_slash
    logical :: exist

    call get_unused_unit (unit)
    ind_slash=index(runname,"/",.True.)
    if (ind_slash.EQ.0) then !No slash in name
       !Original behaviour
       write(*,*) runname," ."//trim(runname)//".scratch"
       open (unit=unit, file="."//trim(runname)//".scratch")
    else
        !General behaviour
        open (unit=unit, file=trim(runname(1:ind_slash))//"."//trim(runname(ind_slash+1:))//".scratch")
    endif

    write (line, *) index_in
    line = nml//"_"//trim(adjustl(line))
    in_file = input_unit_exist(trim(line), exist)

    if (exist) then
       iunit = input_unit(trim(line))
    else
       write(6,*) "get_indexed_namelist: following namelist not found ",trim(line)
       return
    end if

    read (unit=iunit, fmt="(a)") line
    write (unit=unit, fmt="('&',a)") nml

    do
       read (unit=iunit, fmt="(a)", iostat=iostat) line
       if (iostat /= 0 .or. trim(adjustl(line)) == "/") exit
       write (unit=unit, fmt="(a)") trim(line)
    end do
    write (unit=unit, fmt="('/')")
    rewind (unit=unit)
  end subroutine get_indexed_namelist_unit

  function input_unit_exist (nml,exist)
    implicit none
    character(*), intent (in) :: nml
    logical, intent(out) :: exist
    integer :: input_unit_exist, iostat
    character(500) :: line
    intrinsic adjustl, trim
    input_unit_exist = input_unit_no
    exist = .true.
    if (input_unit_no > 0) then
       rewind (unit=input_unit_no)
       do
          read (unit=input_unit_no, fmt="(a)", iostat=iostat) line
          if (iostat /= 0) then
             rewind (unit=input_unit_no)
             exit
          end if
          if (trim(adjustl(line)) == "&"//nml) then
             backspace (unit=input_unit_no)
             return
          end if
       end do
    end if
    exist = .false.
  end function input_unit_exist

  function input_unit (nml)
    implicit none
    character(*), intent (in) :: nml
    integer :: input_unit, iostat
    character(500) :: line
    intrinsic adjustl, trim
    input_unit = input_unit_no
    if (input_unit_no > 0) then
       rewind (unit=input_unit_no)
       do
          read (unit=input_unit_no, fmt="(a)", iostat=iostat) line
          if (iostat /= 0) then
             rewind (unit=input_unit_no)
             exit
          end if
          if (trim(adjustl(line)) == "&"//nml) then
             backspace (unit=input_unit_no)
             return
          end if
       end do
    end if
    write (unit=error_unit_no, fmt="('Couldn''t find namelist: ',a)") nml
    write (unit=*, fmt="('Couldn''t find namelist: ',a)") nml
  end function input_unit

  subroutine get_unused_unit (unit)
    implicit none
    integer, intent (out) :: unit
    logical :: od
    unit = 50
    do
       inquire (unit=unit, opened=od)
      ! write(*,*)"unit= ",unit,"   opened= ",od
       if (.not.od) return
       unit = unit + 1
    end do
  end subroutine get_unused_unit
!-=-=-=-=-=-

!-=-=-=-=

end module vp_data
