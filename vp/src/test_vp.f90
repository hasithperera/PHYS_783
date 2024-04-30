!=============================================================================!
!=============================================================================!
!!*                                 VP                                      *!!
!!                    1D1V Vlasov-Poisson Nonlinear Solver                   !!
!!                                                                           !!
!!                                Greg Howes                                 !!
!!                         gregory-howes@uiowa.edu                           !!
!!                            University of Iowa                             !! 
!!                    Department of Physics and Astronomy                    !!
!!                                                                           !!
!!                                                                           !!
!!*                                                                         *!!
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
program test_vp
  use vp_data, only: nx,nv,x,dx,v,dv,vmax,fenew,fe,dfedx,dfedv,fe0
  use vp_data, only: phi,dphidx,d2phidx2,rho
  use vp_data, only: t,dt,nt,nl,nk1,k1,dn
  use vp_data, only: initialize_variables,pi
  use vp_funcs, only: ddx,ddv,initial_conditions,get_rho,get_phi
  implicit none
  real, dimension(-nv/2:nv/2) :: fv,dfvdv
  real, dimension(1:nx) :: phi0
  
  !Local
  integer :: i,j
  real :: nn1


  !Initialize variables==============================================
  call initialize_variables
  !Compute initial conditions for distribution function

  !TEST d/dx==============================================
    k1=real(nk1)/nl

    !Integrate over velocity at each position to obtain rho
    do i=1,nx
       phi(i)=sin(2*k1*x(i))
    enddo

    call ddx(phi,dx,dphidx)
    call ddx(dphidx,dx,d2phidx2)

    !Write diagnostic output file
    open(19,file='test_dx.dat',status='replace')
    do i=1,nx
       write(19,'(i6,f12.4,3es15.7)') i,x(i),phi(i),nl*dphidx(i),nl*nl*d2phidx2(i)
    enddo
    close(19)

    !TEST d/dv==============================================
    write(*,'(a,i6)')'nv= ',nv
    do i=-nv/2,nv/2
       !Set background
!       fv(i)=exp(-v(i)**2./2.)/sqrt(pi)+ 0.5*sin(pi*v(i))*exp(-v(i)**2./2.)/sqrt(pi)
       fv(i)=exp(-v(i)**2./2.)/sqrt(pi) + 0.5*sin(pi*v(i))*exp(-v(i)**2./2.)/sqrt(pi)
    enddo
    call ddv(fv,dv,dfvdv)
    !Write diagnostic output file
    open(19,file='test_dv.dat',status='replace')
    do i=-nv/2,nv/2
       write(19,'(i6,f12.4,2es15.7)') i,v(i),fv(i),dfvdv(i)
    enddo
    close(19)

    !TEST Green's function===================================
    nn1=4.
    do i=1,nx
       rho(i)=-1.0+sin(nn1/nl*x(i)+pi/3.)+ cos(0.5*x(i)+pi/6.)-0.5*sin(x(i)/nl+pi/4.)
!       rho(i)=-1.0+sin(0.5/nl*x(i))
    enddo
    call get_phi(rho,nl,x,dx,phi)
    call ddx(phi,dx,dphidx)
    call ddx(dphidx,dx,d2phidx2)
    !Write diagnostic output file
    open(19,file='test_phi.dat',status='replace')
    do i=1,nx
       write(19,'(i6,f12.4,4es15.7)') i,x(i),phi(i),dphidx(i),d2phidx2(i),rho(i)
    enddo
    close(19)
    !TEST Green's function2===================================
    do i=1,nx
       phi0(i)=1.0+ cos(0.25*x(i)/nl+pi/3.)
    enddo
    call ddx(phi0,dx,dphidx)
    call ddx(dphidx,dx,d2phidx2)
    rho(:)=-1.*d2phidx2(:)
    call get_phi(rho,nl,x,dx,phi)

    open(19,file='test_phi2.dat',status='replace')
    do i=1,nx
       write(19,'(i6,f12.4,5es15.7)') i,x(i),phi(i),dphidx(i),d2phidx2(i),rho(i),phi0(i)
    enddo
    close(19)

    
!------------------------------------------------------------------------------
!contains
end program test_vp
 
