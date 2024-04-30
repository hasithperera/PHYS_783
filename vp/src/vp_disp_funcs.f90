!=============================================================================!
!=============================================================================!
!                        Dispersion Relation Functions for                    !
!               VP: 1D1V Vlasov-Poisson Nonlinear Solver WAVES:               !
!=============================================================================!
!=============================================================================!
!
!------------------------------------------------------------------------------
!                             Copyright, 2016
!                                Greg Howes
!
!------------------------------------------------------------------------------
!  
!------------------------------------------------------------------------------
module vp_disp_funcs
  private

  logical :: fail

  public :: vp_disp,rtsec
  
 contains
  
!------------------------------------------------------------------------------
!                           Greg Howes, 2016
!------------------------------------------------------------------------------
!  Calculates the Vlasov-Poisson Dispersion relation as a function of frequency om
! NORMALIZATION:  om=omega/omega_pR, v_ts=sqrt[Ts/ms] [no sqrt(2)]
! k = k lambda_dR, tau=Ts/TR, mu=ms/mR, Qs=qs/qR, vS=V_drift/v_ts
! xi_s= (om /(sqrt(2) k)) * (TR/Ts)(ms/mR)-vS, where om and k are normalized as above
   complex function vp_disp(om)
     use vp_data, only: nspec,spec,kle
     implicit none
     complex :: om                                !Complex Frequency
     complex :: zs                                !Ion/electron Plasma Dispersion
     complex :: xis                               !Plasma Dispersion arguments
     integer :: ispec

     vp_disp=0.
     
     do ispec=1,nspec
        xis=om/(sqrt(2.)*kle)*&
             sqrt(spec(ispec)%mu_s/spec(ispec)%tau_S)&
             -spec(ispec)%v_s
        zs=(zetout(xis))
        if (fail) then
           write(*,'(a,i0)')&
                'ERROR: Plasma Dispersion Function Failed:',ispec
           return
        end if
        vp_disp = vp_disp + spec(ispec)%L2*(1.+xis*zs)
     enddo
     vp_disp = vp_disp + kle**2.
     
     return 
   end function vp_disp
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!                           Greg Howes, 2010
!------------------------------------------------------------------------------
!     NOTE: This routine was adapted from f77 routine by Eliot Quataert
   complex function rtsec(func,x1,x2,xacc,iflag)
     integer, parameter :: maxit=75
     complex :: func, x1, xl, x2
     complex :: fl, f, swap, dx 
     real    :: xacc
     integer :: iflag,j


     fl=func(x1)
     f=func(x2)
     write(*,'(a,2es12.4,a,2es12.4)')'Initial  omega1= ',x1,'  disp= ',fl 
     write(*,'(a,2es12.4,a,2es12.4)')'Initial  omega2= ',x2,'  disp= ',f 
     write(19,'(a,2es12.4,a,2es12.4)')'Initial  omega1= ',x1,'  disp= ',fl 
     write(19,'(a,2es12.4,a,2es12.4)')'Initial  omega2= ',x2,'  disp= ',f 
     if(abs(fl).lt.abs(f))then
        rtsec=x1
        xl=x2
        swap=fl
        fl=f
        f=swap
     else
        xl=x1
        rtsec=x2
     endif
     do  j=1,maxit
        iflag = j
!        write(*,'(a,i4)')'Iteration ',j
        if (abs(f-fl) .gt. 1.0E-37) then
           dx=(xl-rtsec)*f/(f-fl)
	else
           dx = (x2-x1)/25.0
	end if        
        xl=rtsec
        fl=f
!	if (Real(rtsec + dx) .gt. 0.075 .or. Real(rtsec+dx) .lt. 0) then
!		rtsec = (x1 + x2)/2.0
!	else
        !NOTE: Reduce jump by 0.5 to improve convergence GH
!        rtsec=rtsec+dx
        rtsec=rtsec+dx/2.
!	end if

        f=func(rtsec)
!        write(*,'(a,i4,a,2es12.4,a,2es12.4)')'Iteration ',j,'  omega= ',rtsec,'  disp= ',f
        write(19,'(a,i4,a,2es12.4,a,2es12.4)')'Iteration ',j,'  omega= ',rtsec,'  disp= ',f 
        if(abs(dx).lt.xacc.or. Abs(f) .eq. 0.)return
     enddo
 !    stop
     return
!     pause 'rtsec exceed maximum iterations'
   end function rtsec
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!                           Greg Howes, 2010
!------------------------------------------------------------------------------
   complex function zetout(zin)
!
! This is the subroutine used by Linsker to evaluate his Z-function.
!
     complex :: z,dzetaz,term,fmult,terme,an1,bn1,zsquar
     complex :: hold,temp1,temp2,ddzeta,dddzet,zin,zetaoz
     real :: imagte,imagmu,imagse,imagsu

     fail = .false.

     error=1.e-7
     z=zin
     zsquar=z*z
     x=real(z)
     y=aimag(z)
     fn=real(zsquar)
     if (y.gt.0.) go to 99
     if (abs(fn) < 174. .and. abs(aimag(zsquar)) < 5.e4) go to 98
     if (fn.gt.0.) go to 97
     write (3,11) z
     write (*,11) z
11    format (' argument wp of subroutine zetvec has too large a negative imaginary part, wp = '/2e14.7)
 
!     zetout=(0.,0.)
!     stop
! GGH: Modification to return with a fail=.true.
     fail=.true.
     return

97   hold=(0.,0.)
     go to 99
98   hold=(0.,1.77245385090551603)*cexp(-zsquar)
99   if (x*x+y*y > 16.) go to 200
     if (abs(y) >= 1.) go to 300
     realte=-2.*x
     imagte=-2.*y
     realmu=.5*(imagte*imagte-realte*realte)
     imagmu=-imagte*realte
     realsu=realte
     imagsu=imagte
     if (x == 0. .and. y == 0.) go to 103
     fn=3.
100  realse=realte
     imagse=imagte
     realte=(realse*realmu-imagse*imagmu)/fn
     imagte=(realse*imagmu+imagse*realmu)/fn
     realse=realsu
     imagse=imagsu
     realsu=realsu+realte
     imagsu=imagsu+imagte
     fn=fn+2.
     if (abs(realse-realsu) > error .or. abs(imagse-imagsu) > error) go to 100
103  x=realsu
     fn=imagsu
     if (y > 0.) hold=(0.,1.77245385090551603)*cexp(-zsquar)
     zetaoz=cmplx(x,fn)+hold
     go to 401
200  fn=5.
     dddzet=6.
     term=dddzet
     fmult=.5/zsquar
201  terme=term
     term=term*fmult*fn*(fn-1.)/(fn-3.)
     zetaoz=term/terme
     if (abs(real(zetaoz))+abs(aimag(zetaoz)) > 1.) go to 250
     zetaoz=dddzet
     dddzet=dddzet+term
     fn=fn+2.
     if (cabs(zetaoz-dddzet) > error) go to 201
250  dddzet=dddzet/(zsquar*zsquar)
     if (y > 0.) go to 260
     fn=1.
     if (y < 0.) fn=2.
     dddzet=dddzet-4.*fn*hold*z*(2.*zsquar-3.)
260  ddzeta=-(4.+(zsquar-.5)*dddzet)/(z*(2.*zsquar-3.))
     dzetaz=(2.-z*ddzeta)/(2.*zsquar-1.)
     zetaoz=-(1.+.5*dzetaz)/z
     go to 401
300  if (y < 0.) z=conjg(z)
     terme=(1.,0.)
     term=(0.,0.)
     dzetaz=term
     fmult=terme
     n=0
     an1=z
     bn1=-z*z+.5
301  temp1=bn1*term+an1*terme
     temp2=bn1*fmult+an1*dzetaz
     zetaoz=temp1/temp2
     dzetaz=(zetaoz-term/fmult)/zetaoz
     if (abs(real(dzetaz)) < error .and. abs(aimag(dzetaz)) < error) go to 302
     bn1=bn1+2.
     n=n+1
     an1=-.5*float(n*(n+n-1))
     terme=term
     dzetaz=fmult
     term=temp1
     fmult=temp2
     if (n < 30) go to 301
302  if (y >= 0.) go to 401
     zetaoz=conjg(zetaoz)+2.*hold
401  zetout=zetaoz
9999 continue

   end function zetout
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
 end module vp_disp_funcs
