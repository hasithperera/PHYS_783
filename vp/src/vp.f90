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
! NOTES:
! 1) The Density must be initialized such that SUM_s int q_s f_s dv =0
!    This means overall neutrality, otherwise Green's function for phi fails.
! 2) Although unnecessary for the evolution, this code splits the ballistic
!    and wave contributions to the perturbed distribution function
!------------------------------------------------------------------------------
!Sept 2016- Kristopher G Klein:
!           Edited to allow for an arbitrary number of 
!           drifting velocity distributions.
!           Parameters for each 'nspec' distribution are included in 
!           'spec_#' namelists in the input file.
!------------------------------------------------------------------------------
! February 2018- Gregory Howes
!           Addition of a (velocity dependent) Krook collision operator
!              This simple Krook operator is not energy conserving, but we can
!              account for the energy change due to the operator (Incomplete)
!------------------------------------------------------------------------------
program vp
  !PARAMETERS
  use vp_data, only: nx,nv,nl,nonlin,vmax,cfl
  use vp_data, only: coll,nuspec,nspec,spec
  !VARIABLES
  use vp_data, only: fs,fs0,fsb,fsw,dfsbdt,dfswdt,fsbnew,fswnew,fsc,dfscdt,fscnew
  use vp_data, only: fswl,fswnb,fswnw,fswnc,dfswldt,dfswnbdt,dfswnwdt,dfswncdt,fswlnew,fswnbnew,fswnwnew,fswncnew
  use vp_data, only: dfsbdx,dfsbdv,dfswdx,dfswdv,dfs0dv,dfscdx,dfscdv
  use vp_data, only: phi,dphidx,d2phidx2,rho
  use vp_data, only: x,dx,v,dv
  use vp_data, only: t,dt,nt,ntout
  use vp_data, only: ntout_phase,phase_out,phase_unit
  use vp_data, only: nk1,k1,dn
  use vp_data, only: initialize_variables,pi,deallocate_arrays
  use vp_data, only: dfsbdt_tmp,dfswdt_tmp,dfswldt_tmp,dfscdt_tmp
  use vp_data, only: dfswnbdt_tmp,dfswnwdt_tmp,dfswncdt_tmp
  !FUNCTIONS
  use vp_funcs, only: ddx,ddv,initial_conditions,get_rho,get_phi
  use vp_funcs, only: adams_bashforth3
  use vp_funcs, only: calc_dfsbdt,calc_dfscdt
  use vp_funcs, only: calc_dfswldt,calc_dfswnbdt,calc_dfswnwdt,calc_dfswncdt
  use vp_funcs, only: open_units, close_units
  !ENERGY DIAGNOSTICS
  !uncertain if these are needed in the main file.
  use vp_data, only:  etot,ephi,efs,efs0,efsb,efsw,efsc
  use vp_data, only:  efswl,efswnb,efswnw,efswnc
  use vp_data, only:  decs,intdecs,qcs
  use vp_data, only:  edecs,eintdecs,eqcs
  use vp_data, only:  ekos,ekos0,ekos1,ekos2
  use vp_data, only:  fsbeven,fsbodd,fsweven,fswodd,fsceven,fscodd
  !OUTPUT FUNCTIONS
  use vp_funcs, only: periodic_output, energy_output,phase_output,probe_output
  use vp_funcs, only: number_output
  !I/O VARIABLES
  use vp_data, only: runname,full_out,get_unused_unit
  use vp_data, only: xout,nx_out,phi_unit,f_unit
  use vp_data, only: probe_phi_unit,probe_f_unit,moving_probes
  use vp_data, only: full_unit,p_unit,e_unit,q_unit,probe_unit
  use vp_data, only: current_unit


  implicit none
  real :: dt_start,dt2

  !Analysis Option
  integer :: i,j,jt,im,ip,il,ispec
  character(150) :: writeName

  !Initialize variables==============================================
  call initialize_variables
  !Compute initial conditions for distribution function==============
	write(*,*) "ahe -1"
  call initial_conditions
	write(*,*) "ahe -2"
  !Get initial charge density========================================
  call get_rho(fs,dv,rho)
	write(*,*) "ahe -3"
  !Get initial potential=============================================
  call get_phi(rho,nl,x,dx,phi)


  !CHECK INITIALIZATION==============================================
  !Compute fs derivatives
  do ispec=1,nspec
     do j=-nv/2,nv/2
        call ddx(fsw(ispec,:,j),dx,dfswdx(ispec,:,j))
        call ddx(fsb(ispec,:,j),dx,dfsbdx(ispec,:,j))
     enddo
     do j=1,nx
        call ddv(fs0(ispec,j,:),dv,dfs0dv(ispec,j,:))
        call ddv(fsw(ispec,j,:),dv,dfswdv(ispec,j,:))
        call ddv(fsb(ispec,j,:),dv,dfsbdv(ispec,j,:))
     enddo
  enddo
  call ddx(phi,dx,dphidx)
  call ddx(dphidx,dx,d2phidx2)

  !Write diagnostic output file
  if (.false.) then
     write(writeName,'(2a)')trim(runname),'.test_phi'
     open(27,file=trim(writeName),status='replace')
     do i=1,nx
        write(27,'(i6,f12.4,4es15.7)') &
             i,x(i),phi(i),dphidx(i),d2phidx2(i),rho(i)
     enddo
     close(27)     
  endif
  
  !SET TIMESTEP FOR 3RD ORDER ADAMS-BASHFORTH==============================
  dt=cfl*dx/vmax
  write (*,'(a,es15.7)') "Timestep dt=",dt
  
  !INITIALIZE DIAGNOSTICS==================================================
  write(*,'(a)')'===---===---===---===---===---==='
  write(*,'(a)')'OPEN OUTPUT UNITS'
  call open_units

  write(*,'(a)')'===---===---===---===---===---==='
  write(*,'(a)')'OUTPUT PERIODIC DIAGNOSTIC'
  call periodic_output

  if (moving_probes) then
     write(*,'(a)')'===---===---===---===---===---==='
     write(*,'(a)')'OUTPUT PROBE DIAGNOSTIC'
     call probe_output
  endif

  if (phase_out) then
     write(*,'(a)')'===---===---===---===---===---==='
     write(*,'(a)')'OUTPUT PHASE PLANE DIAGNOSTIC'
     call phase_output
  endif

  jt=0
  write(*,'(a)')'===---===---===---===---===---==='
  write(*,'(a)')'ENERGY DIAGNOSTIC'
  call energy_output(jt)
   
  write(*,'(a)')'===---===---===---===---===---==='
  write(*,'(a)')'NUMBER DIAGNOSTIC'
  call number_output(jt)
   
   !COMPUTE INITIAL TIMESTEPS==============================================
   !Compute initial timesteps to start higher order evolution--------------
  do ispec=1,nspec
     dfsbdt(ispec,:,:,:)=0.  !Initialize 
     dfswdt(ispec,:,:,:)=0.  !Initialize 
     dfswldt(ispec,:,:,:)=0.  !Initialize 
     dfswnbdt(ispec,:,:,:)=0.  !Initialize 
     dfswnwdt(ispec,:,:,:)=0.  !Initialize 
     dfswncdt(ispec,:,:,:)=0.  !Initialize 
     dfscdt(ispec,:,:,:)=0.  !Initialize 
  enddo


   write (*,'(a)') "Initializing Simulation Timesteps"
   ! First small step======================================================
   jt=1
   dt_start=dt/32.
   !Get Derivatives------------------------------------------------------
   do ispec=1,nspec
      do j=-nv/2,nv/2
         call ddx(fsw(ispec,:,j),dx,dfswdx(ispec,:,j))
         call ddx(fsb(ispec,:,j),dx,dfsbdx(ispec,:,j))
         call ddx(fsc(ispec,:,j),dx,dfscdx(ispec,:,j))
      enddo
   enddo
   call ddx(phi,dx,dphidx)
   do ispec=1,nspec
      do j=1,nx
         call ddv(fs0(ispec,j,:),dv,dfs0dv(ispec,j,:))
         call ddv(fsw(ispec,j,:),dv,dfswdv(ispec,j,:))
         call ddv(fsb(ispec,j,:),dv,dfsbdv(ispec,j,:))
         call ddv(fsc(ispec,j,:),dv,dfscdv(ispec,j,:))
      enddo
   enddo
   do ispec=1,nspec
      do i=1,nx
         do j=-nv/2,nv/2
            !First step: First-order Eulerian
            !Update Ballistic Term
            dfsbdt(ispec,i,j,1)=&
                 calc_dfsbdt(ispec,v(j),dfsbdx(ispec,i,j),dfswdx(ispec,i,j),dfscdx(ispec,i,j))
            fsbnew(ispec,i,j)=fsb(ispec,i,j)+&
                 dt_start*(dfsbdt(ispec,i,j,1))
            !Update Linear Wave Term
            dfswldt(ispec,i,j,1)=&
                 calc_dfswldt(ispec,dphidx(i),dfs0dv(ispec,i,j))
            fswlnew(ispec,i,j)=fswl(ispec,i,j)+&
                 dt_start*(dfswldt(ispec,i,j,1))
            !Update Nonlinear Wave Terms
            if (nonlin) then
               dfswnbdt(ispec,i,j,1)=&
                    calc_dfswldt(ispec,dphidx(i),dfsbdv(ispec,i,j))
               fswnbnew(ispec,i,j)=fswnb(ispec,i,j)+&
                    dt_start*(dfswnbdt(ispec,i,j,1))
               dfswnwdt(ispec,i,j,1)=&
                    calc_dfswnwdt(ispec,dphidx(i),dfswdv(ispec,i,j))
               fswnwnew(ispec,i,j)=fswnw(ispec,i,j)+&
                    dt_start*(dfswnwdt(ispec,i,j,1))
               dfswncdt(ispec,i,j,1)=&
                    calc_dfswncdt(ispec,dphidx(i),dfscdv(ispec,i,j))
               fswncnew(ispec,i,j)=fswnc(ispec,i,j)+&
                    dt_start*(dfswncdt(ispec,i,j,1))
            endif
            if (coll) then
               dfscdt(ispec,i,j,1)=calc_dfscdt(ispec,j,fsb(ispec,i,j),fsw(ispec,i,j),fsc(ispec,i,j))
               fscnew(ispec,i,j)=fsc(ispec,i,j)+&
                    dt_start*(dfscdt(ispec,i,j,1))
                !NOTE: Collisional energy change need not be calculated here.
            endif
         enddo
      enddo
   enddo
   !Accumulate energy change due to collisions over (x,v)phase space
   !NOTE: No need to accumulate yet, only 1 step
   if (coll) intdecs(:,:,:)=decs(:,:,:)



   t=t+dt_start
   !Leapfrog to first full step===========================================
   do jt=1,6
      dt2=dt_start*2.**real(jt)
      !Sum to get total distribution function
      do ispec=1,nspec
         fswnew(ispec,:,:)= fswlnew(ispec,:,:) + &
              fswnbnew(ispec,:,:) +  fswnwnew(ispec,:,:)  +  fswncnew(ispec,:,:)
         fs(ispec,:,:)=fs0(ispec,:,:)+ &
              fsbnew(ispec,:,:)+ fswnew(ispec,:,:) + fscnew(ispec,:,:)
      enddo
      !Get charge density
      call get_rho(fs,dv,rho)
      !Get potential
      call get_phi(rho,nl,x,dx,phi)
      !Get fe Derivatives------------------------------------------------------
      do ispec=1,nspec
         do j=-nv/2,nv/2
            call ddx(fswnew(ispec,:,j),dx,dfswdx(ispec,:,j))
            call ddx(fsbnew(ispec,:,j),dx,dfsbdx(ispec,:,j))
            call ddx(fscnew(ispec,:,j),dx,dfscdx(ispec,:,j))
         enddo
      enddo
      call ddx(phi,dx,dphidx)
      do ispec=1,nspec
         do j=1,nx
            call ddv(fs0(ispec,j,:),dv,dfs0dv(ispec,j,:))
            call ddv(fswnew(ispec,j,:),dv,dfswdv(ispec,j,:))
            call ddv(fsbnew(ispec,j,:),dv,dfsbdv(ispec,j,:))
            call ddv(fscnew(ispec,j,:),dv,dfscdv(ispec,j,:))
         enddo
      enddo
      do ispec=1,nspec
         do i=1,nx
            do j=-nv/2,nv/2
               !First step: Second-order Leapfrog
               !Update Ballistic Term
               dfsbdt_tmp(ispec)=&
                    calc_dfsbdt(ispec,v(j),dfsbdx(ispec,i,j),dfswdx(ispec,i,j),dfscdx(ispec,i,j))
               fsbnew(ispec,i,j)=fsb(ispec,i,j)+dt2*(dfsbdt_tmp(ispec))

               dfswldt_tmp(ispec)=&
                    calc_dfswldt(ispec,dphidx(i),dfs0dv(ispec,i,j))
               fswlnew(ispec,i,j)=fswl(ispec,i,j)+dt2*(dfswldt_tmp(ispec))
               !Update Nonlinear Wave Terms
               if (nonlin) then
                  dfswnbdt_tmp(ispec)=&
                       calc_dfswnbdt(ispec,dphidx(i),dfsbdv(ispec,i,j))
                  fswnbnew(ispec,i,j)=fswnb(ispec,i,j)+&
                       dt2*(dfswnbdt_tmp(ispec))
                  dfswnwdt_tmp(ispec)=&
                       calc_dfswnwdt(ispec,dphidx(i),dfswdv(ispec,i,j))
                  fswnwnew(ispec,i,j)=fswnw(ispec,i,j)+&
                       dt2*(dfswnwdt_tmp(ispec))
                  dfswncdt_tmp(ispec)=&
                       calc_dfswncdt(ispec,dphidx(i),dfswdv(ispec,i,j))
                  fswncnew(ispec,i,j)=fswnc(ispec,i,j)+&
                       dt2*(dfswncdt_tmp(ispec))
               endif
               if (coll) then
                  dfscdt_tmp(ispec)=calc_dfscdt(ispec,j,fsbnew(ispec,i,j),fswnew(ispec,i,j),fscnew(ispec,i,j))
                  fscnew(ispec,i,j)=fsc(ispec,i,j)+&
                       dt_start*(dfscdt_tmp(ispec))
               endif

               if (jt .eq. 6) then
                  dfsbdt(ispec,i,j,2)= dfsbdt_tmp(ispec)
                  dfswldt(ispec,i,j,2)= dfswldt_tmp(ispec)
                  dfswnbdt(ispec,i,j,2)= dfswnbdt_tmp(ispec)
                  dfswnwdt(ispec,i,j,2)= dfswnwdt_tmp(ispec)
                  dfswncdt(ispec,i,j,2)= dfswncdt_tmp(ispec)
                  dfswncdt(ispec,i,j,2)= dfswncdt_tmp(ispec)
                  dfscdt(ispec,i,j,2)= dfscdt_tmp(ispec)
               endif
            enddo
         enddo
      enddo
      !Update time
      t=dt2 !NOTE: Leapfrog means that this is the correct time after update from t=0.
      if (jt .eq. 5) then
         !Call energy diagnostics
         do ispec=1,nspec
            fswnew(ispec,:,:)= fswlnew(ispec,:,:) + &
                 fswnbnew(ispec,:,:) +  fswnwnew(ispec,:,:) +  fswncnew(ispec,:,:)
            fs(ispec,:,:)=fs0(ispec,:,:)+ &
                 fsbnew(ispec,:,:)+ fswnew(ispec,:,:)+ fscnew(ispec,:,:)
            !Compute energy change due to collisions over (x,v)phase space
            if (coll) then
               do j=-nv/2,nv/2
!                  decs(ispec,:,j)=-spec(ispec)%D_s*spec(ispec)%tau_s*v(j)**2.*nuspec(ispec,j)*(fscnew(ispec,:,j)-fsc(ispec,:,j))*dx*dv
                  decs(ispec,:,j)=-spec(ispec)%tau_s*v(j)**2.*(fscnew(ispec,:,j)-fsc(ispec,:,j))*dx*dv
               enddo
            endif
         enddo
         !Accumulate energy change due to collisions over (x,v)phase space
         if (coll) then
            !NOTE: THIS IS ENERGY CHANGE FROM t=0!
            intdecs(:,:,:)= decs(:,:,:)
            !Compute rate of energy change
            qcs(:,:,:)=decs(:,:,:)/dt2
         endif
         call energy_output(1)         
         call number_output(1)         
         write (*,'(i6,a,es15.7)') 1,"     Timestep t=",t
      elseif (jt .eq. 6) then
         !Call energy diagnostics
         do ispec=1,nspec
            fswnew(ispec,:,:)= fswlnew(ispec,:,:) + &
                 fswnbnew(ispec,:,:) +  fswnwnew(ispec,:,:) +  fswncnew(ispec,:,:)
            fs(ispec,:,:)=fs0(ispec,:,:)+ &
                 fsbnew(ispec,:,:)+ fswnew(ispec,:,:)+ fscnew(ispec,:,:)
            !Compute energy change due to collisions over (x,v)phase space
            if (coll) then
               do j=-nv/2,nv/2
!                  decs(ispec,:,j)=-spec(ispec)%D_s*spec(ispec)%tau_s*v(j)**2.*nuspec(ispec,j)*(fscnew(ispec,:,j)-fsc(ispec,:,j))*dx*dv
!                  decs(ispec,:,j)=-spec(ispec)%D_s*spec(ispec)%tau_s*v(j)**2.*(fscnew(ispec,:,j)-fsc(ispec,:,j))*dx*dv
                  ! NOTE: Don't think I need species density ratio or collision freq
                  ! Collision rate is already included in evolution of fsc
                  decs(ispec,:,j)=-spec(ispec)%tau_s*v(j)**2.*(fscnew(ispec,:,j)-fsc(ispec,:,j))*dx*dv
               enddo
            endif
         enddo
         !Accumulate energy change due to collisions over (x,v)phase space
         if (coll) then
            !NOTE: THIS IS ENERGY CHANGE FROM t=0! (No need to accumulate yet)
            intdecs(:,:,:)= decs(:,:,:)
            !Compute rate of energy change
            qcs(:,:,:)=decs(:,:,:)/dt2
         endif
         call energy_output(2)         
         call number_output(2)         
         write (*,'(i6,a,es15.7)') 2,"     Timestep t=",t
      endif
   enddo

   !Finish with update to full nstep=2
   do ispec=1,nspec
      fsb(ispec,:,:)=fsbnew(ispec,:,:)
      fsw(ispec,:,:)=fswnew(ispec,:,:)
      fswl(ispec,:,:)=fswlnew(ispec,:,:)
      fswnb(ispec,:,:)=fswnbnew(ispec,:,:)
      fswnw(ispec,:,:)=fswnwnew(ispec,:,:)
      fswnc(ispec,:,:)=fswncnew(ispec,:,:)
      fsc(ispec,:,:)=fscnew(ispec,:,:)
      fs(ispec,:,:)=fs0(ispec,:,:)+ fsb(ispec,:,:)+ fsw(ispec,:,:)+ fsc(ispec,:,:)
   enddo
   !Get charge density
   call get_rho(fs,dv,rho)
   !Get potential
   call get_phi(rho,nl,x,dx,phi)
   ! END COMPUTE INITIAL TIMESTEPS======================================

   ! MAIN LOOP: EVOLVE THE DISTIRIBUTION FUNCTION IN TIME:
   write (*,'(a)') "Starting Simulation"
   write (*,'(a)') ""
   do jt=3,nt
      if (mod(jt,100)==0) &
           write (*,'(i6,a,es15.7)') jt,"     Timestep t=",t
      !Get derivatives 
      do ispec=1,nspec
         do j=-nv/2,nv/2
            call ddx(fsw(ispec,:,j),dx,dfswdx(ispec,:,j))
            call ddx(fsb(ispec,:,j),dx,dfsbdx(ispec,:,j))
            call ddx(fsc(ispec,:,j),dx,dfscdx(ispec,:,j))
         enddo
      enddo
      call ddx(phi,dx,dphidx)
      do ispec=1,nspec
         do j=1,nx
            call ddv(fs0(ispec,j,:),dv,dfs0dv(ispec,j,:))
            call ddv(fsw(ispec,j,:),dv,dfswdv(ispec,j,:))
            call ddv(fsb(ispec,j,:),dv,dfsbdv(ispec,j,:))
            call ddv(fsc(ispec,j,:),dv,dfscdv(ispec,j,:))
         enddo
      enddo
      ! Reset collisional energy dissipation to zero
      if (coll) decs(:,:,:)=0.
      ! Update to next timestep
      do ispec=1,nspec
         do i=1,nx
            do j=-nv/2,nv/2
               !Update Ballistic Term
               dfsbdt(ispec,i,j,3)=&
                    calc_dfsbdt(ispec,v(j),dfsbdx(ispec,i,j),dfswdx(ispec,i,j),dfscdx(ispec,i,j))
               call adams_bashforth3(fsbnew(ispec,i,j),dt,fsb(ispec,i,j),dfsbdt(ispec,i,j,3),dfsbdt(ispec,i,j,2),dfsbdt(ispec,i,j,1))
               !Update Linear Wave Term
               dfswldt(ispec,i,j,3)=calc_dfswldt(ispec,dphidx(i),dfs0dv(ispec,i,j))
               call adams_bashforth3(fswlnew(ispec,i,j),dt,fswl(ispec,i,j),dfswldt(ispec,i,j,3),dfswldt(ispec,i,j,2),dfswldt(ispec,i,j,1))
               !Update Nonlinear Wave Terms
               if (nonlin) then
                  dfswnbdt(ispec,i,j,3)=calc_dfswnbdt(ispec,dphidx(i),dfsbdv(ispec,i,j))
                  call adams_bashforth3(fswnbnew(ispec,i,j),dt,fswnb(ispec,i,j),dfswnbdt(ispec,i,j,3),dfswnbdt(ispec,i,j,2),dfswnbdt(ispec,i,j,1))
                  dfswnwdt(ispec,i,j,3)=calc_dfswnwdt(ispec,dphidx(i),dfswdv(ispec,i,j))
                  call adams_bashforth3(fswnwnew(ispec,i,j),dt,fswnw(ispec,i,j),dfswnwdt(ispec,i,j,3),dfswnwdt(ispec,i,j,2),dfswnwdt(ispec,i,j,1))
                  dfswncdt(ispec,i,j,3)=calc_dfswncdt(ispec,dphidx(i),dfscdv(ispec,i,j))
                  call adams_bashforth3(fswncnew(ispec,i,j),dt,fswnc(ispec,i,j),dfswncdt(ispec,i,j,3),dfswncdt(ispec,i,j,2),dfswncdt(ispec,i,j,1))
               endif
               if (coll) then
                  dfscdt(ispec,i,j,3)=calc_dfscdt(ispec,j,fsb(ispec,i,j),fsw(ispec,i,j),fsc(ispec,i,j))
                  call adams_bashforth3(fscnew(ispec,i,j),dt,fsc(ispec,i,j),dfscdt(ispec,i,j,3),dfscdt(ispec,i,j,2),dfscdt(ispec,i,j,1))
                  !Compute energy change due to collisions over (x,v)phase space
!                  decs(ispec,i,j)=-spec(ispec)%D_s*spec(ispec)%tau_s*v(j)**2.*nuspec(ispec,j)*(fscnew(ispec,i,j)-fsc(ispec,i,j))*dx*dv
                  decs(ispec,i,j)=-spec(ispec)%tau_s*v(j)**2.*(fscnew(ispec,i,j)-fsc(ispec,i,j))*dx*dv
                  !GGH NOTE: I am not sure that the spec(ispec)%D_s factor is necessary here!
               endif
            enddo
         enddo
      enddo

      !Accumulate energy change due to collisions over (x,v)phase space
      if (coll) then
         intdecs(:,:,:)=intdecs(:,:,:) + decs(:,:,:)
         !Compute rate of energy change
         qcs(:,:,:)=decs(:,:,:)/dt
      endif
         
      do ispec=1,nspec
         !Copy new steps to fsw,fsb:      
         fsb(ispec,:,:)=fsbnew(ispec,:,:)
         fswl(ispec,:,:)=fswlnew(ispec,:,:)
         fswnb(ispec,:,:)=fswnbnew(ispec,:,:)
         fswnw(ispec,:,:)=fswnwnew(ispec,:,:)
         fswnc(ispec,:,:)=fswncnew(ispec,:,:)
         fsc(ispec,:,:)=fscnew(ispec,:,:)
         !Sum to get total distribution function
         fsw(ispec,:,:)= fswl(ispec,:,:) + fswnb(ispec,:,:) + fswnw(ispec,:,:) + fswnc(ispec,:,:) 
         fs(ispec,:,:)=fs0(ispec,:,:)+ fsb(ispec,:,:)+ fsw(ispec,:,:) + fsc(ispec,:,:)
      enddo
      !Get charge density
      call get_rho(fs,dv,rho)
      !Get potential
      call get_phi(rho,nl,x,dx,phi)
      !Update time
      t=t+dt
      !Move each dfebdt and dfewdt back 1 place
      !Ballistic term
      do ispec=1,nspec
         dfsbdt(ispec,:,:,1)=dfsbdt(ispec,:,:,2)
         dfsbdt(ispec,:,:,2)=dfsbdt(ispec,:,:,3)
         dfsbdt(ispec,:,:,3)=0.
         !Wave terms
         dfswldt(ispec,:,:,1)=dfswldt(ispec,:,:,2)
         dfswldt(ispec,:,:,2)=dfswldt(ispec,:,:,3)
         dfswldt(ispec,:,:,3)=0.
         dfswnbdt(ispec,:,:,1)=dfswnbdt(ispec,:,:,2)
         dfswnbdt(ispec,:,:,2)=dfswnbdt(ispec,:,:,3)
         dfswnbdt(ispec,:,:,3)=0.
         dfswnwdt(ispec,:,:,1)=dfswnwdt(ispec,:,:,2)
         dfswnwdt(ispec,:,:,2)=dfswnwdt(ispec,:,:,3)
         dfswnwdt(ispec,:,:,3)=0.
         dfswncdt(ispec,:,:,1)=dfswncdt(ispec,:,:,2)
         dfswncdt(ispec,:,:,2)=dfswncdt(ispec,:,:,3)
         dfswncdt(ispec,:,:,3)=0.
         if (coll) then
            dfscdt(ispec,:,:,1)=dfscdt(ispec,:,:,2)
            dfscdt(ispec,:,:,2)=dfscdt(ispec,:,:,3)
            dfscdt(ispec,:,:,3)=0.
         endif
      enddo

      !DIAGNOSTICS==========================================
      !Call energy diagnostics
      call energy_output(jt)         
      call number_output(jt)         
      !PERIODIC DIAGNOSTICS----------------------------------
      if (mod(jt,ntout) .eq. 0 ) then
         call periodic_output
         if (moving_probes) &
              call probe_output
      endif
      if ((phase_out).and.(mod(jt,ntout_phase) .eq. 0 )) then
         call phase_output
      endif           
   enddo
   !FINALIZE DIAGNOSTICS==================================================
   write(*,'(a)')'===---===---===---===---===---==='
   write(*,'(a)')'CLOSE OUTPUT UNITS'
   call close_units
   write(*,'(a)')'===---===---===---===---===---==='
   write(*,'(a)')'DEALLOCATE ARRAYS'
   call deallocate_arrays
 
!------------------------------------------------------------------------------
!contains
end program vp
 
