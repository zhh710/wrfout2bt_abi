subroutine tpause(ges_prsl,ges_tsen,geop_hgtl,tropprs,lat2,lon2,nsig)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    tpause      locate tropopause
!   prgmmr: treadon          org: np23                date: 2003-09-13
!
! abstract: locate tropopause using one of two methods.  The default
!           method uses the temperature lapse rate to identify the
!           tropopause.  Method 'pvoz' uses a combination of the 
!           potential voriticy and ozone mixing ration to locate
!           the tropopause.
!
! program history log:
!   2003-09-13 treadon - initial routine
!   2003-12-23 kleist  - generalized to use guess pressure
!   2004-06-15 treadon - reformat documentation
!   2004-07-26 treadon - remove call smooth121 (leads to different
!                           results when running code with different
!                           number of mpi tasks)
!   2004-07-27 treadon - add use only; add intent in/out
!   2005-11-29 derber  - remove psfcg and use ges_lnps instead
!   2006-02-02 treadon - rename prsl as ges_prsl
!   2006-07-28 derber  - use r1000 from constants
!                      - use sensible temperature rather than virtual
!   2006-07-31  kleist - change to ges_ps from ln(ps)
!   2008-04-03  safford - rm unused vars and uses
!   2013-10-19  todling - metguess now holds background
!                         revised how method is chosen based on guess fields
!   2021-11-07  zhh - Remove the dependency on gsi modules
!                     remove method using combination of potential vorticity (pv) and ozone.  
!
!   input argument list:
!            ges_prsl  - pressure(Pa)
!            ges_tsen  - temperature (K)
!            geop_hgtl - geopotential height at mid-layers(m)
!            lon2      - number of longitude
!            lat2      - number of latitude
!            nsig      - number of z
!
!   output argument list:
!            tropprs
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use model_precision,only:r_kind=>SP,i_kind=>INT32
  !use constants, only: rd_over_cp,grav,rad2deg,one,r1000,r0_01
  !use guess_grids, only: tropprs,geop_hgtl,&
  !     ntguessig,ges_prsl,ges_tsen
  !use gridmod, only: istart,nlat,rlats,nsig,lat2,lon2
  !use gsi_bundlemod, only : gsi_bundlegetpointer
  !use gsi_metguess_mod, only : gsxi_metguess_bundle
  implicit none

! Declare passed variables
  integer(i_kind),intent(in   ) :: lon2,lat2,nsig
  real(i_kind),dimension(lat2,lon2,nsig),intent(in   )::ges_tsen
  real(i_kind),dimension(lat2,lon2,nsig),intent(in   )::ges_prsl
  real(i_kind),dimension(lat2,lon2,nsig),intent(in   )::geop_hgtl
  real(i_kind),dimension(lat2,lon2),intent(inout)::tropprs

! Declare local parameters
  character(len=*),parameter::myname='tpause'

! Declare local variables
  logical t_method

  integer(i_kind) i,j,k,mm1

  real(r_kind):: ptrop

  real(r_kind),dimension(nsig):: prs,tdry,hgt

  real(r_kind),dimension(lat2,lon2):: trop_t

  real(r_kind),parameter::r0_01=0.01_r_kind


!================================================================================
! Set local constants
  t_method = .true.

! Locate tropopause based on temperature profile (WMO approach)
  if (t_method) then
     do j=1,lon2
        do i=1,lat2
           do k=1,nsig
              prs(k) = ges_prsl(i,j,k)
              tdry(k)= ges_tsen(i,j,k)
              hgt(k) = geop_hgtl(i,j,k)
           end do
           call tpause_t(nsig,prs,tdry,hgt,ptrop)
           trop_t(i,j) = ptrop*r0_01 !hPa
        end do
     end do

!    Load tropopause pressure (hPa) into output array (passed through module guess_grids)
     do j=1,lon2
        do i=1,lat2
           tropprs(i,j) = trop_t(i,j)
        end do
     end do


! Locate tropopause using combination of potential vorticity (pv) and ozone.  
  else

!    Compute latitudes on subdomain

!     Compute pv.  Use for locating tropopause poleward 30S/N
!     do j=1,lon2
!        do i = 1,lat2
!           psi=one/ges_ps_nt(i,j)
!           do k=1,nsig
!              prs(k) = r1000*ges_prsl(i,j,k,ntguessig)
!           end do
        
!          Compute pv         
           !do k = 2,nsig-1
           !  pm1 = prs(k-1)
           !  pp1 = prs(k+1)
           !  thetam1 = ges_tv_nt(i,j,k-1)*(r1e5/pm1)**(rd_over_cp)
           !   thetap1 = ges_tv_nt(i,j,k+1)*(r1e5/pp1)**(rd_over_cp)
           !   pv = grav*ges_vor_nt(i,j,k)*(thetam1-thetap1)/(pm1-pp1)
           !   pvort(k) = abs(pv)
           !end do
           !pvort(1) = pvort(2)
           !pvort(nsig) = pvort(nsig-1)
           
!          Locate tropopause


!          Search upward (decreasing pressure) for tropopause above sigma 0.7

!                Trop at level where pv=2e-6

               
!                Trop at level where ozone greater than 3e-7

           
!          Merge pv and ozone tropopause levels between 20 and 40 deg latitude


!    Load tropopause pressure (hPa) into output array


! End of tropopause location method blocks     
  endif

!  *** NOTE ***
!  The tropopause pressures are used to deflate the
!  moisture sensitivity vectors for satellite radiance
!  data and for IR quality control;
!  here we are setting bounds on the tropopause
!  pressure to make sure we are deflating at the very
!  minimum above 150 mb, and nowhere below 350 mb

  do j=1,lon2
     do i=1,lat2
        tropprs(i,j)=max(150.0_r_kind,min(350.0_r_kind,tropprs(i,j)))
     end do
  end do


! End of routine
  return
end subroutine tpause
