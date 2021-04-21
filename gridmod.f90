! Adapted from GSI gridmod.F90, 4/21/2021
module gridmod
    use model_precision,only:P,INT32
    use parameters_define,only:pi,half,deg2rad,one,zero
    use read_wrf,only::nx,ny,lat,lon
    implicit none
   ! rlambda0 : Calculate in get_xytilde_domain
   ! sign_pole: Calculate in get_xytilde_domain
   ! pihalf   : Calculate in init_general_transform
   ! atilde_x : Calculate in init_general_transform
! The following is for the generalized transform
    integer(INT32)::nlon=nx
    integer(INT32)::nlat=ny
    real(P) pihalf,sign_pole,rlambda0
    real(P) atilde_x,btilde_x,atilde_y,btilde_y
    real(P) btilde_xinv,btilde_yinv
    integer(INT32) nxtilde,nytilde
    real(P),allocatable::xtilde0(:,:),ytilde0(:,:)
    real(P),allocatable::cos_beta_ref(:,:),sin_beta_ref(:,:)
    integer(INT32),allocatable::i0_tilde(:,:),j0_tilde(:,:)
    integer(INT16),allocatable::ip_tilde(:,:),jp_tilde(:,:)
    contains
!
subroutine destory_spec_vars()
    deallocate(xtilde0)
    deallocate(ytilde0)
    deallocate(cos_beta_ref)
    deallocate(sin_beta_ref)
    deallocate(i0_tilde)
    deallocate(j0_tilde)
    deallocate(ip_tilde(:,:))
    deallocate(jp_tilde(:,:))
end subroutine destory_spec_vars
!
!tll2xy --- convert earth lon-lat to x-y grid coordinates
!
subroutine tll2xy(rlon,rlat,x,y,outside)
        implicit none
! !DESCRIPTION: to convert earth lon-lat to x-y grid units of a 
!           general regional rectangular domain.  Also, decide if
!           point is inside this domain.  As a result, there is
!           no restriction on type of horizontal coordinate for
!           a regional run, other than that it not have periodicity
!           or polar singularities.
!           This is done by first converting rlon, rlat to an
!           intermediate coordinate xtilde,ytilde, which has
!           precomputed pointers and constants for final conversion
!           to the desired x,y via 3 point inverse interpolation.
!           All of the information needed is derived from arrays
!           specifying earth latitude and longitude of every point
!           on the input grid.  Currently, the input x-y grid that
!           this is based on must be non-staggered.  This restriction
!           will eventually be lifted so we can run directly from
!           model grids that are staggered without first resorting
!           to interpolation of the guess to a non-staggered grid.
    real(P),intent(in)::rlon ! earth longitude (radians)
    real(P),intent(in)::rlat ! earth latitude  (radians)
    real(P),intent(  out) :: x  ! x-grid coordinate (grid units)
    real(P),intent(  out) :: y  ! y-grid coordinate (grid units)
    logical     ,intent(  out) :: outside     ! .false., then point is inside x-y domain
    !
    real(P):: clon,slon,r_of_lat,xtilde,ytilde
    real(P):: dtilde,etilde
    real(P):: d1tilde,d2tilde,e1tilde,e2tilde,detinv
    integer(INT32):: itilde,jtilde
    integer(INT32):: i0,j0,ip,jp
! if arrays are ready
    if(.not. allocated(i0_tilde)then
        call  init_general_transform(lat,lon)
    endif

!   first compute xtilde, ytilde

    clon=cos(rlon+rlambda0)
    slon=sin(rlon+rlambda0)
    r_of_lat=pihalf+sign_pole*rlat
    xtilde=atilde_x+btilde_x*r_of_lat*clon
    ytilde=atilde_y+btilde_y*r_of_lat*slon

!  next get interpolation information

    itilde=max(1,min(nint(xtilde),nxtilde))
    jtilde=max(1,min(nint(ytilde),nytilde))

    i0     =   i0_tilde(itilde,jtilde)
    j0     =   j0_tilde(itilde,jtilde)
    ip     =i0+ip_tilde(itilde,jtilde)
    jp     =j0+jp_tilde(itilde,jtilde)
    dtilde =xtilde-xtilde0(i0,j0)
    etilde =ytilde-ytilde0(i0,j0)
    d1tilde=(xtilde0(ip,j0)-xtilde0(i0,j0))*(ip-i0)
    d2tilde=(xtilde0(i0,jp)-xtilde0(i0,j0))*(jp-j0)
    e1tilde=(ytilde0(ip,j0)-ytilde0(i0,j0))*(ip-i0)
    e2tilde=(ytilde0(i0,jp)-ytilde0(i0,j0))*(jp-j0)
    detinv =one/(d1tilde*e2tilde-d2tilde*e1tilde)
    x = i0+detinv*(e2tilde*dtilde-d2tilde*etilde)
    y = j0+detinv*(d1tilde*etilde-e1tilde*dtilde)
    if (i0 == ip .and. j0 == jp) then ! ob at center of domain. 
       x = i0; y = j0
    else
       x = i0+detinv*(e2tilde*dtilde-d2tilde*etilde)
       y = j0+detinv*(d1tilde*etilde-e1tilde*dtilde)
    endif

    !print*, 'BOUNDS:', rlon_min_dd, rlon_max_dd, rlat_min_dd, rlat_max_dd

    outside=x < rlon_min_dd .or. x > rlon_max_dd .or. &
            y < rlat_min_dd .or. y > rlat_max_dd

 end subroutine tll2xy


end subroutine tll2xy
!-------------------------------------------------------------------------------
 subroutine nearest_3(ilast,jlast,i0,j0,ip,jp,x,y,nx0,ny0,x0,y0)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    nearest_3
!   prgmmr:
!
! abstract: find closest 3 points to (x,y) on grid defined by x0,y0
!
! program history log:
!   2009-08-04  lueken - added subprogram doc block
!
!   input argument list:
!    ilast,jlast
!    nx0,ny0
!    x,y
!    x0,y0
!
!   output argument list:
!    ilast,jlast
!    i0,j0
!    ip,jp
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

  implicit none

  integer(INT32),intent(inout) :: ilast,jlast
  integer(INT32),intent(  out) :: i0,j0
  integer(INT16),intent(  out) :: ip,jp
  integer(INT32),intent(in   ) :: nx0,ny0
  real(P)   ,intent(in   ) :: x,y
  real(P)   ,intent(in   ) :: x0(nx0,ny0),y0(nx0,ny0)

  real(P) dista,distb,dist2,dist2min
  integer(INT32) i,inext,j,jnext
  do
     i0=ilast
     j0=jlast
     dist2min=huge(dist2min)
     inext=0
     jnext=0
     do j=max(j0-1,1),min(j0+1,ny0)
        do i=max(i0-1,1),min(i0+1,nx0)
           dist2=(x-x0(i,j))**2+(y-y0(i,j))**2
           if(dist2<dist2min) then
              dist2min=dist2
              inext=i
              jnext=j
           end if
        end do
     end do
     if(inext==i0.and.jnext==j0) exit
     ilast=inext
     jlast=jnext
  end do

!  now find which way to go in x for second point
  ip=0
  if(i0==nx0)  ip=-1
  if(i0==1) ip=1
  if(ip==0) then
     dista=(x-x0(i0-1,j0))**2+(y-y0(i0-1,j0))**2
     distb=(x-x0(i0+1,j0))**2+(y-y0(i0+1,j0))**2
     if(distb<dista) then
        ip=1
     else
        ip=-1
     end if
  end if

!  repeat for y for 3rd point

  jp=0
  if(j0==ny0  ) jp=-1
  if(j0==1 ) jp=1
  if(jp==0) then
     dista=(x-x0(i0,j0-1))**2+(y-y0(i0,j0-1))**2
     distb=(x-x0(i0,j0+1))**2+(y-y0(i0,j0+1))**2
     if(distb<dista) then
        jp=1
     else
        jp=-1
     end if
  end if

  ilast=i0
  jlast=j0

 end subroutine nearest_3
!------------------------------------------------------------------------------
 subroutine get_xytilde_domain(nx0,ny0,rlons0,rlats0, &
                                  nx,ny,xminout,xmaxout,yminout,ymaxout)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    get_xytilde_domain
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-08-04  lueken - added subprogram doc block
!
!   input argument list:
!    nx0,ny0
!    rlons0,rlats0
!
!   output argument list:
!    nx,ny
!    xminout,xmaxout,yminout,ymaxout
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

!  define parameters for xy domain which optimally overlays input grid

  implicit none
  integer(INT32),intent(in   ) :: nx0,ny0
  real(P)   ,intent(in   ) :: rlons0(nx0,ny0),rlats0(nx0,ny0)

  integer(INT32),intent(  out) :: nx,ny
  real(P)   ,intent(  out) :: xminout,xmaxout,yminout,ymaxout

  real(P),parameter:: r37=37.0_P
  real(P),parameter:: r10=10.0_P

  real(P) area,areamax,areamin,extra,rlats0max,rlats0min,testlambda
  real(P) xthis,ythis
  integer(INT32) i,ip1,j,jp1,m

  real(P) coslon0(nx0,ny0),sinlon0(nx0,ny0)
  real(P) coslat0(nx0,ny0),sinlat0(nx0,ny0)
  real(P) icount,delbar
  real(P) dx,dy,disti,distj,distmin,distmax
  real(P) xmin,xmax,ymin,ymax

!  get range of lats for input grid

  rlats0max=maxval(rlats0) ; rlats0min=minval(rlats0)

!   assign hemisphere ( parameter sign_pole )

  sign_pole = zero
  if(rlats0min>-r37*deg2rad) sign_pole=-one   !  northern hemisphere xy domain
  if(rlats0max< r37*deg2rad) sign_pole= one   !  southern hemisphere xy domain
  ! if neither condition satisfied (rlat0max > 37N, rlat0min < 37S), try 
  ! this... 
  if (sign_pole == zero) then
     if (abs(rlats0max) > abs(rlats0min)) then
        sign_pole=-one  ! NH domain 
     else
        sign_pole=one   ! SH 
     endif
  endif


!   get optimum rotation angle rlambda0

  areamin= huge(areamin)
  areamax=-huge(areamax)
  do m=0,359
     testlambda=m*deg2rad
     xmax=-huge(xmax)
     xmin= huge(xmin)
     ymax=-huge(ymax)
     ymin= huge(ymin)
     do j=1,ny0,ny0-1
        do i=1,nx0
           xthis=(pihalf+sign_pole*rlats0(i,j))*cos(rlons0(i,j)+testlambda)
           ythis=(pihalf+sign_pole*rlats0(i,j))*sin(rlons0(i,j)+testlambda)
           xmax=max(xmax,xthis)
           ymax=max(ymax,ythis)
           xmin=min(xmin,xthis)
           ymin=min(ymin,ythis)
        end do
     end do
     do j=1,ny0
        do i=1,nx0,nx0-1
           xthis=(pihalf+sign_pole*rlats0(i,j))*cos(rlons0(i,j)+testlambda)
           ythis=(pihalf+sign_pole*rlats0(i,j))*sin(rlons0(i,j)+testlambda)
           xmax=max(xmax,xthis)
           ymax=max(ymax,ythis)
           xmin=min(xmin,xthis)
           ymin=min(ymin,ythis)
        end do
     end do
     area=(xmax-xmin)*(ymax-ymin)
     areamax=max(area,areamax)
     if(area<areamin) then
        areamin =area
        rlambda0=testlambda
        xmaxout =xmax
        xminout =xmin
        ymaxout =ymax
        yminout =ymin
     end if
  end do


!   now determine resolution of input grid and choose nx,ny of xy grid accordingly
!                 (currently hard-wired at 1/2 the average input grid increment)

  do j=1,ny0
     do i=1,nx0
        coslon0(i,j)=cos(one*rlons0(i,j)) ; sinlon0(i,j)=sin(one*rlons0(i,j))
        coslat0(i,j)=cos(one*rlats0(i,j)) ; sinlat0(i,j)=sin(one*rlats0(i,j))
     end do
  end do
  delbar=zero
  count =zero
  do j=1,ny0-1
     jp1=j+1
     do i=1,nx0-1
        ip1=i+1
        disti=acos(sinlat0(i,j)*sinlat0(ip1,j)+coslat0(i,j)*coslat0(ip1,j)* &
                  (sinlon0(i,j)*sinlon0(ip1,j)+coslon0(i,j)*coslon0(ip1,j)))
        distj=acos(sinlat0(i,j)*sinlat0(i,jp1)+coslat0(i,j)*coslat0(i,jp1)* &
                  (sinlon0(i,j)*sinlon0(i,jp1)+coslon0(i,j)*coslon0(i,jp1)))
        distmax=max(disti,distj)
        distmin=min(disti,distj)
        delbar=delbar+distmax
        icount=icount+one
     end do
  end do
  delbar=delbar/icount
  dx=half*delbar
  dy=dx

!   add extra space to computational grid to push any boundary problems away from
!     area of interest

  extra=r10*dx
  xmaxout=xmaxout+extra
  xminout=xminout-extra
  ymaxout=ymaxout+extra
  yminout=yminout-extra
  nx=1+(xmaxout-xminout)/dx
  ny=1+(ymaxout-yminout)/dy

 end subroutine get_xytilde_domain

!-------------------------------------------------------------------------------
!
 subroutine init_general_transform(glats,glons)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    init_general_transform
!   prgmmr:  parrish
!
! abstract:  set up constants to allow conversion between earth lat lon and analysis grid units.
!     There is no need to specify details of the analysis grid projection.  All that is required
!     is the earth latitude and longitude in radians of each analysis grid point.
!
! program history log:
!   2009-08-04  lueken - added subprogram doc block
!   2010-09-08  parrish - replace computation of wind rotation reference angle cos_beta_ref,sin_beta_ref
!                          with new, more accurate and robust version which works for any orientation
!                          of the analysis grid on the sphere (only restriction for now is that
!                          x-y coordinate of analysis grid is right handed).
!
!   input argument list:
!    glons,glats - lons,lats of input grid points of dimesion nlon,nlat
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

  implicit none

  real(P)   ,intent(in   ) :: glats(nlon,nlat),glons(nlon,nlat)

  real(P),parameter:: rbig =1.0e30_P
  real(P) xbar_min,xbar_max,ybar_min,ybar_max
  real(P) clon,slon,r_of_lat,xbar,ybar
  integer(INT32) i,j,istart0,iend,iinc,itemp,ilast,jlast
  real(P),allocatable:: clata(:,:),slata(:,:),clona(:,:),slona(:,:)
  real(P) clat0,slat0,clon0,slon0
  real(P) clat_m1,slat_m1,clon_m1,slon_m1
  real(P) clat_p1,slat_p1,clon_p1,slon_p1
  real(P) x,y,z,xt,yt,zt,xb,yb,zb
  real(P) rlonb_m1,clonb_m1,slonb_m1
  real(P) rlonb_p1,clonb_p1,slonb_p1
  real(P) crot,srot

  pihalf=half*pi

!  define xtilde, ytilde grid, transform

!      glons,glats are lons, lats of input grid points of dimension nlon,nlat
  call get_xytilde_domain(nlon,nlat,glons,glats,nxtilde,nytilde, &
                   xbar_min,xbar_max,ybar_min,ybar_max)
  allocate(i0_tilde(nxtilde,nytilde),j0_tilde(nxtilde,nytilde))
  allocate(ip_tilde(nxtilde,nytilde),jp_tilde(nxtilde,nytilde))
  allocate(xtilde0(nlon,nlat),ytilde0(nlon,nlat))

! define atilde_x, btilde_x, atilde_y, btilde_y

  btilde_x   =(nxtilde -one     )/(xbar_max-xbar_min)
  btilde_xinv=(xbar_max-xbar_min)/(nxtilde -one     )
  atilde_x   =one-btilde_x*xbar_min
  btilde_y   =(nytilde -one     )/(ybar_max-ybar_min)
  btilde_yinv=(ybar_max-ybar_min)/(nytilde -one     )
  atilde_y   =one-btilde_y*ybar_min

! define xtilde0,ytilde0
  do j=1,nlat
     do i=1,nlon
        r_of_lat=pihalf+sign_pole*glats(i,j)
        clon=cos(glons(i,j)+rlambda0)
        slon=sin(glons(i,j)+rlambda0)
        xbar=r_of_lat*clon
        ybar=r_of_lat*slon
        xtilde0(i,j)=atilde_x+btilde_x*xbar
        ytilde0(i,j)=atilde_y+btilde_y*ybar
     end do
  end do
!  now get i0_tilde, j0_tilde, ip_tilde,jp_tilde
  ilast=1 ; jlast=1
  istart0=nxtilde
  iend=1
  iinc=-1
  do j=1,nytilde
     itemp=istart0
     istart0=iend
     iend=itemp
     iinc=-iinc
     ybar=j
     do i=istart0,iend,iinc
        xbar=i
        call nearest_3(ilast,jlast,i0_tilde(i,j),j0_tilde(i,j), &
                       ip_tilde(i,j),jp_tilde(i,j),xbar,ybar,nlon,nlat,xtilde0,ytilde0)
     end do
  end do
!   new, more accurate and robust computation of cos_beta_ref and sin_beta_ref which is independent
!     of sign_pole and works for any orientation of grid on sphere (only restriction for now is that
!     x-y coordinate of analysis grid is right handed).
  allocate(clata(nlon,nlat),slata(nlon,nlat),clona(nlon,nlat),slona(nlon,nlat))
  allocate(cos_beta_ref(nlon,nlat),sin_beta_ref(nlon,nlat))
  do j=1,nlat
     do i=1,nlon
        clata(i,j)=cos(glats(i,j))
        slata(i,j)=sin(glats(i,j))
        clona(i,j)=cos(glons(i,j))
        slona(i,j)=sin(glons(i,j))
     end do
  end do
  do j=1,nlat
     do i=2,nlon-1

!     do all interior lon points to 2nd order accuracy

!   transform so pole is at rlat0,rlon0 and 0 meridian is tangent to earth latitude at rlat0,rlon0.

        clat0=clata(i,j) ; slat0=slata(i,j) ; clon0=clona(i,j) ; slon0=slona(i,j)

!    now obtain new coordinates for m1 and p1 points.

        clat_m1=clata(i-1,j) ; slat_m1=slata(i-1,j) ; clon_m1=clona(i-1,j) ; slon_m1=slona(i-1,j)
        clat_p1=clata(i+1,j) ; slat_p1=slata(i+1,j) ; clon_p1=clona(i+1,j) ; slon_p1=slona(i+1,j)

        x=clat_m1*clon_m1 ; y=clat_m1*slon_m1 ; z=slat_m1
        xt=x*clon0+y*slon0 ; yt=-x*slon0+y*clon0 ; zt=z
        yb=zt*clat0-xt*slat0
        xb=yt
        zb=xt*clat0+zt*slat0

        rlonb_m1=atan2(-yb,-xb)   !  the minus signs here are so line for m1 is directed same
        clonb_m1=cos(rlonb_m1)
        slonb_m1=sin(rlonb_m1)

        x=clat_p1*clon_p1 ; y=clat_p1*slon_p1 ; z=slat_p1
        xt=x*clon0+y*slon0 ; yt=-x*slon0+y*clon0 ; zt=z
        yb=zt*clat0-xt*slat0
        xb=yt
        zb=xt*clat0+zt*slat0
        rlonb_p1=atan2(yb,xb)
        clonb_p1=cos(rlonb_p1)
        slonb_p1=sin(rlonb_p1)
        crot=half*(clonb_m1+clonb_p1)
        srot=half*(slonb_m1+slonb_p1)
        cos_beta_ref(i,j)=crot*clon0-srot*slon0
        sin_beta_ref(i,j)=srot*clon0+crot*slon0
     end do
!               now do i=1 and i=nlon at 1st order accuracy
     i=1

!   transform so pole is at rlat0,rlon0 and 0 meridian is tangent to earth latitude at rlat0,rlon0.

        clat0=clata(i,j) ; slat0=slata(i,j) ; clon0=clona(i,j) ; slon0=slona(i,j)
!    now obtain new coordinates for m1 and p1 points.

        clat_p1=clata(i+1,j) ; slat_p1=slata(i+1,j) ; clon_p1=clona(i+1,j) ; slon_p1=slona(i+1,j)

        x=clat_p1*clon_p1 ; y=clat_p1*slon_p1 ; z=slat_p1
        xt=x*clon0+y*slon0 ; yt=-x*slon0+y*clon0 ; zt=z
        yb=zt*clat0-xt*slat0
        xb=yt
        zb=xt*clat0+zt*slat0
        rlonb_p1=atan2(yb,xb)
        clonb_p1=cos(rlonb_p1)
        slonb_p1=sin(rlonb_p1)
        crot=clonb_p1
        srot=slonb_p1
        cos_beta_ref(i,j)=crot*clon0-srot*slon0
        sin_beta_ref(i,j)=srot*clon0+crot*slon0

     i=nlon

!   transform so pole is at rlat0,rlon0 and 0 meridian is tangent to earth latitude at rlat0,rlon0.

        clat0=clata(i,j) ; slat0=slata(i,j) ; clon0=clona(i,j) ; slon0=slona(i,j)

!    now obtain new coordinates for m1 and p1 points.

        clat_m1=clata(i-1,j) ; slat_m1=slata(i-1,j) ; clon_m1=clona(i-1,j) ; slon_m1=slona(i-1,j)

        x=clat_m1*clon_m1 ; y=clat_m1*slon_m1 ; z=slat_m1
        xt=x*clon0+y*slon0 ; yt=-x*slon0+y*clon0 ; zt=z
        yb=zt*clat0-xt*slat0
        xb=yt
        zb=xt*clat0+zt*slat0

        rlonb_m1=atan2(-yb,-xb)   !  the minus signs here are so line for m1 is directed same
        clonb_m1=cos(rlonb_m1)
        slonb_m1=sin(rlonb_m1)

        crot=clonb_m1
        srot=slonb_m1
        cos_beta_ref(i,j)=crot*clon0-srot*slon0
        sin_beta_ref(i,j)=srot*clon0+crot*slon0
  end do

end subroutine init_general_transform

!-----------------------------
    

END MODULE gridmod

