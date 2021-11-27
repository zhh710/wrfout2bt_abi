!
!  11/11/2021
!  Huanhuan ZHANG
!
!  Adapted from WRFDA/phys/module_ra_cam_support.F subroutine oznini
!
!

subroutine oznini(ozmixm,pin,levsiz,num_months,XLAT, &
           &      lat2,lon2)
!
! This subroutine assumes uniform distribution of ozone concentration.
! It should be replaced by monthly climatology that varies latitudinally
! and vertically
!
    implicit none
    real,external :: lin_interpol2


    integer,intent(in   ):: levsiz
    integer,intent(in   ):: lat2,lon2
    integer,intent(in   ):: num_months

    real,dimension(lon2,lat2),intent(in   ):: XLAT
    real,dimension(levsiz),   intent(out  ):: pin
    real,dimension(lon2,levsiz,lat2,num_months),intent(out):: ozmixm
! Local
    INTEGER, PARAMETER :: latsiz = 64
    INTEGER, PARAMETER :: lonsiz = 1
    !
    real,dimension(:,:,:,:),allocatable:: ozmixin
    real,dimension(:),allocatable::  lat_ozone, plev_ozone
    !
    integer::k,m,i,j
    real::interp_pt
    !
    allocate(plev_ozone(levsiz),lat_ozone(latsiz))
    allocate(ozmixin(lonsiz, levsiz, latsiz, num_months))
    !
    OPEN(37, FILE='ozone_plev.formatted',FORM='FORMATTED',STATUS='OLD')
    do k = 1,levsiz
        READ (37,*)plev_ozone(k)
    end do
    close(37)
!!-- read in ozone pressure data
    do k = 1,levsiz
        plev_ozone(k) = plev_ozone(k) * 100.
    end do

    pin = plev_ozone
!-- read in ozone lat data
    OPEN(38, FILE='ozone_lat.formatted',FORM='FORMATTED',STATUS='OLD')
    do k=1,latsiz
        READ(38,*)lat_ozone(k)
    enddo
    close(38)
!-- read in ozone data
    OPEN(39, FILE='ozone.formatted',FORM='FORMATTED',STATUS='OLD')
    do m=1,12
    do j=1,latsiz ! latsiz=64
    do k=1,levsiz ! levsiz=59
    do i=1,lonsiz ! lonsiz=1
       READ (39,*)ozmixin(i,k,j,m)
    enddo
    enddo
    enddo
    enddo
    close(39)
!-- latitudinally interpolate ozone data (and extend longitudinally)
!-- using function lin_interpol2(x, f, y) result(g)
! Purpose:
!   interpolates f(x) to point y
!   assuming f(x) = f(x0) + a * (x - x0)
!   where a = ( f(x1) - f(x0) ) / (x1 - x0)
!   x0 <= x <= x1
!   assumes x is monotonically increasing
!    real, intent(in), dimension(:) :: x  ! grid points
!    real, intent(in), dimension(:) :: f  ! grid function values
!    real, intent(in) :: y                ! interpolation point
!    real :: g                            ! interpolated function value
!---------------------------------------------------------------------------
    do m=1,12
    do j=1,lat2
    do k=1,levsiz
    do i=1,lon2
     interp_pt=XLAT(i,j)
     ozmixm(i,k,j,m)=lin_interpol2(lat_ozone(:),ozmixin(1,k,:,m),interp_pt,latsiz)
    enddo
    enddo
    enddo
    enddo

     deallocate(lat_ozone)
     deallocate(plev_ozone)
     deallocate(ozmixin)
end subroutine oznini


function lin_interpol2(x, f, y,n) result(g)
! Purpose:
!   interpolates f(x) to point y
!   assuming f(x) = f(x0) + a * (x - x0)
!   where a = ( f(x1) - f(x0) ) / (x1 - x0)
!   x0 <= x <= x1
!   assumes x is monotonically increasing
! Author: D. Fillmore ::  J. Done changed from r8 to r4
implicit none
integer ,intent(in):: n  ! length of x
real, intent(in), dimension(n) :: x  ! grid points
real, intent(in), dimension(n) :: f  ! grid function values
real, intent(in) :: y                ! interpolation point
real :: g                            ! interpolated function value
integer :: k  ! interpolation point index

real    :: a
!n = size(x)
! find k such that x(k) < y =< x(k+1)
! set k = 1 if y <= x(1)  and  k = n-1 if y > x(n)
if (y <= x(1)) then
  k = 1
else if (y >= x(n)) then
  k = n - 1
else
  k = 1
  do while (y > x(k+1) .and. k < n)
    k = k + 1
  end do
end if
! interpolate
a = (  f(k+1) - f(k) ) / ( x(k+1) - x(k) )
g = f(k) + a * (y - x(k))

end function lin_interpol2  

SUBROUTINE ozn_time_int(julday,ozmixm,ozmixt,lon2,lat2,levsiz,num_months )

! adapted from oznint from CAM module
!  input: ozmixm - read from physics_init
! output: ozmixt - time interpolated

!  USE module_ra_cam_support, ONLY : getfactors

   IMPLICIT NONE

   INTEGER,      INTENT(IN   )    ::   levsiz, num_months,lon2,lat2

   REAL,  DIMENSION( lon2, levsiz, lat2, num_months ),      &
          INTENT(IN   ) ::                                  ozmixm
   integer,    INTENT(IN )      ::        JULDAY

   REAL,  DIMENSION( lon2, levsiz, lat2 ), INTENT(OUT  ) :: ozmixt

   !Local
   integer   :: np1,np,nm,m,k,i,j
   integer, dimension(12) ::  date_oz
   data date_oz/16, 45, 75, 105, 136, 166, 197, 228, 258, 289, 319, 350/
   real, parameter :: daysperyear = 365.  ! number of days in a year
   real      :: cdayozp, cdayozm
   real      :: fact1, fact2, deltat
   real:: intjulian
   logical   :: finddate
   logical   :: ozncyc
   CHARACTER(LEN=256) :: msgstr

   ozncyc = .true.
   np1=1

   intjulian = julday
!  do m=1,num_months
   do m=1,12
      if(date_oz(m).lt.julday) then
        np1=m
      else
        exit
      endif
   enddo
   np1=np1+1
   print*,"ozone.f90: month= ",np1,julday
   cdayozp=date_oz(np1)

   if(np1.gt.1) then
      cdayozm=date_oz(np1-1)
      np=np1
      nm=np-1
   else
      cdayozm=date_oz(12)
      np=np1
      nm=12
   endif

!  call getfactors(ozncyc,np1, cdayozm, cdayozp,intjulian, &
!                   fact1, fact2)
!
! Determine time interpolation factors.  Account for December-January
! interpolation if dataset is being cycled yearly.
!
   if (ozncyc .and. np1 == 1) then                      ! Dec-Jan interpolation
      deltat = cdayozp + daysperyear - cdayozm
      if (intjulian > cdayozp) then                     ! We are in December
         fact1 = (cdayozp + daysperyear - intjulian)/deltat
         fact2 = (intjulian - cdayozm)/deltat
      else                                              ! We are in January
         fact1 = (cdayozp - intjulian)/deltat
         fact2 = (intjulian + daysperyear - cdayozm)/deltat
      end if
   else
      deltat = cdayozp - cdayozm
      fact1 = (cdayozp - intjulian)/deltat
      fact2 = (intjulian - cdayozm)/deltat
   end if
!
! Time interpolation.
!
      do j=1,lat2
      do k=1,levsiz
      do i=1,lon2
            ozmixt(i,k,j) = ozmixm(i,k,j,nm)*fact1 + ozmixm(i,k,j,np)*fact2
      end do
      end do
      end do

END SUBROUTINE ozn_time_int

SUBROUTINE ozn_p_int(p ,pin, levsiz, ozmixt, o3vmr,lon2,lat2,nz)

!-----------------------------------------------------------------------
!
! Purpose: Interpolate ozone from current time-interpolated values to model levels
!
! Method: Use pressure values to determine interpolation levels
!
! Author: Bruce Briegleb
! WW: Adapted for general use
!
!--------------------------------------------------------------------------
   implicit none
!--------------------------------------------------------------------------
!
! Arguments
!

   integer, intent(in) :: levsiz              ! number of ozone layers
   integer, intent(in) :: lon2,lat2,nz

   real, intent(in) :: p(lon2,lat2,nz)   ! level pressures (mks, bottom-up)
   real, intent(in) :: pin(levsiz)        ! ozone data level pressures (mks, top-down)
   real, intent(in) :: ozmixt(lon2,levsiz,lat2) ! ozone mixing ratio

   real, intent(out) :: o3vmr(lon2,lat2,nz) ! ozone volume mixing ratio
!
! local storage
!
   real    pmid(lon2,nz)
   integer i,j                 ! longitude index
   integer k, kk, kkstart, kout! level indices
   integer kupper(lon2)     ! Level indices for interpolation
   integer kount               ! Counter
   integer ncol, pver

   real    dpu                 ! upper level pressure difference
   real    dpl                 ! lower level pressure difference

   ncol = lon2
   pver = nz

   do j=1,lat2
!
! Initialize index array
!
!  do i=1, ncol
   do i=1, lon2
      kupper(i) = 1
   end do
!
! Reverse the pressure array, and pin is in Pa, the same as model pmid
!
      do k = 1,nz
         kk = nz - k + 1
      do i = 1,lon2
         pmid(i,kk) = p(i,j,k)
      enddo
      enddo

   do k=1,pver

      kout = pver - k + 1
!     kout = k
!
! Top level we need to start looking is the top level for the previous k
! for all longitude points
!
      kkstart = levsiz
!     do i=1,ncol
      do i=1,lon2
         kkstart = min0(kkstart,kupper(i))
      end do
      kount = 0
!
! Store level indices for interpolation
!
      do kk=kkstart,levsiz-1
!        do i=1,ncol
         do i=1,lon2
            if (pin(kk).lt.pmid(i,k) .and. pmid(i,k).le.pin(kk+1)) then
               kupper(i) = kk
               kount = kount + 1
            end if
         end do
!
! If all indices for this level have been found, do the interpolation and
! go to the next level
!
         if (kount.eq.ncol) then
!           do i=1,ncol
            do i=1,lon2
               dpu = pmid(i,k) - pin(kupper(i))
               dpl = pin(kupper(i)+1) - pmid(i,k)
               o3vmr(i,j,kout) = (ozmixt(i,kupper(i),j)*dpl + &
                             ozmixt(i,kupper(i)+1,j)*dpu)/(dpl + dpu)
            end do
            goto 35
         end if
      end do
!
! If we've fallen through the kk=1,levsiz-1 loop, we cannot interpolate and
! must extrapolate from the bottom or top ozone data level for at least some
! of the longitude points.
!
!     do i=1,ncol
      do i=1,lon2
         if (pmid(i,k) .lt. pin(1)) then
            o3vmr(i,j,kout) = ozmixt(i,1,j)*pmid(i,k)/pin(1)
         else if (pmid(i,k) .gt. pin(levsiz)) then
            o3vmr(i,j,kout) = ozmixt(i,levsiz,j)
         else
            dpu = pmid(i,k) - pin(kupper(i))
            dpl = pin(kupper(i)+1) - pmid(i,k)
            o3vmr(i,j,kout) = (ozmixt(i,kupper(i),j)*dpl + &
                          ozmixt(i,kupper(i)+1,j)*dpu)/(dpl + dpu)
         end if
      end do

      if (kount.gt.ncol) then
!        call endrun ('OZN_P_INT: Bad ozone data: non-monotonicity suspected')
         stop 'OZN_P_INT: Bad ozone data: non-monotonicity suspected'
      end if
35    continue

   end do
   end do

   return
END SUBROUTINE ozn_p_int
