! Adapted from GSI grdcrd.f90
subroutine grdcrd1(d,x,nx,flg)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    grdcrd1
!   prgmmr: parrish          org: np22                date: 2013-01-26
!
! abstract: same as grdcrd but for d a single value.  This is added to avoid
!            type mismatch errors with debug compile on WCOSS.
!
! program history log:
!   2013-01-26  parrish
!
!   input argument list:
!     d      - input point
!     x      - grid values
!     nx     - number of reference grid points
!     flg    - marks order of values in x
!              (1=increasing, -1=decreasing)
!
!   output argument list:
!     d        - point converted to grid units
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,i_kind
  use constants, only: one
  implicit none

  integer(i_kind)           ,intent(in   ) :: nx
  integer(i_kind)           ,intent(in   ) :: flg
  real(r_kind)              ,intent(inout) :: d
  real(r_kind),dimension(nx),intent(in   ) :: x

  integer(i_kind) ix,isrchf

! Treat "normal" case in which nx>1
  if(nx>1) then
     if (flg==1) then

!       Case in which x is in increasing order
        if(d<=x(1)) then
           ix=1
        else
           ix=isrchf(nx-1,x,d,flg)-1
        end if
        if(ix==nx) ix=ix-1

     else if (flg==(-1)) then

!       Case in which x is in decreasing order
        if(d>=x(1)) then
           ix=1
        else
           ix=isrchf(nx-1,x,d,flg)-1
        end if
     end if
     d=float(ix)+(d-x(ix))/(x(ix+1)-x(ix))

! Treat special case of nx=1
  elseif (nx==1) then
     d = one
  endif

  return
end subroutine grdcrd1

