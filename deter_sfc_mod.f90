! Adapted from GSI deter_sfc_mof.f90 4/21/2021
! : NO time Interpolation Now
MODULE deter_sfc_mod
!
! abstract: subroutine used to determine land surface type
  use read_wrf,only:nx,ny
  use read_wrf,only:zs_full,sst_full
  use read_wrf,only:soil_moi_full,soil_temp_full
  use read_wrf,only:sno_full
  use read_wrf,only:isli_full
  use read_wrf,only:sfc_rough_full
  use read_wrf,only:fact10_full
  use read_wrf,only:veg_type_full=>ivgtyp
  use read_wrf,only:soil_type_full=>isltyp
  use read_wrf,only:veg_frac_full=>vegfrac
  use parameters_define,only:zero,one,regional
  use model_precision,only:i_kind,r_kind
  implicit none

! Set default to private
  private
! Set passed variables to public
  public deter_sfc

  contains
subroutine deter_sfc(alat,alon,dlat_earth,dlon_earth,obstime,isflg, &
       idomsfc,sfcpct,ts,tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10,sfcr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    deter_sfc                     determine land surface type
!   prgmmr: derber           org: np2                date: 2005-01-27
!
! abstract:  determines land surface type based on surrounding land
!            surface types
!
! program history log:
!   2005-01-27 derber
!   2005-03-03 treadon - add implicit none, define zero
!   2006-02-01 parrish  - change names of sno,isli,sst
!
!   input argument list:
!     alat
!     alon
!     obstime- observation time relative to analysis time
!     dlat_earth
!     dlon_earth
!
!   output argument list:
!      isflg    - surface flag
!                 0 sea
!                 1 land
!                 2 sea ice
!                 3 snow
!                 4 mixed
!      sfcpct(0:3)- percentage of 4 surface types
!                 (0) - sea percentage
!                 (1) - land percentage
!                 (2) - sea ice percentage
!                 (3) - snow percentage
!      tsavg - sea surface temperature
!      idomsfc
!      ts
!      dfcr
!      vty,vfr,sty,stp,sm,sn,zz,ff10
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

     implicit none
     real(r_kind)               ,intent(in   ) :: dlat_earth,dlon_earth,obstime,alat,alon
     integer(i_kind)            ,intent(  out) :: isflg,idomsfc
     real(r_kind),dimension(0:3),intent(  out) :: sfcpct
     real(r_kind),dimension(0:3),intent(  out) :: ts
     real(r_kind)               ,intent(  out) :: tsavg,sfcr
     real(r_kind)               ,intent(  out) :: vty,vfr,sty,stp,sm,sn,zz,ff10

     real(r_kind),parameter:: minsnow=0.1_r_kind

     integer(i_kind) istyp00,istyp01,istyp10,istyp11
     integer(i_kind):: ix,iy,ixp,iyp,j
     real(r_kind):: dx,dy,dx1,dy1,w00,w10,w01,w11,wgtmin
     real(r_kind):: sno00,sno01,sno10,sno11,dlat,dlon
     real(r_kind):: sst00,sst01,sst10,sst11
     real(r_kind),dimension(0:3)::wgtavg
     !
     integer(i_kind)::nlon,nlat
     integer(i_kind)::nlon_sfc,nlat_sfc
     !
     nlon=nx;nlat=ny
     nlon_sfc=nx;nlat_sfc=ny

!  First do surface field since it is on model grid
     iy=int(alon); ix=int(alat)
     dy  =alon-iy; dx  =alat-ix
     dx1 =one-dx;    dy1 =one-dy
     w00=dx1*dy1; w10=dx*dy1; w01=dx1*dy; w11=dx*dy

     ix=min(max(1,ix),nlat); iy=min(max(0,iy),nlon)
     ixp=min(nlat,ix+1); iyp=iy+1
     if(iy==0) iy=nlon
     if(iyp==nlon+1) iyp=1

!    Interpolate fields which only vary in space (no time component)
!       zz   = surface height
     zz   = zs_full(ix,iy) *w00 + zs_full(ixp,iy) *w10 + &
            zs_full(ix,iyp)*w01 + zs_full(ixp,iyp)*w11

     if(regional)then
        dlat=alat
        dlon=alon
     else
         print*,"Only process regional domain"
         return
     end if
     iy=int(dlon); ix=int(dlat)
     dy  =dlon-iy; dx  =dlat-ix
     dx1 =one-dx;    dy1 =one-dy
     w00=dx1*dy1; w10=dx*dy1; w01=dx1*dy; w11=one-w00-w10-w01
     ix=min(max(1,ix),nlat_sfc); iy=min(max(0,iy),nlon_sfc)
     ixp=min(nlat_sfc,ix+1); iyp=iy+1
     if(iy==0) iy=nlon_sfc
     if(iyp==nlon_sfc+1) iyp=1

!    Set surface type flag.  Begin by assuming obs over ice-free water

     istyp00 = isli_full(ix ,iy )
     istyp10 = isli_full(ixp,iy )
     istyp01 = isli_full(ix ,iyp)
     istyp11 = isli_full(ixp,iyp)

     ! No time Interpolation now
     ! snow water
     sno00= sno_full(ix ,iy )
     sno01= sno_full(ix ,iyp)
     sno10= sno_full(ixp,iy )
     sno11= sno_full(ixp,iyp)
     !sst
     sst00= sst_full(ix ,iy )
     sst01= sst_full(ix ,iyp)
     sst10= sst_full(ixp,iy )
     sst11= sst_full(ixp,iyp)
!    Interpolate sst to obs location

     tsavg=sst00*w00+sst10*w10+sst01*w01+sst11*w11

     if(istyp00 >=1 .and. sno00 > minsnow)istyp00 = 3
     if(istyp01 >=1 .and. sno01 > minsnow)istyp01 = 3
     if(istyp10 >=1 .and. sno10 > minsnow)istyp10 = 3
     if(istyp11 >=1 .and. sno11 > minsnow)istyp11 = 3

     sfcpct = zero
     sfcpct(istyp00)=sfcpct(istyp00)+w00
     sfcpct(istyp01)=sfcpct(istyp01)+w01
     sfcpct(istyp10)=sfcpct(istyp10)+w10
     sfcpct(istyp11)=sfcpct(istyp11)+w11

     isflg = 0
     if(sfcpct(0) > 0.99_r_kind)then
        isflg = 0
     else if(sfcpct(1) > 0.99_r_kind)then
        isflg = 1
     else if(sfcpct(2) > 0.99_r_kind)then
        isflg = 2
     else if(sfcpct(3) > 0.99_r_kind)then
        isflg = 3
     else
        isflg = 4
     end if

!       vty  = vegetation type
!       sty  = soil type

     ts(0:3)=zero
     wgtavg(0:3)=zero
     vfr=zero
     stp=zero
     sty=zero
     vty=zero
     sm=zero
     sn=zero
     idomsfc=isli_full(ix ,iy )
     wgtmin = w00
     if(istyp00 == 1)then
        vty  = veg_type_full(ix ,iy)
        sty  = soil_type_full(ix ,iy)
        wgtavg(1) = wgtavg(1) + w00
        ts(1)=ts(1)+w00*sst00
        vfr  =vfr  +w00*(veg_frac_full(ix ,iy  )    &
                         )
        stp  =stp  +w00*(soil_temp_full(ix ,iy  )   &
                         )
        sm   =sm   +w00*(soil_moi_full(ix ,iy  )   &
                         )
     else if(istyp00 == 2)then
        wgtavg(2) = wgtavg(2) + w00
        ts(2)=ts(2)+w00*sst00
     else if(istyp00 == 3)then
        wgtavg(3) = wgtavg(3) + w00
        ts(3)=ts(3)+w00*sst00
        sn = sn + w00*sno00
     else
        wgtavg(0) = wgtavg(0) + w00
        ts(0)=ts(0)+w00*sst00
     end if
     if(istyp01 == 1)then
        if(wgtmin < w01 .or. (vty == zero .and. sty == zero))then
           vty  = veg_type_full(ix ,iyp)
           sty  = soil_type_full(ix ,iyp)
        end if
        wgtavg(1) = wgtavg(1) + w01
        ts(1)=ts(1)+w01*sst01
        vfr  =vfr  +w01*(veg_frac_full(ix ,iyp )    &
                         )
        stp  =stp  +w01*(soil_temp_full(ix ,iyp )   &
                         )
        sm   =sm   +w01*(soil_moi_full(ix ,iyp )    &
                         )
     else if(istyp01 == 2)then
        wgtavg(2) = wgtavg(2) + w01
        ts(2)=ts(2)+w01*sst01
     else if(istyp01 == 3)then
        wgtavg(3) = wgtavg(3) + w01
        ts(3)=ts(3)+w01*sst01
        sn = sn + w01*sno01
     else
        wgtavg(0) = wgtavg(0) + w01
        ts(0)=ts(0)+w01*sst01
     end if
     if(wgtmin < w01)then
        idomsfc=isli_full(ix ,iyp)
        wgtmin = w01
     end if
     if(istyp10 == 1)then
        if(wgtmin < w10 .or. (vty == zero .and. sty == zero))then
           vty  = veg_type_full(ixp,iy)
           sty  = soil_type_full(ixp,iy)
        end if
        wgtavg(1) = wgtavg(1) + w10
        ts(1)=ts(1)+w10*sst10
        vfr  =vfr  +w10*(veg_frac_full(ixp,iy  )    &
                         )
        stp  =stp  +w10*(soil_temp_full(ixp,iy  )  &
                         )
        sm   =sm   +w10*(soil_moi_full(ixp,iy  )   &
                         )
     else if(istyp10 == 2)then
        wgtavg(2) = wgtavg(2) + w10
        ts(2)=ts(2)+w10*sst10
     else if(istyp10 == 3)then
        wgtavg(3) = wgtavg(3) + w10
        ts(3)=ts(3)+w10*sst10
        sn = sn + w10*sno10
     else
        wgtavg(0) = wgtavg(0) + w10
        ts(0)=ts(0)+w10*sst10
     end if
     if(wgtmin < w10)then
        idomsfc=isli_full(ixp,iy )
        wgtmin = w10
     end if
     if(istyp11 == 1)then
        if(wgtmin < w11 .or. (vty == zero .and. sty == zero))then
           vty  = veg_type_full(ixp,iyp)
           sty  = soil_type_full(ixp,iyp)
        endif
        wgtavg(1) = wgtavg(1) + w11
        ts(1)=ts(1)+w11*sst11
        vfr  =vfr  +w11*(veg_frac_full(ixp,iyp )    &
                         )
        stp  =stp  +w11*(soil_temp_full(ixp,iyp )   &
                         )
        sm   =sm   +w11*(soil_moi_full(ixp,iyp )    &
                         )
     else if(istyp11 == 2)then
        wgtavg(2) = wgtavg(2) + w11
        ts(2)=ts(2)+w11*sst11
     else if(istyp11 == 3)then
        wgtavg(3) = wgtavg(3) + w11
        ts(3)=ts(3)+w11*sst11
        sn = sn + w11*sno11
     else
        wgtavg(0) = wgtavg(0) + w11
        ts(0)=ts(0)+w11*sst11
     end if
     if(wgtmin < w11)then
        idomsfc=isli_full(ixp,iyp)
        wgtmin = w11
     end if
     if(wgtavg(0) > zero)then
        ts(0) = ts(0)/wgtavg(0)
     else
        ts(0) = tsavg
     end if
     if(wgtavg(1) > zero)then
        ts(1) = ts(1)/wgtavg(1)
        sm = sm/wgtavg(1)
        vfr = vfr/wgtavg(1)
        stp = stp/wgtavg(1)
     else
        ts(1) = tsavg
        sm=one
     end if
     if(wgtavg(2) > zero)then
        ts(2) = ts(2)/wgtavg(2)
     else
        ts(2) = tsavg
     end if
     if(wgtavg(3) > zero)then
        ts(3) = ts(3)/wgtavg(3)
        sn = sn/wgtavg(3)
     else
        ts(3) = tsavg
     end if
!    ts(0)=max(ts(0),270._r_kind)
!    ts(2)=min(ts(2),280._r_kind)
!    ts(3)=min(ts(3),280._r_kind)

     ff10=(fact10_full(ix ,iy  )*w00+ &
           fact10_full(ixp,iy  )*w10+ &
           fact10_full(ix ,iyp )*w01+ &
           fact10_full(ixp,iyp )*w11)

     sfcr=(sfc_rough_full(ix ,iy  )*w00+ &
           sfc_rough_full(ixp,iy  )*w10+ &
           sfc_rough_full(ix ,iyp )*w01+ &
           sfc_rough_full(ixp,iyp )*w11)

     return
end subroutine deter_sfc

END MODULE deter_sfc_mod
