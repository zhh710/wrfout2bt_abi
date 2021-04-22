! Adapted from GSI deter_sfc_mof.f90 4/21/2021
MODULE deter_sfc_mod
!
! abstract: subroutine used to determine land surface type
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

     real(r_kind),parameter:: minsnow=one_tenth

     integer(i_kind) istyp00,istyp01,istyp10,istyp11
     integer(i_kind):: itsfc,itsfcp
     integer(i_kind):: ix,iy,ixp,iyp,j
     real(r_kind):: dx,dy,dx1,dy1,w00,w10,w01,w11,dtsfc,dtsfcp,wgtmin
     real(r_kind):: sno00,sno01,sno10,sno11,dlat,dlon
     real(r_kind):: sst00,sst01,sst10,sst11
     real(r_kind),dimension(0:3)::wgtavg

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
!       call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)
        dlat=alat
        dlon=alon
     else
        dlat=dlat_earth
        dlon=dlon_earth
        call grdcrd1(dlat,rlats_sfc,nlat_sfc,1)
        call grdcrd1(dlon,rlons_sfc,nlon_sfc,1)
     end if
     iy=int(dlon); ix=int(dlat)
     dy  =dlon-iy; dx  =dlat-ix
     dx1 =one-dx;    dy1 =one-dy
     w00=dx1*dy1; w10=dx*dy1; w01=dx1*dy; w11=one-w00-w10-w01
     ix=min(max(1,ix),nlat_sfc); iy=min(max(0,iy),nlon_sfc)
     ixp=min(nlat_sfc,ix+1); iyp=iy+1
     if(iy==0) iy=nlon_sfc
     if(iyp==nlon_sfc+1) iyp=1

!    Get time interpolation factors for surface files
     if(obstime > hrdifsfc(1) .and. obstime <= hrdifsfc(nfldsfc))then
        do j=1,nfldsfc-1
           if(obstime > hrdifsfc(j) .and. obstime <= hrdifsfc(j+1))then
              itsfc=j
              itsfcp=j+1
              dtsfc=(hrdifsfc(j+1)-obstime)/(hrdifsfc(j+1)-hrdifsfc(j))
           end if
        end do
     else if(obstime <=hrdifsfc(1))then
        itsfc=1
        itsfcp=1
        dtsfc=one
     else
        itsfc=nfldsfc
        itsfcp=nfldsfc
        dtsfc=one
     end if
     dtsfcp=one-dtsfc

!    Set surface type flag.  Begin by assuming obs over ice-free water

     istyp00 = isli_full(ix ,iy )
     istyp10 = isli_full(ixp,iy )
     istyp01 = isli_full(ix ,iyp)
     istyp11 = isli_full(ixp,iyp)

     sno00= sno_full(ix ,iy ,itsfc)*dtsfc+sno_full(ix ,iy ,itsfcp)*dtsfcp
     sno01= sno_full(ix ,iyp,itsfc)*dtsfc+sno_full(ix ,iyp,itsfcp)*dtsfcp
     sno10= sno_full(ixp,iy ,itsfc)*dtsfc+sno_full(ixp,iy ,itsfcp)*dtsfcp
     sno11= sno_full(ixp,iyp,itsfc)*dtsfc+sno_full(ixp,iyp,itsfcp)*dtsfcp

     sst00= sst_full(ix ,iy ,itsfc)*dtsfc+sst_full(ix ,iy ,itsfcp)*dtsfcp
     sst01= sst_full(ix ,iyp,itsfc)*dtsfc+sst_full(ix ,iyp,itsfcp)*dtsfcp
     sst10= sst_full(ixp,iy ,itsfc)*dtsfc+sst_full(ixp,iy ,itsfcp)*dtsfcp
     sst11= sst_full(ixp,iyp,itsfc)*dtsfc+sst_full(ixp,iyp,itsfcp)*dtsfcp
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
        vfr  =vfr  +w00*(veg_frac_full(ix ,iy ,itsfc ) *dtsfc+   &
                         veg_frac_full(ix ,iy ,itsfcp) *dtsfcp)
        stp  =stp  +w00*(soil_temp_full(ix ,iy ,itsfc )*dtsfc+   &
                         soil_temp_full(ix ,iy ,itsfcp)*dtsfcp)
        sm   =sm   +w00*(soil_moi_full(ix ,iy ,itsfc ) *dtsfc+   &
                         soil_moi_full(ix ,iy ,itsfcp) *dtsfcp)
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
        vfr  =vfr  +w01*(veg_frac_full(ix ,iyp,itsfc ) *dtsfc+   &
                         veg_frac_full(ix ,iyp,itsfcp) *dtsfcp)
        stp  =stp  +w01*(soil_temp_full(ix ,iyp,itsfc )*dtsfc+   &
                         soil_temp_full(ix ,iyp,itsfcp)*dtsfcp)
        sm   =sm   +w01*(soil_moi_full(ix ,iyp,itsfc ) *dtsfc+   &
                         soil_moi_full(ix ,iyp,itsfcp) *dtsfcp)
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
        vfr  =vfr  +w10*(veg_frac_full(ixp,iy ,itsfc ) *dtsfc+   &
                         veg_frac_full(ixp,iy ,itsfcp) *dtsfcp)
        stp  =stp  +w10*(soil_temp_full(ixp,iy ,itsfc )*dtsfc+   &
                         soil_temp_full(ixp,iy ,itsfcp)*dtsfcp)
        sm   =sm   +w10*(soil_moi_full(ixp,iy ,itsfc ) *dtsfc+   &
                         soil_moi_full(ixp,iy ,itsfcp) *dtsfcp)
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
        vfr  =vfr  +w11*(veg_frac_full(ixp,iyp,itsfc ) *dtsfc+   &
                         veg_frac_full(ixp,iyp,itsfcp) *dtsfcp)
        stp  =stp  +w11*(soil_temp_full(ixp,iyp,itsfc )*dtsfc+   &
                         soil_temp_full(ixp,iyp,itsfcp)*dtsfcp)
        sm   =sm   +w11*(soil_moi_full(ixp,iyp,itsfc ) *dtsfc+   &
                         soil_moi_full(ixp,iyp,itsfcp) *dtsfcp)
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
!    Space-time interpolation of fields from surface wind speed

     ff10=(fact10_full(ix ,iy ,itsfc )*w00+ &
           fact10_full(ixp,iy ,itsfc )*w10+ &
           fact10_full(ix ,iyp,itsfc )*w01+ &
           fact10_full(ixp,iyp,itsfc )*w11)*dtsfc + &
          (fact10_full(ix ,iy ,itsfcp)*w00+ &
           fact10_full(ixp,iy ,itsfcp)*w10+ &
           fact10_full(ix ,iyp,itsfcp)*w01+ &
           fact10_full(ixp,iyp,itsfcp)*w11)*dtsfcp

     sfcr=(sfc_rough_full(ix ,iy ,itsfc )*w00+ &
           sfc_rough_full(ixp,iy ,itsfc )*w10+ &
           sfc_rough_full(ix ,iyp,itsfc )*w01+ &
           sfc_rough_full(ixp,iyp,itsfc )*w11)*dtsfc + &
          (sfc_rough_full(ix ,iy ,itsfcp)*w00+ &
           sfc_rough_full(ixp,iy ,itsfcp)*w10+ &
           sfc_rough_full(ix ,iyp,itsfcp)*w01+ &
           sfc_rough_full(ixp,iyp,itsfcp)*w11)*dtsfcp

     return
end subroutine deter_sfc

END MODULE deter_sfc_mod
