!4/24/2021 
!
module module_abi
    use model_precision,only:i_kind=>INT32,r_kind=>P,P
    !
    use crtm_interface,only:init_crtm,nsigradjac,nsigaerojac
    use crtm_interface,only:call_crtm,destroy_crtm 
    use crtm_interface,only:channelinfo
    use crtm_interface,only:iqv,icw
    !
    use parameters_define,only:lcloud_fwd,lallsky,cld_sea_only
    use parameters_define,only:n_actual_clouds,n_clouds_fwd,n_clouds_jac
    use parameters_define,only:cloud_names,cloud_names_fwd,cloud_names_jac
    use parameters_define,only:laerosol_fwd,laerosol,n_actual_aerosols
    use parameters_define,only:n_aerosols_fwd,n_aerosols_jac
    use parameters_define,only:aerosol_names,aerosol_names_fwd
    use parameters_define,only:aerosol_names_jac
    use parameters_define,only:n_ghg,ghg_names
    use parameters_define,only:obstype
    use parameters_define,only:isis,nchanl,subset_start,subset_end
    use parameters_define,only:crtm_coeffs_path
    use parameters_define,only:nvege_type0=>nvege_type,regional0=>regional
    use parameters_define,only:tbfout,datapath
    use parameters_define,only:ldiagout
    !
    use read_wrf,only: tropprs
    use read_wrf,only: lat2=>nx,lon2=>ny,nsig=>nz
    use read_wrf,only: psges=>psfc,tkges=>tk,qges=>qv
    use read_wrf,only: prslges=>pmid,prsiges=>pml
    use read_wrf,only: u,v
    use read_wrf,only: isli2=>isli_full ,sno2=>sno_full
    use read_wrf,only: qrges=>qr,qsges=>qs,qgges=>qg
    use read_wrf,only: qiges=>qi,qhges=>qh,qlges=>qc
    use read_wrf,only: iadate=>idates
    use read_wrf,only: lon,lat
    use read_wrf,only: ozges=>oz
    use read_wrf,only: t,pmid,hgt,qv
    !
    use goesabi_obs,only: data_s=>data_obsabi
    use goesabi_obs,only: read_abiobsarray_from_file
    use goesabi_obs,only: nreal,nabiobs
    use goesabi_obs,only: ilone,ilate
    !

    implicit none

    public::setuprad 
contains

    subroutine setuprad()
    !
    ! 1.calculate clear  sky tb and cloudy tb
    ! 2. write out for further use
    !
        implicit none
        character(len=*),parameter::myname_="setuprad"
        logical  :: init_pass
        integer(i_kind) :: mype_diaghdr,mype
        logical::cold_start 
        character(10)                          :: obstype
        !
        integer(i_kind)::n,i,k
        real(r_kind)::obstime
        ! variables used in call_crtm
        real(r_kind),dimension(lat2,lon2)::dsfct
        real(r_kind)                           :: trop5,tzbgr
        real(r_kind)                           :: sfc_speed,dtsavg
        integer(i_kind)                        :: error_status
        real(r_kind)                           :: clw_guess
        real(r_kind),dimension(:),allocatable :: tsim_clr 
        real(r_kind),dimension(:),allocatable :: h,q,prsl
        real(r_kind),dimension(:),allocatable :: prsi
        real(r_kind),dimension(:),allocatable :: tsim,emissivity,ts,emissivity_k
        real(r_kind),dimension(:,:),allocatable ::temp,ptau5,wmix 
        real(r_kind),dimension(:,:),allocatable ::jacobian
        ! varibales to write out
        real(r_kind),dimension(:,:),allocatable:: tb_cld,tb_clr,tb_obs
        real(r_kind),dimension(:,:),allocatable:: ca2019_obs
        real(r_kind),dimension(:,:),allocatable:: ca2014_obs
        real(r_kind),dimension(:,:),allocatable:: pwmax_obs
        real(r_kind),dimension(:,:,:),allocatable:: pwmax_grid
        real(r_kind),dimension(:,:,:),allocatable:: ca2014_grid,ca2019_grid
        real(r_kind),dimension(:),allocatable:: lonobs,latobs
        !
        real(r_kind)::ptau5deriv
        real(r_kind)::ptau5derivmax
        !
        ! init some variables
        init_pass = .TRUE.
        mype=0
        mype_diaghdr=0
        cold_start = .false.
        print*,myname_,"* INIT CRTM *lcloud_fwd:",lcloud_fwd
        call init_crtm(init_pass,mype_diaghdr,mype,                            &
             & nchanl,isis,obstype,subset_start,subset_end,                    &
             & nsig,                                                           &
             & lcloud_fwd,lallsky,cld_sea_only,n_actual_clouds,n_clouds_fwd,   &
             & n_clouds_jac,cloud_names,cloud_names_fwd,cloud_names_jac,       &
             & laerosol_fwd,laerosol,n_actual_aerosols,n_aerosols_fwd,         &
             & n_aerosols_jac,aerosol_names,aerosol_names_fwd,                 &
             & aerosol_names_jac,                                              &
             & n_ghg,ghg_names,                                                &
             & crtm_coeffs_path,                                               &
             & regional0,nvege_type0,cold_start )
!!!!!!!!!!!!
        print*,myname_,'*channelinfo: WMO_Satellite_ID :',channelinfo(1)%WMO_Satellite_ID
        print*,myname_,'*channelinfo: WMO_Sensor_ID    :',channelinfo(1)%WMO_Sensor_ID
        print*,myname_,'*channelinfo: n_Channels       :',channelinfo(1)%n_Channels
        !
        print*,myname_,'*channelinfo: Process_Channel  :',channelinfo(1)%Process_Channel
        print*,myname_,'*channelinfo: Sensor_Channel   :',channelinfo(1)%Sensor_Channel
        print*,myname_,'*channelinfo: Channel_Index    :',channelinfo(1)%Channel_Index
        !
        ! loop obs
        !
        call read_abiobsarray_from_file()
        dsfct = 0.
        !
        allocate(h(nsig),q(nsig),prsl(nsig))
        allocate(prsi(nsig+1))
        allocate(tsim_clr(nchanl))
        allocate(tsim(nchanl))
        allocate(emissivity(nchanl))
        allocate(ts(nchanl))
        allocate(emissivity_k(nchanl))
        allocate(temp(nsig,nchanl))
        allocate(ptau5(nsig,nchanl))
        allocate(wmix(nsig,nchanl))
        allocate(jacobian(nsigradjac,nchanl))
        allocate(tb_cld(nabiobs,nchanl))
        allocate(tb_clr(nabiobs,nchanl))
        allocate(tb_obs(nabiobs,nchanl))
        allocate(ca2014_obs(nchanl,nabiobs))
        allocate(ca2019_obs(nchanl,nabiobs))
        allocate(pwmax_obs(nchanl,nabiobs))
        allocate(ca2014_grid(lat2,lon2,nchanl))
        allocate(ca2019_grid(lat2,lon2,nchanl))
        allocate(pwmax_grid(lat2,lon2,nchanl))
        allocate(lonobs(nabiobs))
        allocate(latobs(nabiobs))
        tsim=0.
        tsim_clr=0.
!!$omp parallel do  schedule(dynamic,1) private(n,i)
        do n=1,nabiobs
        !do n=1,5
        ! init arrays
            call  call_crtm(obstype,iadate,data_s(:,n),nchanl,nreal, &
             &   mype, lon2,lat2,&            
             &   psges/1000.,tkges,qges ,prslges/1000.,prsiges/1000.,&
             &   u(:,:,1),v(:,:,1),tropprs, &
             &   isli2,sno2, dsfct,    &
             &   h,q,clw_guess,prsl,prsi, &
             &   trop5,tzbgr,dtsavg,sfc_speed,&
             &   tsim,emissivity,ptau5,ts, &
             &   emissivity_k,temp,wmix,jacobian,error_status,tsim_clr=tsim_clr, &
             &   qrges=qrges,qsges=qsges,qgges=qgges,qiges=qiges,                &
             &   qhges=qhges,qlges=qlges,ozges=ozges) 
        !     &   qhges=qhges,qlges=qlges) 
             !print*,'Simulated TB CLEAR: ',tsim_clr
             !print*,'Simulated TB ALL SKY: ',tsim
             !print*,"OBS TB      : ",data_s(nreal+1:nreal+nchanl,n)
             !print*,"Qvapor Sensitive :",minval(wmix),maxval(wmix)
             !print*,"temperature Sensitive :",minval(temp),maxval(temp)
             do i=1,nchanl
                 tb_cld(n,i) = tsim(i)
                 tb_clr(n,i) = tsim_clr(i)
                 tb_obs(n,i) = data_s(nreal+i,n)
             enddo
             !
             lonobs(n)=data_s(ilone,n)
             latobs(n)=data_s(ilate,n)
             ! maximum of weighting function is level at which transmittanc
             ! (ptau5) is changing the fastest.  This is used for the level
             !  assignment (needed for vertical localization).
             do i=1,nchanl
                 ptau5derivmax = -9.9e31_r_kind
                 pwmax_obs(i,n) = 0.0_r_kind
                         !kmax(i) = k
                 do k=2,nsig
                     ptau5deriv = abs( (ptau5(k-1,i)-ptau5(k,i))/ &
                         (log(prsl(k-1))-log(prsl(k))) )
                     if (ptau5deriv > ptau5derivmax) then
                         ptau5derivmax = ptau5deriv
                         !kmax(i) = k
                         pwmax_obs(i,n) = 10.0*prsi(k) ! cb to mb.
                     end if
                 end do
             end do ! nchanl
             
             print*,'--------------------------------'
             print*,n
             print*,'--------------------------------'
         enddo ! n=1,nabiobs
!!$omp barrier
         ! calc cloud effect patameters
         call calc_ca(tb_obs,tb_clr,tb_cld,ca2019_obs,ca2014_obs,nabiobs,nchanl)
         do i=1,nchanl
             if(i == 1)call cressman_lookup(latobs,lonobs,ca2014_obs(i,:),nabiobs,lat,lon,lat2,lon2,30.,-999)
             call interp_cressman(latobs,lonobs,ca2014_obs(i,:),nabiobs,lat,lon,ca2014_grid(:,:,i),lat2,lon2,30.,-999.)
             call interp_cressman(latobs,lonobs,ca2019_obs(i,:),nabiobs,lat,lon,ca2019_grid(:,:,i),lat2,lon2,30.,-999.)
             call interp_cressman(latobs,lonobs,pwmax_obs(i,:),nabiobs,lat,lon,pwmax_grid(:,:,i),lat2,lon2,30.,-999.)
         end do
         ! write out
         call write_tb(tb_cld,tb_clr,tb_obs,lonobs,latobs,lon,lat, & 
              ca2014_grid,ca2019_grid,pwmax_grid,nabiobs,nchanl,lat2,lon2)

         ! destory allocatable arrays
         deallocate(h,q,prsl)
         deallocate(prsi)
         deallocate(tsim_clr)
         deallocate(tsim)
         deallocate(emissivity)
         deallocate(ts)
         deallocate(emissivity_k)
         deallocate(temp)
         deallocate(ptau5)
         deallocate(wmix)
         deallocate(jacobian)
         deallocate(tb_cld)
         deallocate(tb_clr)
         deallocate(ca2014_obs)
         deallocate(ca2019_obs)
         deallocate(ca2014_grid)
         deallocate(ca2019_grid)
         deallocate(lonobs)
         deallocate(latobs)
         deallocate(pwmax_obs)
         deallocate(pwmax_grid)
         call destroy_crtm()

    end subroutine setuprad
    !
    subroutine calc_ca(yo,hclr,h,ca_2019,ca_2014,nobs,nc)
        !
        !  ca_2019 = Co - Cm
        !  ca_2014 = (Co+Cm)/2.
        !
        implicit none
        integer(i_kind),intent(in   ):: nobs,nc
        real(r_kind),dimension(nobs,nc),intent(in   )::yo
        real(r_kind),dimension(nobs,nc),intent(in   )::hclr
        real(r_kind),dimension(nobs,nc),intent(in   )::h
        real(r_kind),dimension(nc,nobs),intent(inout)::ca_2019
        real(r_kind),dimension(nc,nobs),intent(inout)::ca_2014
        !
        !local variables
        integer(i_kind)::n,i
        real(r_kind)::co,cm
!$omp parallel do  schedule(dynamic,1) private(n,i,co,cm)
        do n=1,nobs
            do i=1,nc
                co = abs(yo(n,i)-hclr(n,i))
                cm = abs(h(n,i)-hclr(n,i))
                ca_2019(i,n)=co-cm
                ca_2014(i,n)=(co+cm)*0.5
            enddo
        enddo
    end subroutine calc_ca
    !
    subroutine write_tb(tb_cld,tb_clr,tb_obs,lonobs,latobs,lon,lat, &
            ca2014,ca2019,pwmax,nobs,nc,n1,n2)
        use netcdf, only: nf90_create,nf90_close
        use netcdf, only: nf90_clobber
        use netcdf, only: nf90_put_att
        use netcdf, only: nf90_put_var
        use netcdf, only: nf90_global
        use netcdf, only: nf90_real
        use netcdf, only: nf90_enddef
        use netcdf, only: nf90_noerr
        use netcdf, only: nf90_strerror
        use netcdf, only: nf90_def_dim
        use netcdf, only: nf90_def_var
        implicit none
        integer(i_kind),intent(in   )                :: nobs
        integer(i_kind),intent(in   )                :: nc
        integer(i_kind),intent(in   )                :: n1,n2
        real(r_kind),dimension(nobs,nc),intent(in   ):: tb_cld
        real(r_kind),dimension(nobs,nc),intent(in   ):: tb_clr
        real(r_kind),dimension(nobs,nc),intent(in   ):: tb_obs
        real(r_kind),dimension(nobs   ),intent(in   ):: lonobs
        real(r_kind),dimension(nobs   ),intent(in   ):: latobs
        real(r_kind),dimension(n1,n2  ),intent(in   ):: lon
        real(r_kind),dimension(n1,n2  ),intent(in   ):: lat
        real(r_kind),dimension(n1,n2,nc),intent(in  ):: ca2014
        real(r_kind),dimension(n1,n2,nc),intent(in  ):: ca2019
        real(r_kind),dimension(n1,n2,nc),intent(in  ):: pwmax
        !
        character (len = *), parameter :: EW_NAME = "east_west"
        character (len = *), parameter :: SN_NAME = "south_north"
        character (len = *), parameter :: BTP_NAME = "bottom_top"
        character (len = *), parameter :: NOBS_NAME = "nobs"
        character (len = *), parameter :: CHANNEL_NAME = "channel"
        character (len = *), parameter :: TBCLD_NAME = "tb_cld"
        character (len = *), parameter :: TBCLR_NAME = "tb_clr"
        character (len = *), parameter :: LONOBS_NAME = "lonobs"
        character (len = *), parameter :: LATOBS_NAME = "latobs"
        character (len = *), parameter :: LON_NAME = "longrid"
        character (len = *), parameter :: LAT_NAME = "latgrid"
        character (len = *), parameter :: TBOBS_NAME = "tb_obs"
        character (len = *), parameter :: CA2019_NAME = "ca_2019"
        character (len = *), parameter :: CA2014_NAME = "ca_2014"
        character (len = *), parameter :: PWMAX_NAME = "pwmax"
        character(len=500):: filename
        !
        integer :: ncid
        integer :: longrid_varid, latgrid_varid
        integer :: lonobs_varid, latobs_varid
        integer :: tbcld_varid, tbclr_varid,tbobs_varid
        integer :: ca2014_varid,temp_varid(22)
        integer :: ca2019_varid
        integer :: pwmax_varid
        integer :: ew_dimid, sn_dimid,btp_dimid
        integer :: ch_dimid,nobs_dimid
        !
        integer::e_w,s_n
        !
        integer :: dimids1(2)
        integer :: dimids2(2)
        integer :: dimids3(3)
        integer :: dimids4(3)
        integer :: count0(2),count3(3)
        integer :: start(2),start3(3)
        !
        ! Loop indices
        integer :: i,j,k
        !
        filename = trim(datapath)//'/'//trim(tbfout)
        !
        e_w = n2
        s_n = n1
        ! Create the file.
        call nc_check( nf90_create(trim(filename), nf90_clobber, ncid) )
        ! Define the dimensions. 
        call nc_check( nf90_def_dim(ncid, EW_NAME, e_w, ew_dimid) )
        call nc_check( nf90_def_dim(ncid, SN_NAME, s_n, sn_dimid) )
        call nc_check( nf90_def_dim(ncid, BTP_NAME, nsig, btp_dimid) )
        call nc_check( nf90_def_dim(ncid, NOBS_NAME,nobs, nobs_dimid) )
        call nc_check( nf90_def_dim(ncid, CHANNEL_NAME, nc, ch_dimid) )
        dimids1 = (/ ew_dimid,sn_dimid/)
        dimids2 = (/ nobs_dimid, ch_dimid /)
        dimids3 = (/ ew_dimid,sn_dimid,ch_dimid/)
        if(ldiagout)then
            dimids4 = (/ ew_dimid,sn_dimid,btp_dimid/)
        endif
        !! Define the netCDF variables
        call nc_check( nf90_def_var(ncid, LON_NAME, NF90_REAL, & 
             dimids1,longrid_varid))
        call nc_check( nf90_def_var(ncid, LAT_NAME, NF90_REAL, & 
             dimids1,latgrid_varid))
        call nc_check( nf90_def_var(ncid, LONOBS_NAME, NF90_REAL, & 
             dimids2(1),lonobs_varid))
        call nc_check( nf90_def_var(ncid, LATOBS_NAME, NF90_REAL, & 
             dimids2(1),latobs_varid))
        call nc_check( nf90_def_var(ncid, TBCLD_NAME, NF90_REAL, & 
             dimids2,tbcld_varid))
        call nc_check( nf90_def_var(ncid, TBCLR_NAME, NF90_REAL, & 
             dimids2,tbclr_varid))
        call nc_check( nf90_def_var(ncid, TBOBS_NAME, NF90_REAL, & 
             dimids2,tbobs_varid))
        call nc_check( nf90_def_var(ncid, CA2014_NAME, NF90_REAL, & 
             dimids3,ca2014_varid))
        call nc_check( nf90_def_var(ncid, CA2019_NAME, NF90_REAL, & 
             dimids3,ca2019_varid))
        call nc_check( nf90_def_var(ncid, PWMAX_NAME, NF90_REAL, & 
             dimids3,pwmax_varid))
        if(ldiagout)then
            !1: u,u-wind;
            !2: v,v-wind;
            !3: t,potential temperature;
            !4: qv,qvapor;
            !5: qr, qrain;
            !6: qs, qsnow;
            !7: qg, qgraup;
            !8: qh, qhail;
            !9: qc, qcloud;
            !10: qi, qice;
            !11: pmid, mid layer pressure;
            !12: hgt, geopotential height;
            call nc_check( nf90_def_var(ncid, "u", NF90_REAL, &
             dimids4,temp_varid(1)))
            call nc_check( nf90_put_att(ncid,temp_varid(1), 'units', 'm s-1') )
            call nc_check( nf90_put_att(ncid,temp_varid(1), 'name', 'u-wind') )
            !
            call nc_check( nf90_def_var(ncid, "v", NF90_REAL, &
             dimids4,temp_varid(2)))
            call nc_check( nf90_put_att(ncid,temp_varid(2), 'units', 'm s-1') )
            call nc_check( nf90_put_att(ncid,temp_varid(2), 'name', 'v-wind') )
            !
            call nc_check( nf90_def_var(ncid, "t", NF90_REAL, &
             dimids4,temp_varid(3)))
            call nc_check( nf90_put_att(ncid,temp_varid(3), 'units', 'K') )
            call nc_check( nf90_put_att(ncid,temp_varid(3), 'name', 'potential temperature') )
            !
            call nc_check( nf90_def_var(ncid, "qv", NF90_REAL, &
             dimids4,temp_varid(4)))
            call nc_check( nf90_put_att(ncid,temp_varid(4), 'units', 'kg kg-1') )
            call nc_check( nf90_put_att(ncid,temp_varid(4), 'name', 'Water vapor mixing ratio') )
            !
            call nc_check( nf90_def_var(ncid, "qr", NF90_REAL, &
             dimids4,temp_varid(5)))
            call nc_check( nf90_put_att(ncid,temp_varid(5), 'units', 'kg kg-1') )
            call nc_check( nf90_put_att(ncid,temp_varid(5), 'name', 'Rain water mixing ratio') )
            !
            call nc_check( nf90_def_var(ncid, "qs", NF90_REAL, &
             dimids4,temp_varid(6)))
            call nc_check( nf90_put_att(ncid,temp_varid(6), 'units', 'kg kg-1') )
            call nc_check( nf90_put_att(ncid,temp_varid(6), 'name', 'SNOW mixing ratio') )
            !
            call nc_check( nf90_def_var(ncid, "qg", NF90_REAL, &
             dimids4,temp_varid(7)))
            call nc_check( nf90_put_att(ncid,temp_varid(7), 'units', 'kg kg-1') )
            call nc_check( nf90_put_att(ncid,temp_varid(7), 'name', 'Graupel mixing ratio') )
            !
            call nc_check( nf90_def_var(ncid, "qh", NF90_REAL, &
             dimids4,temp_varid(8)))
            call nc_check( nf90_put_att(ncid,temp_varid(8), 'units', 'kg kg-1') )
            call nc_check( nf90_put_att(ncid,temp_varid(8), 'name', 'HAIL mixing ratio') )
            !
            call nc_check( nf90_def_var(ncid, "qc", NF90_REAL, &
             dimids4,temp_varid(9)))
            call nc_check( nf90_put_att(ncid,temp_varid(9), 'units', 'kg kg-1') )
            call nc_check( nf90_put_att(ncid,temp_varid(9), 'name', ' cloud water mixing ratio') )
            !
            call nc_check( nf90_def_var(ncid, "qi", NF90_REAL, &
             dimids4,temp_varid(10)))
            call nc_check( nf90_put_att(ncid,temp_varid(10), 'units', 'kg kg-1') )
            call nc_check( nf90_put_att(ncid,temp_varid(10), 'name', 'ice mixing ratio') )
            !
            call nc_check( nf90_def_var(ncid, "pmid", NF90_REAL, &
             dimids4,temp_varid(11)))
            call nc_check( nf90_put_att(ncid,temp_varid(11), 'units', 'Pa') )
            call nc_check( nf90_put_att(ncid,temp_varid(11), 'name', 'pressure') )
            !
            call nc_check( nf90_def_var(ncid, "hgt", NF90_REAL, &
             dimids4,temp_varid(12)))
            call nc_check( nf90_put_att(ncid,temp_varid(12), 'units', 'm2 s-2') )
            call nc_check( nf90_put_att(ncid,temp_varid(12), 'name', 'geopotential') )
        endif
        ! Assign units attributes to the netCDF variables.
        call nc_check( nf90_put_att(ncid, latgrid_varid, 'units', 'degrees_north') )
        call nc_check( nf90_put_att(ncid, longrid_varid, 'units', 'degrees_east') )
        call nc_check( nf90_put_att(ncid, latobs_varid, 'units', 'degrees_north') )
        call nc_check( nf90_put_att(ncid, lonobs_varid, 'units', 'degrees_east') )
        call nc_check( nf90_put_att(ncid, tbcld_varid, 'units', 'K') )
        call nc_check( nf90_put_att(ncid, tbclr_varid, 'units', 'K') )
        call nc_check( nf90_put_att(ncid, tbobs_varid, 'units', 'K') )
        call nc_check( nf90_put_att(ncid, ca2014_varid, 'units', 'K') )
        call nc_check( nf90_put_att(ncid, ca2019_varid, 'units', 'K') )
        call nc_check( nf90_put_att(ncid, pwmax_varid, 'units', 'hPa') )
        call nc_check( nf90_put_att(ncid, ca2014_varid, 'desc', 'Co=abs(Yo-hclr);Cm=abs(h-hclr);Ca=Co-Cm') )
        call nc_check( nf90_put_att(ncid, ca2019_varid, 'desc', 'Co=abs(Yo-hclr);Cm=abs(h-hclr);Ca=(Co+Cm)/2.') )
        call nc_check( nf90_put_att(ncid, pwmax_varid, 'desc', 'weight max level') )
        !! End define mode.
        call nc_check( nf90_enddef(ncid) )
        !
        count0 = (/e_w,s_n/)
        start = (/ 1, 1  /)
        ! Write the pretend data.
        call nc_check( nf90_put_var(ncid, longrid_varid, lon, start, count0) )
        call nc_check( nf90_put_var(ncid, latgrid_varid, lat, start, count0) )
        !
        count0 = (/nobs,nc/)
        call nc_check( nf90_put_var(ncid, tbcld_varid, tb_cld, start, count0) )
        call nc_check( nf90_put_var(ncid, tbclr_varid, tb_clr, start, count0) )
        call nc_check( nf90_put_var(ncid, tbobs_varid, tb_obs, start, count0) )
        !
        call nc_check( nf90_put_var(ncid, lonobs_varid, lonobs))
        call nc_check( nf90_put_var(ncid, latobs_varid, latobs ))
        !
        count3 = (/e_w,s_n,nc/)
        start3= (/ 1, 1 ,1 /)
        call nc_check( nf90_put_var(ncid, ca2019_varid, ca2019 ))
        call nc_check( nf90_put_var(ncid, ca2014_varid, ca2014 ))
        call nc_check( nf90_put_var(ncid, pwmax_varid, pwmax ))
        ! if diagout
        if(ldiagout)then
            count3 = (/e_w,s_n,nsig/)
            start3= (/ 1, 1 ,1 /)
            ! Define the netCDF variables
            ! Assign units attributes to the netCDF variables.
            ! Write the pretend data.
            !1: u,u-wind;
            !2: v,v-wind;
            !3: t,potential temperature;
            !4: qv,qvapor;
            !5: qr, qrain;
            !6: qs, qsnow;
            !7: qg, qgraup;
            !8: qh, qhail;
            !9: qc, qcloud;
            !10: qi, qice;
            !11: pmid, mid layer pressure;
            !12: hgt, geopotential height;
            call nc_check( nf90_put_var(ncid, temp_varid(1), u,     start3, count3) )
            call nc_check( nf90_put_var(ncid, temp_varid(2), v,     start3, count3) )
            call nc_check( nf90_put_var(ncid, temp_varid(3), t,     start3, count3) )
            call nc_check( nf90_put_var(ncid, temp_varid(4), qv,    start3, count3) )
            call nc_check( nf90_put_var(ncid, temp_varid(5), qrges, start3, count3) )
            call nc_check( nf90_put_var(ncid, temp_varid(6), qsges, start3, count3) )
            call nc_check( nf90_put_var(ncid, temp_varid(7), qgges, start3, count3) )
            call nc_check( nf90_put_var(ncid, temp_varid(8), qhges, start3, count3) )
            call nc_check( nf90_put_var(ncid, temp_varid(9), qlges, start3, count3) )
            call nc_check( nf90_put_var(ncid, temp_varid(10), qiges,start3, count3) )
            call nc_check( nf90_put_var(ncid, temp_varid(11), pmid, start3, count3) )
            call nc_check( nf90_put_var(ncid, temp_varid(12), hgt,  start3, count3) )
        end if
        !
        call nc_check( nf90_close(ncid) )

    contains
        subroutine nc_check(status)
            integer, intent ( in) :: status
            if(status /= nf90_noerr) then
                print *, trim(nf90_strerror(status))
                stop "Stopped"
            endif
        end subroutine nc_check


    end subroutine write_tb




end module module_abi
