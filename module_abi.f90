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
    !
    use goesabi_obs,only: data_s=>data_obsabi
    use goesabi_obs,only: read_abiobsarray_from_file
    use goesabi_obs,only: nreal,nabiobs
    !

    implicit none

    public::setuprad 
    contains

        subroutine setuprad()
            implicit none
            character(len=*),parameter::myname_="setuprad"
            logical  :: init_pass
            integer(i_kind) :: mype_diaghdr,mype
            logical::cold_start 
            character(10)                          :: obstype
            !
            integer(i_kind)::n
            real(r_kind)::obstime
            !
            real(r_kind),dimension(lat2,lon2)::dsfct
            !
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
            !
            !
            init_pass = .TRUE.
            mype=0
            mype_diaghdr=0
            cold_start = .false.

            print*,myname_,"* INIT CRTM *lcloud_fwd:",lcloud_fwd
            call init_crtm(init_pass,mype_diaghdr,mype,               &
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
             tsim=0.
             tsim_clr=0.
             do n=300,300
                 ! init arrays

                print*,"OBS : ",data_s(:,n)
                print*,iadate,lon2,lat2
                print*,nchanl

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
                &   qhges=qhges,qlges=qlges) 
                    print*,'Simulated TB CLEAR: '
                    print*,'Simulated TB CLEAR: ',tsim_clr
                    print*,'Simulated TB ALL SKY: ',tsim
                    print*,"OBS TB      : ",data_s(nreal+1:nreal+nchanl,n)
                    print*,"Qvapor Sensitive :",minval(wmix),maxval(wmix)
                    print*,"temperature Sensitive :",minval(temp),maxval(temp)

            enddo
            !
            print*,iqv,icw
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
            call destroy_crtm()

        end subroutine setuprad




end module module_abi
