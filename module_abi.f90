!4/24/2021 
!
module module_abi
    use model_precision,only:i_kind=>INT32,r_kind=>P,P
    use crtm_interface,only:init_crtm
    use crtm_interface,only:channelinfo
    use parameters_define,only:lcloud_fwd,lallsky,cld_sea_only
    use parameters_define,only:n_actual_clouds,n_clouds_fwd,n_clouds_jac
    use parameters_define,only:cloud_names,cloud_names_fwd
    use parameters_define,only:laerosol_fwd,laerosol,n_actual_aerosols
    use parameters_define,only:n_aerosols_fwd,n_aerosols_jac
    use parameters_define,only:aerosol_names,aerosol_names_fwd
    use parameters_define,only:n_ghg,ghg_names
    use parameters_define,only:obstype
    use parameters_define,only:isis,nchanl,subset_start,subset_end
    use parameters_define,only:crtm_coeffs_path
    use parameters_define,only:nvege_type0=>nvege_type,regional0=>regional

    implicit none
    real(P),allocatable,dimension(:,:)::tb
    real(P),allocatable,dimension(:,:)::tb_clr

    public::setuprad 
    contains

        subroutine setuprad()
            implicit none
            character(len=*),parameter::myname_="setuprad"
            logical  :: init_pass
            integer(i_kind) :: mype_diaghdr,mype
            integer(i_kind) :: nsig
            integer(i_kind) :: msig
            !
            !
            msig=50
            nsig=msig
            init_pass = .TRUE.
            mype=0
            mype_diaghdr=0

            print*,myname_,"* INIT CRTM *"
            call init_crtm(init_pass,mype_diaghdr,mype,               &
	& nchanl,isis,obstype,subset_start,subset_end,                    &
	& msig,nsig,                                                      &
	& lcloud_fwd,lallsky,cld_sea_only,n_actual_clouds,n_clouds_fwd,   &
	& n_clouds_jac,cloud_names,cloud_names_fwd,                       &
	& laerosol_fwd,laerosol,n_actual_aerosols,n_aerosols_fwd,         &
	& n_aerosols_jac,aerosol_names,aerosol_names_fwd,                 &
	& n_ghg,ghg_names,                                                &
	& crtm_coeffs_path,                                               &
	& regional0,nvege_type0)
            
             print*,myname_,'*channelinfo: WMO_Satellite_ID :',channelinfo(1)%WMO_Satellite_ID
             print*,myname_,'*channelinfo: WMO_Sensor_ID    :',channelinfo(1)%WMO_Sensor_ID
             print*,myname_,'*channelinfo: n_Channels       :',channelinfo(1)%n_Channels
             !
             print*,myname_,'*channelinfo: Process_Channel  :',channelinfo(1)%Process_Channel
             print*,myname_,'*channelinfo: Sensor_Channel   :',channelinfo(1)%Sensor_Channel
             print*,myname_,'*channelinfo: Channel_Index    :',channelinfo(1)%Channel_Index



        end subroutine setuprad




end module module_abi
