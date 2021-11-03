!
! 11/1/2021 Huanhuan ZHANG
!
! Adapted from GSI/src/crtm_interface.f90
! Remove dependence of GSI MODULE
!
!
!$$$ module documentation block
!           .      .    .                                       
! module:   crtm_interface module for setuprad. Calculates profile and calls crtm
!  prgmmr:
!
! abstract: crtm_interface module for setuprad. Initializes CRTM, Calculates profile and 
!         calls CRTM and destroys initialization
!
!
!$$$ end documentation block

module crtm_interface

!use kinds,only: r_kind,i_kind,r_single
use model_precision,only:r_kind=>SP,i_kind=>INT32,r_single=>P
use crtm_module, only: crtm_atmosphere_type,crtm_surface_type,crtm_geometry_type, &
    crtm_options_type,crtm_rtsolution_type,crtm_destroy,crtm_options_destroy, &
    crtm_options_create,crtm_options_associated,success,crtm_atmosphere_create, &
    crtm_surface_create,crtm_k_matrix,crtm_forward, &   
    ssu_input_setvalue, &
    crtm_channelinfo_type, &
    crtm_surface_destroy, crtm_surface_associated, crtm_surface_zero, &
    crtm_atmosphere_associated, &
    crtm_atmosphere_destroy,crtm_atmosphere_zero, &
    crtm_rtsolution_type, crtm_rtsolution_create, &
    crtm_rtsolution_destroy, crtm_rtsolution_associated, &
    crtm_irlandcoeff_classification, &
    crtm_kind => fp, &
    crtm_microwave_sensor => microwave_sensor,&
    RT_SOI
use crtm_aod_module, only: crtm_aod_k

implicit none

private
public init_crtm            ! Subroutine initializes crtm for specified instrument
!public call_crtm            ! Subroutine creates profile for crtm, calls crtm, then adjoint of create
public destroy_crtm         ! Subroutine destroys initialization for crtm
public sensorindex
public surface
public channelinfo
public isatid               ! = 1  index of satellite id
public itime                ! = 2  index of analysis relative obs time
public ilon                 ! = 3  index of grid relative obs location (x)
public ilat                 ! = 4  index of grid relative obs location (y)
public ilzen_ang            ! = 5  index of local (satellite) zenith angle (radians)
public ilazi_ang            ! = 6  index of local (satellite) azimuth angle (radians)
public iscan_ang            ! = 7  index of scan (look) angle (radians)
public iscan_pos            ! = 8  index of integer scan position
public iszen_ang            ! = 9  index of solar zenith angle (degrees)
public isazi_ang            ! = 10 index of solar azimuth angle (degrees)
public ifrac_sea            ! = 11 index of ocean percentage
public ifrac_lnd            ! = 12 index of land percentage
public ifrac_ice            ! = 13 index of ice percentage
public ifrac_sno            ! = 14 index of snow percentage
public its_sea              ! = 15 index of ocean temperature
public its_lnd              ! = 16 index of land temperature
public its_ice              ! = 17 index of ice temperature
public its_sno              ! = 18 index of snow temperature
public itsavg               ! = 19 index of average temperature
public ivty                 ! = 20 index of vegetation type
public ivfr                 ! = 21 index of vegetation fraction
public isty                 ! = 22 index of soil type
public istp                 ! = 23 index of soil temperature
public ism                  ! = 24 index of soil moisture
public isn                  ! = 25 index of snow depth
public izz                  ! = 26 index of surface height
public idomsfc              ! = 27 index of dominate surface type
public isfcr                ! = 28 index of surface roughness
public iff10                ! = 29 index of ten meter wind factor
public ilone                ! = 30 index of earth relative longitude (degrees)
public ilate                ! = 31 index of earth relative latitude (degrees)
public iclr_sky             ! = 7  index of clear sky amount (goes_img, seviri)
public isst_navy            ! = 7  index of navy sst retrieval (K) (avhrr_navy)
public idata_type           ! = 32 index of data type (151=day, 152=night, avhrr_navy)
public iclavr               ! = 32 index of clavr cloud flag (avhrr)
public isst_hires           ! = 33 index of interpolated hires sst
public itref                ! = 34/36 index of Tr
public idtw                 ! = 35/37 index of d(Tw)
public idtc                 ! = 36/38 index of d(Tc)
public itz_tr               ! = 37/39 index of d(Tz)/d(Tr)

!  Note other module variables are only used within this routine

  character(len=*), parameter :: myname='crtm_interface'
  
  ! Indices for the CRTM NPOESS EmisCoeff file
  integer(i_kind), parameter :: INVALID_LAND = 0
  integer(i_kind), parameter :: COMPACTED_SOIL = 1
  integer(i_kind), parameter :: TILLED_SOIL = 2
  integer(i_kind), parameter :: IRRIGATED_LOW_VEGETATION = 5
  integer(i_kind), parameter :: MEADOW_GRASS = 6
  integer(i_kind), parameter :: SCRUB = 7
  integer(i_kind), parameter :: BROADLEAF_FOREST = 8
  integer(i_kind), parameter :: PINE_FOREST = 9
  integer(i_kind), parameter :: TUNDRA = 10
  integer(i_kind), parameter :: GRASS_SOIL = 11
  integer(i_kind), parameter :: BROADLEAF_PINE_FOREST = 12
  integer(i_kind), parameter :: GRASS_SCRUB = 13
  integer(i_kind), parameter :: URBAN_CONCRETE = 15
  integer(i_kind), parameter :: BROADLEAF_BRUSH = 17
  integer(i_kind), parameter :: WET_SOIL = 18
  integer(i_kind), parameter :: SCRUB_SOIL = 19

  real(r_kind)   , save ,allocatable,dimension(:,:) :: aero         ! aerosol (guess) profiles at obs location
  real(r_kind)   , save ,allocatable,dimension(:,:) :: aero_conc    ! aerosol (guess) concentrations at obs location
  real(r_kind)   , save ,allocatable,dimension(:)   :: auxrh        ! temporary array for rh profile as seen by CRTM

  character(len=20),save,allocatable,dimension(:)   :: ghg_names    ! names of green-house gases

  integer(i_kind), save ,allocatable,dimension(:)   :: icloud       ! cloud index for those considered here 
  integer(i_kind), save ,allocatable,dimension(:)   :: jcloud       ! cloud index for those fed to CRTM
  real(r_kind)   , save ,allocatable,dimension(:,:) :: cloud        ! cloud considered here
  real(r_kind)   , save ,allocatable,dimension(:,:) :: cloudefr     ! effective radius of cloud type in CRTM
  real(r_kind)   , save ,allocatable,dimension(:,:) :: cloud_cont   ! cloud content fed into CRTM 
  real(r_kind)   , save ,allocatable,dimension(:,:) :: cloud_efr    ! effective radius of cloud type in CRTM

  real(r_kind)   , save ,allocatable,dimension(:,:,:,:)  :: gesqsat ! qsat to calc rh for aero particle size estimate
  real(r_kind)   , save ,allocatable,dimension(:)  :: lcloud4crtm_wk ! cloud info usage index for each channel

  integer(i_kind),save, allocatable,dimension(:) :: map_to_crtm_ir
  integer(i_kind),save, allocatable,dimension(:) :: map_to_crtm_mwave 
  integer(i_kind),save, allocatable,dimension(:) :: icw
  integer(i_kind),save, allocatable,dimension(:) :: iaero_jac
  integer(i_kind),save :: isatid,itime,ilon,ilat,ilzen_ang,ilazi_ang,iscan_ang
  integer(i_kind),save :: iscan_pos,iszen_ang,isazi_ang,ifrac_sea,ifrac_lnd,ifrac_ice
  integer(i_kind),save :: ifrac_sno,its_sea,its_lnd,its_ice,its_sno,itsavg
  integer(i_kind),save :: ivty,ivfr,isty,istp,ism,isn,izz,idomsfc,isfcr,iff10,ilone,ilate
  integer(i_kind),save :: iclr_sky,isst_navy,idata_type,isst_hires,iclavr
  integer(i_kind),save :: itref,idtw,idtc,itz_tr,istype
  integer(i_kind),save :: sensorindex
  integer(i_kind),save :: ico2,ico24crtm
  integer(i_kind),save :: n_actual_aerosols_wk      ! number of aerosols considered
  integer(i_kind),save :: n_aerosols_fwd_wk         ! number of aerosols considered
  integer(i_kind),save :: n_aerosols_jac_wk         ! number of aerosols considered
  integer(i_kind),save :: n_actual_clouds_wk        ! number of clouds considered
  integer(i_kind),save :: n_clouds_fwd_wk           ! number of clouds considered
  integer(i_kind),save :: n_clouds_jac_wk           ! number of clouds considered
  integer(i_kind),save :: n_ghg              ! number of green-house gases
  integer(i_kind),save :: itv,iqv,ioz,ius,ivs,isst
  integer(i_kind),save :: indx_p25, indx_dust1, indx_dust2
  logical        ,save :: lwind
  logical        ,save :: cld_sea_only_wk
  logical        ,save :: mixed_use

  logical        ,save :: regional
  integer(i_kind),save :: nvege_type

  integer(i_kind), parameter :: min_n_absorbers = 2

  type(crtm_atmosphere_type),save,dimension(1)   :: atmosphere
  type(crtm_surface_type),save,dimension(1)      :: surface
  type(crtm_geometry_type),save,dimension(1)     :: geometryinfo
  type(crtm_options_type),save,dimension(1)      :: options
  type(crtm_channelinfo_type),save,dimension(1)  :: channelinfo



  type(crtm_atmosphere_type),save,allocatable,dimension(:,:):: atmosphere_k
  type(crtm_atmosphere_type),save,allocatable,dimension(:,:):: atmosphere_k_clr
  type(crtm_surface_type),save,allocatable,dimension(:,:):: surface_k
  type(crtm_surface_type),save,allocatable,dimension(:,:):: surface_k_clr
  type(crtm_rtsolution_type),save,allocatable,dimension(:,:):: rtsolution
  type(crtm_rtsolution_type),save,allocatable,dimension(:,:):: rtsolution0              
  type(crtm_rtsolution_type),save,allocatable,dimension(:,:):: rtsolution_clr              
  type(crtm_rtsolution_type),save,allocatable,dimension(:,:):: rtsolution_k
  type(crtm_rtsolution_type),save,allocatable,dimension(:,:):: rtsolution_k_clr

! Mapping land surface type of GFS to CRTM
!  Notes: index 0 is water, and index 13 is ice. The two indices are not
!         used and just assigned to COMPACTED_SOIL. Also, since there
!         is currently one relevant mapping for the global we apply
!         'crtm' in the naming convention.  
  integer(i_kind), parameter, dimension(0:13) :: gfs_to_crtm=(/COMPACTED_SOIL, &
     BROADLEAF_FOREST, BROADLEAF_FOREST, BROADLEAF_PINE_FOREST, PINE_FOREST, &
     PINE_FOREST, BROADLEAF_BRUSH, SCRUB, SCRUB, SCRUB_SOIL, TUNDRA, &
     COMPACTED_SOIL, TILLED_SOIL, COMPACTED_SOIL/)
! Mapping surface classification to CRTM
  integer(i_kind), parameter :: USGS_N_TYPES = 24
  integer(i_kind), parameter :: IGBP_N_TYPES = 20
  integer(i_kind), parameter :: GFS_N_TYPES = 13
  integer(i_kind), parameter :: SOIL_N_TYPES = 16
  integer(i_kind), parameter :: GFS_SOIL_N_TYPES = 9
  integer(i_kind), parameter :: GFS_VEGETATION_N_TYPES = 13
  integer(i_kind), parameter, dimension(1:USGS_N_TYPES) :: usgs_to_npoess=(/URBAN_CONCRETE, &
     COMPACTED_SOIL, IRRIGATED_LOW_VEGETATION, GRASS_SOIL, MEADOW_GRASS, &
     MEADOW_GRASS, MEADOW_GRASS, SCRUB, GRASS_SCRUB, MEADOW_GRASS, &
     BROADLEAF_FOREST, PINE_FOREST, BROADLEAF_FOREST, PINE_FOREST, &
     BROADLEAF_PINE_FOREST, INVALID_LAND, WET_SOIL, WET_SOIL, &
     IRRIGATED_LOW_VEGETATION, TUNDRA, TUNDRA, TUNDRA, TUNDRA, &
     INVALID_LAND/)
  integer(i_kind), parameter, dimension(1:IGBP_N_TYPES) :: igbp_to_npoess=(/PINE_FOREST, &
    BROADLEAF_FOREST, PINE_FOREST, BROADLEAF_FOREST, BROADLEAF_PINE_FOREST, &
    SCRUB, SCRUB_SOIL, BROADLEAF_BRUSH, BROADLEAF_BRUSH, SCRUB, BROADLEAF_BRUSH, &
    TILLED_SOIL, URBAN_CONCRETE, TILLED_SOIL, INVALID_LAND, COMPACTED_SOIL, &
    INVALID_LAND, TUNDRA, TUNDRA, TUNDRA/)
  integer(i_kind), parameter, dimension(1:USGS_N_TYPES) :: usgs_to_usgs=(/1, &
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, &
    20, 21, 22, 23, 24/)
  integer(i_kind), parameter, dimension(1:IGBP_N_TYPES) :: igbp_to_igbp=(/1, &
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, &
    20/)
  integer(i_kind), parameter, dimension(1:IGBP_N_TYPES) :: igbp_to_gfs=(/4, &
    1, 5, 2, 3, 8, 9, 6, 6, 7, 8, 12, 7, 12, 13, 11, 0, 10, 10, 11/)
  integer(i_kind), parameter, dimension(1:USGS_N_TYPES) :: usgs_to_gfs=(/7, &
    12, 12, 12, 12, 12, 7, 9, 8, 6, 2, 5, 1, 4, 3, 0, 8, 8, 11, 10, 10, &
    10, 11, 13/)
 ! Mapping soil types to CRTM
 ! The CRTM soil types for microwave calculations are based on the 
 ! GFS use of the 9 category Zobler dataset. The regional soil types
 ! are based on a 16 category representation of FAO/STATSGO. 
  integer(i_kind), parameter, dimension(1:SOIL_N_TYPES) :: map_soil_to_crtm=(/1, &
    1, 4, 2, 2, 8, 7, 2, 6, 5, 2, 3, 8, 1, 6, 9/)

contains

subroutine init_crtm(init_pass,mype_diaghdr,mype,                     &
	& nchanl,isis,obstype,subset_start,subset_end,                    &
	& msig,nsig,                                                      &
	& lcloud_fwd,lallsky,cld_sea_only,n_actual_clouds,n_clouds_fwd,   &
	& n_clouds_jac,cloud_names,cloud_names_fwd,                       &
	& laerosol_fwd,laerosol,n_actual_aerosols,n_aerosols_fwd,         &
	& n_aerosols_jac,aerosol_names,aerosol_names_fwd,                 &
	& n_ghg,ghg_names,                                                &
	& crtm_coeffs_path,                                               &
	& regional0,nvege_type0)

    !msig: n_layers
    !

  use crtm_module, only: mass_mixing_ratio_units,co2_id,o3_id,crtm_init, &
      crtm_channelinfo_subset, crtm_channelinfo_n_channels, toa_pressure,max_n_layers, &
      volume_mixing_ratio_units,h2o_id,ch4_id,n2o_id,co_id

  implicit none

! argument 
  logical        ,intent(in) :: init_pass
  integer(i_kind),intent(in) :: nchanl,mype_diaghdr,mype
  character(20)  ,intent(in) :: isis
  character(10)  ,intent(in) :: obstype
  integer(i_kind),intent(in) :: subset_start,subset_end
  integer(i_kind),intent(in) :: msig,nsig
  logical        ,intent(in) :: lcloud_fwd,lallsky,cld_sea_only
  integer(i_kind),intent(in) :: n_actual_clouds,n_clouds_fwd
  integer(i_kind),intent(in) :: n_clouds_jac
  character(*)  ,dimension(6),intent(in) :: cloud_names,cloud_names_fwd
  logical        ,intent(in) :: laerosol_fwd,laerosol
  integer(i_kind),intent(in) :: n_actual_aerosols,n_aerosols_fwd
  integer(i_kind),intent(in) :: n_aerosols_jac
  character(*)  ,dimension(3),intent(in) :: aerosol_names,aerosol_names_fwd
  integer(i_kind),intent(in) :: n_ghg
  character(*)  ,dimension(5),intent(in) :: ghg_names
  character(200) ,intent(in) :: crtm_coeffs_path
  logical,        intent(in) :: regional0
  integer(i_kind),intent(in) :: nvege_type0


! local parameters
  character(len=*), parameter :: myname_=myname//'*init_crtm'

  real(r_kind)::zero=0.0

! local variables
  integer(i_kind) :: ier,ii,error_status,iderivative
  integer(i_kind) :: k
  logical :: ice,Load_AerosolCoeff,Load_CloudCoeff
  character(len=20),dimension(1) :: sensorlist
  integer(i_kind) :: indx,iii,icloud4crtm
! ...all "additional absorber" variables
  integer(i_kind) :: j,icount
  integer(i_kind) :: ig
  integer(i_kind) :: n_absorbers
  
  !
  regional = regional0
  nvege_type = nvege_type0

  !
 if (lcloud_fwd) then
  ! 
    allocate(cloud_cont(msig,n_clouds_fwd))
    allocate(cloud_efr(msig,n_clouds_fwd))
    allocate(jcloud(n_clouds_fwd))
    allocate(cloud(nsig,n_clouds_fwd))
    allocate(cloudefr(nsig,n_clouds_fwd))
    allocate(icloud(n_actual_clouds))
    cloud_cont=zero
    cloud_efr =zero
    cloud     =zero
    cloudefr  =zero
  !
    n_actual_clouds_wk = n_actual_clouds
    n_clouds_fwd_wk = n_clouds_fwd
    n_clouds_jac_wk = n_clouds_jac
    cld_sea_only_wk = cld_sea_only
    Load_CloudCoeff = .true.
 else
    n_actual_clouds_wk = 0
    n_clouds_fwd_wk = 0
    n_clouds_jac_wk = 0
    cld_sea_only_wk = .false. 
    Load_CloudCoeff = .false.
 endif
 ! Set up index for input satellite data array

 isatid    = 1  ! index of satellite id
 itime     = 2  ! index of analysis relative obs time
 ilon      = 3  ! index of grid relative obs location (x)
 ilat      = 4  ! index of grid relative obs location (y)
 ilzen_ang = 5  ! index of local (satellite) zenith angle (radians)
 ilazi_ang = 6  ! index of local (satellite) azimuth angle (radians)
 iscan_ang = 7  ! index of scan (look) angle (radians)
 iscan_pos = 8  ! index of integer scan position
 iszen_ang = 9  ! index of solar zenith angle (degrees)
 isazi_ang = 10 ! index of solar azimuth angle (degrees)
 ifrac_sea = 11 ! index of ocean percentage
 ifrac_lnd = 12 ! index of land percentage
 ifrac_ice = 13 ! index of ice percentage
 ifrac_sno = 14 ! index of snow percentage
 its_sea   = 15 ! index of ocean temperature
 its_lnd   = 16 ! index of land temperature
 its_ice   = 17 ! index of ice temperature
 its_sno   = 18 ! index of snow temperature
 itsavg    = 19 ! index of average temperature
 ivty      = 20 ! index of vegetation type
 ivfr      = 21 ! index of vegetation fraction
 isty      = 22 ! index of soil type
 istp      = 23 ! index of soil temperature
 ism       = 24 ! index of soil moisture
 isn       = 25 ! index of snow depth
 izz       = 26 ! index of surface height
 idomsfc   = 27 ! index of dominate surface type
 isfcr     = 28 ! index of surface roughness
 iff10     = 29 ! index of ten meter wind factor
 ilone     = 30 ! index of earth relative longitude (degrees)
 ilate     = 31 ! index of earth relative latitude (degrees)
 icount=ilate
 if (obstype == 'abi') then
    iclr_sky      =  7 ! index of clear sky amount
 endif

  !trace gase
  n_absorbers = min_n_absorbers + n_ghg
  !
! Are there aerosols to affect CRTM?
 if (laerosol_fwd) then 
    if(.not.allocated(aero)) allocate(aero(nsig,n_actual_aerosols))
    if(.not.allocated(aero_conc)) allocate(aero_conc(msig,n_actual_aerosols),auxrh(msig))
    n_actual_aerosols_wk=n_actual_aerosols
    n_aerosols_fwd_wk=n_aerosols_fwd
    n_aerosols_jac_wk=n_aerosols_jac
    Load_AerosolCoeff=.true.
 else
    n_actual_aerosols_wk=0
    n_aerosols_fwd_wk=0
    n_aerosols_jac_wk=0
    Load_AerosolCoeff=.false.
 endif

! Initialize radiative transfer

 sensorlist(1)=isis
 if( crtm_coeffs_path /= "" ) then
    if(init_pass .and. mype==mype_diaghdr) write(6,*)myname_,': crtm_init() on path "'//trim(crtm_coeffs_path)//'"'
    if(init_pass .and. mype==mype_diaghdr) &
        write(6,*)myname_,': crtm_init():sensorlisti=',sensorlist
    if(init_pass .and. mype==mype_diaghdr) &
        write(6,*)myname_,': crtm_init():channelinfo=',channelinfo%sensor_id
    if(init_pass .and. mype==mype_diaghdr) &
        write(6,*)myname_,': crtm_init():Load_CloudCoeff=',Load_CloudCoeff
    if(init_pass .and. mype==mype_diaghdr) &
        write(6,*)myname_,': crtm_init():Load_AerosolCoeff=',Load_AerosolCoeff
    error_status = crtm_init(sensorlist,channelinfo,&
       Process_ID=mype,Output_Process_ID=mype_diaghdr, &
       Load_CloudCoeff=Load_CloudCoeff,Load_AerosolCoeff=Load_AerosolCoeff, &
       File_Path = crtm_coeffs_path )
    if(init_pass .and. mype==mype_diaghdr) &
        write(6,*)myname_,': crtm_init():channelinfo=',channelinfo%sensor_id
 else
    error_status = crtm_init(sensorlist,channelinfo,&
       Process_ID=mype,Output_Process_ID=mype_diaghdr, &
       Load_CloudCoeff=Load_CloudCoeff,Load_AerosolCoeff=Load_AerosolCoeff)
 endif
 if (error_status /= success) then
    write(6,*)myname_,':  ***ERROR*** crtm_init error_status=',error_status,&
       '   TERMINATE PROGRAM EXECUTION'
    call stop2(71)
 endif

! Channel Subsets
 sensorindex = 0
if (channelinfo(1)%sensor_id == "abi_g16") then
    sensorindex =  sensorindex + 1
    error_status = crtm_channelinfo_subset(channelinfo(1), &
           channel_subset = (/(ii,ii=subset_start,subset_end)/))
endif
 if (sensorindex == 0 ) then
    write(6,*)myname_,':  ***WARNING*** problem with sensorindex=',isis,&
       ' --> CAN NOT PROCESS isis=',isis,'   TERMINATE PROGRAM EXECUTION found ',&
       channelinfo(1)%sensor_id
    call stop2(71)
 endif
! Check for consistency between user specified number of channels (nchanl)
! and those defined by CRTM channelinfo structure.   Return to calling
! routine if there is a mismatch.

 if (nchanl /= crtm_channelinfo_n_channels(channelinfo(sensorindex))) then
    write(6,*)myname_,':  ***WARNING*** mismatch between nchanl=',&
       nchanl,' and n_channels=',crtm_channelinfo_n_channels(channelinfo(sensorindex)),&
       ' --> CAN NOT PROCESS isis=',isis,'   TERMINATE PROGRAM EXECUTION'
    call stop2(71)
 endif

! Allocate structures for radiative transfer
 mixed_use=.False.
 if (lcloud_fwd .and. (.not. mixed_use) .and. (.not. allocated(rtsolution0)) ) & 
    allocate(rtsolution0(channelinfo(sensorindex)%n_channels,1))   

    allocate( &
    rtsolution  (channelinfo(sensorindex)%n_channels,1),&
    rtsolution_k(channelinfo(sensorindex)%n_channels,1),&
    atmosphere_k(channelinfo(sensorindex)%n_channels,1),&
    surface_k   (channelinfo(sensorindex)%n_channels,1))
 if (mixed_use)  allocate(&
    rtsolution_clr  (channelinfo(sensorindex)%n_channels,1),&
    rtsolution_k_clr(channelinfo(sensorindex)%n_channels,1),&
    atmosphere_k_clr(channelinfo(sensorindex)%n_channels,1),&
    surface_k_clr   (channelinfo(sensorindex)%n_channels,1))

!  Check to ensure that number of levels requested does not exceed crtm max

 if(msig > max_n_layers)then
    write(6,*) myname_,':  msig > max_n_layers - increase crtm max_n_layers ',&
       msig,max_n_layers
    call stop2(36)
 end if

 !  Create structures for radiative transfer
 ! atmosphere
 call crtm_atmosphere_create(atmosphere(1),msig,n_absorbers,n_clouds_fwd_wk,n_aerosols_fwd_wk)
 !
 if(init_pass .and. mype==mype_diaghdr) write(6,*)myname_,&
        ': crtm_init():crtm_microwave_sensor=',crtm_microwave_sensor
 if(init_pass .and. mype==mype_diaghdr) write(6,*)myname_,&
        ': crtm_init():sensor_type=',channelinfo(sensorindex)%sensor_type
 ! microwave_sensor
 if ( channelinfo(sensorindex)%sensor_type == crtm_microwave_sensor ) then
   call crtm_surface_create(surface(1),channelinfo(sensorindex)%n_channels)
   if (.NOT.(crtm_surface_associated(surface(1)))) then
      write(6,*)myname_,' ***ERROR** creating surface.'
   else
      surface(1)%sensordata%sensor_id        = channelinfo(sensorindex)%sensor_id
      surface(1)%sensordata%wmo_sensor_id    = channelinfo(sensorindex)%wmo_sensor_id
      surface(1)%sensordata%wmo_satellite_id = channelinfo(sensorindex)%wmo_satellite_id
      surface(1)%sensordata%sensor_channel   = channelinfo(sensorindex)%sensor_channel
   end if
 end if
!rtsolution ,options
 if (lcloud_fwd .and. (.not. mixed_use)) &                                       
 call crtm_rtsolution_create(rtsolution0,msig) 
 call crtm_rtsolution_create(rtsolution,msig)
 call crtm_rtsolution_create(rtsolution_k,msig)
 call crtm_options_create(options,nchanl)

 if (mixed_use) then
    call crtm_rtsolution_create(rtsolution_clr,msig)
    call crtm_rtsolution_create(rtsolution_k_clr,msig)
 end if

 if (.NOT.(crtm_atmosphere_associated(atmosphere(1)))) &
    write(6,*)myname_,' ***ERROR** creating atmosphere.'
 if (.NOT.(ANY(crtm_rtsolution_associated(rtsolution)))) &
    write(6,*)myname_,' ***ERROR** creating rtsolution.'
 if (lcloud_fwd .and. (.not. mixed_use)) then 
 if (.NOT.(ANY(crtm_rtsolution_associated(rtsolution0)))) &  
    write(6,*)' ***ERROR** creating rtsolution0.'             
 endif                                                        
 if (.NOT.(ANY(crtm_rtsolution_associated(rtsolution_k)))) &
    write(6,*)myname_,' ***ERROR** creating rtsolution_k.'
 if (.NOT.(ANY(crtm_options_associated(options)))) &
    write(6,*)myname_,' ***ERROR** creating options.'

! Turn off antenna correction
 options(1)%use_antenna_correction = .false. 

! Load surface sensor data structure

 surface(1)%sensordata%n_channels = channelinfo(sensorindex)%n_channels

 atmosphere(1)%n_layers = msig
 atmosphere(1)%absorber_id(1) = H2O_ID
 atmosphere(1)%absorber_id(2) = O3_ID
 atmosphere(1)%absorber_units(1) = MASS_MIXING_RATIO_UNITS
 atmosphere(1)%absorber_units(2) = VOLUME_MIXING_RATIO_UNITS
 atmosphere(1)%level_pressure(0) = TOA_PRESSURE

! Currently all considered trace gases affect CRTM. Load trace gases into CRTM atmosphere
 if(init_pass .and. mype==mype_diaghdr) write(6,*)myname_,&
        ': crtm_init():n_ghg=',n_ghg
 ico2=-1
 if (n_ghg>0) then
    do ig=1,n_ghg
       j = min_n_absorbers + ig
       select case(trim(ghg_names(ig)))
         case('co2'); atmosphere(1)%absorber_id(j) = CO2_ID
         case('ch4'); atmosphere(1)%absorber_id(j) = CH4_ID
         case('n2o'); atmosphere(1)%absorber_id(j) = N2O_ID
         case('co') ; atmosphere(1)%absorber_id(j) = CO_ID
         case default
           write(0,*)myname_,':  invalid absorber  TERMINATE PROGRAM'//trim(ghg_names(ig)),71
           call stop2(71)
       end select
       atmosphere(1)%absorber_units(j) = VOLUME_MIXING_RATIO_UNITS
       if (trim(ghg_names(ig))=='co2') ico2=j
    enddo
 endif
 ico24crtm=-1
  ! GEt CO2 DATA ?

!  Allocate structure for _k arrays (jacobians)

 do ii=1,nchanl
    atmosphere_k(ii,1) = atmosphere(1)
    surface_k(ii,1)   = surface(1)
 end do

 if (mixed_use) then
    do ii=1,nchanl
       atmosphere_k_clr(ii,1) = atmosphere(1)
       surface_k_clr(ii,1)   = surface(1)
    end do
 end if
! Mapping land surface type to CRTM surface fields
 if (regional .or. nvege_type==IGBP_N_TYPES) then
    allocate(map_to_crtm_ir(nvege_type))
    allocate(map_to_crtm_mwave(nvege_type))
    if(nvege_type==USGS_N_TYPES)then
       ! Assign mapping for CRTM microwave calculations
       map_to_crtm_mwave=usgs_to_gfs
       ! map usgs to CRTM
       select case ( TRIM(CRTM_IRlandCoeff_Classification()) ) 
         case('NPOESS'); map_to_crtm_ir=usgs_to_npoess
         case('USGS')  ; map_to_crtm_ir=usgs_to_usgs
       end select
    else if(nvege_type==IGBP_N_TYPES)then
       ! Assign mapping for CRTM microwave calculations
       map_to_crtm_mwave=igbp_to_gfs
       ! nmm igbp to CRTM 
       select case ( TRIM(CRTM_IRlandCoeff_Classification()) )
         case('NPOESS'); map_to_crtm_ir=igbp_to_npoess
         case('IGBP')  ; map_to_crtm_ir=igbp_to_igbp
       end select
    else
       write(6,*)myname_,':  ***ERROR*** invalid vegetation types' &
          //' for the CRTM IRland EmisCoeff file used.', &
          ' (only 20 and 24 are setup)  nvege_type=',nvege_type, &
          '  ***STOP IN SETUPRAD***'
       call stop2(71)
    endif ! nvege_type
 endif ! regional or IGBP
    
! Calculate RH when aerosols are present and/or cloud-fraction used
 if (n_actual_aerosols_wk>0) then
    !allocate(gesqsat(lat2,lon2,nsig,nfldsig))
    !ice=.true.
    !iderivative=0
    !do ii=1,nfldsig
       !call genqsat(gesqsat(1,1,1,ii),ges_tsen(1,1,1,ii),ges_prsl(1,1,1,ii),lat2,lon2,nsig,ice,iderivative)
    !end do
 endif

 return

end subroutine init_crtm
!
subroutine destroy_crtm
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    destroy_crtm  deallocates crtm arrays
!   prgmmr: parrish          org: np22                date: 2005-01-22
!
! abstract: deallocates crtm arrays
!
! program history log:
!   2010-08-17  derber 
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  implicit none

  character(len=*),parameter::myname_ = myname//'*destroy_crtm'
  integer(i_kind) error_status

  error_status = crtm_destroy(channelinfo)
  if (error_status /= success) &
     write(6,*)myname_,':  ***ERROR*** error_status=',error_status
  if (n_actual_aerosols_wk>0) then
     deallocate(gesqsat)
  endif
  call crtm_atmosphere_destroy(atmosphere(1))
  call crtm_surface_destroy(surface(1))
  if (n_clouds_fwd_wk>0 .and. (.not. mixed_use)) &     
  call crtm_rtsolution_destroy(rtsolution0)    
  call crtm_rtsolution_destroy(rtsolution)
  call crtm_rtsolution_destroy(rtsolution_k)
  if (mixed_use) then
     call crtm_rtsolution_destroy(rtsolution_clr)
     call crtm_rtsolution_destroy(rtsolution_k_clr)
  end if
  call crtm_options_destroy(options)
  if (crtm_atmosphere_associated(atmosphere(1))) &
     write(6,*)myname_,' ***ERROR** destroying atmosphere.'
  if (crtm_surface_associated(surface(1))) &
     write(6,*)myname_,' ***ERROR** destroying surface.'
  if (ANY(crtm_rtsolution_associated(rtsolution))) &
     write(6,*)myname_,' ***ERROR** destroying rtsolution.'
  if (n_clouds_fwd_wk>0 .and. (.not. mixed_use)) then            
  if (ANY(crtm_rtsolution_associated(rtsolution0))) &    
     write(6,*)' ***ERROR** destroying rtsolution0.'    
  endif                                                 
  if (ANY(crtm_rtsolution_associated(rtsolution_k))) &
     write(6,*)myname_,' ***ERROR** destroying rtsolution_k.'
  if (ANY(crtm_options_associated(options))) &
     write(6,*)myname_,' ***ERROR** destroying options.'
  deallocate(rtsolution,atmosphere_k,surface_k,rtsolution_k)
  if (mixed_use) deallocate(rtsolution_clr,atmosphere_k_clr, & 
                             surface_k_clr,rtsolution_k_clr)
  if (n_clouds_fwd_wk>0 .and. (.not. mixed_use)) &         
  deallocate(rtsolution0) 
  if(n_actual_aerosols_wk>0)then
     deallocate(aero,aero_conc,auxrh)
     if(n_aerosols_jac_wk>0) deallocate(iaero_jac)
  endif
  if (n_ghg>0) then
     deallocate(ghg_names)
  endif
  if(allocated(icloud)) deallocate(icloud)
  if(allocated(cloud)) deallocate(cloud)
  if(allocated(cloudefr)) deallocate(cloudefr)
  if(allocated(jcloud)) deallocate(jcloud)
  if(allocated(cloud_cont)) deallocate(cloud_cont)
  if(allocated(cloud_efr)) deallocate(cloud_efr)
  if(allocated(icw)) deallocate(icw)
  if(allocated(lcloud4crtm_wk)) deallocate(lcloud4crtm_wk)
  if(regional .or. nvege_type==IGBP_N_TYPES)deallocate(map_to_crtm_ir)
  if(regional .or. nvege_type==IGBP_N_TYPES)deallocate(map_to_crtm_mwave)

  return
end subroutine destroy_crtm
!
 subroutine stop2(ierror_code)
  !use mpimod, only: mpi_comm_world,ierror
  implicit none

  integer(i_kind) ierror_code

  write(6,*)'****STOP2****  ABORTING EXECUTION w/code=',ierror_code
  write(0,*)'****STOP2****  ABORTING EXECUTION w/code=',ierror_code
  !call mpi_abort(mpi_comm_world,ierror_code,ierror)
  stop
  return
end subroutine stop2
!
end module crtm_interface
