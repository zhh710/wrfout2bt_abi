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
use model_precision,only:r_kind=>SP,i_kind=>INT32,r_single=>P,r_double=>DP
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
use crtm_module, only: limit_exp,o3_id
use crtm_module, only: crtm_atmosphere_inspect

implicit none

private
public init_crtm            ! Subroutine initializes crtm for specified instrument
public call_crtm            ! Subroutine creates profile for crtm, calls crtm, then adjoint of create
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
  character(len=10),save,allocatable,dimension(:)::cloud_names
  character(len=10),save,allocatable,dimension(:)::cloud_names_fwd
  character(len=10),save,allocatable,dimension(:)::cloud_names_jac
  character(len=10),save,allocatable,dimension(:)::aerosol_names
  character(len=10),save,allocatable,dimension(:)::aerosol_names_fwd
  character(len=10),save,allocatable,dimension(:)::aerosol_names_jac
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
  logical        ,save :: wrf_mass_regional
  logical        ,save :: cold_start
  integer(i_kind),save :: nvege_type
  integer(i_kind),save :: nsig
  integer(i_kind),save :: msig
  integer(i_kind),save,public:: nsigaerojac
  integer(i_kind),save,public:: nsigradjac
  public:: iqv,icw
  integer(i_kind),dimension(200),save :: nlayers

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
     nchanl,isis,obstype,subset_start,subset_end,                    &
     nsig0,                                                      &
     lcloud_fwd,lallsky,cld_sea_only,n_actual_clouds,n_clouds_fwd,   &
     n_clouds_jac,cloud_names0,cloud_names_fwd0,cloud_names_jac0,    &
     laerosol_fwd,laerosol,n_actual_aerosols,n_aerosols_fwd,         &
     n_aerosols_jac,aerosol_names0,aerosol_names_fwd0,aerosol_names_jac0,&
     n_ghg0,ghg_names0,                                                &
     crtm_coeffs_path,                                               &
     regional0,nvege_type0,cold_start0)

    !msig: n_layers
    ! These parameters can be read from namelist


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
  integer(i_kind),intent(in) :: nsig0
  logical        ,intent(in) :: lcloud_fwd,lallsky,cld_sea_only
  integer(i_kind),intent(in) :: n_actual_clouds,n_clouds_fwd
  integer(i_kind),intent(in) :: n_clouds_jac
  character(*)  ,dimension(6),intent(in) :: cloud_names0,cloud_names_fwd0
  character(*)  ,dimension(6),intent(in) :: cloud_names_jac0
  logical        ,intent(in) :: laerosol_fwd,laerosol
  integer(i_kind),intent(in) :: n_actual_aerosols,n_aerosols_fwd
  integer(i_kind),intent(in) :: n_aerosols_jac
  character(*)  ,dimension(3),intent(in) :: aerosol_names0,aerosol_names_fwd0
  character(*)  ,dimension(3),intent(in) :: aerosol_names_jac0
  integer(i_kind),intent(in) :: n_ghg0
  character(*)  ,dimension(5),intent(in) :: ghg_names0
  character(200) ,intent(in) :: crtm_coeffs_path
  logical,        intent(in) :: regional0
  logical,        intent(in) :: cold_start0
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
! ...
  integer::mvl
  
  !
  regional = regional0
  wrf_mass_regional = regional
  cold_start = cold_start0
  nvege_type = nvege_type0
  nsig=nsig0
  n_ghg = n_ghg0
  !
! Sum total number of vertical layers for RTM call
    msig = 0
    if(size(nlayers)<nsig) call stop2(99) !'insufficient size of nlayers',99)
    do ii=1,size(nlayers)
        nlayers(ii)=1
    end do
    !
    do k=1,nsig
       msig = msig + nlayers(k)
    end do

  !
  isst=-1
  ivs=-1
  ius=-1
  ioz=-1
  iqv=-1
  itv=-1
! Get indexes of variables composing the jacobian
  iqv=0
  mvl = 0
  if(isst>=0)mvl=mvl+1
  if(ivs>=0)mvl=mvl+1
  if(ius>=0)mvl=mvl+1
  if(ioz>=0)mvl=mvl+1
  if(iqv>=0)mvl=mvl+1
  if(itv>=0)mvl=mvl+1

!

 if (lcloud_fwd) then
  ! 
    allocate(cloud_cont(msig,n_clouds_fwd))
    allocate(cloud_efr(msig,n_clouds_fwd))
    allocate(jcloud(n_clouds_fwd))
    allocate(cloud(nsig,n_clouds_fwd))
    allocate(cloudefr(nsig,n_clouds_fwd))
    allocate(icloud(n_actual_clouds))

    !
    allocate(cloud_names(n_actual_clouds))
    allocate(cloud_names_fwd(n_clouds_fwd))
    allocate(cloud_names_jac(n_clouds_jac))


    do ii=1,n_actual_clouds
       cloud_names(ii) = cloud_names0(ii)
    end do

    do ii=1, n_clouds_fwd
       cloud_names_fwd(ii) = cloud_names_fwd0(ii)
       !print*,"cloud_names_fwd(ii)",cloud_names_fwd(ii)
    end do

    do ii=1,n_clouds_jac
       cloud_names_jac(ii) = cloud_names_jac0(ii)
    end do

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
! Get indexes of variables for cloud jacobians
  if (n_clouds_jac>0) then
     allocate(icw(max(n_clouds_jac,1)))
     icw=-1
     icount=0
     do ii=1,n_clouds_jac
        indx=getindex(cloud_names_fwd,trim(cloud_names_jac(ii)))
        if (indx>0) then
           icount=icount+1
           mvl=mvl+1
           icw(icount)=(mvl-1)*nsig 
        end if
     end do
  end if
  nsigradjac = mvl*nsig

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
 if (obstype == 'abi_g16') then
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

    allocate(aerosol_names(n_actual_aerosols))
    allocate(aerosol_names_fwd(n_aerosols_fwd))
    allocate(aerosol_names_jac(n_aerosols_jac))

    do ii=1,n_actual_aerosols
       aerosol_names(ii) = aerosol_names0(ii)
    end do

    do ii=1,n_aerosols_fwd
       aerosol_names_fwd(ii) = aerosol_names_fwd0(ii)
    end do

    do ii=1,n_aerosols_jac
       aerosol_names_jac(ii) = aerosol_names_jac0(ii)
    end do

! Get indexes of variables composing the jacobian_aero
  mvl=0
  if (n_actual_aerosols > 0) then
     indx_p25   = getindex(aerosol_names,'p25')
     indx_dust1 = getindex(aerosol_names,'dust1')
     indx_dust2 = getindex(aerosol_names,'dust2')
     if (n_aerosols_jac >0) then
        allocate(iaero_jac(n_aerosols_jac))
        iaero_jac=-1
        icount=0
        do ii=1,n_actual_aerosols
           indx=getindex(aerosol_names_jac,trim(aerosol_names(ii)))
           if(indx>0) then
              icount=icount+1
              mvl=mvl+1
              iaero_jac(icount)=(mvl-1)*nsig
           endif
        end do
     endif
  endif
  nsigaerojac = mvl*nsig

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
 print*,'msig,n_absorbers,n_clouds_fwd_wk,n_aerosols_fwd_wk: ',&
  & msig,n_absorbers,n_clouds_fwd_wk,n_aerosols_fwd_wk
 call crtm_atmosphere_create(atmosphere(1),msig,n_absorbers,n_clouds_fwd_wk,n_aerosols_fwd_wk)
 !
 if(init_pass .and. mype==mype_diaghdr) write(6,*)myname_,&
        ': crtm_init():crtm_microwave_sensor=',crtm_microwave_sensor
 if(init_pass .and. mype==mype_diaghdr) write(6,*)myname_,&
        ': crtm_init():sensor_type=',channelinfo(sensorindex)%sensor_type
 if(init_pass .and. mype==mype_diaghdr) write(6,*)myname_,'crtm_microwave_sensor,',&
        crtm_microwave_sensor 
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
 if (lcloud_fwd .and. (.not. mixed_use)) call crtm_rtsolution_create(rtsolution0,msig) 
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
! SOI
 options%RT_Algorithm_Id = RT_SOI

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
 allocate(ghg_names(min_n_absorbers+n_ghg))
 ghg_names(1)='h2o'
 ghg_names(2)='o3'
 if (n_ghg>0) then
    do ig=1,n_ghg
       j = min_n_absorbers + ig
       
       ghg_names(j) = ghg_names0(ig)
       !print*,'min_n_absorbers + ig',j,ghg_names(j)

       select case(trim(ghg_names(j)))
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
       !call genqsat(gesqsat(1,1,1,ii),tkges(1,1,1,ii),prslges(1,1,1,ii),lat2,lon2,nsig,ice,iderivative)
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
subroutine call_crtm(obstype,iadate,data_s,nchanl,nreal,&
           &   mype, lon2,lat2,&            
           &   psges,tkges,qges ,prslges,prsiges,&
           &   uu5ges,vv5ges,tropprs, &
           &   isli2,sno2, dsfct,    &
                &   h,q,clw_guess,prsl,prsi, &
                &   trop5,tzbgr,dtsavg,sfc_speed,&
                &   tsim,emissivity,ptau5,ts, &
                &   emissivity_k,temp,wmix,jacobian,error_status,tsim_clr, &
                &   layer_od,jacobian_aero,&
                &   qrges,qsges,qgges,qiges,qhges,qlges,ozges)  
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    call_crtm   creates vertical profile of t,q,oz,p,zs,etc., 
!             calls crtm, and does adjoint of creation (where necessary) for setuprad    
!   prgmmr: parrish          org: np22                date: 1990-10-11
!
! abstract: creates vertical profile of t,q,oz,p,zs,etc., 
!             calls crtm, and does adjoint of creation (where necessary) for setuprad
!
! program history log:
!   2010-08-17  derber - modify from intrppx and add threading
!   2011-02-23  todling/da silva - revisit interface to fill in aerosols
!   2011-03-25  yang   - turn off the drop-off of co2 amount when using climatological CO2
!   2011-05-03  todling - merge with Min-Jeong's MW cloudy radiance; combine w/ metguess
!                         (did not include tendencies since they were calc but not used)
!   2011-05-17  auligne/todling - add handling for hydrometeors
!   2011-06-29  todling - no explict reference to internal bundle arrays
!   2011-07-05  zhu - add cloud_efr & cloudefr; add cloud_efr & jcloud in the interface of Set_CRTM_Cloud
!   2011-07-05  zhu - rewrite cloud_cont & cwj for cloud control variables (lcw4crtm)
!   2012-03-12  veldelst-- add a internal interpolation function (option)
!   2012-04-25  yang - modify to use trace gas chem_bundle. Trace gas variables are 
!                       invoked by the global_anavinfo.ghg.l64.txt
!   2013-02-25  zhu - add cold_start option for regional applications
!   2014-01-31  mkim-- remove 60.0degree boundary for icmask for all-sky MW radiance DA 
!   2014-02-26  zhu - add non zero jacobian so jacobian will be produced for            
!                     clear-sky background or background with small amount of cloud     
!   2014-04-27  eliu - add option to calculate clear-sky Tb under all-sky condition                
!   2015-02-27  eliu-- wind direction fix for using CRTM FASTEM model 
!   2015-03-23  zaizhong ma - add Himawari-8 ahi
!   2015-09-10  zhu  - generalize enabling all-sky and aerosol usage in radiance assimilation,
!                      use n_clouds_fwd_wk,n_aerosols_fwd_wk,cld_sea_only_wk, cld_sea_only_wk,cw_cv,etc
!
!   input argument list:
!     obstype      - type of observations for which to get profile
!     obstime      - time of observations for which to get profile, not used here
!     iadate       - array(5), year,month,day,hour,minute
!     data_s       - array containing input data information
!     nchanl       - number of channels
!     nreal        - number of descriptor information in data_s
!     ich          - channel number array
!   added input : 
!     mype         - mpi rank
!     lon2         - n lon
!     lat2         - n lat
!     psges        - surface pressure         ,  hPa
!     tkges        - temperature              ,  K
!     qges         - specific humidity        ,  kg kg-1
!     prslges      - layer pressure (...,nsig),  kPa
!     prsiges      - level pressure (...,nsig+1),kPa
!     uu5ges       - bottom sigma level zonal wind 
!     vv5ges       - bottom sigma level zonal wind 
!     tropprs      - tropopause pressure      ,  hPa
!     lsli2        - dominate surface type, 0 sea/water,1 land,2 sea ice,3 snow
!     sno2         - snow water depth
!     dsfc         - delta Surface temperatures 
!     
!   optional input argument list:
!     qrges        - QRAIN   crtm: RAIN_CLOUD
!     qsges        - QSNOW   crtm: SNOW_CLOUD
!     qgges        - QGRAUP  crtm: GRAUPEL_CLOUD
!     qiges        - QICE    crtm: ICE_CLOUD
!     qhges        - QHAIL   crtm: HAIL_CLOUD
!     qlges        - QCLOUD  crtm: WATER_CLOUD
!     ozges        - o3
!
!   output argument list:
!     h            - interpolated temperature
!     q            - interpolated specific humidity (max(qsmall,q))
!     prsl         - interpolated layer pressure (nsig),kPa
!     prsi         - interpolated level pressure (nsig+1),kPa
!     trop5        - interpolated tropopause pressure
!     tzbgr        - water surface temperature used in Tz retrieval
!     dtsavg       - delta average skin temperature over surface types
!     uu5          - interpolated bottom sigma level zonal wind    
!     vv5          - interpolated bottom sigma level meridional wind  
!     tsim         - simulated brightness temperatures
!     emissivity   - surface emissivities
!     ptau5        - level transmittances
!     ts           - skin temperature sensitivities
!     emissivity_k - surface emissivity sensitivities             
!     temp         - temperature sensitivities
!     wmix         - humidity sensitivities
!     jacobian     - nsigradjac level jacobians for use in intrad and stprad
!     error_status - error status from crtm
!     layer_od     - layer optical depth
!     jacobian_aero- nsigaerojac level jacobians for use in intaod
!     tsim_clr     - option to output simulated brightness temperatures for clear sky                  
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
!--------
  implicit none
! Declare passed variables
!  real(r_kind)                          ,intent(in   ) :: obstime ! not used
  integer(i_kind) ,dimension(5)         ,intent(in   ) :: iadate
  integer(i_kind)                       ,intent(in   ) :: nchanl,nreal

  integer(i_kind)                       ,intent(in   ) :: lon2,lat2
  integer(i_kind)                       ,intent(in   ) :: mype
  real(r_kind),dimension(lat2,lon2     ),intent(  in ) :: psges
  real(r_kind),dimension(lat2,lon2,nsig),intent(  in ) :: tkges,qges,prslges
  real(r_kind),dimension(lat2,lon2,nsig+1),intent(  in ) :: prsiges
  real(r_kind),dimension(lat2,lon2     ),intent(  in ) :: uu5ges,vv5ges
  real(r_kind),dimension(lat2,lon2     ),intent(  in ) :: tropprs
  real(r_kind),dimension(lat2,lon2     ),intent(  in ) :: isli2
  real(r_kind),dimension(lat2,lon2     ),intent(  in ) :: sno2
  real(r_kind),dimension(lat2,lon2     ),intent(  in ) :: dsfct

  real(r_kind)                          ,intent(  out) :: trop5,tzbgr
  real(r_kind),dimension(nsig)          ,intent(  out) :: h,q,prsl
  real(r_kind),dimension(nsig+1)        ,intent(  out) :: prsi
  real(r_kind)                          ,intent(  out) :: sfc_speed,dtsavg
  real(r_kind),dimension(nchanl+nreal)  ,intent(in   ) :: data_s
  real(r_kind),dimension(nchanl)        ,intent(  out) :: tsim,emissivity,ts,emissivity_k
  character(10)                         ,intent(in   ) :: obstype
  integer(i_kind)                       ,intent(  out) :: error_status
  real(r_kind),dimension(nsig,nchanl)   ,intent(  out) :: temp,ptau5,wmix
  real(r_kind),dimension(nsigradjac,nchanl),intent(out):: jacobian
  real(r_kind)                          ,intent(  out) :: clw_guess
  real(r_kind),dimension(nchanl)        ,intent(  out), optional  :: tsim_clr      
  real(r_kind),dimension(nsigaerojac,nchanl),intent(out),optional :: jacobian_aero
  real(r_kind),dimension(nsig,nchanl)   ,intent(  out)  ,optional :: layer_od

  real(r_kind),dimension(lat2,lon2,nsig),intent( in  )  ,optional :: qrges,qsges
  real(r_kind),dimension(lat2,lon2,nsig),intent( in  )  ,optional :: qgges,qiges
  real(r_kind),dimension(lat2,lon2,nsig),intent( in  )  ,optional :: qhges,qlges
  real(r_kind),dimension(lat2,lon2,nsig),intent( in  )  ,optional :: ozges

! Declare local parameters
  character(len=*),parameter::myname_=myname//'*call_crtm'
  real(r_kind),parameter:: minsnow=0.10
  real(r_kind),parameter:: qsmall  = 1.e-6_r_kind
  real(r_kind),parameter:: ozsmall = 1.e-10_r_kind
  real(r_kind),parameter:: small_wind = 1.e-3_r_kind
  real(r_kind),parameter:: windscale = 999999.0_r_kind
  real(r_kind),parameter:: windlimit = 0.0001_r_kind
  real(r_kind),parameter:: quadcof  (4, 2  ) =      &
      reshape((/0.0_r_kind, 1.0_r_kind, 1.0_r_kind, 2.0_r_kind, 1.0_r_kind, &
               -1.0_r_kind, 1.0_r_kind, -1.0_r_kind/), (/4, 2/))

! Declare local variables  
  integer(i_kind):: iquadrant  
  integer(i_kind):: ier,ii,kk,kk2,i,itype,leap_day,day_of_year
  integer(i_kind):: ig,istatus
  integer(i_kind):: j,k,m1,ix,ix1,ixp,iy,iy1,iyp,m,iii
  integer(i_kind):: itsig,itsigp,itsfc,itsfcp
  integer(i_kind):: istyp00,istyp01,istyp10,istyp11
  integer(i_kind):: iqs,iozs
  integer(i_kind):: error_status_clr
  integer(i_kind),dimension(8)::obs_time,anal_time
  integer(i_kind),dimension(msig) :: klevel

! ****************************** 
! Constrained indexing for lai
! CRTM 2.1 implementation change
! ******************************
  integer(i_kind):: lai_type

  real(r_kind):: wind10,wind10_direction,windratio,windangle 
  real(r_kind):: w00,w01,w10,w11,kgkg_kgm2,f10,panglr,dx,dy
! real(r_kind):: w_weights(4)
  real(r_kind):: delx,dely,delx1,dely1,dtsig,dtsigp,dtsfc,dtsfcp
  real(r_kind):: sst00,sst01,sst10,sst11,total_od,term,uu5,vv5, ps
  real(r_kind):: sno00,sno01,sno10,sno11,secant_term
  real(r_kind),dimension(0:3):: wgtavg
  real(r_kind),dimension(nsig,nchanl):: omix
  real(r_kind),dimension(nsig,nchanl,n_aerosols_jac_wk):: jaero
  real(r_kind),dimension(nchanl) :: uwind_k,vwind_k
  real(r_kind),dimension(msig+1) :: prsi_rtm
  real(r_kind),dimension(msig)  :: prsl_rtm
  real(r_kind),dimension(msig)  :: auxq,auxdp
  real(r_kind),dimension(nsig)  :: poz
  real(r_kind),dimension(nsig)  :: rh,qs
  real(r_kind),dimension(5)     :: tmp_time
  real(r_kind),dimension(0:3)   :: dtskin
  real(r_kind),dimension(msig)  :: c6
  real(r_kind),dimension(nsig)  :: c2,c3,c4,c5
  real(r_kind),dimension(nsig) :: ugkg_kgm2,cwj
  real(r_kind),allocatable,dimension(:,:) :: tgas1d

  logical :: sea,icmask,logaod

  integer(i_kind),parameter,dimension(12):: mday=(/0,31,59,90,&
       120,151,181,212,243,273,304,334/)
  real(r_kind) ::   lai
  ! constants
  real(r_kind),parameter::zero=0.0_r_kind
  real(r_kind),parameter::r10=10.0_r_kind
  real(r_kind),parameter::r100=100.0_r_kind
  real(r_kind),parameter::r1000=1000.0_r_kind
  real(r_kind),parameter::one =1.0_r_kind
  real(r_kind),parameter::one_tenth =0.10_r_kind
  real(r_kind),parameter::rd     = 2.8705e+2_r_kind
  real(r_kind),parameter::rv     = 4.6150e+2_r_kind
  real(r_kind),parameter::constoz= 604229.0_r_kind
  real(r_kind),parameter::grav=9.81_r_kind
  real(r_kind),parameter::t0c=273.15_r_kind
  real(r_kind)::pi
  real(r_kind)::deg2rad
  real(r_kind)::rad2deg
  real(r_kind)::fv
  integer(i_kind),parameter::nst_gsi=0
  logical::cw_cv
  real(r_kind)::sqrt_tiny_r_kind
  real(r_kind)::tiny_r_kind
  !
  pi = acos(-one)
  deg2rad = pi/180.0_r_kind
  rad2deg = one/deg2rad
  fv = rv/rd - one

  tiny_r_kind = tiny(zero)
  sqrt_tiny_r_kind = r10*sqrt(tiny_r_kind)
  !
  logaod=(index(obstype,'aod') > 0)

  m1=mype+1

  dx  = data_s(ilat)                 ! grid relative latitude
  dy  = data_s(ilon)                 ! grid relative longitude

! Set spatial interpolation indices and weights
  ix1=dx
  ix1=max(1,min(ix1,lat2))
  delx=dx-ix1
  delx=max(zero,min(delx,one))
  ix=ix1
  ixp=ix+1
  if(ix1==lat2) then
     ixp=ix
  end if
  delx1=one-delx

  iy1=dy
  dely=dy-iy1
  iy=iy1
  if(iy<1) then
     iy1=iy1+lon2
     iy=iy1
  end if
  if(iy>lon2+1) then
     iy1=iy1-lon2
     iy=iy1
  end if
  iyp=iy+1
  dely1=one-dely

  w00=delx1*dely1; w10=delx*dely1; w01=delx1*dely; w11=delx*dely
 ! w_weights = (/w00,w10,w01,w11/)
  if(m1 .eq. 1) write(6,*)'SETUPRAD*CRTM_INTERFACE*CALL_CRTM: ',&
      'w00,w10,w01,w11=',w00,w10,w01,w11,&
      'dx,dy= ',dx,dy,&
      'ix1,delx,ix,ixp,delx1',ix1,delx,ix,ixp,delx1,&
      'iy1,dely,iy,iyp,dely1',iy1,dely,iy,iyp,dely1,&
      'nlat,dx-ix1',lat2,dx-ix1,&
      'nlon',lon2

! Space interpolation of temperature (tk) and q fields from sigma files
!$omp parallel do  schedule(dynamic,1) private(k,ii,iii)
  do k=1,nsig
      if(k == 1)then
        jacobian=zero

!    Set surface type flag.  (Same logic as in subroutine deter_sfc)
        istyp00 = isli2(ix ,iy )
        istyp10 = isli2(ixp,iy )
        istyp01 = isli2(ix ,iyp)
        istyp11 = isli2(ixp,iyp)

        sno00= sno2(ix ,iy )
        sno01= sno2(ix ,iyp)
        sno10= sno2(ixp,iy )
        sno11= sno2(ixp,iyp)
        if(istyp00 >= 1 .and. sno00 > minsnow)istyp00 = 3
        if(istyp01 >= 1 .and. sno01 > minsnow)istyp01 = 3
        if(istyp10 >= 1 .and. sno10 > minsnow)istyp10 = 3
        if(istyp11 >= 1 .and. sno11 > minsnow)istyp11 = 3
        !
!        Find delta Surface temperatures for all surface types

        sst00= dsfct(ix ,iy) ; sst01= dsfct(ix ,iyp)
        sst10= dsfct(ixp,iy ); sst11= dsfct(ixp,iyp) 
        dtsavg=sst00*w00+sst10*w10+sst01*w01+sst11*w11

        dtskin(0:3)=zero
        wgtavg(0:3)=zero
!print*,"Find delta Surface temperatures for all surface types",sno00,istyp00,sst00
        if(istyp00 == 1)then
           wgtavg(1) = wgtavg(1) + w00
           dtskin(1)=dtskin(1)+w00*sst00
        else if(istyp00 == 2)then
           wgtavg(2) = wgtavg(2) + w00
           dtskin(2)=dtskin(2)+w00*sst00
        else if(istyp00 == 3)then
           wgtavg(3) = wgtavg(3) + w00
           dtskin(3)=dtskin(3)+w00*sst00
        else
           wgtavg(0) = wgtavg(0) + w00
           dtskin(0)=dtskin(0)+w00*sst00
        end if

        if(istyp01 == 1)then
           wgtavg(1) = wgtavg(1) + w01
           dtskin(1)=dtskin(1)+w01*sst01
        else if(istyp01 == 2)then
           wgtavg(2) = wgtavg(2) + w01
           dtskin(2)=dtskin(2)+w01*sst01
        else if(istyp01 == 3)then
           wgtavg(3) = wgtavg(3) + w01
           dtskin(3)=dtskin(3)+w01*sst01
        else
           wgtavg(0) = wgtavg(0) + w01
           dtskin(0)=dtskin(0)+w01*sst01
        end if

        if(istyp10 == 1)then
           wgtavg(1) = wgtavg(1) + w10
           dtskin(1)=dtskin(1)+w10*sst10
        else if(istyp10 == 2)then
           wgtavg(2) = wgtavg(2) + w10
           dtskin(2)=dtskin(2)+w10*sst10
        else if(istyp10 == 3)then
           wgtavg(3) = wgtavg(3) + w10
           dtskin(3)=dtskin(3)+w10*sst10
        else
           wgtavg(0) = wgtavg(0) + w10
           dtskin(0)=dtskin(0)+w10*sst10
        end if

        if(istyp11 == 1)then
           wgtavg(1) = wgtavg(1) + w11
           dtskin(1)=dtskin(1)+w11*sst11
        else if(istyp11 == 2)then
           wgtavg(2) = wgtavg(2) + w11
           dtskin(2)=dtskin(2)+w11*sst11
        else if(istyp11 == 3)then
           wgtavg(3) = wgtavg(3) + w11
           dtskin(3)=dtskin(3)+w11*sst11
        else
           wgtavg(0) = wgtavg(0) + w11
           dtskin(0)=dtskin(0)+w11*sst11
        end if

        if(wgtavg(0) > zero)then
           dtskin(0) = dtskin(0)/wgtavg(0)
        else
           dtskin(0) = dtsavg
        end if
        if(wgtavg(1) > zero)then
           dtskin(1) = dtskin(1)/wgtavg(1)
        else
           dtskin(1) = dtsavg
        end if
        if(wgtavg(2) > zero)then
           dtskin(2) = dtskin(2)/wgtavg(2)
        else
           dtskin(2) = dtsavg
        end if
        if(wgtavg(3) > zero)then
           dtskin(3) = dtskin(3)/wgtavg(3)
        else
           dtskin(3) = dtsavg
        end if

        !

        if (n_clouds_fwd_wk>0) then
            ps=psges (ix,iy )*w00+psges (ixp,iy )*w10+ &
                psges (ix,iyp)*w01+psges (ixp,iyp)*w11

        endif
!       skip loading surface structure if obstype is modis_aod/viirs_aod (logaod)
!print*,"skip loading surface structure if obstype is modis_aod/viirs_aod (logaod)"
        if (.not. logaod) then

!       Load surface structure

! **NOTE:  The model surface type --> CRTM surface type
!          mapping below is specific to the versions NCEP
!          GFS and NNM as of Summer 2016
           itype  = nint(data_s(ivty))
           istype = nint(data_s(isty))
           if (regional .or. nvege_type==IGBP_N_TYPES) then
              itype  = min(max(1,itype),nvege_type)
              istype = min(max(1,istype),SOIL_N_TYPES)
              if (ChannelInfo(sensorindex)%sensor_type == crtm_microwave_sensor)then
                 surface(1)%land_type = max(1,map_to_crtm_mwave(itype))
              else
                 surface(1)%land_type = max(1,map_to_crtm_ir(itype))
              end if
              surface(1)%Vegetation_Type = max(1,map_to_crtm_mwave(itype))
              surface(1)%Soil_Type = map_soil_to_crtm(istype)
              lai_type = map_to_crtm_mwave(itype)
           elseif (nvege_type==GFS_N_TYPES) then
              itype  = min(max(0,itype),GFS_VEGETATION_N_TYPES)
              istype = min(max(1,istype),GFS_SOIL_N_TYPES)
              surface(1)%land_type = gfs_to_crtm(itype)
              surface(1)%Vegetation_Type = max(1,itype)
              surface(1)%Soil_Type = istype
              lai_type = itype
           else
              write(6,*)myname_,':  ***ERROR*** invalid vegetation types' &
                 //' the information does not match any currenctly.', &
                 ' supported surface type maps to the CRTM,', &
                 '  ***STOP IN SETUPRAD***'
                 call stop2(71)
           end if

           if (lwind) then
!        Interpolate lowest level winds to observation location 

             uu5=uu5ges (ix,iy )*w00+uu5ges (ixp,iy )*w10+ &
                  uu5ges (ix,iyp)*w01+uu5ges (ixp,iyp)*w11
               
             vv5=vv5ges (ix,iy )*w00+vv5ges (ixp,iy )*w10+ &
                  vv5ges (ix,iyp)*w01+vv5ges (ixp,iyp)*w11
               
             f10=data_s(iff10)
             sfc_speed = f10*sqrt(uu5*uu5+vv5*vv5)
             wind10    = sfc_speed 
             if (uu5*f10 >= 0.0_r_kind .and. vv5*f10 >= 0.0_r_kind) iquadrant = 1 
             if (uu5*f10 >= 0.0_r_kind .and. vv5*f10 <  0.0_r_kind) iquadrant = 2 
             if (uu5*f10 <  0.0_r_kind .and. vv5*f10 >= 0.0_r_kind) iquadrant = 4 
             if (uu5*f10 <  0.0_r_kind .and. vv5*f10 <  0.0_r_kind) iquadrant = 3 
             if (abs(vv5*f10) >= windlimit) then 
                 windratio = (uu5*f10) / (vv5*f10) 
             else 
                 windratio = 0.0_r_kind 
                 if (abs(uu5*f10) > windlimit) then 
                     windratio = windscale * uu5*f10 
                 endif 
             endif 
             windangle        = atan(abs(windratio))   ! wind azimuth is in radians 
             wind10_direction = quadcof(iquadrant, 1) * pi + windangle * quadcof(iquadrant, 2)   
             surface(1)%wind_speed           = sfc_speed
             surface(1)%wind_direction       = rad2deg*wind10_direction
           else !RTodling: not sure the following option makes any sense
             surface(1)%wind_speed           = zero
             surface(1)%wind_direction       = zero
           endif
!       CRTM will reject surface coverages if greater than one and it is possible for
!       these values to be larger due to round off.

           !print*, 'SEA', data_s(ifrac_sea)
           !print*, 'LAND', data_s(ifrac_lnd)
           !print*, 'ICE', data_s(ifrac_ice)
           !print*, 'SNOW', data_s(ifrac_sno)

           surface(1)%water_coverage        = min(max(zero,data_s(ifrac_sea)),one)
           surface(1)%land_coverage         = min(max(zero,data_s(ifrac_lnd)),one)
           surface(1)%ice_coverage          = min(max(zero,data_s(ifrac_ice)),one)
           surface(1)%snow_coverage         = min(max(zero,data_s(ifrac_sno)),one)

!
!       get vegetation lai from summer and winter values.
!

           surface(1)%Lai  = zero
           if (surface(1)%land_coverage>zero) then
              if(lai_type>0)then
                call get_lai(data_s,nchanl,nreal,itime,ilate,lai_type,iadate,lai)
                surface(1)%Lai  = lai   ! LAI  
              endif     
     
              ! for Glacial land ice soil type and vegetation type
              if(surface(1)%Soil_Type == 9 .OR. surface(1)%Vegetation_Type == 13) then
                 surface(1)%ice_coverage = min(real(surface(1)%ice_coverage + surface(1)%land_coverage), one)
                 surface(1)%land_coverage = zero
              endif
           endif

           surface(1)%water_temperature     = max(data_s(its_sea)+dtskin(0),270._r_kind)
           if(nst_gsi>1 .and. surface(1)%water_coverage>zero) then
              surface(1)%water_temperature  = max(data_s(itref)+data_s(idtw)-data_s(idtc)+dtskin(0),271._r_kind)
           endif
           surface(1)%land_temperature      = data_s(its_lnd)+dtskin(1)
           surface(1)%ice_temperature       = min(data_s(its_ice)+dtskin(2),280._r_kind)
           surface(1)%snow_temperature      = min(data_s(its_sno)+dtskin(3),280._r_kind)
           surface(1)%soil_moisture_content = data_s(ism)
           surface(1)%vegetation_fraction   = data_s(ivfr)
           surface(1)%soil_temperature      = data_s(istp)
           surface(1)%snow_depth            = data_s(isn)

           sea = min(max(zero,data_s(ifrac_sea)),one)  >= 0.99_r_kind 
           icmask = (sea .and. cld_sea_only_wk) .or. (.not. cld_sea_only_wk) 

!       assign tzbgr for Tz retrieval when necessary
           tzbgr = surface(1)%water_temperature

        endif ! end of loading surface structure

!       Load geometry structure

!       skip loading geometry structure if obstype is modis_aod/viirs_aod(logaod)   
!       iscan_ang,ilzen_ang,ilazi_ang are not available in the modis aod bufr file
!       also, geometryinfo is not needed in crtm aod calculation
        if ( .not. logaod ) then
           panglr = data_s(iscan_ang)
           if(obstype == 'goes_img' .or. obstype == 'seviri' .or. obstype == 'abi_g16')panglr = zero
           geometryinfo(1)%sensor_zenith_angle = data_s(ilzen_ang)*rad2deg  ! local zenith angle
           geometryinfo(1)%source_zenith_angle = data_s(iszen_ang)          ! solar zenith angle
           geometryinfo(1)%sensor_azimuth_angle = data_s(ilazi_ang)         ! local zenith angle
           geometryinfo(1)%source_azimuth_angle = data_s(isazi_ang)         ! solar zenith angle
           geometryinfo(1)%sensor_scan_angle   = panglr*rad2deg             ! scan angle
           geometryinfo(1)%ifov                = nint(data_s(iscan_pos))    ! field of view position
           !
           !ifov is used in  antenna correction process,antenna correction was turned off in init_crtm
           !

!        For some microwave instruments the solar and sensor azimuth angles can be
!        missing  (given a value of 10^11).  Set these to zero to get past CRTM QC.

           if (geometryinfo(1)%source_azimuth_angle > 360.0_r_kind .OR. &
               geometryinfo(1)%source_azimuth_angle < zero ) &
               geometryinfo(1)%source_azimuth_angle = zero
           if (geometryinfo(1)%sensor_azimuth_angle > 360.0_r_kind .OR. &
               geometryinfo(1)%sensor_azimuth_angle < zero ) &
               geometryinfo(1)%sensor_azimuth_angle = zero

        endif ! end of loading geometry structure
!       Special block for SSU cell pressure leakage correction.   Need to compute
!       observation time and load into Time component of geometryinfo structure.
!       geometryinfo%time is only defined in CFSRR CRTM.
        if (obstype == 'ssu') then
!          Compute absolute observation time
           print*,myname_,"*add ssu?"
           call stop2(71)
        endif

!       Load surface sensor data structure

        do i=1,nchanl


!        Set-up to return Tb jacobians.                                         

           rtsolution_k(i,1)%radiance = zero
           rtsolution_k(i,1)%brightness_temperature = one
           if (mixed_use) then 
              rtsolution_k_clr(i,1)%radiance = zero
              rtsolution_k_clr(i,1)%brightness_temperature = one
           end if

           if (.not. logaod)then

!        Pass CRTM array of tb for surface emissiviy calculations
           if ( channelinfo(1)%sensor_type == crtm_microwave_sensor .and. & 
                crtm_surface_associated(surface(1)) ) & 
                surface(1)%sensordata%tb(i) = data_s(nreal+i) 

!       set up to return layer_optical_depth jacobians
              rtsolution_k(i,1)%layer_optical_depth = one
              if (mixed_use) rtsolution_k_clr(i,1)%layer_optical_depth = one
           endif

        end do

     end if
     !end if k==1
     h(k)  = tkges(ix ,iy ,k)*w00+ &
             tkges(ixp,iy ,k)*w10+ &
             tkges(ix ,iyp,k)*w01+ &
             tkges(ixp,iyp,k)*w11
! Interpolate layer pressure to observation point
     prsl(k)= prslges(ix ,iy ,k)*w00+ &
              prslges(ixp,iy ,k)*w10+ &
              prslges(ix ,iyp,k)*w01+ &
              prslges(ixp,iyp,k)*w11
! Interpolate level pressure to observation point
     prsi(k)= prsiges(ix ,iy ,k)*w00+ &
              prsiges(ixp,iy ,k)*w10+ &
              prsiges(ix ,iyp,k)*w01+ &
              prsiges(ixp,iyp,k)*w11


     q(k)  =qges(ix ,iy ,k)*w00+ &
            qges(ixp,iy ,k)*w10+ &
            qges(ix ,iyp,k)*w01+ &
            qges(ixp,iyp,k)*w11
     q(k)=max(qsmall,q(k))
     !
     c2(k)=one/(one+fv*q(k))
     c3(k)=one/(one-q(k))
     c4(k)=fv*h(k)*c2(k)
     c5(k)=r1000*c3(k)*c3(k)
! Space interpolation of ozone(poz)
     if (present(ozges)) then
         poz(k)= (ozges (ix ,iy ,k)*w00+ &
                  ozges (ixp,iy ,k)*w10+ &
                  ozges (ix ,iyp,k)*w01+ &
                  ozges (ixp,iyp,k)*w11)*constoz

!        Ensure ozone is greater than ozsmall

         poz(k)=max(ozsmall,poz(k))
     endif ! oz
! Quantities required for MW cloudy radiance calculations
     if (n_clouds_fwd_wk>0) then

        do ii=1,n_clouds_fwd_wk
           ! water cloud
           if(trim(cloud_names_fwd(ii)) == 'ql')then
               if(present(qlges))then
                   cloud(k,ii) = qlges(ix ,iy ,k)*w00 + qlges(ixp,iy ,k)*w10 + &
                                qlges(ix ,iyp,k)*w01 + qlges(ixp,iyp,k)*w11
                end if
           end if
           ! ice cloud
           if(trim(cloud_names_fwd(ii)) == 'qi')then
               if(present(qiges))then
                   cloud(k,ii) = qiges(ix ,iy ,k)*w00 + qiges(ixp,iy ,k)*w10 + &
                                qiges(ix ,iyp,k)*w01 + qiges(ixp,iyp,k)*w11
                end if
           end if
           ! snow cloud
           if(trim(cloud_names_fwd(ii)) == 'qs')then
               if(present(qsges))then
                   cloud(k,ii) = qsges(ix ,iy ,k)*w00 + qsges(ixp,iy ,k)*w10 + &
                                qsges(ix ,iyp,k)*w01 + qsges(ixp,iyp,k)*w11
                end if
           end if
           ! graup cloud
           if(trim(cloud_names_fwd(ii)) == 'qg')then
               if(present(qgges))then
                   cloud(k,ii) = qgges(ix ,iy ,k)*w00 + qgges(ixp,iy ,k)*w10 + &
                                qgges(ix ,iyp,k)*w01 + qgges(ixp,iyp,k)*w11
                end if
           end if
           ! hail cloud
           if(trim(cloud_names_fwd(ii)) == 'qh')then
               if(present(qhges))then
                   cloud(k,ii) = qhges(ix ,iy ,k)*w00 + qhges(ixp,iy ,k)*w10 + &
                                qhges(ix ,iyp,k)*w01 + qhges(ixp,iyp,k)*w11
                end if
           end if
           ! rain cloud
           if(trim(cloud_names_fwd(ii)) == 'qr')then
               if(present(qrges))then
                   cloud(k,ii) = qrges(ix ,iy ,k)*w00 + qrges(ixp,iy ,k)*w10 + &
                                qrges(ix ,iyp,k)*w01 + qrges(ixp,iyp,k)*w11
                end if
           end if
           cloud(k,ii)=max(cloud(k,ii),zero)

           !if (regional .and. (.not. wrf_mass_regional)) then
           if (regional .and. .FALSE.) then
           !
           !

            end if

        end do  
     endif ! <n_clouds_fwd_wk>
  end do ! < do k=1,nsig>

! Interpolate level pressure to observation point for top interface
  prsi(nsig+1)=prsiges(ix ,iy ,nsig+1)*w00+ &
               prsiges(ixp,iy ,nsig+1)*w10+ &
               prsiges(ix ,iyp,nsig+1)*w01+ &
               prsiges(ixp,iyp,nsig+1)*w11
! Add additional crtm levels/layers to profile       

  call add_rtm_layers(prsi,prsl,prsi_rtm,prsl_rtm,klevel)



!       check trace gases
  if (n_ghg>0) then
    !allocate (tgas1d(nsig,n_ghg))
    !do ig=1,n_ghg

     !     do k=1,nsig
    !         tgas1d(k,ig) = tgasges(ix ,iy ,k)*w00+ &
    !                        tgasges(ixp,iy ,k)*w10+ &
    !                        tgasges(ix ,iyp,k)*w01+ &
    !                        tgasges(ixp,iyp,k)*w11
    !      enddo
    !enddo

  endif

! Space-time interpolation of aerosol fields from sigma files

  if(n_actual_aerosols_wk>0)then

      ! do ii=1,n_actual_aerosols_wk
      !      do k=1,nsig
      !        aero(k,ii) = aeroges(ix ,iy ,k)*w00+ &
      !                     aeroges(ixp,iy ,k)*w10+ &
      !                     aeroges(ix ,iyp,k)*w01+ &
      !                     aeroges(ixp,iyp,k)*w11
      !     end do
      ! enddo

    !do k=1,nsig
    !!    rh(k) = q(k)/qs(k)
    !end do
  endif

! Find tropopause height at observation

  trop5= one_tenth*(tropprs(ix,iy )*w00+tropprs(ixp,iy )*w10+ &
                    tropprs(ix,iyp)*w01+tropprs(ixp,iyp)*w11)

!  Zero atmosphere jacobian structures

  call crtm_atmosphere_zero(atmosphere_k(:,:))
  call crtm_surface_zero(surface_k(:,:))
  if (mixed_use) then
     call crtm_atmosphere_zero(atmosphere_k_clr(:,:))
     call crtm_surface_zero(surface_k_clr(:,:))
  end if

  clw_guess = zero

  if (n_actual_aerosols_wk>0) then
     do k = 1, nsig
!       Convert mixing-ratio to concentration
        ugkg_kgm2(k)=1.0e-9_r_kind*(prsi(k)-prsi(k+1))*r1000/grav
        aero(k,:)=aero(k,:)*ugkg_kgm2(k)
     enddo
  endif

  sea = min(max(zero,data_s(ifrac_sea)),one)  >= 0.99_r_kind
  icmask = (sea .and. cld_sea_only_wk) .or. (.not. cld_sea_only_wk)

  if(m1 .eq. 1) write(6,*)'SETUPRAD*CRTM_INTERFACE*CALL_CRTM: ',&
      'icmask(mask determining where to consider cloud): ',icmask,&
      'mixed_use : ',mixed_use


  do k = 1,msig

! Load profiles into extended RTM model layers

     kk = msig - k + 1
     atmosphere(1)%level_pressure(k) = r10*prsi_rtm(kk)
     atmosphere(1)%pressure(k)       = r10*prsl_rtm(kk)

     kk2 = klevel(kk)
     atmosphere(1)%temperature(k) = h(kk2)
     atmosphere(1)%absorber(k,1)  = r1000*q(kk2)*c3(kk2)
     if(present(ozges)) then
        atmosphere(1)%absorber(k,2)  = poz(kk2)
     else
        atmosphere(1)%absorber(k,2)  = O3_ID
     endif

     if (n_ghg > 0) then
        do ig=1,n_ghg
           j=min_n_absorbers+ ig
           ! default value for CO2
           if(j==3)atmosphere(1)%absorber(k,j) =380.0
           
        enddo
     endif

     if (n_actual_aerosols_wk>0) then
        !aero_conc(k,:)=aero(kk2,:)
        !auxrh(k)      =rh(kk2)
     endif

! Include cloud guess profiles in mw radiance computation

    cw_cv = .true.

     if (n_clouds_fwd_wk>0) then
        kgkg_kgm2=(atmosphere(1)%level_pressure(k)-atmosphere(1)%level_pressure(k-1))*r100/grav
        if (cw_cv) then
          if (icmask) then 
              c6(k) = kgkg_kgm2
              auxdp(k)=abs(prsi_rtm(kk+1)-prsi_rtm(kk))*r10
              auxq (k)=q(kk2)

              if (.false.) then
                 do ii=1,n_clouds_fwd_wk
                    cloud_cont(k,ii)=cloud(kk2,ii)*c6(k)
                    cloud_efr (k,ii)=cloudefr(kk2,ii)
                 end do
              else
                 do ii=1,n_clouds_fwd_wk
                    cloud_cont(k,ii)=cloud(kk2,ii)*c6(k)
                 end do
              end if

              clw_guess = clw_guess +  cloud_cont(k,1)
              do ii=1,n_clouds_fwd_wk
                 if (ii==1 .and. atmosphere(1)%temperature(k)-t0c>-20.0_r_kind) &
                    cloud_cont(k,1)=max(1.001_r_kind*1.0E-6_r_kind, cloud_cont(k,1))
                 if (ii==2 .and. atmosphere(1)%temperature(k)<t0c) &
                    cloud_cont(k,2)=max(1.001_r_kind*1.0E-6_r_kind, cloud_cont(k,2))
              end do
          endif   
        else 
           if (icmask) then
              do ii=1,n_clouds_fwd_wk
                 cloud_cont(k,ii)=cloud(kk2,ii)*kgkg_kgm2
                 if (trim(cloud_names_fwd(ii))=='ql' .and.  atmosphere(1)%temperature(k)-t0c>-20.0_r_kind) &
                     cloud_cont(k,ii)=max(1.001_r_kind*1.0E-6_r_kind, cloud_cont(k,ii))
                 if (trim(cloud_names_fwd(ii))=='qi' .and.  atmosphere(1)%temperature(k)<t0c) &
                     cloud_cont(k,ii)=max(1.001_r_kind*1.0E-6_r_kind, cloud_cont(k,ii))
              end do
           end if
        endif
     endif

!    Add in a drop-off to absorber amount in the stratosphere to be in more
!    agreement with ECMWF profiles.  The drop-off is removed when climatological CO2 fields
!    are used.
     if(n_ghg == 0)then
        if (atmosphere(1)%level_pressure(k) < 200.0_r_kind) &
            atmosphere(1)%absorber(k,ico2) = atmosphere(1)%absorber(k,ico2) * &
           (0.977_r_kind + 0.000115_r_kind * atmosphere(1)%pressure(k))
     endif
  end do  !< do k=1,msig>

! Set clouds for CRTM
  if(m1 .eq. 1) write(6,*)'SETUPRAD*CRTM_INTERFACE*CALL_CRTM: ',&
      'n_clouds_fwd_wk: ',n_clouds_fwd_wk
  if(n_clouds_fwd_wk>0) then
     atmosphere(1)%n_clouds = n_clouds_fwd_wk  
     call Set_CRTM_Cloud (msig,n_actual_clouds_wk,cloud_names(1:n_actual_clouds_wk), &
          & icmask,n_clouds_fwd_wk,cloud_cont,cloud_efr,auxdp, &
          & real(atmosphere(1)%temperature),real(atmosphere(1)%pressure),auxq,atmosphere(1)%cloud)
  endif

! Set aerosols for CRTM
  if(m1 .eq. 1) write(6,*)'SETUPRAD*CRTM_INTERFACE*CALL_CRTM: ',&
      'n_actual_aerosols_wk: ',n_actual_aerosols_wk
  if(n_actual_aerosols_wk>0) then
     !call Set_CRTM_Aerosol ( msig, n_actual_aerosols_wk, n_aerosols_fwd_wk, aerosol_names, aero_conc, auxrh, &
     !                        atmosphere(1)%aerosol )
  endif
! Call CRTM K Matrix model
  print*,myname_,"*q:",minval(q),maxval(q)
  print*,myname_,"*prsi:",minval(prsi),maxval(prsi)
  print*,myname_,"*prsi_rtm:",minval(prsi_rtm),maxval(prsi_rtm)
  do ig=1,n_clouds_fwd_wk
      print*,myname_,"*cloud_name :",cloud_names(ig),minval(cloud_cont(:,ig)),maxval(cloud_cont(:,ig))
  end do

  !call crtm_atmosphere_inspect(atmosphere)

  error_status = 0
  if ( .not. logaod  ) then
     print*,"call crtm_k_matrix: rtsolution"
     error_status = crtm_k_matrix(atmosphere,surface,rtsolution_k,&
        geometryinfo,channelinfo(sensorindex:sensorindex),atmosphere_k,&
        surface_k,rtsolution,options=options)

     if (mixed_use) then 
        ! Zero out data array in cloud structure
        atmosphere(1)%n_clouds = 0
        print*,"call crtm_k_matrix: n_clouds = 0"
        error_status_clr = crtm_k_matrix(atmosphere,surface,rtsolution_k_clr,&
           geometryinfo,channelinfo(sensorindex:sensorindex),atmosphere_k_clr,&
           surface_k_clr,rtsolution_clr,options=options)
     end if
  else

     do i=1,nchanl
        rtsolution_k(i,1)%layer_optical_depth(:) = one
     enddo

     error_status = crtm_aod_k(atmosphere,rtsolution_k,&
        channelinfo(sensorindex:sensorindex),rtsolution,atmosphere_k)
  end if

! If the CRTM returns an error flag, do not assimilate any channels for this ob
! and set the QC flag to 10 (done in setuprad).

  if (error_status /=0) then
     write(6,*)myname_,':  ***ERROR*** during crtm_k_matrix call ',&
        error_status
  end if
! Calculate clear-sky Tb for AMSU-A over sea when allsky condition is on
  if (n_clouds_fwd_wk>0 .and. present(tsim_clr) .and. (.not. mixed_use)) then
     ! Zero out data array in cloud structure: water content, effective
     ! radius and variance

     atmosphere(1)%n_clouds = 0
!    call crtm_cloud_zero(atmosphere(1)%cloud)

     ! call crtm forward model for clear-sky calculation
     !call crtm_atmosphere_inspect(atmosphere)
     print*,"call crtm_forward: n_clouds = 0"
     error_status = crtm_forward(atmosphere,surface,&
                                 geometryinfo,channelinfo(sensorindex:sensorindex),&
                                 rtsolution0,options=options)
     !print*,"call crtm_forward: n_clouds = 0, finished"
     ! If the CRTM returns an error flag, do not assimilate any channels for this ob
     ! and set the QC flag to 10 (done in setuprad).
     if (error_status /=0) then
        write(6,*)'CRTM_FORWARD  ***ERROR*** during crtm_forward call ',&
        error_status
     end if
  endif 

  if (.not. logaod ) then
! Secant of satellite zenith angle

    secant_term = one/cos(data_s(ilzen_ang))

    if (mixed_use) then
       do i=1,nchanl
          !if (lcloud4crtm_wk(i)<0) then
          if (.false.) then
             rtsolution(i,1) = rtsolution_clr(i,1)
             rtsolution_k(i,1) = rtsolution_k_clr(i,1)
             atmosphere_k(i,1) = atmosphere_k_clr(i,1)
             surface_k(i,1) = surface_k_clr(i,1)
          end if
       end do
    end if

!$omp parallel do  schedule(dynamic,1) private(i) &
!$omp private(total_od,k,kk,m,term,ii,cwj)
    do i=1,nchanl
!   Zero jacobian and transmittance arrays
       do k=1,nsig
         omix(k,i)=zero
         temp(k,i)=zero
         ptau5(k,i)=zero
         wmix(k,i)=zero
       end do
       !!print*,"Zero jacobian and transmittance arrays"
!  Simulated brightness temperatures
       tsim(i)=rtsolution(i,1)%brightness_temperature

       if (n_clouds_fwd_wk>0 .and. present(tsim_clr)) then
          if (mixed_use) then 
             tsim_clr(i)=rtsolution_clr(i,1)%brightness_temperature  
          else
             tsim_clr(i)=rtsolution0(i,1)%brightness_temperature  
          end if
       end if
       !print*,"Simulated brightness temperatures"
!  Estimated emissivity
       emissivity(i)   = rtsolution(i,1)%surface_emissivity

!  Emissivity sensitivities
       emissivity_k(i) = rtsolution_k(i,1)%surface_emissivity
       !print*,"Emissivity sensitivities"

!  Surface temperature sensitivity
!        nst_gsi = 0

       if(nst_gsi>1  ) then
          if (data_s(itz_tr) > zero .and. data_s(itz_tr) <= one) &
          ts(i)   = surface_k(i,1)%water_temperature*data_s(itz_tr) + &
                    surface_k(i,1)%land_temperature + &
                    surface_k(i,1)%ice_temperature + &
                    surface_k(i,1)%snow_temperature
       else
          ts(i)   = surface_k(i,1)%water_temperature + &
                    surface_k(i,1)%land_temperature + &
                    surface_k(i,1)%ice_temperature + &
                    surface_k(i,1)%snow_temperature
       endif
 

       if (abs(ts(i))<sqrt_tiny_r_kind) ts(i) = sign(sqrt_tiny_r_kind,ts(i))
       !print*,"Surface temperature sensitivity"


!  Surface wind sensitivities
       if (surface(1)%wind_speed>small_wind) then
          term = surface_k(i,1)%wind_speed * f10*f10 / surface(1)%wind_speed
          uwind_k(i) = term * uu5
          vwind_k(i) = term * vv5
       else
          uwind_k(i)    = zero
          vwind_k(i)    = zero
       endif
       !print*,"Surface wind sensitivities"


       total_od = zero
!   Accumulate values from extended into model layers
!   temp  - temperature sensitivity
!   wmix  - moisture sensitivity
!   omix  - ozone sensitivity
!   ptau5 - layer transmittance
       do k=1,msig
          kk = klevel(msig-k+1)
          temp(kk,i) = temp(kk,i) + atmosphere_k(i,1)%temperature(k)
          wmix(kk,i) = wmix(kk,i) + atmosphere_k(i,1)%absorber(k,1)
          omix(kk,i) = omix(kk,i) + atmosphere_k(i,1)%absorber(k,2)
          total_od   = total_od + rtsolution(i,1)%layer_optical_depth(k)
          ptau5(kk,i) = exp(-min(real(limit_exp),real(total_od*secant_term)))
       end do
       !print*,"Accumulate values from extended into model layers"

!  Load jacobian array
       do k=1,nsig

!  Small sensitivities for temp
          if (abs(temp(k,i))<sqrt_tiny_r_kind) temp(k,i)=sign(sqrt_tiny_r_kind,temp(k,i))
       end do ! <nsig>
       !print*,"Small sensitivities for temp"
!  Deflate moisture jacobian above the tropopause.
       if (itv>=0) then
          do k=1,nsig
             jacobian(itv+k,i)=temp(k,i)*c2(k)               ! virtual temperature sensitivity
          end do ! <nsig>
       endif
       if (iqv>=0) then

          do k=1,nsig
             jacobian(iqv+k,i)=c5(k)*wmix(k,i)-c4(k)*temp(k,i)        ! moisture sensitivity
             if (prsi(k) < trop5) then
                term = (prsi(k)-trop5)/(trop5-prsi(nsig))
                !jacobian(iqv+k,i) = exp(ifactq(m)*term)*jacobian(iqv+k,i)
                jacobian(iqv+k,i) = exp(15.*term)*jacobian(iqv+k,i)
             endif
          end do ! <nsig>
       endif
       if (ioz>=0) then
!        if (.not. regional .or. use_gfs_ozone)then
          do k=1,nsig
             jacobian(ioz+k,i)=omix(k,i)*constoz       ! ozone sensitivity
          end do ! <nsig>
!        end if
       endif
       !print*,"Deflate moisture jacobian above the tropopause."
       !
       if (n_clouds_fwd_wk>0 .and. n_clouds_jac_wk>0) then
          !if (lcloud4crtm_wk(i)<=0) then 
          if (.false.) then 
             do ii=1,n_clouds_jac_wk
                do k=1,nsig
                   jacobian(icw(ii)+k,i) = zero
                end do 
             end do
          else
             if (icmask) then
                do ii=1,n_clouds_jac_wk
                   do k=1,nsig
                      cwj(k)=zero 
                   end do
                   do k=1,msig
                      kk = klevel(msig-k+1)
                      cwj(kk) = cwj(kk) + atmosphere_k(i,1)%cloud(ii)%water_content(k)*c6(k)
                   end do
                   do k=1,nsig
                      jacobian(icw(ii)+k,i) = cwj(k)
                   end do ! <nsig>
                 end do
             else
                do ii=1,n_clouds_jac_wk
                   do k=1,nsig
                      jacobian(icw(ii)+k,i) = zero
                   end do ! <nsig>
                end do
             endif
          endif
       endif

       if (ius>=0) then
           jacobian(ius+1,i)=uwind_k(i)         ! surface u wind sensitivity
       endif
       if (ivs>=0) then
           jacobian(ivs+1,i)=vwind_k(i)         ! surface v wind sensitivity
       endif
       if (isst>=0) then
           jacobian(isst+1,i)=ts(i)             ! surface skin temperature sensitivity
       endif
       !print*,"ius,ivs,isst,k:",ius,ivs,isst,i
    end do
!$omp barrier

  else                                    !       obstype == 'modis_aod'/viirs_aod
     ! initialize intent(out) variables that are not available with modis_aod
     tzbgr        = zero
     sfc_speed    = zero
     tsim         = zero
     emissivity   = zero
     ts           = zero
     emissivity_k = zero
     ptau5        = zero
     temp         = zero
     wmix         = zero
     jaero        = zero
     if(present(layer_od)) layer_od = zero
     if(present(jacobian_aero)) jacobian_aero = zero
     do i=1,nchanl
        do k=1,msig
           kk = klevel(msig-k+1)
           if(present(layer_od)) then
              layer_od(kk,i) = layer_od(kk,i) + rtsolution(i,1)%layer_optical_depth(k)
           endif
           do ii=1,n_aerosols_jac_wk
              if ( n_aerosols_jac_wk > n_aerosols_fwd_wk .and. ii == indx_p25 ) then
                 jaero(kk,i,ii) = jaero(kk,i,ii) + &
                                  (0.5_r_kind*(0.78_r_kind*atmosphere_k(i,1)%aerosol(indx_dust1)%concentration(k) + &
                                               0.22_r_kind*atmosphere_k(i,1)%aerosol(indx_dust2)%concentration(k)) )
              else
                 jaero(kk,i,ii) = jaero(kk,i,ii) + atmosphere_k(i,1)%aerosol(ii)%concentration(k)
              endif
           enddo
        enddo
        if (present(jacobian_aero)) then
           do k=1,nsig
              do ii=1,n_aerosols_jac_wk
                 jacobian_aero(iaero_jac(ii)+k,i) = jaero(k,i,ii)*ugkg_kgm2(k)
              end do
           enddo
        endif
     enddo
  endif

if(allocated(tgas1d))deallocate (tgas1d)

return
end subroutine call_crtm
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
!-------------------------------------------------------------------------
!    NOAA/NCEP, National Centers for Environmental Prediction GSI        !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: add_rtm_layers --- Add pressure layers for RTM use
!
! !INTERFACE:
!
  subroutine add_rtm_layers(prsitmp,prsltmp,prsitmp_ext,prsltmp_ext,klevel)

! !USES:

    use crtm_module, only: toa_pressure

    implicit none

! !INPUT PARAMETERS:
    integer(i_kind),dimension(msig)  ,intent(  out) :: klevel

    real(r_kind)   ,dimension(nsig+1),intent(in   ) :: prsitmp
    real(r_kind)   ,dimension(nsig)  ,intent(in   ) :: prsltmp

    real(r_kind)   ,dimension(msig+1),intent(  out) :: prsitmp_ext
    real(r_kind)   ,dimension(msig)  ,intent(  out) :: prsltmp_ext


! !DESCRIPTION:  Add pressure layers for use in RTM
!
! !REVISION HISTORY:
!   2005-06-01  treadon
!   2006-05-10  derber modify how levels are added above model top
!   2013-03-27  rancic fix for toa units: crtm(hPa); prsitmp(kPa)
!
! !REMARKS:
!   language: f90
!   machine:  ibm rs/6000 sp; SGI Origin 2000; Compaq/HP
!
! !AUTHOR:
!   treadon          org: w/nmc20      date: 2005-06-01
!
!EOP
!-------------------------------------------------------------------------

!   Declare local variables
    integer(i_kind) k,kk,l
    real(r_kind) dprs,toa_prs_kpa
    real(r_kind),parameter::half=0.5
    real(r_kind),parameter::ten=10.
    real(r_kind),parameter::one=1.0
    real(r_kind),parameter::one_tenth=0.10

!   Convert toa_pressure to kPa
!   ---------------------------
    toa_prs_kpa = toa_pressure*one_tenth

!   Check if model top pressure above rtm top pressure, where prsitmp
!   is in kPa and toa_pressure is in hPa.
    if (prsitmp(nsig) < toa_prs_kpa)then
       write(6,*)'ADD_RTM_LAYERS:  model top pressure(hPa)=', &
            ten*prsitmp(nsig),&
            ' above rtm top pressure(hPa)=',toa_pressure
       call stop2(35)
    end if

!   Linear in pressure sub-divsions
    kk=0
    do k = 1,nsig
       if ((k)<=1) then
          kk = kk + 1
          prsltmp_ext(kk) = prsltmp(k)
          prsitmp_ext(kk) = prsitmp(k)
          klevel(kk) = k
       else
          if (k/=nsig) then
             dprs = (prsitmp(k+1)-prsitmp(k))/nlayers(k)
          else
             dprs = (toa_prs_kpa -prsitmp(k))/nlayers(k)
          end if
          prsitmp_ext(kk+1) = prsitmp(k)
          do l=1,nlayers(k)
             kk=kk + 1
             prsitmp_ext(kk+1) = prsitmp(k) + dprs*l
             prsltmp_ext(kk) = half*(prsitmp_ext(kk+1)+prsitmp_ext(kk))
             klevel(kk) = k
          end do
       endif
    end do

!   Set top of atmosphere pressure
    prsitmp_ext(msig+1) = toa_prs_kpa

   end subroutine add_rtm_layers
!
!
integer(i_kind) function getindex(varnames,usrname)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    getid_
!   prgmmr: todling         org: gmao                date:
!
! abstract:
!
! program history log:
!   2010-05-28  todling - initial code
!
!   input argument list:
!    varnames - array w/ variable names (e.g., cvars3d, or nrf_var)
!    usrname  - name of desired control variable
!
!   output argument list:
!     getindex_ - variable index in varnames (control variable name array)
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

  implicit none
  character(len=*),intent(in) :: varnames(:)
  character(len=*),intent(in) :: usrname
  integer(i_kind) i
  getindex=-1
  do i=1,size(varnames)
     if(trim(usrname)==trim(varnames(i))) then
        getindex=i
        exit
     endif
  enddo
end function getindex
!
!
subroutine Set_CRTM_Cloud ( km, nac, cloud_name, icmask, nc, cloud_cont, cloud_efr, dp, tp, pr, qh, cloud)
  use CRTM_Cloud_Define, only: CRTM_Cloud_type
  implicit none

  integer(i_kind) , intent(in)    :: km                ! number of levels
  integer(i_kind) , intent(in)    :: nac               ! number of actual clouds
  character(len=*), intent(in)    :: cloud_name(nac)   ! [nac]   Model cloud names: qi, ql, etc.
  logical,          intent(in)    :: icmask            ! mask determining where to consider clouds
  integer(i_kind),  intent(in)    :: nc                ! number of clouds
 ! integer(i_kind),  intent(in)    :: jcloud(nc)        ! cloud index
  real(r_kind),     intent(in)    :: cloud_cont(km,nc) ! cloud content 
  real(r_kind),     intent(in)    :: cloud_efr (km,nc) ! cloud effective radius
  real(r_kind),     intent(in)    :: dp(km)            ! [km]    
  real(r_kind),     intent(in)    :: tp(km)            ! [km]   atmospheric temperature (K)
  real(r_kind),     intent(in)    :: pr(km)            ! [km]   atmospheric pressure  
  real(r_kind),     intent(in)    :: qh(km)            ! [km]   specific humidity

  type(CRTM_Cloud_type), intent(inout) :: cloud(nc)    ! [nc]   CRTM Cloud object

  call setCloud (cloud_name, icmask, cloud_cont, cloud_efr, dp, tp, pr, qh, cloud)

end subroutine Set_CRTM_Cloud
  !
  !
subroutine setCloud (cloud_name, icmask, cloud_cont, cloud_efr,dp, tp, pr, qh, cloud)
  use CRTM_Cloud_Define, only: CRTM_Cloud_type

  implicit none

! !ARGUMENTS:

  character(len=*), intent(in)    :: cloud_name(:)     ! [nc]    Model cloud names: Water, Ice, etc.
  logical,          intent(in)    :: icmask            !         mask for where to consider clouds  
 ! integer(i_kind),  intent(in)    :: jcloud(:)         !         cloud order
  real(r_kind),     intent(in)    :: cloud_cont(:,:)   ! [km,nc] cloud contents  (kg/m2)
  real(r_kind),     intent(in)    :: cloud_efr (:,:)   ! [km,nc] cloud effective radius (microns)
  real(r_kind),     intent(in)    :: dp(:)             ! [km]    layer thickness   
  real(r_kind),     intent(in)    :: tp(:)             ! [km]    atmospheric temperature (K)
  real(r_kind),     intent(in)    :: pr(:)             ! [km]    atmospheric pressure (??)
  real(r_kind),     intent(in)    :: qh(:)             ! [km]    atmospheric specific humidity (??)

  type(CRTM_Cloud_type), intent(inout) :: cloud(:)     ! [nc]   CRTM Cloud object

  real(r_kind),parameter::zero=0.0_r_kind
  real(r_kind),parameter::r10=10.0_r_kind
  real(r_kind),parameter::r100=100.0_r_kind
  real(r_kind),parameter::r1000=1000.0_r_kind
  real(r_kind),parameter::one =1.0_r_kind
  real(r_kind),parameter::two =2.0_r_kind
  real(r_kind),parameter::five =5.0_r_kind
  real(r_kind),parameter::r0_05 =0.050_r_kind
  real(r_kind),parameter::one_tenth =0.10_r_kind

  real(r_kind),parameter::rd     = 2.8705e+2_r_kind
  real(r_kind),parameter::rv     = 4.6150e+2_r_kind
  real(r_kind),parameter::constoz= 604229.0_r_kind
  real(r_kind),parameter::grav=9.81_r_kind
  real(r_kind),parameter::t0c=273.15_r_kind
  real(r_kind)::fv

! !DESCRIPTION: Set the CRTM Cloud object given Model cloud properties.
!
! !REVISION HISTORY:
!
! 03May2011  Min-Jeong  Initial version.
! 14May2011  Todling    Encapsulate Min-Jeong's code in present module.
! 01July2011 Zhu        Add jcloud and cloud_efr; add codes for the regional 
! 19Feb2013  Zhu        Add cold_start for the regional
!
!EOP
!-----------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'setCloud'
  integer(i_kind) :: na, nc, km, n, k
  real(r_kind)    :: tem1,tem2,tem3,tem4

  km = size(cloud_cont,1)
  nc = size(cloud_cont,2)
  na = size(cloud_name)

  fv = rv/rd - one

! Handle hand-split case as particular case
! -----------------------------------------
  if (cold_start .or. (na /= nc .and. (.not. regional))) then

!    Initialize Loop over clouds ...
     do n = 1, nc
        Cloud(n)%Type = CloudType_(cloud_name(n))
        Cloud(n)%water_content(:) = zero
        cloud(n)%Effective_Radius(:) = zero
        cloud(n)%effective_variance(:) = two
     enddo

     if(icmask) then
        Cloud(1)%water_content(:) = cloud_cont(:,1)
        Cloud(2)%water_content(:) = cloud_cont(:,2)
     else
        Cloud(1)%water_content(:) = zero
        Cloud(2)%water_content(:) = zero
     endif
!    Calculate effective radius for each cloud component (wired to 2)
!    ----------------------------------------------------------------
     if(icmask) then
        do k=1,km
           ! liquid water cloud drop size
           tem4=max(zero,(t0c-tp(k))*r0_05)
           cloud(1)%effective_radius(k) = five + five * min(one, tem4)

           ! ice water cloud particle size
           tem2 = tp(k) - t0c
           tem1 = grav/rd
           tem3 = tem1 * cloud(2)%water_content(k) * (pr(k)/dp(k)) &
                 /tp(k) * (one + fv * qh(k))

           if (tem2 < -50.0_r_kind ) then
              cloud(2)%effective_radius(k) =  (1250._r_kind/9.917_r_kind)*tem3**0.109_r_kind
           elseif (tem2 < -40.0_r_kind ) then
              cloud(2)%effective_radius(k) =  (1250._r_kind/9.337_r_kind)*tem3**0.08_r_kind
           elseif (tem2 < -30.0_r_kind ) then
              cloud(2)%effective_radius(k) =  (1250._r_kind/9.208_r_kind)*tem3**0.055_r_kind
           else
              cloud(2)%effective_radius(k) =  (1250._r_kind/9.387_r_kind)*tem3**0.031_r_kind
           endif

           cloud(1)%effective_radius(k)=max(0.0_r_double, cloud(1)%effective_radius(k))
           cloud(2)%effective_radius(k)=max(0.0_r_double, cloud(2)%effective_radius(k))

        end do
        cloud(1)%effective_variance(:) = two
        cloud(2)%effective_variance(:) = two

     else
        cloud(1)%effective_radius  (:) = zero
        cloud(2)%effective_radius  (:) = zero
        cloud(1)%effective_variance(:) = two
        cloud(2)%effective_variance(:) = two
     endif

  else ! Handle general case with arbitray number of clouds
       ! --------------------------------------------------

!    Loop over clouds ...
!    --------------------
     do n = 1, nc

!       Map Model cloud names into CRTM Cloud indices
!       ---------------------------------------------
        Cloud(n)%Type = CloudType_(cloud_name(n))
        !print*,myname,"*CloudType_:",CloudType_(cloud_name(n)),cloud_name(n)

        if(icmask) then
           Cloud(n)%water_content(:) = cloud_cont(:,n)
        else
           Cloud(n)%water_content(:) = zero
        endif

!       Calculate effective radius of given cloud type
!       ----------------------------------------------
        if(icmask) then
           if (regional .and. (.not. wrf_mass_regional)) then
              cloud(n)%Effective_Radius(:) = cloud_efr(:,n)
           else
              cloud(n)%Effective_Radius(:) = EftSize_(cloud_name(n))
              !print*,myname,"*EftSize_(cloud_name(n)):",EftSize_(cloud_name(n)),cloud_name(n)

           end if
        else
           cloud(n)%Effective_Radius(:) = zero
        endif
        cloud(n)%effective_variance(:) = two

     enddo

  endif

end subroutine setCloud
!
function CloudType_(name) Result(ctype)
    use CRTM_Cloud_Define, only: WATER_CLOUD,ICE_CLOUD,RAIN_CLOUD, &
      SNOW_CLOUD,GRAUPEL_CLOUD,HAIL_CLOUD 

    character(len=*), parameter :: myname = 'CloudType_'
    character(len=*) :: name  ! Model cloud name
    integer(i_kind)  :: ctype ! CRTM cloud type

    if ( trim(name) == 'ql' ) then
       ctype = WATER_CLOUD
    else if ( trim(name) == 'qi' ) then
       ctype = ICE_CLOUD
    else if ( trim(name) == 'qh' ) then
       ctype = HAIL_CLOUD
    else if ( trim(name) == 'qg' ) then
       ctype = GRAUPEL_CLOUD
    else if ( trim(name) == 'qr' ) then
       ctype = RAIN_CLOUD
    else if ( trim(name) == 'qs' ) then
       ctype = SNOW_CLOUD

    else
       !call die2_(myname,"cannot recognize cloud name <"//trim(name)//">")
       call stop2(71)
    end if

end function CloudType_

function EftSize_(name) Result(csize)
    implicit none
    character(len=*), parameter :: myname_ = 'EftSize_'
    character(len=*) :: name  ! Model cloud name
    real(r_kind)     :: csize ! CRTM cloud type

! Note: Values below from Tom Auligne
    if ( trim(name) == 'ql' ) then
       csize = 10.0_r_kind  ! in micron
    else if ( trim(name) == 'qi' ) then
       csize = 30.0_r_kind
    else if ( trim(name) == 'qh' ) then
!       csize = zero ! RT: can somebody fill this in?
       csize = 1000.0_r_kind
    else if ( trim(name) == 'qg' ) then
       csize = 600.0_r_kind
    else if ( trim(name) == 'qr' ) then
       csize = 300.0_r_kind
    else if ( trim(name) == 'qs' ) then
       csize = 600.0_r_kind

    else
       !call die2_(myname,"cannot recognize cloud name <"//trim(name)//">")
       call stop2(71)
    end if

end function EftSize_
!
subroutine get_lai(data_s,nchanl,nreal,itime,ilate,lai_type,iadate,lai)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    get_lai   interpolate vegetation LAI data for call_crtm
!
!   prgmmr:
!
! abstract:
!
! program history log:
!
!   input argument list:
!     data_s       - array containing input data information
!     nchanl       - number of channels
!     nreal        - number of descriptor information in data_s
!     itime        - index of analysis relative obs time
!     ilate        - index of earth relative latitude (degrees)
!
!   output argument list:
!     lai          - interpolated vegetation leaf-area-index for various types (13)
!
!   language: f90
!   machine:  ibm RS/6000 SP
!   
!$$$
!--------
  !use kinds, only: r_kind,i_kind
  !!use constants, only: zero
  !use obsmod, only: iadate
  use w3nco,only:w3movdat,w3doxdat
  implicit none

! Declare passed variables
  integer(i_kind)                       ,intent(in   ) :: nchanl,nreal
  real(r_kind),dimension(nchanl+nreal)  ,intent(in   ) :: data_s
  integer(i_kind)                       ,intent(in   ) :: itime, ilate,lai_type
  integer(i_kind) ,dimension(5)         ,intent(in   ) :: iadate
  real(r_kind)                          ,intent(  out) :: lai

! Declare local variables
  integer(i_kind),dimension(8)::obs_time,anal_time
  real(r_kind),dimension(5)     :: tmp_time

  
  integer(i_kind) jdow, jdoy, jday
  real(r_kind)    rjday
  real(r_kind),dimension(3):: dayhf
  data dayhf/15.5_r_kind, 196.5_r_kind, 380.5_r_kind/
  real(r_kind),dimension(13):: lai_min, lai_max
  data lai_min/3.08_r_kind, 1.85_r_kind, 2.80_r_kind, 5.00_r_kind, 1.00_r_kind, &
               0.50_r_kind, 0.52_r_kind, 0.60_r_kind, 0.50_r_kind, 0.60_r_kind, &
               0.10_r_kind, 1.56_r_kind, 0.01_r_kind            /
  data lai_max/6.48_r_kind, 3.31_r_kind, 5.50_r_kind, 6.40_r_kind, 5.16_r_kind, &
               3.66_r_kind, 2.90_r_kind, 2.60_r_kind, 3.66_r_kind, 2.60_r_kind, &
               0.75_r_kind, 5.68_r_kind, 0.01_r_kind            /
  real(r_kind),dimension(2):: lai_season
  real(r_kind)    wei1s, wei2s
  integer(i_kind) n1, n2, mm, mmm, mmp
  real(r_kind),parameter::zero=0.0_r_kind
!
      anal_time=0
      obs_time=0
      tmp_time=zero
      tmp_time(2)=data_s(itime)
      anal_time(1)=iadate(1)
      anal_time(2)=iadate(2)
      anal_time(3)=iadate(3)
      anal_time(5)=iadate(4)
      call w3movdat(tmp_time,anal_time,obs_time)

      jdow = 0
      jdoy = 0
      jday = 0
      call w3doxdat(obs_time,jdow,jdoy,jday)
      rjday=jdoy+obs_time(5)/24.0_r_kind
      if(rjday.lt.dayhf(1)) rjday=rjday+365.0

      DO MM=1,2
        MMM=MM
        MMP=MM+1
        IF(RJDAY.GE.DAYHF(MMM).AND.RJDAY.LT.DAYHF(MMP)) THEN
            N1=MMM
            N2=MMP
            GO TO 10
        ENDIF
      ENDDO
      PRINT *,'WRONG RJDAY',RJDAY
   10 CONTINUE
      WEI1S = (DAYHF(N2)-RJDAY)/(DAYHF(N2)-DAYHF(N1))
      WEI2S = (RJDAY-DAYHF(N1))/(DAYHF(N2)-DAYHF(N1))
      IF(N2.EQ.3) N2=1

      lai_season(1) = lai_min(lai_type)
      lai_season(2) = lai_max(lai_type)
      if(data_s(ilate) < 0.0_r_kind) then
         lai = wei1s * lai_season(n2) + wei2s * lai_season(n1)
      else
         lai = wei1s * lai_season(n1) + wei2s * lai_season(n2)
      endif

  return
  end subroutine get_lai
end module crtm_interface
