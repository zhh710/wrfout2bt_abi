module parameters_define
    use model_precision,only:P,INT32
    implicit none
    !
    real(P),parameter::PI=3.141592653589793
    real(P),parameter::DTR=3.1415926/180.
    real(P),parameter::deg2rad=3.1415926/180.
    real(P),parameter::RTD=1./DTR
    real(P),parameter::rad2deg=1./DTR
    real(P),parameter::P1000=1000.E2
    real(P),parameter :: CAPA=0.28589641E0
    real(P), parameter :: H1000=1000. 
    real(P), parameter :: H1=1.
    real(P), parameter :: G=9.81
    real(P), parameter :: RD=287.04
    real(P), parameter :: D608=0.608
    real(P), parameter :: SMALL=1.E-6
    real(P), parameter :: Qconv=0.1E-3
    real(P), parameter :: TFRZ=273.15
    real(P), parameter :: T_ICE=-30.
    real(P), parameter :: TRAD_ice=0.5*T_ICE+TFRZ
    real(P), parameter :: half=0.5
    real(P), parameter :: zero=0.0

    INTEGER(INT32), parameter :: one=1
    !
    !variables in namelist
    character(LEN=40)::wrf_file
    character(LEN=40)::abi_file
    character(LEN=40)::abiobs_mid_file
    character(LEN=40)::tbfout
    character(LEN=200)::datapath
    INTEGER(INT32)::nchanl_abi
    logical,ALLOCATABLE,DIMENSION(:,:):: iuseabi !(nchannel,npass)
    REAL(P),ALLOCATABLE,DIMENSION(:,:):: cld_abi_err !(nchannel,npass)
    REAL(P),ALLOCATABLE,DIMENSION(:,:):: clr_abi_err !(nchannel,npass)
    NAMELIST /adas_abi/ abi_file,abiobs_mid_file, &
        &               nchanl_abi,iuseabi,wrf_file,tbfout,datapath

    ! variables related to crtm
    ! cloud
    logical::lcloud_fwd
    logical::lallsky
    logical::cld_sea_only
    integer::n_actual_clouds
    integer::n_clouds_fwd
    integer::n_clouds_jac
    character(len=10),dimension(6)::cloud_names
    character(len=10),dimension(6)::cloud_names_fwd
    character(len=10),dimension(6)::cloud_names_jac
    NAMELIST /crtm_cloud/lcloud_fwd,lallsky,cld_sea_only ,&
    & n_actual_clouds,n_clouds_fwd,n_clouds_jac,          &
    & cloud_names,cloud_names_fwd,cloud_names_jac
    ! aerosol
    logical::laerosol_fwd
    logical::laerosol
    integer::n_actual_aerosols
    integer::n_aerosols_fwd
    integer::n_aerosols_jac
    character(len=10),dimension(3)::aerosol_names
    character(len=10),dimension(3)::aerosol_names_fwd
    character(len=10),dimension(3)::aerosol_names_jac
    NAMELIST/crtm_aerosol/laerosol_fwd,laerosol,        &
    & n_actual_aerosols,n_aerosols_fwd,n_aerosols_jac,  &
    & aerosol_names,aerosol_names_fwd,aerosol_names_jac
    ! trace gase
    integer::n_ghg
    character(len=8),dimension(6)::ghg_names
    NAMELIST/crtm_tracegase/n_ghg,ghg_names
    !
    character(len=20)::isis
    character(len=10)::obstype
    integer::nchanl
    integer::subset_start
    integer::subset_end
    NAMELIST/goesabi/isis,obstype,nchanl,subset_start,subset_end
    !
    character(len=200)::crtm_coeffs_path
    NAMELIST/crtm_coeffs/crtm_coeffs_path
    !
    integer:: nvege_type
    logical::regional
    NAMELIST/gridmod/nvege_type,regional

    contains
    subroutine read_nml()
        implicit none
        integer(INT32)::istatus
        character(LEN=*),parameter::parameter_file='wrf2abi.namelist'
        !set default value
        wrf_file="wrfinput_d01.2"
        tbfout="tb.nc"
        datapath="./"
        abi_file="goes.nc"
        abiobs_mid_file="abiobs.dat"
        nchanl_abi=10
        allocate(iuseabi(nchanl_abi,1),stat=istatus)
        iuseabi(:,:) = .False.
        iuseabi(4,1) = .True.
        allocate(cld_abi_err(nchanl_abi,1),stat=istatus)
        cld_abi_err(:,:)=3.0
        allocate(clr_abi_err(nchanl_abi,1),stat=istatus)
        clr_abi_err(:,:)=1.75
        !
        lcloud_fwd=.TRUE.
        lallsky=.TRUE.
        cld_sea_only=.FALSE.
        n_actual_clouds=6
        n_clouds_fwd=6
        n_clouds_jac=0
        cloud_names=(/"ql","qi","qr","qs","qg","qh"/)
        cloud_names_fwd=(/"ql","qi","qr","qs","qg","qh"/)
        !
        laerosol_fwd=.FALSE.
        laerosol=.FALSE.
        n_actual_aerosols=0
        n_aerosols_fwd=0
        n_aerosols_jac=0
        aerosol_names=(/"p25  ","dust1","dust2"/)
        aerosol_names_fwd=(/"p25  ","dust1","dust2"/)
        ! trace gase
        n_ghg=0
        !
        isis='abi'
        nchanl=10
        subset_start=7
        subset_end=16
        !
        crtm_coeffs_path=""
        !
        nvege_type=20
        regional=.TRUE.
        !
        open(88,file=parameter_file,status='old')
        read(88,NML=adas_abi)
        read(88,NML=crtm_cloud)
        read(88,NML=crtm_aerosol)
        read(88,NML=crtm_tracegase)
        read(88,NML=goesabi)
        read(88,NML=crtm_coeffs)
        read(88,NML=gridmod)
        close(88)

    end subroutine


end module parameters_define
