module parameters_define
    use model_precision,only:P,INT32
    implicit none
    !
    real(P),parameter::PI=3.141592653589793
    real(P),parameter::DTR=3.1415926/180.
    real(P),parameter::RTD=1./DTR
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
    !variables in namelist
    character(LEN=40)::wrf_file
    character(LEN=40)::abi_file
    character(LEN=40)::abiobs_mid_file
    INTEGER(INT32)::nchanl_abi
    INTEGER(INT32),ALLOCATABLE,DIMENSION(:,:):: iuseabi !(nchannel,npass)


    contains
    subroutine read_nml(parameter_file)
        implicit none
        integer(INT32)::istatus
        character(LEN=*)::parameter_file
        !set default value
        wrf_file="wrfinput_d01.2"
        abi_file="goes.nc"
        abiobs_mid_file="abiobs.dat"
        nchanl_abi=10
        allocate(iuseabi(nchanl_abi,1),stat=istatus)

    end subroutine


end module parameters_define
