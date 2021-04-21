module read_wrf
    USE model_precision
    USE netcdf
    use parameters_define,only:wrf_file
    implicit none
    !3d arrays allocated in init_wrfinput_array
    real(P),allocatable,dimension(:,:,:)::pmid! pressure at mass levels, Pa
    real(P),allocatable,dimension(:,:,:)::pml ! pressure at full levels, Pa
    real(P),allocatable,dimension(:,:,:)::t   ! full poten TEMP, K
    real(P),allocatable,dimension(:,:,:)::tk  ! full TEMP, K
    real(P),allocatable,dimension(:,:,:)::u,v ! ms-1
    real(P),allocatable,dimension(:,:,:)::qv  ! kg kg-1
    real(P),allocatable,dimension(:,:,:)::qi  ! QICE, kg kg-1
    real(P),allocatable,dimension(:,:,:)::qc  ! QCLOUD, kg kg-1
    real(P),allocatable,dimension(:,:,:)::qnc ! QNCCN, kg-1
    real(P),allocatable,dimension(:,:,:)::qr  ! QRAIN, kg kg-1
    real(P),allocatable,dimension(:,:,:)::qnr ! QNRAIN, kg-1
    real(P),allocatable,dimension(:,:,:)::qg  ! QGRAUP, kg kg-1
    real(P),allocatable,dimension(:,:,:)::qng ! QNGRAUPEL, kg-1
    real(P),allocatable,dimension(:,:,:)::qh  ! QHAIL, kg kg-1
    real(P),allocatable,dimension(:,:,:)::qnh ! QNHAIL, kg-1
    real(P),allocatable,dimension(:,:,:)::qs  ! QSNOW, kg kg-1
    real(P),allocatable,dimension(:,:,:)::qns ! QSNOW, kg-1
    real(P),allocatable,dimension(:,:,:)::oz  ! ppmv
    !2d arrays allocated in init_wrfinput_array
    real(P),allocatable,dimension(:,:)::lon,lat ! degree
    real(P),allocatable,dimension(:,:)::psfc    ! SFC PRESSURE,Pa,PSFC
    real(P),allocatable,dimension(:,:)::czen    ! cos(solar zenith angle),COSCEN
    integer(INT32),allocatable,dimension(:,:)::ivgtyp  ! VEGETATION CATEGORY,IVGTYP
    real(P),allocatable,dimension(:,:):: vegfrac!  VEGETATION FRACTION,VRGFRAC
    real(P),allocatable,dimension(:,:):: sno!  SNOW WATER EQUIVALENT, SNOW
    real(P),allocatable,dimension(:,:):: si !  PHYSICAL SNOW DEPTH, SNOWH
    real(P),allocatable,dimension(:,:):: pctsno !  FLAG INDICATING SNOW COVERAGE (1 FOR SNOW COVER),SNOWC
    real(P),allocatable,dimension(:,:):: ths !  SURFACE SKIN TEMPERATURE,TSK
    real(P),allocatable,dimension(:,:):: u10,v10 !  U,V at 10 m, U10,V10
    real(P),allocatable,dimension(:,:):: sice    !  "SEA ICE FLAG, SEAICE


    ! 1d array allocated in init_wrfinput_array
    real(P),allocatable,dimension(:)::eta1 !eta values on full(w) levels, ZNW,nz_stagger
    real(P),allocatable,dimension(:)::eta2 !eta values on mass levels, ZNU,nz
    ! 
    !parameters read in get_dims()
    integer(INT32),private::ncid
    integer(INT32)::nx,ny,nz
    integer(INT32)::nx_stagger,ny_stagger,nz_stagger
    real(P)::ctrlat, ctrlon, trulat1, trulat2, trulon
    integer(INT32)::iproj
    !dim id of wrfout file
    integer(INT32),allocatable,dimension(:)::dim_
    !constant from get_base_value
    real(P)::t00 !BASE STATE TEMPERATURE
    real(P)::t_base ! 300 K ?
    real(P)::p00 ! 100000Pa
    real(P)::g,r_d,r_v,cp_mass
    real(P)::cp,rd_over_cp_mass
    real(P)::p_top ! PRESSURE TOP OF THE MODEL,Pa
    real(P)::grid_ratio!  PARENT_GRID_RATIO
    integer(INT32)::imp_physics !MP_PHYSICS
    !
    public:: load_data,destory_wrfinput_array
    private:: dim_,get_dims,init_wrfinput_array
    !
    contains
        ! get dims
        subroutine get_dims(wrf_file)
            implicit none
            character(len=*),intent(in)::wrf_file
            character(len=1)::filemode='r'
            integer(INT32)::istatus
            integer(INT32)::nzsoil,nstyps,nscalar
            integer(INT32)::P_QC,P_QR,P_QI,P_QS,P_QG,P_QH
            integer(INT32)::P_NC,P_NR,P_NI,P_NS,P_NG,P_NH
            integer(INT32)::P_ZR,P_ZI,P_ZS,P_ZG,P_ZH,P_NCCN
            CHARACTER(LEN=40)::qnames(30)
            real(P)::dx,dy
            call open_ncd_wrf_file ( wrf_file, filemode, ncid, istatus )
            call get_wrf_dimensions(ncid,nx,ny,nz,nzsoil,nstyps,nscalar,    &
                                    P_QC,P_QR,P_QI,P_QS,P_QG,P_QH,          &
                                    P_NC,P_NR,P_NI,P_NS,P_NG,P_NH,          &
                                    P_ZR,P_ZI,P_ZS,P_ZG,P_ZH,               &
                                    P_NCCN,qnames,                          &
                                    iproj,ctrlat,ctrlon,trulat1,trulat2,trulon,&
                                    dx,dy, istatus )
            nx_stagger = nx;nx = nx_stagger -1
            ny_stagger = ny;ny = ny_stagger -1
            nz_stagger = nz;nz = nz_stagger -1
            !
            !write(6,*)qnames
            !write(6,*)iproj,ctrlat,ctrlon,trulat1,trulat2,trulon
            !write(6,*)nx,ny,nz 
            !write(6,*)P_QC,P_QR,P_QI,P_QS,P_QG,P_QH
        end subroutine get_dims
        !
        ! allocate arrays
        !
        subroutine init_wrfinput_array()
            implicit none
            integer(INT32)::istatus
            call get_dims(wrf_file)
            call get_base_value()
            !1d arrays
            allocate(eta1(nz_stagger),stat=istatus)
            if(istatus/=0)write(6,*)"ALLOCATE eta1(nz_stagger) error"
            allocate(eta2(nz),stat=istatus)
            if(istatus/=0)write(6,*)"ALLOCATE eta2(nz) error"
            !2d arrays
            allocate(lat(nx,ny),stat=istatus)
            if(istatus/=0)write(6,*)"ALLOCATE lat(nx,ny) error"
            allocate(lon(nx,ny),stat=istatus)
            if(istatus/=0)write(6,*)"ALLOCATE lon(nx,ny) error"
            allocate(psfc(nx,ny),stat=istatus)
            if(istatus/=0)write(6,*)"ALLOCATE psfc(nx,ny) error"
            allocate(czen(nx,ny),stat=istatus)
            if(istatus/=0)write(6,*)"ALLOCATE czen(nx,ny) error"
            allocate(ivgtyp(nx,ny),stat=istatus)
            if(istatus/=0)write(6,*)"ALLOCATE ivgtyp(nx,ny) error"
            allocate(vegfrac(nx,ny),stat=istatus)
            if(istatus/=0)write(6,*)"ALLOCATE vegfrac(nx,ny) error"
            allocate(sno(nx,ny),stat=istatus)
            if(istatus/=0)write(6,*)"ALLOCATE sno(nx,ny) error"
            allocate(si(nx,ny),stat=istatus)
            allocate(pctsno(nx,ny),stat=istatus)
            if(istatus/=0)write(6,*)"ALLOCATE pctsno(nx,ny) error"
            allocate(ths(nx,ny),stat=istatus)
            if(istatus/=0)write(6,*)"ALLOCATE ths(nx,ny) error"
            allocate(u10(nx,ny),stat=istatus)
            if(istatus/=0)write(6,*)"ALLOCATE u10(nx,ny) error"
            allocate(v10(nx,ny),stat=istatus)
            if(istatus/=0)write(6,*)"ALLOCATE v10(nx,ny) error"
            allocate(sice(nx,ny),stat=istatus)
            if(istatus/=0)write(6,*)"ALLOCATE sice(nx,ny) error"
            ! 3d arrays
            !Pressure, Pa
            allocate(pmid(nx,ny,nz),stat=istatus)
            if(istatus/=0)write(6,*)"ALLOCATE PMID(nx,ny,nz) error"
            !
            allocate(pml(nx,ny,nz_stagger),stat=istatus)
            if(istatus/=0)write(6,*)"ALLOCATE PML(nx,ny,nz) error"
            !Potential temperature
            allocate(t(nx,ny,nz),stat=istatus)
            if(istatus/=0)write(6,*)"ALLOCATE t(nx,ny,nz) error"
            !U wind, m/s
            allocate(u(nx,ny,nz),stat = istatus)
            if(istatus/=0)write(6,*)"ALLOCATE u(nx,ny,nz) error"
            !V wind, m/s
            allocate(v(nx,ny,nz),stat = istatus)
            if(istatus/=0)write(6,*)"ALLOCATE v(nx,ny,nz) error"
            !O3, ppmv
            allocate(oz(nx,ny,nz),stat = istatus)
            if(istatus/=0)write(6,*)"ALLOCATE oz(nx,ny,nz) error"
            !water vapor
            allocate(qv(nx,ny,nz),stat = istatus)
            if(istatus/=0)write(6,*)"ALLOCATE qv(nx,ny,nz) error"


        end subroutine init_wrfinput_array
        !
        !deallocate arrays
        !
        subroutine destory_wrfinput_array()
            implicit none
            integer(INT32)::istatus
            deallocate(eta1,stat=istatus)
            deallocate(eta2,stat=istatus)
            !
            deallocate(lat,stat=istatus)
            deallocate(lon,stat=istatus)
            deallocate(psfc,stat=istatus)
            deallocate(czen,stat=istatus)
            deallocate(vegfrac,stat=istatus)
            deallocate(ivgtyp,stat=istatus)
            deallocate(sno,stat=istatus)
            deallocate(si,stat=istatus)
            deallocate(pctsno,stat=istatus)
            deallocate(ths,stat=istatus)
            deallocate(u10,stat=istatus)
            deallocate(v10,stat=istatus)
            deallocate(sice,stat=istatus)
            !
            deallocate(pmid,stat=istatus)
            deallocate(pml,stat=istatus)
            deallocate(u,stat=istatus)
            deallocate(v,stat=istatus)
            deallocate(oz,stat=istatus)
            deallocate(qv,stat=istatus)
        end subroutine destory_wrfinput_array
        !
        !load data
        !
        subroutine load_data()
            implicit none
            integer(INT32)::istatus
            !local  array 
            real(P),allocatable,dimension(:,:,:)::tmp1_3d,tmp2_3d
            call init_wrfinput_array()
            !1d array
            call get_ncd_1d(ncid,1,"ZNW",nz_stagger,eta1,istatus)
            call get_ncd_1d(ncid,1,"ZNU",nz,eta2,istatus)
            !2d array
            call get_ncd_2d(ncid,1,"XLAT",nx,ny,lat,istatus)
            call get_ncd_2d(ncid,1,"XLONG",nx,ny,lon,istatus)
            call get_ncd_2d(ncid,1,"PSFC",nx,ny,psfc,istatus)
            call get_ncd_2d(ncid,1,"COSZEN",nx,ny,czen,istatus)
            call get_ncd_2d(ncid,1,"VEGFRA",nx,ny,vegfrac,istatus)
            call get_ncd_2d_int(ncid,1,"IVGTYP",nx,ny,IVGTYP,istatus)
            call get_ncd_2d(ncid,1,"SNOW",nx,ny,sno,istatus)
            call get_ncd_2d(ncid,1,"SNOWH",nx,ny,si,istatus)
            call get_ncd_2d(ncid,1,"SNOWC",nx,ny,pctsno,istatus)
            call get_ncd_2d(ncid,1,"TSK",nx,ny,ths,istatus)
            call get_ncd_2d(ncid,1,"U10",nx,ny,u10,istatus)
            call get_ncd_2d(ncid,1,"V10",nx,ny,v10,istatus)
            call get_ncd_2d(ncid,1,"SEAICE",nx,ny,sice,istatus)
            ! 3d array
            ! READ PRESSURE , Pa
            allocate(tmp1_3d(nx,ny,nz),stat=istatus)
            allocate(tmp2_3d(nx,ny,nz),stat=istatus) 
            call get_ncd_3d(ncid,1,"PB",nx,ny,nz,tmp1_3d,istatus)
            call get_ncd_3d(ncid,1,"P",nx,ny,nz,tmp2_3d,istatus)
            pmid = tmp1_3d +tmp2_3d
            deallocate(tmp1_3d,stat=istatus)
            deallocate(tmp2_3d,stat =istatus)
            ! Calculate Pressure on full levels
            call get_full_pressure()
            ! READ U, m/s
            ! READ V, m/s
            CALL READ_UV()
            ! READ Potential temperature
            call get_ncd_3d(ncid,1,"T",nx,ny,nz,t,istatus)
            t = t + t_base
            tk = t*(pmid/p00)**rd_over_cp_mass
            ! READ O3
            call get_ncd_3d(ncid,1,"O3RAD",nx,ny,nz,oz,istatus)
            ! READ QVAPOR
            call get_ncd_3d(ncid,1,"QVAPOR",nx,ny,nz,qv,istatus)


        end subroutine load_data
        !
        !base value
        !
        subroutine get_base_value()
            implicit none
            integer(INT32)::istatus
            p00 = 100000. ! Pa
            t_base = 300. ! K
            call get_ncd_scalar(ncid,1,"T00",t00,istatus)
            call get_ncd_scalar(ncid,1,"P_TOP",p_top,istatus)
            !
            g = 9.81
            r_d = 287.04
            r_v = 461.6
            cp_mass = 1004.67
            cp = 1.0046e+3
            rd_over_cp_mass = r_d/cp_mass
            !
            istatus = nf90_get_att(ncid,nf90_global,"PARENT_GRID_RATIO",grid_ratio)
            istatus = nf90_get_att(ncid,nf90_global,"MP_PHYSICS",imp_physics)


        end subroutine get_base_value
        !
        !READ U , m/s
        !
        subroutine read_UV()
            implicit none
            integer(INT32)::i,j,k
            integer(INT32)::istatus
            integer(INT32)::var_id
            integer(INT32)::ndim
            integer(INT32),allocatable,dimension(:)::dim_id
            real(P),allocatable,dimension(:,:,:)::temp_3d
            !
            if ( .NOT. allocated(dim_) )then
                call get_dim_id()
            endif
            ! U
            istatus =  nf90_inq_varid(ncid,'U',var_id)
            istatus =  nf90_inquire_variable(ncid,var_id,ndims=ndim)
            allocate(dim_id(ndim),stat=istatus)
            istatus =  nf90_inquire_variable(ncid,var_id,dimids=dim_id)
            allocate(temp_3d(dim_(dim_id(1)),dim_(dim_id(2)),dim_(dim_id(3))))
            !
            istatus =  nf90_get_var(ncid,var_id,temp_3d)
            !
            ! INTERPOLATE TO MASS GRID
            do k=1,dim_(dim_id(3))
               do j=1,dim_(dim_id(2))
                  do i=1,dim_(dim_id(1))-1
                     u(j,i,k)=.5*(temp_3d(i,j,k)+temp_3d(i+1,j,k))
                  enddo
               enddo
            enddo
            print*,dim_(dim_id(:))
            deallocate(temp_3d)
            deallocate(dim_id)
            ! V
            istatus =  nf90_inq_varid(ncid,'V',var_id)
            istatus =  nf90_inquire_variable(ncid,var_id,ndims=ndim)
            allocate(dim_id(ndim),stat=istatus)
            istatus =  nf90_inquire_variable(ncid,var_id,dimids=dim_id)
            allocate(temp_3d(dim_(dim_id(1)),dim_(dim_id(2)),dim_(dim_id(3))))
            !
            istatus =  nf90_get_var(ncid,var_id,temp_3d)
            !
            ! INTERPOLATE TO MASS GRID
            do k=1,dim_(dim_id(3))
               do j=1,dim_(dim_id(2))-1
                  do i=1,dim_(dim_id(1))
                     v(j,i,k)=.5*(temp_3d(i,j,k)+temp_3d(i,j+1,k))
                  enddo
               enddo
            enddo
            deallocate(temp_3d)
            deallocate(dim_id)
        end subroutine read_UV
        !
        !get full pressure
        !
        subroutine get_full_pressure()
            implicit none
            integer(INT32)::i,j,k
            real(P)::mu
            pml(:,:,1) = psfc(:,:)
            do i=1,nx
                do j=1,ny
                    mu = psfc(i,j)-p_top
                    do k=1,nz
                        pml(i,j,k+1)=mu*eta1(k+1)+p_top
                    enddo
                enddo
            enddo

        end subroutine get_full_pressure
        !
        ! dim_,dim_id
        !
        subroutine get_dim_id()
            implicit none
            integer(INT32)::Time_id
            integer(INT32)::s_n_id,w_e_id
            integer(INT32)::b_t_id
            integer(INT32)::s_n_stag_id,w_e_stag_id
            integer(INT32)::b_t_stag_id
            integer(INT32)::Time_len
            integer(INT32)::s_n_len,w_e_len
            integer(INT32)::b_t_len
            integer(INT32)::s_n_stag_len,w_e_stag_len
            integer(INT32)::b_t_stag_len
            integer(INT32)::istatus
            integer(INT32)::d_max
            istatus = nf90_inq_dimid(ncid,'Time',Time_id)
            IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_dim_id:nf90_inq_dimid(Times)',.TRUE.)
            istatus = nf90_inq_dimid(ncid,'south_north',s_n_id)
            IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_dim_id:nf90_inq_dimid(south_north)',.TRUE.)
            istatus = nf90_inq_dimid(ncid,'west_east',w_e_id)
            IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_dim_id:nf90_inq_dimid(west_east)',.TRUE.)
            istatus = nf90_inq_dimid(ncid,'bottom_top',b_t_id)
            IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_dim_id:nf90_inq_dimid(bottom_top)',.TRUE.)
            istatus = nf90_inq_dimid(ncid,'south_north_stag',s_n_stag_id)
            IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_dim_id:nf90_inq_dimid(south_north_stag)',.TRUE.)
            istatus = nf90_inq_dimid(ncid,'west_east_stag',w_e_stag_id)
            IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_dim_id:nf90_inq_dimid(west_east_stag)',.TRUE.)
            istatus = nf90_inq_dimid(ncid,'bottom_top_stag',b_t_stag_id)
            IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_dim_id:nf90_inq_dimid(bottom_top)',.TRUE.)
            d_max=max(Time_id, s_n_id, w_e_id, b_t_id, s_n_stag_id, w_e_stag_id, b_t_stag_id)
            allocate(dim_(d_max))
            dim_(:)  = -999
            !
            istatus = nf90_inquire_dimension(ncid,Time_id,len=Time_len)
            IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_dim_id:nf90_inquire_dimension(Time_id)',.TRUE.)
            istatus = nf90_inquire_dimension(ncid,s_n_id,len=s_n_len)
            IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_dim_id:nf90_inquire_dimension(s_n_id)',.TRUE.)
            istatus = nf90_inquire_dimension(ncid,w_e_id,len=w_e_len)
            IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_dim_id:nf90_inquire_dimension(w_e_id)',.TRUE.)
            istatus = nf90_inquire_dimension(ncid,b_t_id,len=b_t_len)
            IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_dim_id:nf90_inquire_dimension(b_t_id)',.TRUE.)
            istatus = nf90_inquire_dimension(ncid,s_n_stag_id,len=s_n_stag_len)
            IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_dim_id:nf90_inquire_dimension(n_s_stag)',.TRUE.)
            istatus = nf90_inquire_dimension(ncid,w_e_stag_id,len=w_e_stag_len)
            IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_dim_id:nf90_inquire_dimension(w_e_stag)',.TRUE.)
            istatus = nf90_inquire_dimension(ncid,b_t_stag_id,len=b_t_stag_len)
            IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_dim_id:nf90_inquire_dimension(b_t_stag)',.TRUE.)

            dim_(Time_id)=Time_len
            dim_(s_n_id)=s_n_len
            dim_(w_e_id)=w_e_len
            dim_(b_t_id)=b_t_len
            dim_(s_n_stag_id)=s_n_stag_len
            dim_(w_e_stag_id)=w_e_stag_len
            dim_(b_t_stag_id)=b_t_stag_len

        end subroutine get_dim_id
        !
        !From here on, the code is adapted from arps/wrfapi/wrf_ncd_io.f90
        !
        SUBROUTINE open_ncd_wrf_file ( filename, filemode, ncid, istatus )
            IMPLICIT NONE

            CHARACTER(LEN=*), INTENT(IN)  :: filename
            CHARACTER(LEN=1), INTENT(IN)  :: filemode
            INTEGER(INT32),          INTENT(OUT) :: ncid, istatus

            !------------------------------------------------------------------

            LOGICAL :: fexists
            INTEGER :: file_mode
            istatus = 0

            IF (filemode == 'r') THEN
                INQUIRE(FILE = filename, EXIST = fexists)
                IF (fexists) THEN
                   file_mode = NF90_NOWRITE
                ELSE
                   istatus = -1
                   WRITE(6,'(2a)') 'ERROR: File not found: ', filename
                   RETURN
                ENDIF

                istatus = NF90_OPEN(TRIM(filename),file_mode,ncid)
                CALL handle_ncd_error(istatus,'NF90_OPEN("'//TRIM(filename)//'") in open_wrf_file.',.FALSE.)

            ELSE IF (filemode == 'w') THEN
           !file_mode = IOR(NF90_CLOBBER,NF90_64BIT_OFFSET)
               file_mode =  IOR(IOR(NF90_CLOBBER,NF90_NETCDF4),NF90_CLASSIC_MODEL)
               istatus = NF90_CREATE(TRIM(filename),file_mode,ncid)
               CALL handle_ncd_error(istatus,'NF90_CREATE('//TRIM(filename)//') in open_wrf_file.',.TRUE.)

            ELSE IF (filemode == 'm') THEN
               INQUIRE(FILE = filename, EXIST = fexists)
               IF (fexists) THEN
                  file_mode = NF90_WRITE !NF90_SHARE  !NF90_WRITE
               ELSE
                  istatus = -1
                  WRITE(6,'(2a)') 'ERROR: File not found: ', filename
                  RETURN
               ENDIF

               istatus = NF90_OPEN(TRIM(filename),file_mode,ncid)
               CALL handle_ncd_error(istatus,'NF90_OPEN('//TRIM(filename)//' in open_wrf_file.',.TRUE.)

           ELSE
              istatus = -2
              WRITE(6,'(2a)') 'ERROR: Unsupported file opening mode : ', filemode
           END IF


           RETURN
       END SUBROUTINE open_ncd_wrf_file
       !
!###################### Close an opened WRF file ######################
SUBROUTINE close_ncd_wrf_file ( ncid, istatus )
  IMPLICIT NONE

  INTEGER(INT32), INTENT(IN)  :: ncid
  INTEGER(INT32), INTENT(OUT) :: istatus

!------------------------------------------------------------------


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  istatus = NF90_CLOSE(ncid)
  CALL handle_ncd_error(istatus,'NF90_CLOSE in close_wrf_file',.TRUE.)

  RETURN
END SUBROUTINE close_ncd_wrf_file
!###################### Close an opened WRF file ######################

SUBROUTINE get_wrf_dimensions( ncid, nx,ny,nz,nzsoil,nstyps,nscalar,    &
                             P_QC,P_QR,P_QI,P_QS,P_QG,P_QH,             &
                             P_NC,P_NR,P_NI,P_NS,P_NG,P_NH,             &
                                  P_ZR,P_ZI,P_ZS,P_ZG,P_ZH,             &
                             P_NCCN,qnames,                             &
                             iproj,ctrlat,ctrlon,trulat1,trulat2,trulon,&
                             dx,dy, istatus )
  IMPLICIT NONE

  INTEGER(INT32),           INTENT(IN)  :: ncid
  INTEGER(INT32),           INTENT(OUT) :: nx, ny, nz, nzsoil, nstyps, nscalar
  INTEGER(INT32),           INTENT(OUT) :: P_QC,P_QR,P_QI,P_QS,P_QG,P_QH,    &
                                    P_NC,P_NR,P_NI,P_NS,P_NG,P_NH,     &
                                    P_ZR,P_ZI,P_ZS,P_ZG,P_ZH
  INTEGER(INT32),           INTENT(OUT) :: P_NCCN
  CHARACTER(LEN=40), INTENT(OUT) :: qnames(30)

  INTEGER(INT32),           INTENT(OUT) :: iproj
  REAL(P),              INTENT(OUT) :: ctrlat, ctrlon, trulat1, trulat2, trulon
  REAL(P),              INTENT(OUT) :: dx, dy
  INTEGER(INT32),           INTENT(OUT) :: istatus

!------------------------------------------------------------------
  INTEGER(INT32), PARAMETER :: ntotal = 13
  INTEGER(INT32), PARAMETER :: I_QC = 1, I_QR = 2,  I_QI = 3,  I_QS = 4,  I_QG = 5,  I_QH = 6,  &
                        I_NC = 7, I_NR = 8,  I_NI = 9,  I_NS = 10, I_NG = 11, I_NH = 12, &
                        I_NCCN = 13,                                                     &
                                  I_ZR = 14, I_ZI = 15, I_ZS = 16, I_ZG = 17, I_ZH = 18
                ! denotes the order of names in qnames_wrf

  CHARACTER(LEN=40), PARAMETER :: qnames_wrf(ntotal) = (/ 'QCLOUD    ', &
                'QRAIN     ', 'QICE      ', 'QSNOW     ', 'QGRAUP    ', &
                'QHAIL     ', 'QNCLOUD   ', 'QNRAIN    ', 'QNICE     ', &
                'QNSNOW    ', 'QNGRAUPEL ', 'QNHAIL    ', 'QNCCN     '  /)

  INTEGER(INT32) :: P_ARR(20)   ! remember the order of qnames_wrf and map to the variable P_??.

  INTEGER(INT32) :: dimid, varid
  INTEGER(INT32) :: nq

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0
  istatus = NF90_inq_dimid(ncid,'west_east_stag',dimid)
  IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_wrf_dimensions:NF90_inq_dimid(west_east_stag)',.TRUE.)
  istatus = nf90_inquire_dimension(ncid,dimid,len=nx)
  IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_wrf_dimensions:nf90_inquire_dimension(nx)',.TRUE.)

  istatus = NF90_inq_dimid(ncid,'south_north_stag',dimid)
  IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_wrf_dimensions:NF90_inq_dimid(south_north_stag)',.TRUE.)
  istatus = nf90_inquire_dimension(ncid,dimid,len=ny)
  IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_wrf_dimensions:nf90_inquire_dimension(ny)',.TRUE.)

  istatus = NF90_inq_dimid(ncid,'bottom_top_stag',dimid)
  IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_wrf_dimensions:NF90_inq_dimid(bottom_top_stag)',.TRUE.)
  istatus = nf90_inquire_dimension(ncid,dimid,len=nz)
  IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_wrf_dimensions:nf90_inquire_dimension(nz)',.TRUE.)

  istatus = NF90_inq_dimid(ncid,'soil_layers_stag',dimid)
  IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_wrf_dimensions:NF90_inq_dimid(soil_layers_stag)',.TRUE.)
  istatus = nf90_inquire_dimension(ncid,dimid,len=nzsoil)
  IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_wrf_dimensions:nf90_inquire_dimension(nzsoil)',.TRUE.)

  nstyps = 1
  istatus = NF90_GET_ATT(ncid,NF90_GLOBAL,'MAP_PROJ',iproj)
  IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_wrf_dimensions:NF90_GET_ATT(MAP_PROJ)',.TRUE.)

  istatus = NF90_GET_ATT(ncid,NF90_GLOBAL,'CEN_LAT',ctrlat)
  IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_wrf_dimensions:NF90_GET_ATT(CEN_LAT)',.TRUE.)

  istatus = NF90_GET_ATT(ncid,NF90_GLOBAL,'CEN_LON',ctrlon)
  IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_wrf_dimensions:NF90_GET_ATT(CEN_LON)',.TRUE.)

  istatus = NF90_GET_ATT(ncid,NF90_GLOBAL,'TRUELAT1',trulat1)
  IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_wrf_dimensions:NF90_GET_ATT(TRUELAT1)',.TRUE.)

  istatus = NF90_GET_ATT(ncid,NF90_GLOBAL,'TRUELAT2',trulat2)
  IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_wrf_dimensions:NF90_GET_ATT(TRUELAT2)',.TRUE.)

  istatus = NF90_GET_ATT(ncid,NF90_GLOBAL,'STAND_LON',trulon)
  IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_wrf_dimensions:NF90_GET_ATT(STAND_LON)',.TRUE.)

  istatus = NF90_GET_ATT(ncid,NF90_GLOBAL,'DX',dx)
  IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_wrf_dimensions:NF90_GET_ATT(DX)',.TRUE.)

  istatus = NF90_GET_ATT(ncid,NF90_GLOBAL,'DY',dy)
  IF (istatus /= NF90_NOERR) CALL handle_ncd_error(istatus,'get_wrf_dimensions:NF90_GET_ATT(DY)',.TRUE.)


  nscalar  = 0
  P_ARR(:) = 0
  DO nq = 1, ntotal
    istatus = NF90_INQ_VARID(ncid,trim(qnames_wrf(nq)),varid)
    IF (istatus /= NF90_NOERR) THEN
      
    ELSE
      nscalar = nscalar + 1
      P_ARR(nq) = nscalar
      qnames(nscalar) = qnames_wrf(nq)
    END IF
  END DO

  P_QC = P_ARR(I_QC);  P_NC = P_ARR(I_NC)
  P_QR = P_ARR(I_QR);  P_NR = P_ARR(I_NR);  P_ZR = P_ARR(I_ZR)
  P_QI = P_ARR(I_QI);  P_NI = P_ARR(I_NI);  P_ZI = P_ARR(I_ZI)
  P_QS = P_ARR(I_QS);  P_NS = P_ARR(I_NS);  P_ZS = P_ARR(I_ZS)
  P_QG = P_ARR(I_QG);  P_NG = P_ARR(I_NG);  P_ZG = P_ARR(I_ZG)
  P_QH = P_ARR(I_QH);  P_NH = P_ARR(I_NH);  P_ZH = P_ARR(I_ZH)
  P_NCCN = P_ARR(I_NCCN)

  RETURN
END SUBROUTINE get_wrf_dimensions

SUBROUTINE handle_ncd_error(ierr,sub_name,progstop)
  IMPLICIT NONE
  INTEGER(INT32),          INTENT(IN) :: ierr
  CHARACTER(LEN=*), INTENT(IN) :: sub_name
  LOGICAL,          INTENT(IN) :: progstop
  CHARACTER(LEN=80) :: errmsg
  IF(ierr /= NF90_NOERR) THEN
    errmsg = NF90_STRERROR(ierr)
    WRITE(6,'(1x,2a)') 'NetCDF error: ',errmsg
    WRITE(6,'(1x,3a)') 'ERROR while calling <', sub_name,'>.'
    WRITE(6,*) progstop
    STOP
  END IF

  RETURN
END SUBROUTINE handle_ncd_error

SUBROUTINE get_ncd_2d(ncid,itime,varname,nx,ny,var2d,istatus)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Read in a 2D array from the WRF NetCDF file.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER(INT32),          INTENT(IN)  :: ncid
  INTEGER(INT32),          INTENT(IN)  :: itime
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  INTEGER(INT32),          INTENT(IN)  :: nx
  INTEGER(INT32),          INTENT(IN)  :: ny
  REAL(P),         INTENT(OUT) :: var2d(nx,ny)
  INTEGER(INT32),          INTENT(OUT) :: istatus
!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------


  INTEGER(INT32)                    :: varid
  CHARACTER(LEN=NF90_MAX_NAME) :: namein
  INTEGER(INT32)                    :: vartype, ndims,natts,dimlen
  INTEGER(INT32)                    :: dimids(NF90_MAX_VAR_DIMS)

  INTEGER(INT32), PARAMETER         :: VAR_NOTEXIST = -1

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  istatus = NF90_INQ_VARID(ncid,varname,varid)
  IF(istatus == NF90_ENOTVAR) THEN
     WRITE(6,'(3a)') ' WARNING: variable ',TRIM(varname),' does not exist.'
     var2d(:,:) = -9999.0
     istatus = VAR_NOTEXIST
     RETURN
  END IF
  CALL handle_ncd_error(istatus,'NF90_INQ_VARID in get_ncd_2d with '//TRIM(varname),.FALSE.)

  istatus = NF90_INQUIRE_VARIABLE(ncid,varid,namein,vartype,ndims,dimids,natts)
  CALL handle_ncd_error(istatus,'NF90_INQ_VAR in get_ncd_2d',.FALSE.)

  IF(vartype /= NF90_FLOAT) THEN
    WRITE(6,'(3a)') 'Variable ',varname, ' is not REAL.'
    STOP 'WRONG_VAR_TYPE'
  END IF

  IF(ndims /= 3) THEN
    WRITE(6,'(3a)') 'Variable ', varname, ' is not a 2D array.'
    STOP 'WRONG_VAR_DIMENSIONS'
  END IF

  istatus = NF90_inquire_dimension(ncid,dimids(1),len=dimlen)
  CALL handle_ncd_error(istatus,'NF90_INQ_DIMLEN in get_wrf_2d',.FALSE.)
  IF(dimlen /= nx) THEN
    WRITE(6,'(3a,I3,a,I0)') 'First dimension of variable ', varname,    &
                    ' is ',dimlen, ' and it should be ',nx
    STOP 'WRONG_DIM_length'
  END IF

  istatus = NF90_INQUIRE_DIMENSION(ncid,dimids(2),len=dimlen)
  CALL handle_ncd_error(istatus,'NF90_INQ_DIMLEN in get_wrf_2d',.FALSE.)
  IF(dimlen /= ny) THEN
    WRITE(6,'(3a,I3,a,I0)') 'Second dimension of variable ', varname,   &
                    ' is ',dimlen, ' and it should be ',ny
    STOP 'WRONG_DIM_length'
  END IF

  istatus = NF90_INQUIRE_DIMENSION(ncid,dimids(3),len=dimlen)
  CALL handle_ncd_error(istatus,'NF90_INQ_DIMLEN in get_wrf_2d',.FALSE.)
  IF(dimlen < itime) THEN
    WRITE(6,'(a,I3,a,I3)') 'The total records number is ', dimlen,     &
                    ' however, the required time level is ',itime
    STOP 'itime_tool_large'
  END IF
  istatus = NF90_GET_VAR(ncid,varid,var2d,(/1,1,itime/),(/nx,ny,1/))
  CALL handle_ncd_error(istatus,'NF90_GET_VAR in get_wrf_2d',.FALSE.)

  RETURN
  END SUBROUTINE get_ncd_2d
!
SUBROUTINE get_ncd_2d_int(ncid,itime,varname,nx,ny,var2d,istatus)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Read in a 2D array from the WRF NetCDF file.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER(INT32),          INTENT(IN)  :: ncid
  INTEGER(INT32),          INTENT(IN)  :: itime
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  INTEGER(INT32),          INTENT(IN)  :: nx
  INTEGER(INT32),          INTENT(IN)  :: ny
  INTEGER(INT32),         INTENT(OUT) :: var2d(nx,ny)
  INTEGER(INT32),          INTENT(OUT) :: istatus
!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------


  INTEGER(INT32)                    :: varid
  CHARACTER(LEN=NF90_MAX_NAME) :: namein
  INTEGER(INT32)                    :: vartype, ndims,natts,dimlen
  INTEGER(INT32)                    :: dimids(NF90_MAX_VAR_DIMS)

  INTEGER(INT32), PARAMETER         :: VAR_NOTEXIST = -1

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  istatus = NF90_INQ_VARID(ncid,varname,varid)
  IF(istatus == NF90_ENOTVAR) THEN
     WRITE(6,'(3a)') ' WARNING: variable ',TRIM(varname),' does not exist.'
     var2d(:,:) = -9999.0
     istatus = VAR_NOTEXIST
     RETURN
  END IF
  CALL handle_ncd_error(istatus,'NF90_INQ_VARID in get_ncd_2d_int with '//TRIM(varname),.FALSE.)

  istatus = NF90_INQUIRE_VARIABLE(ncid,varid,namein,vartype,ndims,dimids,natts)
  CALL handle_ncd_error(istatus,'NF90_INQ_VAR in get_ncd_2d_int',.FALSE.)

  IF(vartype /= NF90_INT) THEN
    WRITE(6,'(3a)') 'Variable ',varname, ' is not INTEGER.'
    STOP 'WRONG_VAR_TYPE'
  END IF

  IF(ndims /= 3) THEN
    WRITE(6,'(3a)') 'Variable ', varname, ' is not a 2D array.'
    STOP 'WRONG_VAR_DIMENSIONS'
  END IF

  istatus = NF90_inquire_dimension(ncid,dimids(1),len=dimlen)
  CALL handle_ncd_error(istatus,'NF90_INQ_DIMLEN in get_wrf_2d',.FALSE.)
  IF(dimlen /= nx) THEN
    WRITE(6,'(3a,I3,a,I0)') 'First dimension of variable ', varname,    &
                    ' is ',dimlen, ' and it should be ',nx
    STOP 'WRONG_DIM_length'
  END IF

  istatus = NF90_INQUIRE_DIMENSION(ncid,dimids(2),len=dimlen)
  CALL handle_ncd_error(istatus,'NF90_INQ_DIMLEN in get_wrf_2d',.FALSE.)
  IF(dimlen /= ny) THEN
    WRITE(6,'(3a,I3,a,I0)') 'Second dimension of variable ', varname,   &
                    ' is ',dimlen, ' and it should be ',ny
    STOP 'WRONG_DIM_length'
  END IF

  istatus = NF90_INQUIRE_DIMENSION(ncid,dimids(3),len=dimlen)
  CALL handle_ncd_error(istatus,'NF90_INQ_DIMLEN in get_wrf_2d',.FALSE.)
  IF(dimlen < itime) THEN
    WRITE(6,'(a,I3,a,I3)') 'The total records number is ', dimlen,     &
                    ' however, the required time level is ',itime
    STOP 'itime_tool_large'
  END IF
  istatus = NF90_GET_VAR(ncid,varid,var2d,(/1,1,itime/),(/nx,ny,1/))
  CALL handle_ncd_error(istatus,'NF90_GET_VAR in get_wrf_2d',.FALSE.)

  RETURN
  END SUBROUTINE get_ncd_2d_int
!

SUBROUTINE get_ncd_3d(ncid,itime,varname,nx,ny,nz,var3d,istatus)
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read in a 3D array from the WRF NetCDF file
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER(INT32),          INTENT(IN)  :: ncid
  INTEGER(INT32),          INTENT(IN)  :: itime
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  INTEGER(INT32),          INTENT(IN)  :: nx
  INTEGER(INT32),          INTENT(IN)  :: ny
  INTEGER(INT32),          INTENT(IN)  :: nz
  REAL(P),         INTENT(OUT) :: var3d(nx,ny,nz)
  INTEGER(INT32),          INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  INTEGER(INT32)                    :: varid
  CHARACTER(LEN=NF90_MAX_NAME) :: namein
  INTEGER(INT32)                    :: vartype, ndims,natts,dimlen
  INTEGER(INT32)                    :: dimids(NF90_MAX_VAR_DIMS)

  INTEGER(INT32), PARAMETER         :: VAR_NOTEXIST = -1

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@END SUBROUTINE get_ncd_2d

  istatus = NF90_INQ_VARID(ncid,varname,varid)
  IF(istatus == NF90_ENOTVAR) THEN
     WRITE(6,'(3a)') ' WARNING: variable ',TRIM(varname),' does not exist.'
     var3d(:,:,:) = -999.0
     istatus = VAR_NOTEXIST
     RETURN
  END IF
  CALL handle_ncd_error(istatus,'NF90_INQ_VARID in get_ncd_3d',.FALSE.)

  istatus = NF90_INQUIRE_VARIABLE(ncid,varid,namein,vartype,ndims,dimids,natts)
  CALL handle_ncd_error(istatus,'NF90_INQ_VAR in get_ncd_3d',.FALSE.)

  IF(vartype /= NF90_FLOAT) THEN
    WRITE(6,'(3a)') 'Variable ',TRIM(varname), ' is not REAL.'
    STOP 'WRONG_VAR_TYPE'
  END IF

  IF(ndims /= 4) THEN
    WRITE(6,'(3a)') 'Variable ', TRIM(varname), ' is not a 3D array.'
    STOP 'WRONG_VAR_DIMENSIONS'
  END IF

  istatus = NF90_INQUIRE_DIMENSION(ncid,dimids(1),len=dimlen)
  CALL handle_ncd_error(istatus,'NF90_INQ_DIMLEN in get_ncd_3d',.FALSE.)
  IF(dimlen /= nx) THEN
    WRITE(6,'(3a,I3,a,I3)') 'First dimension of variable ', TRIM(varname),    &
                    ' is ',dimlen, ' and it should be ',nx
    STOP 'WRONG_DIM_length'
  END IF

  istatus = NF90_INQUIRE_DIMENSION(ncid,dimids(2),len=dimlen)
  CALL handle_ncd_error(istatus,'NF90_INQ_DIMLEN in get_ncd_3d',.FALSE.)
  IF(dimlen /= ny) THEN
    WRITE(6,'(3a,I3,a,I3)') 'Second dimension of variable ', TRIM(varname),   &
                    ' is ',dimlen, ' and it should be ',ny
    STOP 'WRONG_DIM_length'
  END IF

  istatus = NF90_INQUIRE_DIMENSION(ncid,dimids(3),len=dimlen)
  CALL handle_ncd_error(istatus,'NF90_INQ_DIMLEN in get_ncd_3d',.FALSE.)
  IF(dimlen /= nz) THEN
    WRITE(6,'(3a,I3,a,I3)') 'Third dimension of variable ', TRIM(varname),   &
                    ' is ',dimlen, ' and it should be ',nz
    STOP 'WRONG_DIM_length'
  END IF
  istatus = NF90_INQUIRE_DIMENSION(ncid,dimids(4),len=dimlen)
  CALL handle_ncd_error(istatus,'NF90_INQ_DIMLEN in get_ncd_3d',.FALSE.)
  IF(dimlen < itime) THEN
    WRITE(6,'(a,I3,a,I3)') 'The total records number is ', dimlen,     &
                    ' however, the required time level is ',itime
    STOP 'itime_tool_large'
  END IF

  istatus = NF90_GET_VAR(ncid,varid,var3d,(/1,1,1,itime/),               &
                             (/nx,ny,nz,1/))
  IF(istatus /= NF90_NOERR .OR. istatus /= NF90_EEXIST) THEN
    CALL handle_ncd_error(istatus,'NF90_GET_VARA_REAL in get_ncd_3d:'//varname,.FALSE.)
  END IF

  RETURN
END SUBROUTINE get_ncd_3d
!
SUBROUTINE get_ncd_scalar(ncid,itime,varname,var,istatus)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Read in a scalar from the WRF NetCDF file.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER(INT32),          INTENT(IN)  :: ncid
  INTEGER(INT32),          INTENT(IN)  :: itime
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  REAL(P),         INTENT(OUT) :: var
  INTEGER(INT32),          INTENT(OUT) :: istatus
!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------


  INTEGER(INT32)                    :: varid
  CHARACTER(LEN=NF90_MAX_NAME) :: namein
  INTEGER(INT32)                    :: vartype, ndims,natts,dimlen
  INTEGER(INT32)                    :: dimids(NF90_MAX_VAR_DIMS)

  !REAL :: varin(1)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  istatus = NF90_INQ_VARID(ncid,varname,varid)
  CALL handle_ncd_error(istatus,'NF90_INQ_VARID ('//TRIM(varname)//') in get_ncd_scalar',.FALSE.)

  istatus = NF90_INQUIRE_VARIABLE(ncid,varid,namein,vartype,ndims,dimids,natts)
  CALL handle_ncd_error(istatus,'NF90_INQ_VAR ('//TRIM(varname)//') in get_ncd_scalar',.FALSE.)

  IF(vartype /= NF90_FLOAT) THEN
    WRITE(6,'(3a)') 'Variable ',varname, ' is not REAL.'
    STOP 'WRONG_VAR_TYPE'
  END IF

  IF(ndims /= 1) THEN
    WRITE(6,'(3a)') 'Variable ', varname, ' is not a scalar.'
    STOP 'WRONG_VAR_DIMENSIONS'
  END IF

  istatus = NF90_INQUIRE_DIMENSION(ncid,dimids(1),len=dimlen)
  CALL handle_ncd_error(istatus,'NF90_INQ_DIMLEN in get_wrf_1d',.FALSE.)
  IF(dimlen < itime) THEN
    WRITE(6,'(a,I3,a,I3)') 'The total records number is ', dimlen,     &
                    ' however, the required time level is ',itime
    STOP 'itime_tool_large'
  END IF

  istatus = NF90_GET_VAR(ncid,varid,var,start=(/itime/))
  !istatus = NF90_GET_VAR(ncid,varid,var,(/itime/),(/1/))
  CALL handle_ncd_error(istatus,'NF90_GET_VARin get_wrf_scalar',.FALSE.)
  RETURN
END SUBROUTINE get_ncd_scalar

SUBROUTINE get_ncd_1d(ncid,itime,varname,nx,var1d,istatus)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Read in a 1D array from the WRF NetCDF file.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER(INT32),          INTENT(IN)  :: ncid
  INTEGER(INT32),          INTENT(IN)  :: itime
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  INTEGER(INT32),          INTENT(IN)  :: nx
  REAL(P),         INTENT(OUT) :: var1d(nx)
  INTEGER(INT32),          INTENT(OUT) :: istatus
!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------


  INTEGER(INT32)                    :: varid
  CHARACTER(LEN=NF90_MAX_NAME) :: namein
  INTEGER(INT32)                    :: vartype, ndims,natts,dimlen
  INTEGER(INT32)                    :: dimids(NF90_MAX_VAR_DIMS)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  istatus = NF90_INQ_VARID(ncid,varname,varid)
  CALL handle_ncd_error(istatus,'NF90_INQ_VARID ('//TRIM(varname)//') in get_wrf_1d',.FALSE.)

  istatus = NF90_INQUIRE_VARIABLE(ncid,varid,namein,vartype,ndims,dimids,natts)
  CALL handle_ncd_error(istatus,'NF90_INQURE_VARIABLE ('//TRIM(varname)//') in get_wrf_1d',.FALSE.)

  IF(vartype /= NF90_FLOAT) THEN
    WRITE(6,'(3a)') 'Variable ',varname, ' is not REAL.'
    STOP 'WRONG_VAR_TYPE'
  END IF

  IF(ndims /= 2) THEN
    WRITE(6,'(3a)') 'Variable ', varname, ' is not a 1D array.'
    STOP 'WRONG_VAR_DIMENSIONS'
  END IF

  istatus = NF90_INQUIRE_DIMENSION(ncid,dimids(1),len=dimlen)
  CALL handle_ncd_error(istatus,'NF90_INQ_DIMLEN ('//TRIM(varname)//') in get_wrf_1d',.FALSE.)
  IF(dimlen /= nx) THEN
    WRITE(6,'(3a,I3,a,I0)') 'The dimension of variable ', varname,    &
                    ' is ',dimlen, ' and it should be ',nx
    STOP 'WRONG_DIM_length'
  END IF

  istatus = NF90_INQUIRE_DIMENSION(ncid,dimids(2),len=dimlen)
  CALL handle_ncd_error(istatus,'NF90_INQ_DIMLEN ('//TRIM(varname)//') in get_wrf_1d',.FALSE.)
  IF(dimlen < itime) THEN
    WRITE(6,'(a,I3,a,I3)') 'The total records number is ', dimlen,     &
                    ' however, the required time level is ',itime
    STOP 'itime_tool_large'
  END IF

  istatus = NF90_GET_VAR(ncid,varid,var1d,(/1,itime/),(/nx,1/))
  CALL handle_ncd_error(istatus,'NF90_GET_VAR ('//TRIM(varname)//') in get_wrf_1d',.FALSE.)

  RETURN
END SUBROUTINE get_ncd_1d


end module read_wrf

