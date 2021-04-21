Module goesabi_obs
    ! ADAPTED FROM read_goesabi_netcdf.f90 of GSI
! isatid    = 1     ! index of satellite id
! itime     = 2     ! index of analysis relative obs time 
! ilon      = 3     ! index of grid relative obs location (x)
! ilat      = 4     ! index of grid relative obs location (y)
! ilzen_ang = 5     ! index of local (satellite) zenith angle (radians)
! ilazi_ang = 6     ! index of local (satellite) azimuth angle (radians)
! icloud_frc= 7     ! index of clear sky amount,1.0=>clear
! iscan_pos = 8     ! index of scan position
! iszen_ang = 9     ! index of solar zenith angle (degrees)
! isazi_ang = 10    ! index of solar azimuth angle (degrees)
! ifrac_sea = 11    ! index of ocean percentage
! ifrac_lnd = 12    ! index of land percentage
! ifrac_ice = 13    ! index of ice percentage
! ifrac_sno = 14    ! index of snow percentage
! its_sea   = 15    ! index of ocean temperature
! its_lnd   = 16    ! index of land temperature
! its_ice   = 17    ! index of ice temperature
! its_sno   = 18    ! index of snow temperature
! itsavg    = 19    ! index of average skin temperature
! ivty      = 20    ! index of vegetation type
! ivfr      = 21    ! index of vegetation fraction
! isty      = 22    ! index of soil type
! istp      = 23    ! index of soil temperature
! ism       = 24    ! index of soil moisture
! isn       = 25    ! index of snow depth
! izz       = 26    ! index of surface height
! idomsfc   = 27    ! index of dominate surface type
! isfcr     = 28    ! index of surface roughness
! iff10     = 29    ! index of ten meter wind factor
! ilone     = 30    ! index of earth relative longitude (degrees)
! ilate     = 31    ! index of earth relative latitude (degrees)
! iassim    = 32    ! 0:not assimilate

    use model_precision,only:P,INT32,DP
    use parameters_define,only:abi_file,abiobs_mid_file,iuseabi,nchanl_abi
    use parameters_define,only:regional
    use parameters_define,only:deg2rad
    use netcdf
    implicit none
    INTEGER(INT32):: isatid    = 1     ! index of satellite id
    INTEGER(INT32):: itime     = 2     ! index of analysis relative obs time 
    INTEGER(INT32):: ilon      = 3     ! index of grid relative obs location (x)
    INTEGER(INT32):: ilat      = 4     ! index of grid relative obs location (y)
    INTEGER(INT32):: ilzen_ang = 5     ! index of local (satellite) zenith angle (radians)
    INTEGER(INT32):: ilazi_ang = 6     ! index of local (satellite) azimuth angle (radians)
    INTEGER(INT32):: icloud_frc= 7     ! index of clear sky amount,1.0=>clear
    INTEGER(INT32):: iscan_pos = 8     ! index of scan position
    INTEGER(INT32):: iszen_ang = 9     ! index of solar zenith angle (degrees)
    INTEGER(INT32):: isazi_ang = 10    ! index of solar azimuth angle (degrees)
    INTEGER(INT32):: ifrac_sea = 11    ! index of ocean percentage
    INTEGER(INT32):: ifrac_lnd = 12    ! index of land percentage
    INTEGER(INT32):: ifrac_ice = 13    ! index of ice percentage
    INTEGER(INT32):: ifrac_sno = 14    ! index of snow percentage
    INTEGER(INT32):: its_sea   = 15    ! index of ocean temperature
    INTEGER(INT32):: its_lnd   = 16    ! index of land temperature
    INTEGER(INT32):: its_ice   = 17    ! index of ice temperature
    INTEGER(INT32):: its_sno   = 18    ! index of snow temperature
    INTEGER(INT32):: itsavg    = 19    ! index of average skin temperature
    INTEGER(INT32):: ivty      = 20    ! index of vegetation type
    INTEGER(INT32):: ivfr      = 21    ! index of vegetation fraction
    INTEGER(INT32):: isty      = 22    ! index of soil type
    INTEGER(INT32):: istp      = 23    ! index of soil temperature
    INTEGER(INT32):: ism       = 24    ! index of soil moisture
    INTEGER(INT32):: isn       = 25    ! index of snow depth
    INTEGER(INT32):: izz       = 26    ! index of surface height
    INTEGER(INT32):: idomsfc   = 27    ! index of dominate surface type
    INTEGER(INT32):: isfcr     = 28    ! index of surface roughness
    INTEGER(INT32):: iff10     = 29    ! index of ten meter wind factor
    INTEGER(INT32):: ilone     = 30    ! index of earth relative longitude (degrees)
    INTEGER(INT32):: ilate     = 31    ! index of earth relative latitude (degrees)
    INTEGER(INT32):: iassim    = 32    ! 0:not assimilate
    !
    real(P),allocatable,dimension(:,:)::abi_data
    contains
    !
    subroutine read_goesabi_netcdf()
        implicit none
        integer(INT32),parameter:: maxinfo=35
        real(P),parameter:: r360=360.0_P
        real(P),parameter:: tbmin=50.0_P
        real(P),parameter:: tbmax=550.0_P
        !
        character(4)  idate5s(5)
        INTEGER(INT32)::idate5(5),mins_an
        !
        real(DP),dimension(12) :: dataabi
        !------------------
        !  NETCDF-RELATED
        !------------------
        INTEGER(INT32)   :: ncdfID
        INTEGER(INT32)   :: istatus
        INTEGER(INT32)   :: datestlen, varID, DimID
        INTEGER(INT32)   :: vardim, natts
        INTEGER(INT32)   :: vartype, id_time
        INTEGER(INT32)   :: numdim, numvars, numatt
        INTEGER(INT32)   :: nn
        !------------------
        !  SATELLITE DATA
        !------------------
        REAL(P), ALLOCATABLE, DIMENSION( :, : )    :: tb
        REAL(P), ALLOCATABLE, DIMENSION( : )       :: sataz
        REAL(P), ALLOCATABLE, DIMENSION( : )       :: solaz
        REAL(P), ALLOCATABLE, DIMENSION( : )       :: vza
        REAL(P), ALLOCATABLE, DIMENSION( : )       :: sza
        REAL(P), ALLOCATABLE, DIMENSION( : )       :: lon
        REAL(P), ALLOCATABLE, DIMENSION( : )       :: lat
        INTEGER(P), ALLOCATABLE, DIMENSION( : )    :: cmask
        CHARACTER(LEN=10)::Satellite
        INTEGER(INT32)   :: kidsat
        !
        INTEGER(INT32)::nreal,nele
        INTEGER(INT32)::ii,itx,cc,i ! loop index
        INTEGER(INT32)::ndata !number of profiles retained for further processing
        !
        real(P)::thislon,thislat
        real(P)::dlon_earth,dlat_earth
        real(P)::dlon,dlat
        !
        LOGICAL::outside
        !!!!!!!!!!!!
        ! OPEN NETCDF FILE
        istatus = nf90_open(TRIM(abi_file), NF90_NOWRITE, ncdfID)
        print*, '*** OPENING GOES-ABI RADIANCE NETCDF FILE', istatus
        ! Get date information
        istatus = nf90_get_att( ncdfID, nf90_global, 'year', idate5s(1) )
        istatus = nf90_get_att( ncdfID, nf90_global, 'month', idate5s(2) )
        istatus = nf90_get_att( ncdfID, nf90_global, 'day', idate5s(3) )
        istatus = nf90_get_att( ncdfID, nf90_global, 'hour', idate5s(4) )
        istatus = nf90_get_att( ncdfID, nf90_global, 'minute', idate5s(5) )
        istatus = nf90_get_att( ncdfID, nf90_global, 'satellite', Satellite )
        read(idate5s(:) , *) idate5(:)
        print*,idate5
        print*,"Satellite:",satellite
        !
        ! Get Dimension Info (1-D)
        istatus = nf90_inq_varid( ncdfID, 'numobs', varID )
        istatus = nf90_get_var( ncdfID, varID, nn )
        !! Allocate data arrays
        ALLOCATE( lat( nn ) ) 
        ALLOCATE( lon( nn ) )
        ALLOCATE( sza( nn ) )
        ALLOCATE( vza( nn ) )
        ALLOCATE( solaz( nn ) )
        ALLOCATE( sataz( nn ) )
        ALLOCATE( cmask( nn ) )
        ALLOCATE( tb( nn, nchanl_abi ) )
        ! Get useful data arrays
        ! LAT
        istatus = nf90_inq_varid( ncdfID, 'lat', varID )
        istatus = nf90_get_var( ncdfID, varID, lat )
        ! LON
        istatus = nf90_inq_varid( ncdfID, 'lon', varID )
        istatus = nf90_get_var( ncdfID, varID, lon )
        ! VZA
        istatus = nf90_inq_varid( ncdfID, 'vza', varID )
        istatus = nf90_get_var( ncdfID, varID, vza )
        ! SZA
        istatus = nf90_inq_varid( ncdfID, 'sza', varID )
        istatus = nf90_get_var( ncdfID, varID, sza )
        ! Satellite azimuth
        istatus = nf90_inq_varid( ncdfID, 'solaz', varID )
        istatus = nf90_get_var( ncdfID, varID, solaz )
        ! Solar azimuth
        istatus = nf90_inq_varid( ncdfID, 'sataz', varID )
        istatus = nf90_get_var( ncdfID, varID, sataz )
        ! Brightess temperature datea
        istatus = nf90_inq_varid( ncdfID, 'value', varID )
        istatus = nf90_get_var( ncdfID, varID, tb )
        ! Cloud Mask (0 = clear, 1 = cloud)
        istatus = nf90_inq_varid( ncdfID, 'cloud_mask', varID )
        istatus = nf90_get_var( ncdfID, varID, cmask )
        ! CLOSE NETCDF FILE
        istatus = nf90_close( ncdfID )

        !SAT ID
        if (trim(Satellite) == 'G16)') kidsat = 270
        if (trim(Satellite) == 'G17)') kidsat = 271
        !! Allocate arrays to hold all data for given satellite 
        nreal = maxinfo
        nele = nreal + nchanl_abi
        allocate(abi_data(nele,nn),stat=istatus)
        !
        itx=1
        ndata=0
        ! CHECK TIME WINDOW: ALL OBSERVATIONS HAVE THE SAME TIME
        ! call w3fs21(idate5,mins_an) !mins_an -integer number of mins snce 01/01/1978
        !
        !LOOP OVER DATA
        outside = .TRUE.
        DO i=1,nn
            !Check if there is any missing obs. Skip all channels if this is the case
            do cc=1,nchanl_abi
                if (tb(i,cc) /= tb(i,cc)) then
                    cycle
                endif
            enddo
            ! Convert obs location from degrees to radians
            thislon=lon(i)
            thislat=lat(i)
            if (thislon >= 360.) thislon=thislon-360.
            if (thislon < 0.)  thislon=thislon+360.
            dlon_earth=thislon*deg2rad
            dlat_earth=thislat*deg2rad
            if (regional)then
               ! call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)
                if (outside)  cycle
            else

            endif
            !
            !!!!!!!!!!!!!!!!!!!!!!!!SFC
            ! Locate the observation on the analysis grid.  Get sst and land/sea/ice  mask.  
            !            call deter_sfc(dlat,dlon,dlat_earth,dlon_earth,t4dv,isflg,idomsfc,sfcpct, &
            !    ts,tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10,sfcr)

        ENDDO



    end subroutine read_goesabi_netcdf


END Module goesabi_obs
