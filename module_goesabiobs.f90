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

    use model_precision,only:P,INT32,DP,r_kind,i_kind
    use parameters_define,only:abi_file,abiobs_mid_file,iuseabi,nchanl_abi
    use parameters_define,only:regional
    use parameters_define,only:deg2rad,rad2deg
    use gridmod,only:tll2xy
    use read_wrf,only:iwinbgn
    use w3nco,only:W3FS21
    use deter_sfc_mod,only:deter_sfc
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
    integer(INT32)::nabiobs
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
        real(P),allocatable,dimension(:,:)::data_abi
        !
         character(20)  :: isis="abi_g16"
         character(10)  :: obstype="abi"
        !
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
        INTEGER(INT32)::ii,itx,cc,i,n,k ! loop index
        INTEGER(INT32)::ndata !number of profiles retained for further processing
        !
        real(P)::thislon,thislat
        real(P)::dlon_earth,dlat_earth
        real(P)::dlon,dlat
        !
        LOGICAL::outside
        !sfc
        integer(INT32)::isflg,idomsfc
        real(P) :: tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10
        real(P) :: t4dv ! obstime observation time relative to analysis time
        real(P) :: sfcr !
        real(r_kind),dimension(0:3):: sfcpct
        real(r_kind),dimension(0:3):: ts

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
        allocate(data_abi(nele,nn),stat=istatus)
        !
        itx=1
        ndata=0
        ! CHECK TIME WINDOW: ALL OBSERVATIONS HAVE THE SAME TIME
        call w3fs21(idate5,mins_an) !mins_an -integer number of mins snce 01/01/1978
        t4dv = real(mins_an - iwinbgn)
        print*,"OBS time - ANA time (mins):",t4dv
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
                call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)
                if (outside)  cycle
            else

            endif
            !
            !!!!!!!!!!!!!!!!!!!!!!!!SFC
            ! Locate the observation on the analysis grid.  Get sst and land/sea/ice  mask.  
                        call deter_sfc(dlat,dlon,dlat_earth,dlon_earth,t4dv,isflg,idomsfc,sfcpct, &
                ts,tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10,sfcr)
            !       Transfer information to work array
            data_abi( 1,itx) = kidsat                    ! satellite id
            data_abi( 2,itx) = t4dv                       ! analysis relative time
            data_abi( 3,itx) = dlon                       ! grid relative longitude
            data_abi( 4,itx) = dlat                       ! grid relative latitude
            data_abi( 5,itx) = vza(i)*deg2rad             ! satellite zenith angle (radians)
            data_abi( 6,itx) = sataz(i)*deg2rad           ! satellite azimuth angle (radians)
            if ( cmask(i) == 0 ) data_abi( 7,itx) = 1.0    !clear
            if ( cmask(i) == 1 ) data_abi( 7,itx) = 0.0    !cloud
            data_abi( 8,itx) =  0.                     ! integer scan position,not use
            data_abi( 9,itx) = sza(i)                     ! solar zenith angle
            data_abi(10,itx) = solaz(i)                   ! solar azimuth angle
            data_abi(11,itx) = sfcpct(0)                  ! sea percentage of
            data_abi(12,itx) = sfcpct(1)                  ! land percentage
            data_abi(13,itx) = sfcpct(2)                  ! sea ice percentage
            data_abi(14,itx) = sfcpct(3)                  ! snow percentage
            data_abi(15,itx)= ts(0)                       ! ocean skin temperature
            data_abi(16,itx)= ts(1)                       ! land skin temperature
            data_abi(17,itx)= ts(2)                       ! ice skin temperature
            data_abi(18,itx)= ts(3)                       ! snow skin temperature
            data_abi(19,itx)= tsavg                       ! average skin temperature
            data_abi(20,itx)= vty                         ! vegetation type
            data_abi(21,itx)= vfr                         ! vegetation fraction
            data_abi(22,itx)= sty                         ! soil type
            data_abi(23,itx)= stp                         ! soil temperature
            data_abi(24,itx)= sm                          ! soil moisture
            data_abi(25,itx)= sn                          ! snow depth
            data_abi(26,itx)= zz                          ! surface height
            data_abi(27,itx)= idomsfc + 0.001_r_kind      ! dominate surface type
            data_abi(28,itx)= sfcr                        ! surface roughness
            data_abi(29,itx)= ff10                        ! ten meter wind factor
            data_abi(30,itx)= dlon_earth*rad2deg          ! earth relative longitude (degrees)
            data_abi(31,itx)= dlat_earth*rad2deg          ! earth relative latitude (degrees)

            data_abi(32,itx)=  1.         ! 
            !       Transfer observation location and other data to local arrays
            do k=1,nchanl_abi
                data_abi(k+nreal,itx) = tb(i,k)
            enddo
            ndata = ndata + 1
            itx = itx + 1

        ENDDO
        ! store in file
        if(ndata >0 )then
            open(998,file=abiobs_mid_file,form="unformatted",status="replace")
            write(998)obstype,isis,nreal,nchanl_abi,ilat,ilon,ndata
            write(998)((data_abi(k,n),k=1,nele),n=1,ndata)
            close(998)
        endif
        !
        nabiobs = ndata
        !
        DEALLOCATE( lat )
        DEALLOCATE( lon )
        DEALLOCATE( sza )
        DEALLOCATE( vza )
        DEALLOCATE( solaz )
        DEALLOCATE( sataz)
        DEALLOCATE( cmask)
        DEALLOCATE( tb )
        if(allocated(data_abi))deallocate(data_abi)

    end subroutine read_goesabi_netcdf
    !
    ! read abi obs array from file "abiobs_mid_file"
    subroutine read_abiobsarray_from_file(dy_array)
        implicit none
        character(10)::obstype
        character(20)::isis
        integer(INT32)::ndata,nchanl,nreal,ilat,ilon
        integer(INT32)::k,n,nele
        integer(INT32)::istatus
        real(r_kind),allocatable,dimension(:,:)::dy_array
        !
        open(998,file=abiobs_mid_file,form="unformatted",status="old")
        read(998)obstype,isis,nreal,nchanl,ilat,ilon,ndata
        nele = nreal + nchanl
        allocate(dy_array(nele,ndata),stat=istatus)
        read(998)((dy_array(k,n),k=1,nele),n=1,ndata)
        close(998)
        
    end subroutine


END Module goesabi_obs
