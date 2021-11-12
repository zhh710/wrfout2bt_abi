use model_precision,only: print_precision
use parameters_define,only:read_nml
use read_wrf,only:lat,lon,rlats,rlons
use read_wrf,only:load_data,destory_wrfinput_array
use read_wrf,only:oz,u,v,qv,psfc,pmid,pml
use read_wrf,only:eta1,eta2
use read_wrf,only:czen,vegfrac,ivgtyp
use read_wrf,only:isltyp
use read_wrf,only:sno_full,si,pctsno
use read_wrf,only:soil_moi_full,soil_temp_full
use read_wrf,only:sfc_rough_full
use read_wrf,only:isli_full
use read_wrf,only:fact10_full
use read_wrf,only:ths
use read_wrf,only:sice
use read_wrf,only:zs_full
use read_wrf,only:u10,v10
use read_wrf,only:p_top
use read_wrf,only:grid_ratio,imp_physics
!
use  goesabi_obs,only:read_goesabi_netcdf
use  goesabi_obs,only:nabiobs,data_obsabi
use  goesabi_obs,only:read_abiobsarray_from_file
!
use module_abi,only:setuprad
!
implicit none

call print_precision()
call read_nml()
call load_data()
print*,"min/max of XLAT: ",minval(lat),maxval(lat)
print*,"min/max of XLONG: ",minval(lon),maxval(lon)
print*,"min/max of radian XLAT: ",minval(rlats),maxval(rlats)
print*,"min/max of radian XLONG: ",minval(rlons),maxval(rlons)
print*,"min/max of O3RAD: ",minval(oz),maxval(oz)
print*,"min/max of u: ",minval(u),maxval(u)
print*,"min/max of v: ",minval(v),maxval(v)
print*,"min/max of qv: ",minval(qv),maxval(qv)
print*,"min/max of psfc: ",minval(psfc),maxval(psfc)
print*,"min/max of pmid: ",minval(pmid),maxval(pmid)
print*,"min/max of pml: ",minval(pml),maxval(pml)
print*,"min/max of COSZEN: ",minval(czen),maxval(czen)
print*,"min/max of ZEN: ",acos(minval(czen))*180/3.14,acos(maxval(czen))*180/3.14
print*,"min/max of VEGFRAC: ",minval(vegfrac),maxval(vegfrac)
print*,"min/max of IVGTYP: ",minval(ivgtyp),maxval(ivgtyp)
print*,"min/max of islTYP: ",minval(isltyp),maxval(isltyp)
print*,"min/max of SNOW: ",minval(sno_full),maxval(sno_full)
print*,"min/max of SNOWH: ",minval(si),maxval(si)
print*,"min/max of SNOWC: ",minval(pctsno),maxval(pctsno)
print*,"min/max of TSK: ",minval(ths),maxval(ths)
print*,"min/max of u10: ",minval(u10),maxval(u10)
print*,"min/max of v10: ",minval(v10),maxval(v10)
print*,"min/max of SEAICE: ",minval(sice),maxval(sice)
print*,"min/max of HGT: ",minval(zs_full),maxval(zs_full)
print*,"min/max of SIOL MOISTURE: ",minval(soil_moi_full),maxval(soil_moi_full)
print*,"min/max of SIOL TEMPERATURE: ",minval(soil_temp_full),maxval(soil_temp_full)
print*,"min/max of SURFACE roughness: ",minval(sfc_rough_full),maxval(sfc_rough_full)
print*,"min/max of 10m wind factor: ",minval(fact10_full),maxval(fact10_full)
print*,"min/max of water,land ,sea ice,snow: ",minval(isli_full),maxval(isli_full)
print*,"P top:",p_top
print*,"GRID RATIO:",grid_ratio
print*,"MP_PHYSICS:",imp_physics
!
call read_goesabi_netcdf()
!call read_abiobsarray_from_file()
!print*,"min/max of tb channel 7 ",minval(data_obsabi(36,:)),maxval(data_obsabi(36,:))
!print*,"min/max of tb channel 8 ",minval(data_obsabi(37,:)),maxval(data_obsabi(37,:))
!print*,"min/max of tb channel 10 ",minval(data_obsabi(39,:)),maxval(data_obsabi(39,:))
!print*,"min/max of tb obs  lon ",minval(data_obsabi(30,:)),maxval(data_obsabi(30,:))
!print*,"min/max of tb obs  lat ",minval(data_obsabi(31,:)),maxval(data_obsabi(31,:))

!
call setuprad()
!
print*,"destory_wrfinput_array()"
call destory_wrfinput_array()

end program
