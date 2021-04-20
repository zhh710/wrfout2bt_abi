use model_precision,only: print_precision
use parameters_define,only:read_nml
use read_wrf,only:load_data,destory_wrfinput_array
use read_wrf,only:oz,u,v,qv,psfc,pmid,pml
use read_wrf,only:eta1,eta2
use read_wrf,only:czen,vegfrac,ivgtyp
use read_wrf,only:sno,si
use read_wrf,only:p_top
call print_precision()
call read_nml("input.namelist")
call load_data()
print*,"MODEL TOP PRESSURE: ",p_top
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
print*,"min/max of SNOW: ",minval(sno),maxval(sno)
print*,"min/max of SNOWH: ",minval(si),maxval(si)

print*,"P top:",p_top
call destory_wrfinput_array()

end program
