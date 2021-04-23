! Adapted from WRF/var/da/da_physics/da_roughness_from_lanu.inc 4/22/2021
! : Get surface roughness
! : it can be got in wrfoutput file if edit Registry/Registry.EM_COMMON
subroutine da_roughness_from_lanu(ltbl, mminlu, date, lanu, rough)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer             ,   intent(in)    :: ltbl
   character (len=*)   ,   intent(in)    :: mminlu
   character (len=19)  ,   intent(in)    :: date
   real, dimension(ims:ime,jms:jme),   intent(in)    :: lanu 
   real, dimension(ims:ime,jms:jme),   intent(out)   :: rough 

   integer                               :: LS, LC, LI, LUCATS, LuseAS, &
                                           LUMATCH, year, month, day,  &
                                           julday, Isn, io_error, &
                                           m1, m2, n1, n2 
   real                                  :: albd, slmo, sfem
   real(kind=4), dimension(50,2)         :: sfz0
   character (len=256)                   :: LUtype
   logical                               :: iexist


   read(unit=date,fmt='(I4,1x,I2,1X,I2)') year, month, day
   call da_julian_day(year,month,day,Julday, 1)
   Isn = 1
   if (JULDAY < 105 .OR. JULDAY > 288) Isn=2

   inquire (file = 'LANDUSE.TBL', exist = iexist)

   if (exist) then
      open (unit = ltbl, file = 'LANDUSE.TBL', form='formatted', &
                     action = 'read', iostat = io_error)
   else
         print*,"Cannot open file LANDUSE.TBL for conversion of roughness"
         rough = 0.
         return
   end if

   lumatch=0  

   do
      read (unit=ltbl,fmt='(A)', iostat=io_error) lutype
      if (io_error /= 0) exit
      read (unit=ltbl,fmt=*, iostat=io_error) lucats,luseas

      if (trim(lutype) == trim(mminlu)) lumatch=1 

      do LS=1,LuseAS 
         read (unit=ltbl,fmt=*)  
         do lc=1,lucats 
            if (trim(lutype) == trim(mminlu)) then 
               read (unit=ltbl,fmt=*) li, albd, slmo, sfem, sfz0(lc,ls)
               ! prevent compiler whinge
               if (albd == 0.0 .or. sfem == 0.0 .or. slmo == 0.0) then
               end if
               if (LC /= LI) then
                 call da_error(__FILE__,__LINE__, &
                   (/"Missing landuse: lc"/))
               end if
            else 
               read (unit=ltbl,fmt=*) 
            end if 
         end do 
      end do
   end do

   close (unit=ltbl)

   if (lumatch == 0)then
    call da_error(__FILE__,__LINE__,&
      (/"landuse in input file does not match lutable"/))
   end if   

   m1 = 1
   m2 = nx
   n1 = 1
   n2 = ny

   do lc = m1,m2
      do ls = n1,n2
         Li = int(lanu(lc,ls)+0.001)
         ! cm => M?
         rough(lc,ls) =  sfz0(Li,Isn)/100.0
      end do
   end do


end subroutine da_roughness_from_lanu


