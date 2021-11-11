!
! 11/10/2021
! huanhuan ZHANG
! adapted from ?? 
!
!  interpolates fields to another grid as specified 
!  by lat., lon. values.
!  uses a Cressman interpolation with a specified 
!  search radius 
!  uses output from cressman_lookup.f90 to determine search points
!
! input :
!    lati  -  lat of source data
!    loni  -  lon of source data
!    datai -  source data
!    ni    -  length of source data
!    lato  -  lat of target data
!    lono  -  lon of target data
!    datao -  target data
!    nxo   -  first dimension of target data
!    nyo   -  2nd dimension of target data
!    r     -  search radius
!    xmiss -  value not used
!
subroutine interp_cressman(lati,loni,datai,ni,lato,lono,datao,nxo,nyo,r,xmiss)
    implicit none
    integer,parameter:: maxnum = 500
    real,   parameter:: pi = 3.1415926
    !
    integer,                intent(in   ):: ni,nxo,nyo
    real,dimension(ni),     intent(in   ):: lati,loni,datai
    real,dimension(nxo,nyo),intent(in   ):: lato,lono
    real,dimension(nxo,nyo),intent(inout):: datao
    real,      intent(in   ):: xmiss
    real,      intent(in   ):: r
    ! LOCAL VARIABLES
    integer:: ilevel,i,j,n,kk
    integer,dimension(maxnum)::indx
    real:: sum1,sum2,sum0
    real:: dim0,dam,dam1,dist
    real:: guess,xmiss0
    real:: zpiggy,zpiggy1
    real::scale0
    !
    guess=100!km
    xmiss0=xmiss
    scale0=1.
    guess=r

    open(30,file='lookup.tab1',status='unknown')
    open(31,file='lookup.tab2',status='unknown')
! reinterpolate to NMC Octagonal grid using cressman weights with a specified
! search radius
    do j=1,nyo
        do i=1,nxo
            sum1=0.0
            sum2=0.0
            !read in indices from lookup table
            read(30,*)ilevel
            if(ilevel.ne.0.)then
                read(31,*) (indx(n),n=1,ilevel)
            else
                cycle
            end if
            do kk=1,ilevel
                n = indx(kk)


            ! find distance between NMC point to be interpolated to
            ! and original point. 
                dim0=(sin(pi*lati(n)/180.)*sin(pi*lato(i,j)/180.))
                dam=(cos(pi*lati(n)/180.)*cos(pi*lato(i,j)/180.))
                dam1=lono(i,j)-loni(n)
                dam1=cos(pi*dam1/180.)
                dist=dim0+dam*dam1
                dist=acos(dist)*110.949*(180./pi)

                if(dist.le.guess .and. datai(n).ne.xmiss0) then
                    zpiggy=(guess**2)-(dist**2)
                    zpiggy1=(guess**2)+(dist**2)
                    sum0=zpiggy/zpiggy1

                    sum1=sum1+(sum0*datai(n))
                    sum2=sum2+sum0
                end if
            end do
            !
            if(ilevel.ge.1) then
                datao(i,j)=(sum1/sum2)*scale0
            else
                write(*,*) 'Insufficient search radius ...exiting'
                datao(i,j)=xmiss
            end if
        end do
    end do

    close(30)
    close(31)

    return

end subroutine interp_cressman
