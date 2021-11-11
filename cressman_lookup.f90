!
! 11/10/2021
! huanhuan ZHANG
! adapted from ?? 
!
!creates lookup table for cressman interpolation routine
!
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
subroutine cressman_lookup(lati,loni,datai,ni,lato,lono,nxo,nyo,r,xmiss)
    implicit none
    integer,parameter:: maxnum = 500
    real,   parameter:: pi = 3.1415926
    !
    integer,                intent(in   ):: ni,nxo,nyo
    real,dimension(ni),     intent(in   ):: lati,loni,datai
    real,dimension(nxo,nyo),intent(in   ):: lato,lono
    real,   intent(in   ):: xmiss
    real,   intent(in   ):: r
    ! LOCAL VARIABLES
    integer:: ilevel,i,j,n,kk
    integer,dimension(maxnum)::indx
    real:: sum1,sum2
    real:: dim0,dam,dam1,dist
    real:: guess,xmiss0
    character(len=*),parameter:: myname_="cressman_lookup"
    !
    print*,myname_,'*OPEN FILE to WRITE'
    xmiss0=xmiss
    guess=r

    open(30,file='lookup.tab1',status='unknown')
    open(31,file='lookup.tab2',status='unknown')
! reinterpolate to NMC Octagonal grid using cressman weights with a specified
! search radius
    do j=1,nyo
        do i=1,nxo
            sum1=0.0
            sum2=0.0
            ilevel=0
            !loop of source data
            do n=1,ni
            ! find distance between NMC point to be interpolated to
            ! and original point. 
                dim0=(sin(pi*lati(n)/180.)*sin(pi*lato(i,j)/180.))
                dam=(cos(pi*lati(n)/180.)*cos(pi*lato(i,j)/180.))
                dam1=lono(i,j)-loni(n)
                dam1=cos(pi*dam1/180.)
                dist=dim0+dam*dam1
                dist=acos(dist)*110.949*(180./pi)

                if(dist.le.guess .and. datai(n).ne.xmiss0) then
                    ilevel=ilevel+1
                    indx(ilevel)=n
                end if
            end do
            write(30,*) ilevel
            if(ilevel .ne. 0)write(31,*) (indx(kk),kk=1,ilevel)
        end do
    end do

    close(30)
    close(31)

    return

end subroutine cressman_lookup
