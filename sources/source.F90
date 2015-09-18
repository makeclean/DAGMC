!+ $Id: source.F90,v 1.2 2006/10/03 02:10:16 mashnik Exp $
! Copyright LANS/LANL/DOE - see file COPYRIGHT_INFO

subroutine source
  ! dummy subroutine.  aborts job if source subroutine is missing.
  ! if nsr=0, subroutine source must be furnished by the user.
  ! at entrance, a random set of uuu,vvv,www has been defined.  the
  ! following variables must be defined within the subroutine:
  ! xxx,yyy,zzz,icl,jsu,erg,wgt,tme and possibly ipt,uuu,vvv,www.
  ! subroutine srcdx may also be needed.
  use mcnp_global
  use mcnp_debug

  implicit real(dknd) (a-h,o-z)
  logical, save :: first_run = .true.
  integer,save :: icl_tmp
  integer :: j
  integer :: ec ! error code
  real(dknd) :: randoms(6) ! randoms to pass to functions
  real(dknd) :: x,y,z,w,e
  
  if(first_run.eqv..true.)then
  ! sets up the plasma source
     if((idum(1).eq.1).or.(idum(1).eq.0))then
        write(*,*) 'Using the parametric source'
        write(*,*) 'rdum(1) = ion density in the pedistal region'
        write(*,*) 'rdum(2) = ion density in the seperatrix region'
        write(*,*) 'rdum(3) = ion density at the origin region'
        write(*,*) 'rdum(4) = ion temperature in the pedistal region'
        write(*,*) 'rdum(5) = ion temperature in the seperatrix'
        write(*,*) 'rdum(6) = ion temperature at the origin '
        write(*,*) 'rdum(7) = pedistal radius'
        write(*,*) 'rdum(8) = ion density peaking factor'
        write(*,*) 'rdum(9) = ion temp peaking factor'
        write(*,*) 'rdum(10) = minor radius'
        write(*,*) 'rdum(11) = major radius'
        write(*,*) 'rdum(12) = elongation'
        write(*,*) 'rdum(13) = triangularity'
        write(*,*) 'rdum(14) = shafranov shift'
        write(*,*) 'rdum(15) = start angle'
        write(*,*) 'rdum(16) = end angle'
        write(*,*) 'idum(2) = number of bins'        
        write(*,*) 'idum(3) = the start cell containing the plasma'        
     else if (idum(1).eq.2)then
        write(*,*) 'Using the RZ source'
        write(*,*) 'rdum(1) = start angle'
        write(*,*) 'rdum(2) = end angle'
        write(*,*) 'idum(3) = the start cell containing the plasma'        
     else if (idum(1).eq.3)then
        write(*,*) 'Using the R profile source'
        write(*,*) 'rdum(1) = minor radius'
        write(*,*) 'rdum(2) = major radius'
        write(*,*) 'rdum(3) = elongation'
        write(*,*) 'rdum(4) = triangularity'
        write(*,*) 'rdum(5) = shafranov shift'
        write(*,*) 'rdum(6) = start angle'
        write(*,*) 'rdum(7) = end angle'
        write(*,*) 'idum(3) = the start cell containing the plasma'        
     else
        write(*,*) 'Uknown source idum/rdum combination'
        stop
     endif

     ! call the setup function
     call setup(idum,rdum,ec)
     if(ec.eq.0)then
        write(*,*) 'Failed in setup'
     endif
     
     ! set the first run is no longer true
     first_run = .false.
     
     do m=1,mxa
        if(idum(3).eq.ncl(m))then
           icl_tmp = m
           cycle
        endif
     enddo

  endif

! resample source point
200 continue     
  ! sample source pos
  
  randoms(1) = rang()
  randoms(2) = rang()
  randoms(3) = rang()
  randoms(4) = rang()
  randoms(5) = rang()
  randoms(6) = rang()

  call sample(randoms,x,y,z,w,e,ec)
  jsu=0
  tme=0.0

  if(ec.eq.0)then
     write(*,*) 'Failed to sample from plasma_source'
     stop
  endif

  ! turn the coordinates into cm
  ipt = 1
  xxx = x*100
  yyy = y*100
  zzz = z*100
  erg = e
  wgt = w

  

  ! if the source cell was set make sure the 
  call chkcel(icl_tmp,2,j)
  if(j.ne.0)then
     icl_tmp=i_find_cell()
  endif

  icl = icl_tmp
  !write(*,*) icl,sqrt(xxx**2+yyy**2+zzz**2),xxx,yyy,zzz,erg,wgt
  ! 
  return
end subroutine source

integer function i_find_cell()
! This function to determines the current MCNP cell index location and exits if
! no valid cell is found. This only works if there are no repeated geometries or
! universes present in the model.
    use mcnp_global
    use mcnp_debug
    implicit none
    ! xxx,yyy,zzz are global variables
    ! mxa is global
    integer :: i ! iterator variable
    integer :: j ! tempory cell test
    integer :: icl_tmp ! temporary cell variable

    icl_tmp = -1

    do i = 1, mxa
      call chkcel(i, 2, j)
      if (j .eq. 0) then
         ! valid cel set
         icl_tmp = i
         exit
      endif
    enddo
    ! icl is now set

    if(icl_tmp .le. 0) then
      write(*,*) 'Nonsense cell number stopping'
      stop
    endif

    i_find_cell = icl_tmp
end function i_find_cell
