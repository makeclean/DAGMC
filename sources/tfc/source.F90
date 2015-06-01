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
  double precision :: dummy
  logical, save :: first_time = .true.

  double precision, save :: y_min,y_max
  double precision, save :: z_min,z_max

  if(first_time .eqv. .true.)then
    call setup()
    y_min = -40.
    y_max = 40.
    z_min = -399.
    z_max = 399.
    first_time = .false.
  endif

  jsu = 0
  wgt = 1.0
  tme = 0.0
  xxx = 350.0
  
  call linear_sample(y_max,y_min,rang(),yyy)
  call linear_sample(z_max,z_min,rang(),zzz)
  call sample(rang(),rang(),erg,dummy)
  ipt = 1 !neutrons
  icl = idum(1) ! get start icl from idum

  if(wgt.le.0.0)then
     write(*,*) xxx, yyy ,zzz ,erg ,ipt , icl
  endif

  return
end subroutine source
