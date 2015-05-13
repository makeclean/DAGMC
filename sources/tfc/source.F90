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

  logical first_time = .true.

  implicit real(dknd) (a-h,o-z)

  if(first_time.eqv..true.)
    call setup()
  else
    jsu = 0
    wgt = 1.0
    tme = 0.0
    xxx = 350.0
    call sample_linear(40.,-40.,rand(),yyy)
    call sample_linear(400.,-400.,rand(),zzz)
    call sample(rand(),rand(),erg,dummy)
    ipt = 1 !neutrons
    icl = idum(1) ! get start icl from idum
  endif

  return
end subroutine source
