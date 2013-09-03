#ifdef _WMAP7_

module wmapwrap

use wmap_likelihood_7yr
use wmap_options
use wmap_util
logical::initswitch=.true.
contains

subroutine wmaplikeness_actual(cltt,clte,clee,clbb,likeout)

	integer i
	integer :: tt_npix, teeebb_npix
	real(8) cltt(3000),clte(3000),clee(3000),clbb(3000),ell(3000)
	real(8) like(num_WMAP)
	real(8) likeout
	
	!these settings will cause the likelihood function to consider only
	!the TT correlation spectrum, and to work directly in C_l space,
	!rather than pixel space.
	
	use_TT               = .true.
	use_TE               = .false. 
	use_lowl_TT          = .false.
	use_lowl_pol         = .false.
	use_TT_beam_ptsrc    = .true.
	
	use_gibbs=.false.
	!lowl_max=30
	
	if(initswitch) then
	call wmap_likelihood_init
	initswitch=.false.
	end if
	call wmap_likelihood_dof( tt_npix, teeebb_npix )
	call wmap_likelihood_compute(cltt,clte,clee,clbb,like)
	
	likeout=0.0
	
	likeout=2.0*like(1)+2.0*like(4)

	!do i=1,num_WMAP
	!write(*,*)'i ',i,' likepart ',like(i)
	!end do
	
	

end subroutine wmaplikeness_actual

end module

subroutine wmaplikeness(cltt,clte,clee,clbb,likeout)
use wmapwrap
real(8) cltt(3000),clte(3000),clee(3000),clbb(3000),ell(3000)
real(8) likeout
call wmaplikeness_actual(cltt,clte,clee,clbb,likeout)

end subroutine wmaplikeness

#endif
