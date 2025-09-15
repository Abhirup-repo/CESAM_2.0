subroutine numbmod(n)
  use pumamod
!notafmod!!  use tafmod
  implicit none
  integer,intent(out) :: n
!SB    set the number of control variables;
!SB    the value can be computed using plasim routines or modules
!Si:   from subroutine tafini:      nx_taf=NHOR
!notafmod!!     nx_taf=NHOR               ! for module tafmod
  n=NHOR
end subroutine numbmod
!----------------------------------------------------------
subroutine initmod( n, x)
  use pumamod
  implicit none
  integer,intent(in) :: n
  real,intent(inout) :: x(n)
  integer :: i
!SB initialise the model; for PlaSim, this would include a call to "prolog" and
!SB a copying from the control variables to a single array x. Ideally use
!SB something like subroutine x2mod and mod2x for this.

!SB   do i=1,n
!SB      x(i) = float(i)
!SB   enddo

  call mpstart(-1)       ! -1: Start MPI   >=0 arg = MPI_COMM_WORLD
  call setfilenames
  call opendiag
  if (mrnum == 2) then
    call mrdimensions
  endif
  call allocate_arrays
  call tafsave
  call prolog
  call tafini(n,x)
!SiS!	  x(:)=x_taf(:)
end subroutine initmod
!----------------------------------------------------------
subroutine model( n, x, fc )
  use pumamod
  use tafmod 
  use observation
  implicit none
  integer,intent(in) :: n
  real,intent(in)    :: x(n)
  real,intent(out)   :: fc
!SB make sure the state variables are reset to their initial values
!SB in order to have the same result on repeated calls to model
!SB For plasim "internal restart" was used for this.

!SB unpack control variables from control vector x, e.g. by a call to something
!SB like subroutine x2mod
!SB   a(:) = x(:)
!SB run the model
!SB   a(:) = a(:)**2
!SB compute the cost function
!SB   fc = sum(a)
#define STATIC
#ifndef STATIC
! file tapes:  
!$TAF INIT euler     = "taftapes/euler"
!$TAF INIT main      = "taftapes/main"
!$TAF INIT main_nlev = "taftapes/mainn"
!$TAF INIT obsn      = "taftapes/obsn"
!$TAF INIT obs1      = "taftapes/obs1"
#else
! static memory tapes
!$TAF INIT euler      = static, nkits
!$TAF INIT main       = static,  n_run_months*31+n_run_days*ntspd+n_run_steps
!$TAF INIT main_nlev  = static, (n_run_months*31+n_run_days*ntspd+n_run_steps) * nlev 
!$TAF INIT obs1       = static, (nstep-nstep0+1)
!$TAF INIT obsn       = static, nobs*(nstep-nstep0+1)
#endif
! FO: variable size ("dynamic") tape in memory for do while counters (last tape defined is automatically used):
!$TAF INIT dowtape    = dynamic

! to initialize model for multiple runs, e.g. for gradient check
  call dealloc_fields
  call tafreset
  call prolog
  nstep2=0
!! only in initmod!  call tafini(n,x)
  call xtaf2mod(n,x) 
  call reinitobs
!$TAF STORE cost=obs1
  call master
  call mod2ytaf(fc) 
end subroutine model
!----------------------------------------------------------
subroutine postmod(n,x,fc)
  integer,intent(in) :: n
  real,intent(inout)    :: x(n)
  real,intent(in)   :: fc

  call tafstop(n,x) 
  call epilog
  call mpstop
  print*,"The value of the cost function is",fc
!SB call epilog or whatever shall happen after a cost function evaluation
end subroutine postmod
!----------------------------------------------------------
subroutine postadm(n,x,fc,adx)
  integer,intent(in) :: n
  real,intent(in)    :: x(n), adx(n)
  real,intent(in)   :: fc
  integer :: i

integer, parameter :: NLAT=32                ! nr of latitudes  for t21
integer, parameter :: NLON=NLAT*2            ! nr of longitudes for t21
integer :: ih(8)=(/100,0,0,0,NLON,NLAT,0,0/)  ! srv header
real  :: zmask       ! mask field 
open(31,file="adx.srv",form='unformatted')      ! output
write(31) ih (:)
write(31) adx(:) 
close(31)

open(32,file="adx.sra",form='formatted')      ! output
write(32,*) ih(:)
write(32,'(8E12.4)') adx(:) 
close(32)

! call epilog or whatever shall happen after a cost function evaluation
  print*,"The norm of the gradient is",sqrt(dot_product(adx,adx))
end subroutine postadm




