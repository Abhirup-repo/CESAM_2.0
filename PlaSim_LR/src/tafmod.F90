module tafmod
  implicit none
!
!notafmod!!      integer :: nx_taf
!notafmod!!      integer :: ny_taf
!
!SiS!      real,allocatable :: x_taf(:)
!SiS!      real,allocatable :: y_taf(:)
!notafmod!!      real,allocatable :: x(:)
!!always 1 value:      real,allocatable :: fc(:)
!notafmod!!      real             :: fc
  real,allocatable :: z_mask(:,:)
!!
end module tafmod
!----------------------------------------------------------
subroutine tafini(n,x)
  use pumamod
  use oceanmod,only : ysst
  use radmod, only: th2oc
  use tafmod
  use observation

  implicit none
  integer,intent(in) :: n
  real,intent(inout)   :: x(n)
  integer :: ierr
  integer,dimension(8) :: ih
  real :: z_mask4(NLON,NLAT) ! for compilation with -r8
  if (ncost>0) then
    call initobs 
  endif

  allocate(z_mask(NLON,NLAT))
  open(21,file='mask.srv',form='unformatted',status='old',iostat=ierr)
  if (ierr /= 0) stop 'no inputfile for y_mysk available'
  read(21,iostat=ierr)ih(:)
  if (ierr /= 0) stop 'inputfile for mask is empty'
  write(*,*)'lon,lat=',ih(5),ih(6)
  if(ih(5)*ih(6) /= NHOR) stop 'y_mask input file has wrong resolution'
  read(21)z_mask4
  z_mask = z_mask4
  close(21)

  x(:)=ysst(:,1)       ! real :: ysst(NHOR,NLEV_OCE) = 0.  ! temperature (K)
  
!test SiS	  write(*,*)'x_taf(:)=ysst(:,1),ysst(:NLON,1)',ysst(1:NLON,1)
!
  return
end subroutine tafini
!-------------------------------------------------------
subroutine xtaf2mod(n,x) 
  use pumamod
  use oceanmod,only : ysst
  use tafmod
  implicit none
  integer,intent(in) :: n
  real,intent(in)    :: x(n)
!
!     relate input (x) and model variables 
!
  ysst(:,1)=x(:) !
  return
end subroutine xtaf2mod
!-------------------------------------------------------
subroutine tafstep
  use pumamod
  use tafmod
  implicit none
!
  return
end subroutine tafstep
!-------------------------------------------------------
subroutine mod2ytaf(fc)
  use pumamod
  use tafmod
  use observation 
  implicit none
  real    :: fc, dvar_mean(NHOR)
  integer :: inr_mask,idx
  logical :: lmask(NLON,NLAT)


  fc=cost

  lmask = z_mask(:,:)>0
   write(*,*)'subroutine mod2ytaf, inr_mask=',count(lmask) !SiS

  if (step_flag>=ntspd) then
    print *, "The cost is computed as the mean of the last day"
    dvar_mean=(dvar/(ntspd))
  else 
    print *, "The cost is computed as the mean of the last time step"
    dvar_mean=dtsa(:)
  end if  

  ! fc=dvar_mean(323)
  !  dvar_mean=(dvar/ntspd)
  ! fc=sum( ( real(reshape(dvar_mean,(/NLON,NLAT/))*z_mask,kind=8))) / sum(z_mask)


   print *, "the value of the cost function in mod2ytaf",fc
      
  !print *, "The value of the step_flag", (step_flag-ntspd)
  !print *, "The NTSPD from TAFMOD", ntspd
  !print *, "The value of the mean depended variable", dvar_mean 

  return
end subroutine mod2ytaf
!-------------------------------------------------------
subroutine tafstop(n,x)
  use pumamod
  use tafmod
  implicit none
  integer,intent(in) :: n
  real,intent(inout) :: x(n)
!
!     subroutine to work with tl results
!     (to be completed after applying taf)
!
  open(44,file='xout_dtsa_dt.srv',form='unformatted')
  write(44) 169,0,-1,0,NLON,NLAT,0,0
!SiS!      write(44) x_taf(1:NLON*NLAT)
  write(44) x(1:NLON*NLAT)
  close(44)

  return
end subroutine tafstop
!-------------------------------------------------------
module tafsmod
  use pumamod,only: nkits,t0,tfrc,tdissq,tdisst,tdissz,tdissd,psurf &
       &,dtns,dtep,dttrp,mypid,nroot,NLEV, dtsa, dprc, dprl, devap &
       &,NLEP, du, dv &
       &,gudt,gvdt,gtdt,gqdt,dudt,dvdt,dtdt,dqdt
  use radmod, only: dftu
  use miscmod,only: tnudget
  use landmod,only: tau_veg,tau_soil
  use oceanmod,only: ysst,nhor,nlev_oce, ydsst, yqhd, yicec ! SiS ,ysst_ad
  use surfmod,only: nsurnum
  use observation, only:cost
  implicit none
!
!     module containing variables which need to be saved & reseted 
!     to enable multible prolog, master loops (e.g. for gradcheck) 
!
  integer :: nkits_s
  integer :: nsurnum_s
  real :: psurf_s,dtns_s,dtep_s,dttrp_s
  real :: tnudget_s,tau_veg_s,tau_soil_s
  real :: t0_s(NLEV),tfrc_s(NLEV)
  real :: tdissq_s(NLEV),tdisst_s(NLEV),tdissz_s(NLEV),tdissd_s(NLEV)
  real :: ysst_s(nhor,nlev_oce)
  real :: cost_s


end module tafsmod
!-------------------------------------------------------
subroutine tafsave
  use tafsmod
  implicit none
!
!     subroutine to save variables which need to be saved & reseted
!     to enable multible *prolog, master* loops (e.g. for gradcheck).
!     to be called before 1st call to *prolog*
!
  if(mypid==NROOT) then
    nkits_s     = nkits
    t0_s(:)     = t0(:)
    tfrc_s(:)   = tfrc(:)
    tdissq_s(:) = tdissq(:)
    tdisst_s(:) = tdisst(:)
    tdissz_s(:) = tdissz(:)
    tdissd_s(:) = tdissd(:)
    psurf_s     = psurf
    dtns_s      = dtns
    dtep_s      = dtep
    dttrp_s     = dttrp
    tnudget_s   = tnudget
    tau_veg_s   = tau_veg
    tau_soil_s  = tau_soil
    ysst_s(:,:) = ysst(:,:)  !SiS 
    nsurnum_s   = nsurnum    !SiS
    cost_s      =cost 
!!SiS  ysst_ad_s(:,:) = ysst_ad(:,:)
  endif

  return
end subroutine tafsave
!-------------------------------------------------------

subroutine tafreset
  use tafsmod
  implicit none
!
!     subroutine to reset variables which need to be saved & reset
!     to enable multiple *prolog, master* loops (e.g. for gradcheck).
!     to be called before all calls to *prolog* except the first,
!     where tafsave needs to be called
!
  if(mypid==NROOT) then
    nkits     = nkits_s
    t0(:)     = t0_s(:)
    tfrc(:)   = tfrc_s(:)
    tdissq(:) = tdissq_s(:)
    tdisst(:) = tdisst_s(:)
    tdissz(:) = tdissz_s(:)
    tdissd(:) = tdissd_s(:)
    psurf     = psurf_s
    dtns      = dtns_s
    dtep      = dtep_s
    dttrp     = dttrp_s
    tnudget   = tnudget_s
    tau_veg   = tau_veg_s
    tau_soil  = tau_soil_s
    ysst(:,:) = ysst_s(:,:) !SiS
    nsurnum   = nsurnum_s    !SiS
    cost=cost_s
!!SiS  ysst_ad(:,:) = ysst_ad_s(:,:)
    dtsa = 0. ! used in diag output
    dprc = 0. ! used in diag output
    dprl = 0. ! used in diag output
    devap = 0.! used in diag output
    dftu = 0. ! used in diag output
    ydsst = 0. !!
    yqhd  = 0. !!
    yicec = 0.
    du(:,1:NLEP) = 0.
    dv(:,1:NLEP) = 0.
    gudt(:,:) = 0.
    gvdt(:,:) = 0.
    gtdt(:,:) = 0.
    gqdt(:,:) = 0.
    dudt(:,:) = 0.
    dvdt(:,:) = 0.
    dtdt(:,:) = 0.
    dqdt(:,:) = 0.
  endif

  return
end subroutine tafreset
!-------------------------------------------------------
subroutine dealloc_fields
  use pumamod
  use tafmod
  implicit none
!
!     subroutine to deallocate fields
!     to enable multible *prolog, master* loops (e.g. for gradcheck).
!     to be called before all calls to *prolog* except the first.
!

  if (allocated(meed)) then       !allocated in initrandom
    write(*,*)'deallocate(meed)'
    deallocate(meed)
  endif
!      if (allocated(z_mask)) then
!	    write(*,*)'deallocate(z_mask)'
!	    deallocate(z_mask)
!      endif
end subroutine dealloc_fields
!-------------------------------------------------------


