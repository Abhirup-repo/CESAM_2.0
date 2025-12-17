!-*-f90-*-

! This file is cpp-preprocessed. Be sure to edit mo_mapping.F90

! The following two macros are mutually exclusive. First one wins, usually.
#define PARAMS
#undef VECSTSRF
!#define NX_OCEAN 118800
! 15 level ocn
#define NX_OCEAN 222631200   
!#define NX_OCEAN 222757920
!#define NX_OCEAN 2

MODULE mo_mapping 
! author: FastOpt GmbH
! purpose:
!------------------------------------------------------------------
! map vector of control parameters to physical model variables
!------------------------------------------------------------------
! contains:
!     subroutine cont2mod(n,xarg)
!     subroutine mod2cont(n,par)
!     subroutine addprior(n,x,m,y)
!     subroutine mod2y(m,y)
! public: everything
!
  use substitutes
  use pumamod
  use pumamod, only : gt ,st, stp, stm, nroot, mypid,  nlev,nlep,  ntru
  use pumamod, only : sz,szm,szp,sd,sdm,sdp,sq,sqm,sqp,sp,spp,spm,psurf
  use pumamod, only :gaussw,ct,t0,nugp,nlon,nlat,nesp,npackgp,spnorm,nhor,nlpp,nspp
  use pumamod, only :tdissd,tdissz,tdisst,tdissq,tfrc, dt,dls,ww,nud
  use pumamod, only : genobs  ! added SiS 03.2024
  
!  use miscmod
  use fluxmod, only : vdiff_lamm,vdiff_d,vdiff_c,vdiff_b
  use radmod, only : tswr1,tswr2,tswr3,rcl1,rcl2,acllwr,acl2,tpofmt,clgray,th2oc
  use rainmod, only : rcrit
  use surfmod , only :
  use landmod, only : albsmin,albsmax,albsminf,albsmaxf,albland,albsmax
  use landmod, only : albgmin,albgmax,dz0land,dalbclim,drhsland,rlue,tau_soil
  use landmod, only : nlandt,nlandw,dtcl,dwcl,dts
  use vegmod, only : q10,cveg_a,cveg_b,ct_crit
  use seamod, only: dz0sea, dz0ice, charnock
  use icemod, only : xmind
  use oceanmod, only : vdiffk
  IMPLICIT NONE
  
  integer, parameter :: nx=NX !2 !3 ! dimensions of full control vector
  integer, parameter :: nx_ocean=NX_OCEAN ! dimensions of ocean control vector
  integer, parameter :: ny=NY !dimension of result vector
  integer :: iounit = 15
#if defined VECSTSRF 
  real    :: spatternf(NESP) ! vector for search direction in parameter space
  real    :: spatternp(NSPP) ! vector for search direction in parameter space
  real    :: fpattern(NHOR) ! vector for search direction in parameter space
  real    :: gpattern(NUGP) ! vector for search direction in parameter space
#endif
  real    :: prior(nx),std_prior(nx)
  integer :: ix ! index of x, marks last mapped control variable
  real    :: xsave(NX)
CONTAINS

SUBROUTINE cont2mod(n,xarg)
! purpose: cont2mod maps natural params xarg onto model params
! input:
!    integer :: n     ! dimension of control vector
!    real :: xarg(NX) ! control vector
!  contains
!      subroutine setmod(pmod,kix)
!      subroutine setmodarray(kmod,pmod,kix)
!
  IMPLICIT NONE
  integer :: n
  real :: xarg(NX)
#ifdef PARAMS
  real :: zp
#elif defined (VECSTSRF)
  real :: xghp(NUGP),xsh(NESP),xshp(NSPP)
#endif
  real :: dtmp(NUGP)
  INTEGER k,lock_i,lock_j,lock_k

  
  xsave(:) = xarg(:) ! make a copy to module variable for later use in mo_mapping_ocean
  ix = 0 ! assumes that this is the first mapping routine

#ifdef PARAMS
  if (mypid==nroot) then
#if (NX-NX_OCEAN)==2
!!!     call setmod(zp,1); tdissz(:) = zp   ! diffusion time scale for vorticity 1./(TWOPI*[days])
!!!     call setmod(tpofmt,2)        ! tuning of point of mean (lwr) transmissivity in layer
     call setmod(vdiff_lamm,1)   ! const. for vertical diffusion and surface fluxes
#elif (NX-NX_OCEAN)==6
     call setmodarray(2,tfrc(1:2),1) ! time scale for rayleigh friction 1./(TWOPI*[days]) - interdependence with sztfac & sdtfac!
     call setmod(zp,3); tdissd(:) = zp   ! diffusion time scale for divergence 1./(TWOPI*[days])
     call setmod(zp,4); tdissz(:) = zp   ! diffusion time scale for vorticity 1./(TWOPI*[days])
     call setmod(zp,5); tdisst(:) = zp   ! diffusion time scale for temperature 1./(TWOPI*[days])
     call setmod(tpofmt,6)        ! tuning of point of mean (lwr) transmissivity in layer
!     call setmod(zp,7); dalbclim(:) = min(1.,max(0.,dalbclim(:) + zp)) ! clim. background albedo
#elif (NX-NX_OCEAN)==10

       print*, "=============================="
     print*, "============Entering cont2mod=================="

     call setmodarray(2,tfrc(1:2),1) ! time scale for rayleigh friction 1./(TWOPI*[days]) - interdependence with sztfac & sdtfac!
     call setmod(zp,3); tdissd(:) = zp   ! diffusion time scale for divergence 1./(TWOPI*[days])
     call setmod(zp,4); tdissz(:) = zp   ! diffusion time scale for vorticity 1./(TWOPI*[days])
     call setmod(zp,5); tdisst(:) = zp   ! diffusion time scale for temperature 1./(TWOPI*[days])
     call setmod(tpofmt,6)        ! tuning of point of mean (lwr) transmissivity in layer
!     call setmod(zp,7); dalbclim(:) = min(1.,max(0.,dalbclim(:) + zp)) ! clim. background albedo
     call setmod(vdiff_lamm,7)   ! const. for vertical diffusion and surface fluxes
     call setmod(vdiff_b,8)      ! const. for vertical diffusion and surface fluxes
     call setmod(vdiff_c,9)      ! const. for vertical diffusion and surface fluxes
     call setmod(vdiff_d,10)      ! const. for vertical diffusion and surface fluxes
     vdiff_d = maxx(0.,vdiff_d,0.05)
 
  
#elif (NX-NX_OCEAN)==28672

     print*, "=============================="
     print*, "============Entering cont2mod=================="

     call setmodarray_d01(NUGP,th2ocg(:),1,0.3)
     call setmodarray_d01(NUGP,acllwrg(:),NUGP+1,0.25)
     call setmodarray_d01(NUGP,tswr1g(:),2*NUGP+1,0.2)
     call setmodarray_d01(NUGP,tswr2g(:),3*NUGP+1,0.5)
     call setmodarray_d01(NUGP,tswr3g(:),4*NUGP+1,0.2)
     call setmodarray_d01(NUGP,tpofmtg(:),5*NUGP+1,1.1)
     call setmodarray_d01(NUGP,albgmaxg(:),6*NUGP+1,1.0)
     call setmodarray_d01(NUGP,rcritg(:,7),7*NUGP+1,1.0)
     call setmodarray_d01(NUGP,rcritg(:,8),8*NUGP+1,1.0)
     call setmodarray_d01(NUGP,rcritg(:,9),9*NUGP+1,1.0)
     call setmodarray_d01(NUGP,rcritg(:,10),10*NUGP+1,1.0)
     call setmodarray_d01(NUGP,alpha1g(:),11*NUGP+1,5.0)
     call setmodarray_d01(NUGP,beta1g(:),12*NUGP+1,5.0)
     call setmodarray_d01(NUGP,gamma1g(:),13*NUGP+1,5.0)
  
! smooth paramters to make them smooth
     call smooth_data( th2ocg(:), 1 )
     call smooth_data( acllwrg(:), 1 )
     call smooth_data( tswr1g(:), 1 )
     call smooth_data( tswr2g(:), 1 )
     call smooth_data( tswr3g(:), 1 )
     call smooth_data( tpofmtg(:), 1 )
     call smooth_data( albgmaxg(:), 1 )
     call smooth_data( th2ocg(:), 1 )
     call smooth_data( acllwrg(:), 1 )
     call smooth_data( tswr1g(:), 1 )
     call smooth_data( tswr2g(:), 1 )
     call smooth_data( tswr3g(:), 1 )
     call smooth_data( tpofmtg(:), 1 )
     call smooth_data( albgmaxg(:), 1 )
     dtmp(:)= rcritg(:,7)
     call smooth_data( dtmp(:), 1 )
     rcritg(:,7) =dtmp(:)
     dtmp(:)= rcritg(:,8)
     call smooth_data( dtmp(:), 1 )
     rcritg(:,8) =dtmp(:)
     dtmp(:)= rcritg(:,9)
     call smooth_data( dtmp(:), 1 )
     rcritg(:,9) =dtmp(:)
     dtmp(:)= rcritg(:,10)
     call smooth_data( dtmp(:), 1 )
     rcritg(:,10) =dtmp(:)
     call smooth_data( alpha1g(:), 1 )
     call smooth_data( beta1g(:), 1 )
     call smooth_data( gamma1g(:), 1 )
     call smooth_data( alpha1g(:), 1 )
     call smooth_data( beta1g(:), 1 )
     call smooth_data( gamma1g(:), 1 )

#else
# error 'Control vector inconsistent, check NX, NX_OCEAN.'
#endif
end if
!
  ! call mpbcrn(tfrc,NLEV)
  ! call mpbcrn(tdissd,NLEV)
  ! call mpbcrn(tdissz,NLEV)
  ! call mpbcrn(tdisst,NLEV)
  ! call mpbcr(tpofmt)
  ! call mpbcr(vdiff_lamm)
  ! call mpbcr(vdiff_b)
  ! call mpbcr(vdiff_c)
  ! call mpbcr(vdiff_d)
  ! call mpbci(ix)

! ++++++++++++++++++++ New for 36 parameters 
  call mpbci(ix)
  call mpscgp(th2ocg,th2ocm,1)
  call mpscgp(acllwrg,acllwrm,1)
  call mpscgp(tswr1g,tswr1m,1)
  call mpscgp(tswr2g,tswr2m,1)
  call mpscgp(tswr3g,tswr3m,1)
  call mpscgp(tpofmtg,tpofmtm,1)
  call mpscgp(albgmaxg,albgmaxm,1)
  call mpscgp(rcritg,rcritm,NLEV)
  call mpscgp(alpha1g,alpha1m,1)
  call mpscgp(beta1g,beta1m,1)
  call mpscgp(gamma1g,gamma1m,1)


#elif defined (VECSTSRF)
     select case (n-nx_ocean)
     case (1)
        if (mypid==nroot) print*,"perturbing lowest model level with: PATTERN * ",xsave(ix+1)," K"
        stm(:,NLEV) = stm(:,NLEV) + xsave(1)/ct ! temperature, old timestep
        stp(:,NLEV) = stp(:,NLEV) + xsave(1)/ct ! temperature, perturbation in Kelvin, x being average pert. per gridpt.
     case (NESP)
        if (mypid==nroot) print'(1x,a,g22.15,a)',"perturbing lowest model level with: ",sum(abs(xsave(ix+1:ix+NESP)))/NESP, &
             " K on average at all spectral coefficients."
        call mpscsp(xsave(ix+1:ix+NESP),xshp,1)
        stm(:,NLEV) = stm(:,NLEV) + xshp(:)/ct ! temperature, old timestep
        stp(:,NLEV) = stp(:,NLEV) + xshp(ix+1:ix+n)/ct ! temperature, current timestep
     case (NUGP)
        if (mypid==nroot) print'(1x,a,g22.15,a)',"perturbing lowest model level with: ",sum(abs(xsave(ix+1:ix+NUGP)))/NUGP, &
             " K on average at all gridpoints."
        if (mypid==nroot) print*,"x(1:6)=",xsave(1:6)
        call mpscgp(xsave(ix+1:ix+NUGP)/ct,xghp,1) ! perturbation field, dimensionless
        call gp2fc(xghp(:),NLON,NLPP) ! convert to Fourier coefficients
        call fc2sp(xghp(:),xsh(:))    ! convert to spectral coefficients
        call mpsumsc(xsh(:),xshp(:),1)
        stm(:,NLEV) = stm(:,NLEV) + xshp(:) ! temperature, old timestep
        stp(:,NLEV) = stp(:,NLEV) + xshp(:) ! temperature, current timestep
     case default
        if (mypid==nroot) then
           print*,'No control variable identified in cont2mod.'
           print*,'n = ',n,'nx_ocean =',nx_ocean
        end if
        call mpabort ('No control variable identified in cont2mod.')
     end select
     call mpgallsp( st(:,NLEV),stp(:,NLEV),1 )
     ix = ix + n-nx_ocean
#endif
  contains

    subroutine setmod(pmod,kix)
! purpose: perturb model parameter
! input:
!      integer :: kix ! index
! output:
!      real :: pmod   ! perturbed model parameter value
!
      real :: pmod
      integer :: kix
      if (ix+1 .ne. kix) then
         print*,ix+1,kix
         stop "wrong index in setmod!"
      end if
      pmod = prior(kix) + std_prior(kix) * xsave(kix)
      print '(a,i3,a,g30.20)','setting parameter ',ix+1,' to ',pmod
!	write(nud,'(a,i3,a,g30.20)') 'Setting parameter ',ix+1,' to ',pmod
      ix = kix
    end subroutine setmod

    subroutine setmodarray(kmod,pmod,kix)
! purpose: perturb model parameter array
! input:
!      integer :: kix     ! index
!      integer :: kmod    ! dimension of pmod
! output:
!      real :: pmod(pmod) ! perturbed model parameter value
!
      integer :: kmod,kix
      real :: pmod(kmod)
      if (ix+1 .ne. kix) then
         print*,ix+1,kix
         stop "wrong index in setmodarray!"
      end if
      pmod(:) = prior(kix:kix+kmod-1) + std_prior(kix:kix+kmod-1) * xsave(kix:kix+kmod-1)
      print '(a,i3,a,2g30.20)','setting parameter ',ix+1,' to (max/min)',maxval(pmod(:)),minval(pmod(:))
      ix = kix + kmod - 1
    end subroutine setmodarray

    subroutine setmodarray_d01(kmod,pmod,kix,parr)
  ! purpose: perturb model parameter array
  ! input:
  !      integer :: kix     ! index
  !      integer :: kmod    ! dimension of pmod
  ! output:
  !      real :: pmod(pmod) ! perturbed model parameter value
  !
    integer :: kmod,kix
    real :: pmod(kmod)
    real :: parr
    if (ix+1 .ne. kix) then
        print*,ix+1,kix
        stop "wrong index in setmodarray!"
    end if
    pmod(:) = parr/(1+exp(-(prior(kix:kix+kmod-1) + std_prior(kix:kix+kmod-1) * xsave(kix:kix+kmod-1)))) 
    print '(a,i6,a,2g30.20)','setting parameter ',ix+1,' to (max/min)',maxval(pmod(:)),minval(pmod(:))
    ix = kix + kmod - 1
  end subroutine setmodarray_d01
END SUBROUTINE cont2mod

      subroutine enlarger( adj_var, a, size_z )

!      use pumamod, only : nhor,nlat, nlon
    implicit none
    integer :: size_z
    double precision :: a(nlon+2,nlat+2,size_z)
    double precision :: adj_var(nugp,size_z)
    double precision :: aux(nlon,nlat,size_z)
    integer :: loc_i
    integer :: loc_j
    integer :: loc_k

    loc_k = 0
    do loc_j = 1, nlat
      do loc_i = 1, nlon
        loc_k = loc_k+1
        aux(loc_i,loc_j,:) = adj_var(loc_k,:)
      end do
    end do
    a(:,:,:) = 0.
    a(2:nlon+1,2:nlat+1,:) = aux(1:nlon,1:nlat,:)
    a(1,2:nlat+1,:) = a(nlon+1,2:nlat+1,:)
    a(nlon+2,2:nlat+1,:) = a(2,2:nlat+1,:)
    a(:,1,:) = a(:,2,:)
    a(:,nlat+2,:) = a(:,nlat+1,:)
    end subroutine enlarger


    subroutine smooth_data( adj_var, size_z )
    use pumamod, only : nhor,nlat,nlon

    implicit none

    integer :: size_z
    double precision :: a(nlon+2,nlat+2,size_z)
    double precision :: as(nlon,nlat,size_z)
    double precision :: adj_var(nugp,size_z)


    double precision :: aux(nlon,nlat,size_z)
    integer :: loc_i
    integer :: loc_j
    integer :: loc_k
    call enlarger( adj_var, a, size_z )
    print '(a,2g30.20)','smooth_data: ', maxval(a),minval(a)
    do loc_k= 1, size_z
    do loc_j = 2, nlat+1
      do loc_i = 2, nlon+1
        as(loc_i-1,loc_j-1,loc_k) = sum(a(loc_i-1:loc_i+1,loc_j-1:loc_j+1,loc_k))/9.
      end do
    end do
    end do
    loc_k = 0
    do loc_j = 1, nlat
      do loc_i = 1, nlon
        loc_k = loc_k+1
        adj_var(loc_k,:) = as(loc_i,loc_j,:)
      end do
    end do

    end subroutine smooth_data





SUBROUTINE mod2cont(n,par)
! purpose: mod2cont  maps physical model variables to their natural control parameters 'par'
! input:
!    integer :: n     ! dimension of control vector
! output:
!    real :: xarg(NX) ! control vector
!  contains
!        subroutine setcontrol(pmod)
!        subroutine setcontrolarray(kmod,pmod)
!
  IMPLICIT NONE
  integer :: n
  REAL, DIMENSION(n), INTENT(out) :: par ! control variables

  ix = 0

#ifdef PARAMS
#if (NX-NX_OCEAN)==2
  call setcontrol(vdiff_lamm) ! const. for vertical diffusion and surface fluxes
!!!  call setcontrol(tdissz(1))  ! diffusion time scale for vorticity [days]
!!!  call setcontrol(tpofmt)     ! tuning of point of mean (lwr) transmissivity in layer
#elif (NX-NX_OCEAN)==6
  call setcontrolarray(2,tfrc(1:2))  ! time scale for rayleigh friction (1.0 / (TWOPI * tfrc))
  call setcontrol(tdissd(1))  ! diffusion time scale for divergence [days]
  call setcontrol(tdissz(1))  ! diffusion time scale for vorticity [days]
  call setcontrol(tdisst(1))  ! diffusion time scale for temperature [days]
  call setcontrol(tpofmt)     ! tuning of point of mean (lwr) transmissivity in layer
!  call setcontrol(0.)   ! clim. background albedo - spatially uniform
!  std_prior(i-1) = sum(dalbclim(:))/NUGP
#elif (NX-NX_OCEAN)==10
 print*, "=============================="
 print*, "============Entering mod2cont=================="

  call setcontrolarray(2,tfrc(1:2))  ! time scale for rayleigh friction (1.0 / (TWOPI * tfrc))
  call setcontrol(tdissd(1))  ! diffusion time scale for divergence [days]
  call setcontrol(tdissz(1))  ! diffusion time scale for vorticity [days]
  call setcontrol(tdisst(1))  ! diffusion time scale for temperature [days]
  call setcontrol(tpofmt)     ! tuning of point of mean (lwr) transmissivity in layer
!  call setcontrol(0.)   ! clim. background albedo - spatially uniform
!  std_prior(i-1) = sum(dalbclim(:))/NUGP
  call setcontrol(vdiff_lamm) ! const. for vertical diffusion and surface fluxes
  call setcontrol(vdiff_b)    ! const. for vertical diffusion and surface fluxes
  call setcontrol(vdiff_c)    ! const. for vertical diffusion and surface fluxes
  call setcontrol(vdiff_d)    ! const. for vertical diffusion and surface fluxes

#elif (NX-NX_OCEAN)==28672

  print*, "=============================="
  print*, "============Entering mod2cont=================="

  call setcontrolarray_d01(NUGP,th2ocg,0.3)
  call setcontrolarray_d01(NUGP,acllwrg,0.25)
  call setcontrolarray_d01(NUGP,tswr1g,0.2)
  call setcontrolarray_d01(NUGP,tswr2g,0.5)
  call setcontrolarray_d01(NUGP,tswr3g,0.2)
  call setcontrolarray_d01(NUGP,tpofmtg,1.1)
  call setcontrolarray_d01(NUGP,albgmaxg,1.0)
  call setcontrolarray_d01(NUGP,rcritg(:,7),1.0)
  call setcontrolarray_d01(NUGP,rcritg(:,8),1.0)
  call setcontrolarray_d01(NUGP,rcritg(:,9),1.0)
  call setcontrolarray_d01(NUGP,rcritg(:,10),1.0)
  call setcontrolarray_d01(NUGP,alpha1g,5.0)
  call setcontrolarray_d01(NUGP,beta1g,5.0)
  call setcontrolarray_d01(NUGP,gamma1g,5.0)

#endif 

! change switch genobs to namelist parameter SiS
if (genobs == 0) then
  par(:) = 0 ! de-tune parameters for id-twin experiment
!  only for forward run
else
  par(:) = 0.
endif
#else
if (genobs == 0) then
#if defined (VECSTSRF)
  par(:) = 0. !- .1
!  par(:) = spattern(:) ! perturbation field, dimensionless
!  par(:) = gpattern(:) ! perturbation field, dimensionless
#endif
else
#if defined (VECSTSRF) 
     par(:) = 0 ! unperturbed perturbation, st
#endif
endif
#endif
contains

  subroutine setcontrol(pmod)
! purpose: set first guess of control parameter and its error
! input:
!    real :: pmod ! model quantity
!
    real :: pmod
    prior(ix+1) = pmod
    std_prior(ix+1) = prior(ix+1) ! makes x to be fractions of pmod
    ix = ix + 1
  end subroutine setcontrol

  subroutine setcontrolarray(kmod,pmod)
! purpose: set first guess of control parameter vector section and its error
! input:
!    integer :: kmod    ! dimension of pmod
!    real :: pmod(kmod) ! model quantity
!
    integer :: kmod
    real :: pmod(kmod)
    prior(ix+1:ix+kmod) = pmod(:)
    std_prior(ix+1:ix+kmod) = prior(ix+1:ix+kmod) ! makes x to be fractions of pmod
    ix = ix + kmod
  end subroutine setcontrolarray

  subroutine setcontrolarray_d01(kmod,pmod,parr)
! purpose: set first guess of control parameter vector section and its error
! input:
!    integer :: kmod    ! dimension of pmod
!    real :: pmod(kmod) ! model quantity
!
    integer :: kmod
    real :: pmod(kmod)
    real :: parr  ! real value (0 1] to tune upper range
    prior(ix+1:ix+kmod) = -log(parr/pmod(:)-1.)
    !prior(ix+1:ix+kmod) = pmod(:)
    std_prior(ix+1:ix+kmod) = prior(ix+1:ix+kmod) ! makes x to be fractions of pmod
    ix = ix + kmod
  end subroutine setcontrolarray_d01


END SUBROUTINE mod2cont



subroutine addprior(n,x,m,y)
! purpose: add background term to y
! input:
!    integer :: n ! dimension of x
!    integer :: m ! dimension of y
!    real :: x(n) ! control vector
!    real :: y(m) ! cost function
! output:
!    real :: y(m) ! cost function
!
  integer :: n,m,i
  real :: x(n)
  real :: y(m)
  do i=1,n
      y(1) = y(1) + 0.5*(x(i)/100.)**2 ! x(i) must be normalised with a priori standard deviation
      print '(a,g20.12,a,i6)','adding prior',0.5*(x(i)/100.)**2,' for component ',i
  end do
end subroutine addprior

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine Gaussian_filter(flt,radius)
   implicit none
   real :: flt(1:5,1:5)
   real:: radius,sumflt
   integer :: i,j
    do i=1,5
     do j=1,5
       flt(i,j)=(1./2./3.14/radius**2)*exp(-(i-3.)**2/2./radius**2-(j-3.)**2/2./radius**2);
     enddo
    enddo
  
    sumflt=sum(flt)
     do i=1,5
     do j=1,5
       flt(i,j)=flt(i,j)/sumflt;
     enddo
    enddo
  
end subroutine Gaussian_filter

      subroutine enlarger_data( adj_var, par_mtxe )
      use pumamod, only : nugp
      implicit none
      double precision :: par_mtxe(68,36)
      double precision :: adj_var(nhor)


      double precision :: aux(64,32)
      integer :: loc_i
      integer :: loc_j
      integer :: loc_k

      loc_k = 0
      do loc_j = 1, 32
        do loc_i = 1, 64
          loc_k = loc_k+1
          aux(loc_i,loc_j) = adj_var(loc_k)
        end do
      end do
      
      par_mtxe(:,:) = 0.
      par_mtxe(3:66,3:34) = aux(1:64,1:32)
      par_mtxe(1:2,3:34) = aux(63:64,1:32)
      par_mtxe(67:68,3:34) = aux(1:2,1:32)
      par_mtxe(:,1:2) = 0.
      par_mtxe(:,35:36) = 0.

      end subroutine enlarger_data





SUBROUTINE mod2y(m,y)
! purpose: map model result to y
! input:
!    integer :: m ! dimension of y
! output:
!    real :: y(m) ! cost function
!
  use observation, only : cost,yname
  IMPLICIT NONE
  real :: y(m)
  integer :: m

if (genobs == 1) then   !SiS
  print*,'Wrote pseudo-observations to ',trim(yname)
endif

  y(1) = cost

END SUBROUTINE mod2y

END MODULE mo_mapping

SUBROUTINE init_global_mapping
! author: FastOpt GmbH
! purpose: initialise mapping of model and control variables
!
  USE mo_mapping
  use legmod
  use fftmod
  implicit none
  integer :: get_free_unit
#if defined VECSTSRF
  integer :: ihead(8)
  real    :: fnrm ! vector for search direction in parameter space
  integer :: jlat
#endif
! defaults
#if defined VECSTSRF
  if (mypid==nroot) then
     iounit = get_free_unit(iounit)
! reading horizontal field in gripoint space, SERVICE-format
! header: integer :: ihead(8) ! code level (yy)yymmdd hhmm(ss) nx ny free free
! field:  real(kind=?) :: f(nx*ny)
     open(iounit,file='mapping_pattern.srv',form='unformatted')
     read(iounit) ihead(:)
     read(iounit) gpattern(:)
     close(iounit)
! normalise ?
     fnrm = 0.
     do jlat=1,NLAT
        fnrm = fnrm + gaussw(jlat) * sum((gpattern((jlat-1)*NLON+1:jlat*NLON))**2)
     end do
     gpattern(:) = gpattern(:) / sqrt(fnrm /  NLON / sum(gaussw(1:NLAT)))
     print '(1x,a,3g16.10)',"gpattern (max,min,ave): ",maxval(gpattern(:)),minval(gpattern(:)),sum(gpattern(:))/NUGP
  end if
! spectral and legendre transform
  call mpscgp(gpattern(:),fpattern(:),1) ! perturbation field, dimensionless
  call gp2fc(fpattern(:),NLON,NLPP)
  call fc2sp(fpattern(:),spatternf(:))
  call mpsumsc(spatternf(:),spatternp(:),1)
#endif
END SUBROUTINE init_global_mapping

integer function get_free_unit(kunit)
! author: FastOpt
! purpose: return an unused unit number
! input:
!     integer :: kunit ! unit number before trying higher ones
!
  implicit none
  integer :: kunit
  integer :: iunit
  logical :: l_opened
  l_opened=.true.
  iunit=max(6,kunit-1)
  do while (l_opened) ! find a free unit number
     iunit=iunit+1
     inquire(UNIT=iunit,OPENED=l_opened)
  end do
  get_free_unit=iunit
end function get_free_unit
