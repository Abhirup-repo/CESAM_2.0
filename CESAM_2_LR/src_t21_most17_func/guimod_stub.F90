! Dummy module replacement for guimod.f90
! Use this for environments with no X11

subroutine guistart
return
end subroutine guistart

subroutine guistep_puma
return
end subroutine guistep_puma

subroutine guistep_plasim
return
end subroutine guistep_plasim

subroutine guistop
return
end subroutine guistop

subroutine guihor(yn,f,klev,pm,pa)
use pumamod
character (len=*) :: yn
integer :: klev
real    :: pm,pa
real    :: f(NLON,NLPP,klev)
return
end subroutine guihor

subroutine guigv(yn,f)
use pumamod
character (len=*) :: yn
real :: f(NLON,NLPP,NLEV)
return
end subroutine guigv

subroutine guigvcol(yname,f,klon)
use pumamod
character (len=*) :: yname
integer :: klon
real    :: f(NLON,NLPP,NLEV)
return
end subroutine guigvcol


subroutine guigtcol(f,klon)
use pumamod
integer :: klon
real :: f(NLON,NLPP,NLEV)
return
end subroutine guigtcol


subroutine guid3dcol(yname,f,klon,klev,pm,pa)
use pumamod
character (len=*) :: yname
integer :: klon,klev
real    :: pm,pa
real    :: f(NLON,NLPP,klev)
return
end subroutine guid3dcol

subroutine guips(fp)
use pumamod
real :: fp(NLON,NLAT)
return
end subroutine guips

subroutine guiput(yn,f,k1,k2,k3)
character (len=*) :: yn
integer :: k1,k2,k3
real(kind=4)  :: f(k1)
return
end subroutine guiput

subroutine guihorlsg(yn,f,mask,klev,pm,pa)
character (len=*) :: yn 
return
end subroutine guihorlsg

