!---------------------------------------------------------------------------------------
!The interface function
Function varphi(xx,yy,zz)
use comdata
use pbedata
use molecule

Implicit Double Precision(A-H, O-Z)
r2=xx**2+yy**2+zz**2
varphi=r2-rds**2
Return
End
!------------------------------------------------------------------------------
Function varphix(x,y,z)
Implicit Double Precision(A-H, O-Z)
r2=x**2+y**2+z**2
varphix=2.d0*x
Return
End
!------------------------------------------------------------------------------
Function varphiy(x,y,z)
Implicit Double Precision(A-H, O-Z)
r2=x**2+y**2+z**2
varphiy=2.d0*y
Return
End
!------------------------------------------------------------------------------
Function varphiz(x,y,z)
Implicit Double Precision(A-H, O-Z)
r2=x**2+y**2+z**2
varphiz=2.d0*z
Return
End
