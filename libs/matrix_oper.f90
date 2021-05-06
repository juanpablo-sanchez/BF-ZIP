module matrix_oper

use kinds

implicit none

interface mymatmul
 module procedure mymatmulmmr8,&
                  mymatmulmvr8,&
                  mymatmulvmr8
end interface

integer, save :: num_threads_pregs = 0   ! 0 means all
character :: version*5= "1.000"


contains

function mymatmulmmr8(a,b) result(c)
! matrix x matrix
real(r8) :: a(:,:),b(:,:),c(size(a,1),size(b,2))
integer :: n,m,k
m=size(a,1)
n=size(b,2)
k=size(a,2)
call mymatmulsubr8(m,n,k,a,b,c)
end function

function mymatmulmvr8(a,b) result(c)
! matrix x vector
real(r8) :: a(:,:),b(:),c(size(a,1)),xc(size(a,1),1),xb(size(b),1)
integer :: n,m,k
m=size(a,1)
n=1
k=size(a,2)
call mymatmulsubr8(m,n,k,a,b,xc)
c=xc(:,1)
end function

function mymatmulvmr8(a,b) result(c)
! vector x matrix
real(r8) :: a(:),b(:,:),c(size(b,2)),xa(1,size(a)),xc(1,size(b,2))
integer :: n,m,k
m=1
n=size(b,2)
k=size(a)
call mymatmulsubr8(m,n,k,a,b,c)
end function


end module

