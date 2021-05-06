
subroutine mymatmulsubr8(m,n,k,a,b,c)
use kinds
integer :: m,n,k
real(r8) :: a(m,k), b(k,n), c(m,n)
real(r8) :: myone=1.d0,myzero=0.d0

#if ( _MYMKL ==1)
call dgemm('n','n',m,n,k,myone,a,m,b,k,myzero,c,m)
#else
c=matmul(a,b)
#endif
end subroutine


!#################
!#################
subroutine BlasInverseSym_logdet(x,n,log_det) 
! inteface to call Inverse subroutines from BLAS-MKL libraries
use kinds
!$ use omp_lib
implicit none
integer :: n,i,j
real(r8) :: x(n,n),t11,t22,log_det
real ::t1,t2
logical, save :: first=.true.
integer, save :: evals=1,nthr=1,nthrg
integer :: mkl_get_max_threads
!
#if(_MYMKL==1)
integer :: info
#ifdef _OPENMP
  if (evals == 1) then
     !$omp parallel
     !$omp master
     nthrg=omp_get_num_threads()
     !$omp end master
     !$omp end parallel
     nthr=mkl_get_max_threads()
  endif 

#endif
#endif


#if(_MYMKL==1)
#ifdef _OPENMP
  if (first) t11=omp_get_wtime()
 !              t11=omp_get_wtime()
#else
   nthr=1
   if (first) print*,'Inverse LAPACK ATLAS'
#endif
   call dpotrf('L',n,x,n,info)
   if (info/=0) then   
       print*,'Matrix not positive-definte 1 - info',info
       stop
   endif
 
   log_det=0
   do i=1,n
    log_det=log_det+log(x(i,i)*x(i,i))
   enddo

   call dpotri( 'L', n, x, n, info )
   if (info==0) then
       ! with huge matrices this create a insuficiente virtual memory
       !forall ( i=1:n, j=1:n, j>i ) x(i,j)=x(j,i)
!omp parallel
!omp do private (i,j) shared(n,x)
       do i=1,n
          do j=i,n
             x(i,j) =  x(j,i)
          enddo
       enddo
!omp end do
!omp end parallel

   else
       print*,'Matrix not positive-definte 2 - info',info
       stop
   endif
#ifdef _OPENMP
   t22=omp_get_wtime()
   if (first) print '(a,2i5,a,f10.4)','    Inverse LAPACK MKL dpotrf/i  #threads=',nthr,nthrg,' Elapsed omp_get_time: ',t22-t11
! print '(a,2i5,a,f10.4)','    Inverse LAPACK MKL dpotrf/i  #threads=',nthr,nthrg,' Elapsed omp_get_time: ',t22-t11
#endif
#else
   print*,'MKL inverse specified but Not linked with library compile with _MYMKL = .TRUE. '
   stop
#endif
!call cpu_time(t2)
!print*,'Elapsed time',(t2-t1)/nthr
if (evals>1) first=.false.
evals=evals+1
end subroutine
