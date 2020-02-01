! Contains Electromagnetic fields related code

MODULE fields

   USE mUMFPACK

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE ASSEMBLE_POISSON
      ! UMFPACK via subroutine calls
      integer,parameter :: n=5,nnz=12
      integer :: Ap(0:n)=[0,2,5,9,10,12]
      integer :: Ai(nnz)=[0, 1, 0,  2, 4, 1,  2, 3, 4, 2, 1, 4]
      real(8) :: Ax(nnz)=[2.,3.,3.,-1.,4.,4.,-3.,1.,2.,2.,6.,1.]
      real(8) :: b(n)=[8.,45.,-3.,3.,19.],x(n)
      integer :: status,i
      
      ! Notice these don't need b
      call s_umfpack_symbolic(n,n,Ap,Ai,Ax,status=status)
      call s_umfpack_numeric(Ap,Ai,Ax)
      call s_umfpack_free_symbolic
      ! The solution requires b of course
      call s_umfpack_solve(UMFPACK_A,Ap,Ai,Ax,x,b)
      call s_umfpack_free_numeric
      
      do i=lbound(x,1),ubound(x,1)
        WRITE(*,'(a,i0,a,f0.10)') "x(",i,") = ",x(i)
      enddo
   END SUBROUTINE ASSEMBLE_POISSON

END MODULE fields