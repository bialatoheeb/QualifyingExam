! This subroutine computes the inverse of a matrix A using the LU factorization of the matrix

subroutine  inverse(Ainverse, A, num)
  implicit none
  integer, parameter:: rkind = selected_real_kind(15, 307)
  integer, intent(in):: num
  real (kind = rkind), intent(in):: A(num, num)
  real(kind = rkind), intent(out):: Ainverse(num, num)
  real (kind = rkind), dimension(:), allocatable:: work
  integer, dimension(:), allocatable:: IPIV
  integer :: Lwork, info, LDA

  Lwork = num*num
  LDA = num

  allocate(work(Lwork))
  allocate(IPIV(num))

  ! DGETRF computes the LU factorization of M X N matix using partial pivoting
  ! This subroutine calculates only for N X N matrices
  ! Copy A and use the value of the copied matrix, otherwise LAPACK overwrites the matrix with the inverse

  Ainverse = A  
  call DGETRF(num, num, Ainverse, LDA, IPIV, info)

 ! Check if the LU decomposition is carried out successfully
  if (info .eq. 0) then
     write(*,'(A)') '"Matrix LU decomposition: Successful"'
     write(*, *)
  elseif (info .gt. 0) then
     write(*, '(A,I6, I6)') '"Matrix LU decomposition: Failed with U(', info, ',', info, ') = 0"'
     write(*, *)
  else
     write(*,'(A)') '"Matrix LU decomposition: Failed (Check the values of the matrix)"'
     write(*, *)
  end if


 ! DGETRI computes the inverse of matrix A using the LU factorization obtained for DGETRF
  call DGETRI(num, Ainverse, num, IPIV, work, Lwork, info)

  !Check if inverse was obtained successfully
  if (info .eq. 0) then
     write(*, '(A)') '"Matrix Inversion: Successful"'
     write(*, *)
  else
     write(*,'(A)') '"Matrix Inversion: Failed"'
     write(*, *)
  end if

  deallocate(work)
  deallocate(IPIV)
end subroutine inverse

