! This subroutine uses the matrix Transfer technique to build
! the matrix A^(alpha/2) arising from the semidiscretization of the space-fractional equation.
! with Neumann boundary conditions. The matrix is constructed via the diagonalization 
! of the standard Laplacian. In particular A^(\alpha/2) = H^{-1}Lambda H.

subroutine  MTT1D(matrixA, numGridPoints, diffusionCoefficient, spaceStepSize, alpha)
  implicit none 
  integer:: i, j, numGridPoints
  integer, parameter:: rkind = selected_real_kind(15, 307)
  real, parameter:: Pi = 3.1415927
  real (kind = rkind), intent(out):: matrixA(numGridPoints, numGridPoints)
  real (kind = rkind), dimension(:,:), allocatable:: EigenValues, EigenVectors, EigInv
  real (kind = rkind), intent(in):: alpha, spaceStepSize, diffusionCoefficient
  real (kind = rkind):: Coeff

  ! Allocate arrays for later use
  allocate(EigenVectors(numGridPoints, numGridPoints))
  allocate(EigenValues(numGridPoints, numGridPoints))
  allocate(EigInv(numGridPoints, numGridPoints))
    
  Coeff =  diffusionCoefficient/(spaceStepSize**alpha)

  ! Diagonalization of the matrix
  EigenValues = 0
  do j = 0, numGridPoints-1
     do i = 0, numGridPoints-1
       EigenVectors(i+1, j+1) = (-1)**(j+1)*cos(i*j*Pi/(numGridPoints-1))
    end do
    EigenValues(j+1, j+1) = Coeff*(2* cos(j*Pi/(2*(numGridPoints-1))))**alpha
 end do
 EigenValues(numGridPoints, numGridPoints) = 0
 
 EigenVectors = transpose(EigenVectors)
 
 ! Inverse of matrix is obtained and the matrix is created on the second line below
 call inverse(EigInv, EigenVectors, numGridPoints)

! Compute A = H^{-1}Lambda H.
 matrixA = matmul(EigenVectors, matmul(EigenValues, EigInv))
 
 deallocate(EigenVectors)
 deallocate(Eigenvalues)
 deallocate(EigInv)
end subroutine MTT1D



  
