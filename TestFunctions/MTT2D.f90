! This subroutine uses the matrix Transfer technique to build
! the matrix A^(alpha/2) arising from the semidiscretization of the two-dimensional space-fractional 
! equation with Neumann boundary conditions. The matrix is constructed via the diagonalization 
! of the standard Laplacian. In particular A^(\alpha/2) = H^{-1}Lambda H.

subroutine  MTT2D(matrixA, numGridPoints, diffusionCoefficient, spaceStepSize, alpha)
  implicit none
  integer:: i, j, k, l, numGridPoints, matrixSize, NRHS
  integer, parameter:: rkind = selected_real_kind(15, 307)
  real, parameter:: Pi = 3.1415927
  real (kind = rkind), intent(out):: matrixA(numGridPoints**2, numGridPoints**2)
  real (kind = rkind), dimension(:,:), allocatable:: EigenValues, EigenVectors, EigInv, EigenVals, EigVals
  real (kind = rkind), dimension(:,:,:,:), allocatable:: EigenVects
  real (kind = rkind), intent(in):: alpha, spaceStepSize, diffusionCoefficient
  real (kind = rkind):: Coeff

  NRHS = 1
  matrixSize = numGridPoints**2
  ! Allocate arrays for later use
  allocate(EigenVectors(matrixSize, matrixSize))
  allocate(EigenValues(matrixSize, matrixSize))
  allocate(EigenVals(numGridPoints, numGridPoints))
  allocate(EigInv(matrixSize, matrixSize))
  allocate(EigenVects(numGridPoints, numGridPoints, numGridPoints, numGridPoints))
  allocate(EigVals(matrixSize, NRHS))
  
  Coeff =  diffusionCoefficient/(spaceStepSize**alpha)

  !Diagonalization of the matrix
  do i = 0, numGridPoints-1
     do j = 0, numGridPoints-1
        do k = 0, numGridPoints-1
           do l = 0, numGridPoints-1
              EigenVects(l+1,k+1,j+1,i+1) = (-1)**(l+k+i+j+1)*cos(i*k*Pi/(numGridPoints-1))*cos(j*l*Pi/(numGridPoints-1))
           end do
        enddo
     enddo
  enddo
 
  ! Reshape the matrix to matrixSize
  EigenVectors = reshape(EigenVects, (/ matrixSize, matrixSize /))
  
  ! Eigenvalues of the fractional Laplcian in two dimensions
  do j=0, numGridPoints - 1
     do i= 0, numGridPoints-1
        EigenVals(i+1, j+1) = (4*(cos(i*Pi/(2*(numGridPoints-1)))**2.0 +  cos(j*Pi/(2*(numGridPoints-1)))**2.0))**(alpha/2.0)
     end do
  enddo

  ! Introduce the diffusion coefficient and the space step-size in the m atrix
  EigenVals = Coeff*EigenVals
  EigVals = reshape(EigenVals, (/ matrixSize, NRHS/))

  ! Diagonalize the eigenvalues
  EigenValues = 0
  forall (i=1:matrixSize) EigenValues(i,i) = EigVals(i,NRHS)
 
  ! Inverse of matrix is obtained and the matrix is created on the second line below
  call inverse(EigInv, EigenVectors, matrixSize)

  ! Compute A = H^{-1}Lambda H.
  matrixA = matmul(EigInv, matmul(EigenValues, EigenVectors))
  matrixA = transpose(matrixA)

  deallocate(EigenVectors, EigenVals, EigVals)
 deallocate(Eigenvalues, EigenVects)
 deallocate(EigInv)
end subroutine MTT2D



  
