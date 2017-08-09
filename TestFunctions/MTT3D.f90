! This subroutine uses the matrix Transfer technique to build
! the matrix A^(alpha/2) arising from the semidiscretization of the three-dimensional space-fractional 
! equation with Neumann boundary conditions. The matrix is constructed via the diagonalization 
! of the standard Laplacian. In particular A^(\alpha/2) = H^{-1}Lambda H.

subroutine  MTT3D(matrixA, numGridPoints, diffusionCoefficient, spaceStepSize, alpha)
  implicit none
  integer:: i, j, k, l, m, n, numGridPoints, matrixSize, NRHS
  integer, parameter:: rkind = selected_real_kind(15, 307)
  real, parameter:: Pi = 3.1415927
  real (kind = rkind), intent(out):: matrixA(numGridPoints**3, numGridPoints**3)
  real (kind = rkind), dimension(:,:), allocatable:: EigenValues, EigenVectors, EigInv, EigVals
  real (kind = rkind), dimension(:,:,:,:,:,:), allocatable:: EigenVects
  real (kind = rkind), dimension(:,:,:), allocatable:: EigenVals
  real (kind = rkind), intent(in):: alpha, spaceStepSize, diffusionCoefficient
  real (kind = rkind):: Coeff

  NRHS = 1
  matrixSize = numGridPoints**3
  
  ! Allocate arrays for later use
  allocate(EigenVectors(matrixSize, matrixSize))
  allocate(EigenValues(matrixSize, matrixSize))
  allocate(EigenVals(numGridPoints, numGridPoints, numGridPoints))
  allocate(EigInv(matrixSize, matrixSize))
  allocate(EigenVects(numGridPoints, numGridPoints, numGridPoints, numGridPoints, numGridPoints, numGridPoints))
  allocate(EigVals(matrixSize, NRHS))
  
  Coeff =  diffusionCoefficient/(spaceStepSize**alpha)

  ! Constructing the matrix  P form which H = P^T is obtained 
  do n = 0, numGridPoints - 1
     do m = 0, numGridPoints - 1
        do i = 0, numGridPoints - 1
           do j = 0, numGridPoints - 1
              do k = 0, numGridPoints - 1
                 do l=0, numGridPoints - 1
                    EigenVects(l+1,k+1,j+1,i+1,m+1,n+1) = (-1)**(i+j+m+n+1)*cos(n*j*Pi/(numGridPoints-1))&
                                                           *cos(m*k*Pi/(numGridPoints-1))*cos(i*l*Pi/(numGridPoints-1))
                 enddo
              enddo              
           end do
        enddo
     enddo
  enddo
 
  ! Reshape the matrix to matrixSize
  EigenVectors = reshape(EigenVects, (/ matrixSize, matrixSize /))

 ! Eigenvalues of the fractional Laplcian in three dimensions
  do j=0, numGridPoints - 1
     do i= 0, numGridPoints - 1
        do k=0, numGridPoints - 1
           EigenVals(k+1, i+1, j+1) = (4*(cos(i*Pi/(2*(numGridPoints-1)))**2.0 +  cos(j*Pi/(2*(numGridPoints-1)))**2.0&
                                      + cos(k*Pi/(2*(numGridPoints-1)))**2.0))**(alpha/2.0)
        enddo        
     end do
  enddo

  ! Introduce the diffusion coefficient and the space step-size in the m atrix
  EigenVals = Coeff*EigenVals
  EigVals = reshape(EigenVals, (/ matrixSize, NRHS/))

  ! Diagonalize the eigenvalues
  EigenValues = 0
  forall (i=1:matrixSize) EigenValues(i,i) = EigVals(i,NRHS)
  
  ! Inverse of matrix is obtained 
  call inverse(EigInv, EigenVectors, matrixSize)
  
  ! Compute A = H^{-1}Lambda H.
  matrixA = matmul(EigInv, matmul(EigenValues, EigenVectors))
  matrixA = transpose(matrixA)

 deallocate(EigenVectors, EigenVals, EigVals)
 deallocate(Eigenvalues, EigenVects)
 deallocate(EigInv)
 
end subroutine MTT3D



  
