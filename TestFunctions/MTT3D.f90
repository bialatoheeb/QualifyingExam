! This subroutine uses the matrix Transfer technique to build
! the matrix A^(alpha/2) arising from the semidiscretization of the space-fractional equation.
! The matrix is constructed from the eigenvaliues and eigenvectors of the matrix A
! This is with Neumann boundary conditions
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

  !Eigenvalues and Eigenvectors of the matrix
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
 
  EigenVectors = reshape(EigenVects, (/ matrixSize, matrixSize /))

  do j=0, numGridPoints - 1
     do i= 0, numGridPoints - 1
        do k=0, numGridPoints - 1
           EigenVals(k+1, i+1, j+1) = (4*(cos(i*Pi/(2*(numGridPoints-1)))**2.0 +  cos(j*Pi/(2*(numGridPoints-1)))**2.0&
                                      + cos(k*Pi/(2*(numGridPoints-1)))**2.0))**(alpha/2.0)
        enddo        
     end do
  enddo

  EigenVals = Coeff*EigenVals
  EigVals = reshape(EigenVals, (/ matrixSize, NRHS/))

  
  EigenValues = 0
  forall (i=1:matrixSize) EigenValues(i,i) = EigVals(i,NRHS)
  
 ! Inverse of matrix is obtained and the matrix is created on the second line below
  call inverse(EigInv, EigenVectors, matrixSize)
  
  matrixA = matmul(EigInv, matmul(EigenValues, EigenVectors))
  matrixA = transpose(matrixA)

 deallocate(EigenVectors, EigenVals, EigVals)
 deallocate(Eigenvalues, EigenVects)
 deallocate(EigInv)
 
end subroutine MTT3D



  
