! This program uses the IPCCN scheme to obtain solutions to a two-dimensional
! space-fractional Brusselator equation with Neumann boundary comditions.
! The discretization in space is done
! via the matrix Transfer Technique and the linearly implicit predictor-corrector
! schemes are used in time. 

program Brusselator3D
  implicit none
  character (len = 20):: fileName
  integer, parameter:: rkind = selected_real_kind(15, 307), numDimension = 3
  real (kind = rkind):: alpha, spaceStepSize, x0, x_end
  real (kind = rkind):: timeStepSize, finalTime, A, B
  real (kind = rkind):: kappa1, kappa2
  real (kind = rkind), dimension(:,:), allocatable:: matrixAu, matrixAv
  real (kind = rkind), dimension(:,:), allocatable:: reshapeUvec, reshapeVvec
  real (kind = rkind), dimension(:,:,:,:), allocatable:: Uvec, Vvec
  real (kind = rkind), dimension(:,:,:), allocatable:: Xmesh, Ymesh, Zmesh
  real (kind = rkind), dimension(:), allocatable:: Xvec, Tvec
  integer:: i, j, k, numGridPoints, numTimeSteps
  integer:: LDA, LDB, NRHS, info, matrixSize
  real (kind = rkind), parameter:: Pi = 3.14159265


  x0 = 0                                       ! The problem is solved in the domain [0, 10]^3                           
  x_end = 10.0

  ! File for writing solution values                                          
  fileName = 'Brusselator3D.txt'

  write(*,'(A)') "Please enter the Time of Integration: "
  read(*,*) finalTime

  write(*,'(A)') "Please enter the parameter A: "
  read(*,*) A

  write(*,'(A)') "Please enter the parameter B: "
  read(*,*) B

  write(*,'(A)') "Please enter the parameter kappa1: "
  read(*,*) kappa1

  write(*,'(A)') "Please enter the parameter kappa2: "
  read(*,*) kappa2
  
  
  write(*,'(A)') "Please enter the space step size: "
  read(*,*) spaceStepSize
  
  write(*,'(A)') "Please enter the time step size: "
  read(*,*) timeStepSize
  
  write(*,'(A)') "Please enter the fractional order of the derivative (alpha): "
  read(*,*) alpha


  !Number of grid points
  numGridPoints = (x_end - x0)/spaceStepSize + 1                      ! Number of grid points in one spatial direction
  numTimeSteps = finalTime/timeStepSize                               ! Number of time steps
  LDA = max(1, numGridPoints)
  LDB = max(1, numGridPoints)
  NRHS = 1
  matrixSize = numGridPoints**numDimension
  ! Allocation of arrays for use in the sequel
  allocate(Xvec(numGridPoints))
  allocate(Tvec(numTimeSteps + 2))
  allocate(Uvec(numGridPoints, numGridPoints, numGridPoints, NRHS))
  allocate(Vvec(numGridPoints, numGridPoints, numGridPoints, NRHS))  
  allocate(Xmesh(numGridPoints, numGridPoints, numGridPoints))
  allocate(Ymesh(numGridPoints, numGridPoints, numGridPoints))
  allocate(Zmesh(numGridPoints, numGridPoints, numGridPoints))
  allocate(matrixAu(matrixSize, matrixSize))
  allocate(matrixAv(matrixSize, matrixSize))
  allocate(reshapeUvec(matrixSize, NRHS))
  allocate(reshapeVvec(matrixSize, NRHS))
  
  forall (j=0:numGridPoints-1) Xvec(j+1) = x0 + j*spaceStepSize       ! Stencil in one spatial direction 
  forall (j=0:numTimeSteps+1) Tvec(j+1) = j*timeStepSize              ! Stencil for the time 

  ! This  function builds the matrix A^(alpha/2) of size numGridPoints
  ! depending on the choice of alpha 
  call MTT3D(matrixAu, numGridPoints, kappa1, spaceStepSize, alpha)
  call MTT3D(matrixAv, numGridPoints, kappa2, spaceStepSize, alpha)
  write(*,*)matrixAu(1:20,1:20)
  ! Constructs the three dimensional meshgrid points in each direction X, Y and Z
  call meshgrid3(Xmesh, Ymesh, Zmesh, Xvec, Xvec, Xvec, numGridPoints, numGridPoints, numGridPoints)

  ! Inital  conditions
  Uvec(:, :, :, NRHS) =  A + cos(2*Pi*Xmesh)*cos(2*Pi*Ymesh)*cos(2*Pi*Zmesh)
  Vvec(:, :, :, NRHS) = B/A + cos(2*Pi*Xmesh)*cos(2*Pi*Ymesh)*cos(2*Pi*Zmesh)

  ! Reshape the three dimensional array to a vector for easy implementation
  reshapeUvec = reshape(Uvec, (/ matrixSize, NRHS /))
  reshapeVvec = reshape(Vvec, (/ matrixSize, NRHS /))
  
  ! Solve the semi-discretized problem using the IPC scheme 
  call IPC(reshapeUvec, reshapeVvec, matrixAu, matrixAv, numGridPoints, timeStepSize, numTimeSteps, NRHS, A, B, numDimension)

  ! Reshape the vector back to the three dimensional array 
  Uvec = reshape(reshapeUvec, (/ numGridPoints, numGridPoints, numGridPoints, NRHS /))
  Vvec = reshape(reshapeVvec, (/ numGridPoints, numGridPoints, numGridPoints, NRHS /))

  ! Write solution and the meshgrid to a file for graphics
  open(unit = 40, file = fileName, status = 'replace', action = 'write')
  do k =1, numGridPoints
     do j = 1, numGridPoints
        do i = 1, numGridPoints
           write(40,*) Xmesh(i, j, k), ' ',  Ymesh(i, j, k),' ', Zmesh(i, j, k), ' ', Uvec(i, j, k, NRHS), ' ',  Vvec(i, j, k, NRHS)
        enddo
     end do
  enddo
  
  close(40)
  deallocate(matrixAu, matrixAv,  Uvec, Tvec, Xvec)

end program Brusselator3D

