! This program uses the IPCCN scheme to obtain solutions to a two-dimensional
! space-fractional Brusselator equation with Neumann boundary comditions.
! The discretization in space is done
! via the matrix Transfer Technique and the linearly implicit predictor-corrector
! schemes are used in time. 

program Brusselator2D
  implicit none
  character (len = 20):: fileName
  integer, parameter:: rkind = selected_real_kind(15, 307), numDimension = 1
  real (kind = rkind):: alpha, spaceStepSize, x0, x_end
  real (kind = rkind):: timeStepSize, finalTime, A, B
  real (kind = rkind):: kappa1, kappa2, diffusionCoefficient
  real (kind = rkind), dimension(:,:),   allocatable:: matrixAu, matrixAv
  real (kind = rkind), dimension(:,:), allocatable:: Uvec, Vvec
  real (kind = rkind), dimension(:),     allocatable:: Xvec, Tvec
  integer:: i, j, numGridPoints,  numTimeSteps, info
  integer:: LDA, LDB, NRHS, matrixSize
  real (kind = rkind), parameter:: Pi = 3.14159265


!  diffusionCoefficient = 1.0

  x0 = 0                                       ! The problem is solvedin the domain [0, 1] X [0,10]                          
  x_end = 10.0

  ! File for writing solution values                                          
  fileName = 'Brusselator1D.txt'

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
  allocate(Uvec(numGridPoints, NRHS))
  allocate(Vvec(numGridPoints, NRHS))
  allocate(matrixAu(matrixSize, matrixSize))
  allocate(matrixAv(matrixSize, matrixSize))

  forall (j=0:numGridPoints-1) Xvec(j+1) = x0 + j*spaceStepSize       ! Stencil in one spatial direction 
  forall (j=0:numTimeSteps+1) Tvec(j+1) = j*timeStepSize              ! Stencil for the time 
  

  ! This  function builds the matrix A^(alpha/2) of size matrixSize
  ! depending on the choice of alpha and kappa 
  call MTT1D(matrixAu, numGridPoints, kappa1, spaceStepSize, alpha)

  call MTT1D(matrixAv, numGridPoints, kappa2, spaceStepSize, alpha)

  ! Inital  conditions
  Uvec(:,NRHS) =  A + cos(2*pi*Xvec)
  Vvec(:,NRHS) =  B/A + cos(2*pi*Xvec)

 !Solve the semi-discretized problem using the IPC scheme 
  call IPC(Uvec, Vvec, matrixAu, matrixAv, numGridPoints, timeStepSize, numTimeSteps, A, B, numDimension)

  ! Write solution and the meshgrid to a file for graphics
  open(unit = 40, file = fileName, status = 'replace', action = 'write')
  
  do i = 1, numGridPoints
     write(40,*) Xvec(i), ' ',  Uvec(i,NRHS), ' ',  Vvec(i, NRHS)
  enddo
     

  close(40)
  deallocate(matrixAu, matrixAv, Uvec, Vvec, Tvec, Xvec)

end program Brusselator2D

