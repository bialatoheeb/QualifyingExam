  ! This subroutine solves a nonlinear system of dfferential equations with
  ! the linearly implicit predictor-corrector Crank-Nicolson scheme  

subroutine  IPC(reshapeUvec, reshapeVvec, matrixAu, matrixAv, num, dt, numTimeSteps, A, B, numDimension)
    implicit none
    character (len = 1):: TRANS = 'N'
    integer:: i, j, num, numTimeSteps, matrixSize
    integer:: LDA, LDB, info, int, numDimension
    integer, parameter:: rkind = selected_real_kind(15, 307), NRHS = 1
    real (kind = rkind), intent(in):: matrixAu(num**numDimension, num**numDimension), matrixAv(num**numDimension, num**numDimension)
    real (kind = rkind), intent(inout):: reshapeUvec(num**numDimension, NRHS), reshapeVvec(num**numDimension, NRHS)
    complex*16, allocatable, dimension(:,:):: AuMatrix, AvMatrix, YUvec, YVvec, ASolve
    real (kind = rkind), allocatable, dimension(:,:):: FUvec, F2Uvec, AUvec
    real (kind = rkind), allocatable, dimension(:,:):: FVvec, F2Vvec, AVvec, IdentityMatrix
    integer, allocatable, dimension(:):: IPIV
    real (kind = rkind):: dt, A, B
    complex*16:: c, w1, w2, w3

    c = (-1, 1)
    w1= (0, -1)
    w2 = (0.5, 0)
    w3 = (0.5, -0.5)

    matrixSize = num**numDimension
    LDA = max(1, matrixSize)
    LDB = max(1, matrixSize)
    allocate(AuMatrix(matrixSize, matrixSize))
    allocate(AvMatrix(matrixSize, matrixSize))
    allocate(ASolve(matrixSize, matrixSize))
    allocate(IdentityMatrix(matrixSize, matrixSize))
    allocate(FUvec(matrixSize, NRHS))
    allocate(FVvec(matrixSize, NRHS))
    allocate(F2Uvec(matrixSize,NRHS))
    allocate(F2Vvec(matrixSize,NRHS))
    allocate(AUvec(matrixSize, NRHS))
    allocate(AVvec(matrixSize, NRHS))
    allocate(YUvec(matrixSize, NRHS))
    allocate(YVvec(matrixSize, NRHS))
    allocate(IPIV(matrixSize))

    IdentityMatrix = 0
    forall(j=1:matrixSize) IdentityMatrix(j,j) = 1

    ! ZGETRF computes the LU factorization of a general
    ! complex matrix with double precison    
    AuMatrix = cmplx(dt*matrixAu - c*IdentityMatrix)
    AvMatrix = cmplx(dt*matrixAv - c*IdentityMatrix)

    ! ZGETRF obtains the LU factorization of a general complex matrix
    ! with double precision
    call ZGETRF(matrixSize, matrixSize, AuMatrix, LDA, IPIV, info)
    call ZGETRF(matrixSize, matrixSize, AvMatrix, LDA, IPIV, info)

    ! Check if the LU factorization was successful    
    if (info .eq. 0) then
     !  write(*, '(A)') "LU factorization: Successful"
    elseif (info .gt. 0) then
       write(*,'(A,I3,A,I3,A)') '"LU factorization: Failed with U(', info, ',', info, ') = 0"'
       stop
    else
       write(*,'(A)') "LU factorization:: Failed with illegal ", &
            info , "-th value"
       stop
    end if
    

    ! Implementation of the IPC-CN Scheme
    j = 0
    int = 20

    do i = 1, numTimeSteps
       ! Solve the system by forward and backward substitution
       FUvec = A + reshapeUvec*reshapeUvec*reshapeVvec - (B + 1)*reshapeUvec
       FVvec = B*reshapeUvec  -  reshapeUvec*reshapeUvec*reshapeVvec

       YUvec = w1*reshapeUvec + w3*dt*FUvec
       YVvec = w1*reshapeVvec + w3*dt*FVvec
       
       ! ZGETRS solves the system AX = b using the LU factorization
       ! by the ZGETRF (double precision)
       ASolve = AuMatrix
       call ZGETRS(TRANS, matrixSize, NRHS, ASolve, LDA, IPIV, YUvec, LDB, info)       

       ASolve = AvMatrix
       call ZGETRS(TRANS, matrixSize, NRHS, ASolve, LDA, IPIV, YVvec, LDB, info)

       ! Check on the solution for the Predictor
       if (info .eq. 0) then
       !   write(*, '(A)') "Matrix Solver for Predictor: Succesful"
       else
          write(*,'(A)') "Matrix Solver for Predictor: Failed with illegal ", &
               info , "-th value"
          stop
       end if

       ! Predicted Value 
       AUvec = 2*real(YUvec)
       AVvec = 2*real(YVvec)
       
       
       ! Vectors to be used for computing Corrector
       F2Uvec = A + AUvec*AUvec*AVvec - (B + 1)*AUvec
       F2Vvec = B*AUvec - AUvec*AUvec*AVvec
       
       YUvec = w2*dt*(F2Uvec - FUvec)
       YVvec = w2*dt*(F2Vvec - FVvec)
       
       ! ZGETRS solves the system AX = b using the LU factorization
       ! by the ZGETRF
       ASolve = AuMatrix
       call ZGETRS(TRANS, matrixSize, NRHS, ASolve, LDA, IPIV, YUvec, LDB, info)

       ASolve = AvMatrix
       call ZGETRS(TRANS, matrixSize, NRHS, ASolve, LDA, IPIV, YVvec, LDB, info)

       ! Another check on the solution for the Corrector
       if (info .eq. 0) then
       !   write(*, '(A)') "Matrix Solver for Corrector: Successful"
       else
          write(*,'(A)') "Matrix Solver for Corrector: Failed with illegal ", info , "-th value"
          stop
       end if       

       ! Corrected and Final Solution at time i*dt
       reshapeUvec = AUvec + 2*real(YUvec)
       reshapeVvec = AVvec + 2*real(YVvec)
       
       j = j + 1       
       if (mod(j, int) .eq. 0) then
          write(*,'(A,I5,A,I5,A)') 'Done with  ', j, ' out of ', numTimeSteps, ' iterations'
       end if
       
    end do    
    
    deallocate(FUvec, F2Uvec, AUvec, YUvec, IPIV)
    deallocate(FVvec, F2Vvec, AVvec, YVvec)
    deallocate(IdentityMatrix, ASolve, AuMatrix, AvMatrix)
  end subroutine IPC
  

