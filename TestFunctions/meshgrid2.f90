! This subroutine generates 2D grid coordinates returnX and returnY based
! on the coordinates of th vectors X and Y. returnX and returnY are matrices 
! sizes length(Y) by length(X).

subroutine  meshgrid2(returnX, returnY, X, Y, sizeX, sizeY)
  implicit none
  integer:: i, sizeX, sizeY
  integer, parameter:: rkind = selected_real_kind(15, 307)
  real (kind = rkind), dimension(sizeX), intent(in):: X
  real (kind = rkind), dimension(sizeY), intent(in):: Y
  real (kind = rkind), dimension(sizeY, sizeX), intent(out):: returnX, returnY  
  
  do i=1, sizeY
     returnX(i,:) = X
  end do
  do i=1, sizeX
     returnY(:,i) = Y
  end do
  
end subroutine meshgrid2


  
