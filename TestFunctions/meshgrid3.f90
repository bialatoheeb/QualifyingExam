! This subroutine generates 3D grid coordinates returnX, returnY and returnZ
! based on the coordinates of the vectors X, Y and Z. returnX, returnY and
! returnZ are arrays of size length(Y) by length(X) by length(Z)

subroutine  meshgrid3(returnX, returnY, returnZ, X, Y, Z, sizeX, sizeY, sizeZ)
  implicit none
  integer:: i, sizeX, sizeY, sizeZ
  integer, parameter:: rkind = selected_real_kind(15, 307)
  real (kind = rkind), dimension(sizeX), intent(in):: X
  real (kind = rkind), dimension(sizeY), intent(in):: Y
  real (kind = rkind), dimension(sizeZ), intent(in):: Z
  real (kind = rkind), dimension(sizeY, sizeX, sizeZ),intent(out):: returnX, returnY, returnZ
  real (kind = rkind), dimension(:,:), allocatable:: returnXTemp, returnYTemp, returnZTemp
 
  
  
  allocate(returnXTemp(sizeX, sizeX))
  allocate(returnYTemp(sizeX, sizeX))
  allocate(returnZTemp(sizeX, sizeX))  
  
  call meshgrid2(returnXtemp, returnYTemp, X, Y, sizeX, sizeY)
  
  do i=1, sizeZ
     returnX(:,:,i) = returnXTemp
     returnY(:,:,i) = returnYTemp
     returnZ(:,:,i) = Z(i) 
  end do
  
  deallocate(returnXTemp, returnYTemp, returnZTemp)
end subroutine meshgrid3

  
