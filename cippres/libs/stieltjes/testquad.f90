

integer, parameter :: QR_K = selected_real_kind (32)
real (kind=QR_K) :: MyReal

 MyReal = 1q0
 write(*,*)QR_K,MyReal

end
