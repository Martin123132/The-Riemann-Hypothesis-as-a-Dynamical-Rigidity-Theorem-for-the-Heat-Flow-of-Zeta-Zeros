program rh_binary128_rounding_mode_probe
  use, intrinsic :: ieee_arithmetic
  implicit none

  integer, parameter :: dp = selected_real_kind(33)
  integer, parameter :: hex_kind = selected_int_kind(32)
  real(dp) :: one
  type(ieee_round_type) :: mode

  one = 1.0_dp
  call ieee_get_rounding_mode(mode)

  write(*,'(A,I0)') 'selected_real_kind=', dp
  write(*,'(A,I0)') 'radix=', radix(one)
  write(*,'(A,I0)') 'digits=', digits(one)
  write(*,'(A,I0)') 'min_exponent=', minexponent(one)
  write(*,'(A,I0)') 'max_exponent=', maxexponent(one)
  write(*,'(A,L1)') 'ieee_datatype=', ieee_support_datatype(one)
  write(*,'(A,L1)') 'ieee_nearest_supported=', ieee_support_rounding(ieee_nearest, one)
  if (mode == ieee_nearest) then
     write(*,'(A)') 'rounding_mode=nearest'
  else if (mode == ieee_to_zero) then
     write(*,'(A)') 'rounding_mode=to_zero'
  else if (mode == ieee_up) then
     write(*,'(A)') 'rounding_mode=up'
  else if (mode == ieee_down) then
     write(*,'(A)') 'rounding_mode=down'
  else
     write(*,'(A)') 'rounding_mode=other'
  endif
  write(*,'(A,Z32.32)') 'one_hex=', transfer(one, 0_hex_kind)
  write(*,'(A,Z32.32)') 'epsilon_hex=', transfer(epsilon(one), 0_hex_kind)
  write(*,'(A,Z32.32)') 'tiny_hex=', transfer(tiny(one), 0_hex_kind)
end program rh_binary128_rounding_mode_probe
