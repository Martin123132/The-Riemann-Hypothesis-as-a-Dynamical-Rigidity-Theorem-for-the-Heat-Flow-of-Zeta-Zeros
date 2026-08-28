module rh_hardy_chain_telemetry

  implicit none

  private

  integer, parameter :: dp = selected_real_kind(33)
  integer, parameter :: dp1 = selected_int_kind(16)
  integer, parameter :: hex_kind = selected_int_kind(32)

  logical, save :: telemetry_enabled = .false.
  logical, save :: chain_active = .false.
  integer, save :: telemetry_unit = -1
  integer, save :: telemetry_rank = -1
  integer, save :: telemetry_max_chains = 0
  integer, save :: telemetry_chain = 0
  integer, save :: telemetry_block = -1
  integer, save :: telemetry_branch = 0
  integer, save :: telemetry_initial_degree = 0
  integer(dp1), save :: telemetry_sum = -1_dp1
  integer(dp1), save :: telemetry_max_direct_terms = 1000000_dp1
  real(dp), save :: telemetry_initial_coefficients(3) = 0.0_dp

  public :: rh_tel_init, rh_tel_close, rh_tel_set_position
  public :: rh_tel_set_branch, rh_tel_set_initial_coefficients
  public :: rh_tel_begin_chain, rh_tel_q
  public :: rh_tel_step, rh_tel_end_chain

contains

  subroutine rh_tel_init(rank)
    integer, intent(in) :: rank

    character(len=1024) :: path
    character(len=64) :: value
    character(len=:), allocatable :: line
    integer :: env_length, env_status, ios

    telemetry_rank = rank
    telemetry_enabled = .false.
    chain_active = .false.
    telemetry_chain = 0
    telemetry_initial_degree = 0
    telemetry_initial_coefficients = 0.0_dp

    if (rank.ne.0) return

    call get_environment_variable('RH_HARDY_CHAIN_TELEMETRY', path, &
         length=env_length, status=env_status, trim_name=.true.)
    if (env_status.ne.0.or.env_length.le.0) return

    value = ''
    call get_environment_variable('RH_HARDY_TELEMETRY_MAX_CHAINS', value, &
         length=env_length, status=env_status, trim_name=.true.)
    if (env_status.ne.0.or.env_length.le.0) return
    read(value(1:env_length),*,iostat=ios) telemetry_max_chains
    if (ios.ne.0.or.telemetry_max_chains.le.0) return

    value = ''
    call get_environment_variable('RH_HARDY_TELEMETRY_MAX_DIRECT_TERMS', value, &
         length=env_length, status=env_status, trim_name=.true.)
    if (env_status.eq.0.and.env_length.gt.0) then
       read(value(1:env_length),*,iostat=ios) telemetry_max_direct_terms
       if (ios.ne.0.or.telemetry_max_direct_terms.le.0_dp1) then
          telemetry_max_direct_terms = 1000000_dp1
       endif
    endif

    open(newunit=telemetry_unit, file=trim(path), status='unknown', &
         position='append', action='write', iostat=ios)
    if (ios.ne.0) then
       telemetry_unit = -1
       return
    endif

    telemetry_enabled = .true.
    line = '{"type":"header","schema":"rh_hardy_chain_telemetry_v1",' // &
         '"rank":' // trim(i_text(rank)) // &
         ',"max_chains":' // trim(i_text(telemetry_max_chains)) // &
         ',"max_direct_terms":"' // trim(i8_text(telemetry_max_direct_terms)) // '"}'
    call write_record(line)
  end subroutine rh_tel_init


  subroutine rh_tel_close()
    character(len=:), allocatable :: line

    if (.not.telemetry_enabled) return
    line = '{"type":"summary","chains_recorded":' // &
         trim(i_text(telemetry_chain)) // '}'
    call write_record(line)
    close(telemetry_unit)
    telemetry_unit = -1
    telemetry_enabled = .false.
    chain_active = .false.
  end subroutine rh_tel_close


  subroutine rh_tel_set_position(block, sum_index)
    integer, intent(in) :: block
    integer(dp1), intent(in) :: sum_index

    if (.not.telemetry_enabled) return
    telemetry_block = block
    telemetry_sum = sum_index
  end subroutine rh_tel_set_position


  subroutine rh_tel_set_branch(branch)
    integer, intent(in) :: branch

    if (.not.telemetry_enabled) return
    telemetry_branch = branch
  end subroutine rh_tel_set_branch


  subroutine rh_tel_set_initial_coefficients(degree, coefficients)
    integer, intent(in) :: degree
    real(dp), intent(in) :: coefficients(3)

    if (.not.telemetry_enabled) return
    telemetry_initial_degree = degree
    telemetry_initial_coefficients = coefficients
  end subroutine rh_tel_set_initial_coefficients


  subroutine rh_tel_begin_chain(m, ip, mit, mmax, lengths, degrees, &
       phicoeff, frac_length, xr, conjugate_flag, tpm)
    integer, intent(in) :: m, ip, mit, mmax
    integer, intent(in) :: degrees(50), conjugate_flag(0:50)
    integer(dp1), intent(in) :: lengths(0:50)
    real(dp), intent(in) :: phicoeff(0:10,50), frac_length(0:50), xr(50), tpm

    character(len=:), allocatable :: line
    integer :: level
    logical :: subtract_one

    chain_active = .false.
    if (.not.telemetry_enabled) return
    if (telemetry_chain.ge.telemetry_max_chains) return
    if (telemetry_initial_degree.ne.m) return

    telemetry_chain = telemetry_chain + 1
    chain_active = .true.

    line = '{"type":"chain","chain":' // trim(i_text(telemetry_chain)) // &
         ',"block":' // trim(i_text(telemetry_block)) // &
         ',"sum_index":"' // trim(i8_text(telemetry_sum)) // '"' // &
         ',"branch":' // trim(i_text(telemetry_branch)) // &
         ',"base_degree":' // trim(i_text(m)) // &
         ',"initial_coefficients":' // initial_coeff_text(m) // &
         ',"initial_coefficients_hex":' // initial_coeff_hex_text(m) // &
         ',"ip":' // trim(i_text(ip)) // &
         ',"mit":' // trim(i_text(mit)) // &
         ',"mmax":' // trim(i_text(mmax)) // &
         ',"tpm":"' // trim(r_text(tpm)) // '"' // &
         ',"tpm_hex":"' // real_hex(tpm) // '"' // &
         ',"initial_length":"' // trim(i8_text(lengths(0))) // '"' // &
         ',"kernel_length":"' // trim(i8_text(lengths(mit-1))) // '"}'
    call write_record(line)

    do level=1,mit
       subtract_one = level.gt.1.and.phicoeff(1,level-1).gt.0.0_dp
       line = '{"type":"level","chain":' // trim(i_text(telemetry_chain)) // &
            ',"level":' // trim(i_text(level)) // &
            ',"length":"' // trim(i8_text(lengths(level-1))) // '"' // &
            ',"degree":' // trim(i_text(degrees(level))) // &
            ',"phi0":"' // trim(r_text(phicoeff(0,level))) // '"' // &
            ',"phi0_hex":"' // real_hex(phicoeff(0,level)) // '"' // &
            ',"coefficients":' // coeff_text(phicoeff, level, degrees(level)) // &
            ',"coefficients_hex":' // coeff_hex_text(phicoeff, level, degrees(level)) // &
            ',"xr":"' // trim(r_text(xr(level))) // '"' // &
            ',"xr_hex":"' // real_hex(xr(level)) // '"' // &
            ',"frac_length":"' // trim(r_text(frac_length(level))) // '"' // &
            ',"frac_length_hex":"' // real_hex(frac_length(level)) // '"' // &
            ',"conjugate":' // trim(logical_text(conjugate_flag(level).eq.-1)) // &
            ',"subtract_one":' // trim(logical_text(subtract_one)) // '}'
       call write_record(line)
    enddo
  end subroutine rh_tel_begin_chain


  subroutine rh_tel_q(nit, t1, t2, t3, t4, t5, endpoint, qq)
    integer, intent(in) :: nit
    complex(kind=16), intent(in) :: t1, t2, t3, t4, t5, endpoint, qq

    character(len=:), allocatable :: line

    if (.not.chain_active) return
    line = '{"type":"q_terms","chain":' // trim(i_text(telemetry_chain)) // &
         ',"nit":' // trim(i_text(nit)) // &
         ',"t1":' // complex_text(t1) // &
         ',"t1_hex":' // complex_hex_text(t1) // &
         ',"t2":' // complex_text(t2) // &
         ',"t2_hex":' // complex_hex_text(t2) // &
         ',"t3":' // complex_text(t3) // &
         ',"t3_hex":' // complex_hex_text(t3) // &
         ',"t4":' // complex_text(t4) // &
         ',"t4_hex":' // complex_hex_text(t4) // &
         ',"t5":' // complex_text(t5) // &
         ',"t5_hex":' // complex_hex_text(t5) // &
         ',"endpoint":' // complex_text(endpoint) // &
         ',"endpoint_hex":' // complex_hex_text(endpoint) // &
         ',"qq":' // complex_text(qq) // &
         ',"qq_hex":' // complex_hex_text(qq) // '}'
    call write_record(line)
  end subroutine rh_tel_q


  subroutine rh_tel_step(nit, state_before, multiplier, qq, state_pre_transform, &
       state_after, lengths, degrees, phicoeff, conjugate_flag, tpm)
    integer, intent(in) :: nit
    integer, intent(in) :: degrees(50), conjugate_flag(0:50)
    integer(dp1), intent(in) :: lengths(0:50)
    real(dp), intent(in) :: phicoeff(0:10,50), tpm
    complex(kind=16), intent(in) :: state_before, multiplier, qq
    complex(kind=16), intent(in) :: state_pre_transform, state_after

    character(len=:), allocatable :: line
    complex(kind=16) :: raw_parent, raw_child, adapted_parent, adapted_child
    complex(kind=16) :: model_pre, model_after, local_defect
    complex(kind=16) :: child_state_defect, parent_state_defect
    logical :: parent_ok, child_ok, subtract_one

    if (.not.chain_active) return

    call direct_level_sum(nit, lengths, degrees, phicoeff, tpm, raw_parent, parent_ok)
    call direct_level_sum(nit+1, lengths, degrees, phicoeff, tpm, raw_child, child_ok)
    subtract_one = nit.gt.1.and.phicoeff(1,nit-1).gt.0.0_dp

    if (parent_ok.and.child_ok) then
       adapted_parent = adapt_level_sum(nit, raw_parent, conjugate_flag, phicoeff)
       adapted_child = adapt_level_sum(nit+1, raw_child, conjugate_flag, phicoeff)
       model_pre = multiplier*adapted_child + qq
       model_after = model_pre
       if (conjugate_flag(nit).eq.-1) model_after = conjg(model_after)
       if (subtract_one) model_after = model_after - 1.0_dp
       local_defect = adapted_parent - model_after
       child_state_defect = state_before - adapted_child
       parent_state_defect = state_after - adapted_parent

       line = '{"type":"recurrence","chain":' // trim(i_text(telemetry_chain)) // &
            ',"nit":' // trim(i_text(nit)) // &
            ',"direct_ok":true' // &
            ',"conjugate":' // trim(logical_text(conjugate_flag(nit).eq.-1)) // &
            ',"subtract_one":' // trim(logical_text(subtract_one)) // &
            ',"raw_parent":' // complex_text(raw_parent) // &
            ',"raw_parent_hex":' // complex_hex_text(raw_parent) // &
            ',"raw_child":' // complex_text(raw_child) // &
            ',"raw_child_hex":' // complex_hex_text(raw_child) // &
            ',"adapted_parent":' // complex_text(adapted_parent) // &
            ',"adapted_parent_hex":' // complex_hex_text(adapted_parent) // &
            ',"adapted_child":' // complex_text(adapted_child) // &
            ',"adapted_child_hex":' // complex_hex_text(adapted_child) // &
            ',"state_before":' // complex_text(state_before) // &
            ',"state_before_hex":' // complex_hex_text(state_before) // &
            ',"multiplier":' // complex_text(multiplier) // &
            ',"multiplier_hex":' // complex_hex_text(multiplier) // &
            ',"qq":' // complex_text(qq) // &
            ',"qq_hex":' // complex_hex_text(qq) // &
            ',"state_pre_transform":' // complex_text(state_pre_transform) // &
            ',"state_pre_transform_hex":' // complex_hex_text(state_pre_transform) // &
            ',"state_after":' // complex_text(state_after) // &
            ',"state_after_hex":' // complex_hex_text(state_after) // &
            ',"model_pre_transform":' // complex_text(model_pre) // &
            ',"model_pre_transform_hex":' // complex_hex_text(model_pre) // &
            ',"model_after":' // complex_text(model_after) // &
            ',"model_after_hex":' // complex_hex_text(model_after) // &
            ',"local_defect":' // complex_text(local_defect) // &
            ',"local_defect_hex":' // complex_hex_text(local_defect) // &
            ',"child_state_defect":' // complex_text(child_state_defect) // &
            ',"child_state_defect_hex":' // complex_hex_text(child_state_defect) // &
            ',"parent_state_defect":' // complex_text(parent_state_defect) // &
            ',"parent_state_defect_hex":' // complex_hex_text(parent_state_defect) // '}'
    else
       line = '{"type":"recurrence","chain":' // trim(i_text(telemetry_chain)) // &
            ',"nit":' // trim(i_text(nit)) // &
            ',"direct_ok":false' // &
            ',"conjugate":' // trim(logical_text(conjugate_flag(nit).eq.-1)) // &
            ',"subtract_one":' // trim(logical_text(subtract_one)) // &
            ',"state_before":' // complex_text(state_before) // &
            ',"state_before_hex":' // complex_hex_text(state_before) // &
            ',"multiplier":' // complex_text(multiplier) // &
            ',"multiplier_hex":' // complex_hex_text(multiplier) // &
            ',"qq":' // complex_text(qq) // &
            ',"qq_hex":' // complex_hex_text(qq) // &
            ',"state_pre_transform":' // complex_text(state_pre_transform) // &
            ',"state_pre_transform_hex":' // complex_hex_text(state_pre_transform) // &
            ',"state_after":' // complex_text(state_after) // &
            ',"state_after_hex":' // complex_hex_text(state_after) // '}'
    endif
    call write_record(line)
  end subroutine rh_tel_step


  subroutine rh_tel_end_chain(final_state)
    complex(kind=16), intent(in) :: final_state

    character(len=:), allocatable :: line

    if (.not.chain_active) return
    line = '{"type":"chain_end","chain":' // trim(i_text(telemetry_chain)) // &
         ',"final_state":' // complex_text(final_state) // '}'
    call write_record(line)
    chain_active = .false.
  end subroutine rh_tel_end_chain


  subroutine direct_level_sum(level, lengths, degrees, phicoeff, tpm, value, ok)
    integer, intent(in) :: level
    integer, intent(in) :: degrees(50)
    integer(dp1), intent(in) :: lengths(0:50)
    real(dp), intent(in) :: phicoeff(0:10,50), tpm
    complex(kind=16), intent(out) :: value
    logical, intent(out) :: ok

    integer :: degree, j
    integer(dp1) :: index
    real(dp) :: phase, y

    value = (0.0_dp,0.0_dp)
    ok = .false.
    if (level.lt.1.or.level.gt.50) return
    if (lengths(level-1).lt.0_dp1) return
    if (lengths(level-1).gt.telemetry_max_direct_terms) return

    degree = degrees(level)
    if (degree.lt.1.or.degree.gt.10) return

    do index=0_dp1,lengths(level-1)
       y = real(index,dp)
       phase = phicoeff(degree,level)*y
       do j=degree-1,1,-1
          phase = (phase+phicoeff(j,level))*y
       enddo
       value = value + exp(cmplx(0.0_dp,tpm*phase,kind=16))
    enddo
    ok = .true.
  end subroutine direct_level_sum


  function adapt_level_sum(level, raw_value, conjugate_flag, phicoeff) result(value)
    integer, intent(in) :: level
    integer, intent(in) :: conjugate_flag(0:50)
    real(dp), intent(in) :: phicoeff(0:10,50)
    complex(kind=16), intent(in) :: raw_value
    complex(kind=16) :: value

    value = raw_value
    if (conjugate_flag(level).eq.-1) value = conjg(value)
    if (level.gt.1) then
       if (phicoeff(1,level-1).gt.0.0_dp) value = value - 1.0_dp
    endif
  end function adapt_level_sum


  subroutine write_record(line)
    character(len=*), intent(in) :: line

    if (.not.telemetry_enabled) return
    if (telemetry_unit.eq.-1) return
    write(telemetry_unit,'(A)') trim(line)
    flush(telemetry_unit)
  end subroutine write_record


  function i_text(value) result(text)
    integer, intent(in) :: value
    character(len=32) :: text

    write(text,'(I0)') value
    text = adjustl(text)
  end function i_text


  function i8_text(value) result(text)
    integer(dp1), intent(in) :: value
    character(len=32) :: text

    write(text,'(I0)') value
    text = adjustl(text)
  end function i8_text


  function r_text(value) result(text)
    real(dp), intent(in) :: value
    character(len=64) :: text

    write(text,'(ES50.40E4)') value
    text = adjustl(text)
  end function r_text


  function logical_text(value) result(text)
    logical, intent(in) :: value
    character(len=5) :: text

    if (value) then
       text = 'true '
    else
       text = 'false'
    endif
  end function logical_text


  function complex_text(value) result(text)
    complex(kind=16), intent(in) :: value
    character(len=144) :: text

    text = '["' // trim(r_text(real(value,dp))) // '","' // &
         trim(r_text(aimag(value))) // '"]'
  end function complex_text


  function real_hex(value) result(text)
    real(dp), intent(in) :: value
    integer(hex_kind) :: bits
    character(len=32) :: text

    bits = transfer(value,bits)
    write(text,'(Z32.32)') bits
  end function real_hex


  function complex_hex_text(value) result(text)
    complex(kind=16), intent(in) :: value
    character(len=71) :: text

    text = '["' // real_hex(real(value,dp)) // '","' // &
         real_hex(aimag(value)) // '"]'
  end function complex_hex_text


  function coeff_text(phicoeff, level, degree) result(text)
    real(dp), intent(in) :: phicoeff(0:10,50)
    integer, intent(in) :: level, degree
    character(len=:), allocatable :: text
    integer :: j

    text = '['
    do j=1,degree
       if (j.gt.1) text = text // ','
       text = text // '"' // trim(r_text(phicoeff(j,level))) // '"'
    enddo
    text = text // ']'
  end function coeff_text


  function coeff_hex_text(phicoeff, level, degree) result(text)
    real(dp), intent(in) :: phicoeff(0:10,50)
    integer, intent(in) :: level, degree
    character(len=:), allocatable :: text
    integer :: j

    text = '['
    do j=1,degree
       if (j.gt.1) text = text // ','
       text = text // '"' // real_hex(phicoeff(j,level)) // '"'
    enddo
    text = text // ']'
  end function coeff_hex_text


  function initial_coeff_text(degree) result(text)
    integer, intent(in) :: degree
    character(len=:), allocatable :: text
    integer :: j

    text = '['
    do j=1,degree
       if (j.gt.1) text = text // ','
       text = text // '"' // trim(r_text(telemetry_initial_coefficients(j))) // '"'
    enddo
    text = text // ']'
  end function initial_coeff_text


  function initial_coeff_hex_text(degree) result(text)
    integer, intent(in) :: degree
    character(len=:), allocatable :: text
    integer :: j

    text = '['
    do j=1,degree
       if (j.gt.1) text = text // ','
       text = text // '"' // real_hex(telemetry_initial_coefficients(j)) // '"'
    enddo
    text = text // ']'
  end function initial_coeff_hex_text

end module rh_hardy_chain_telemetry
