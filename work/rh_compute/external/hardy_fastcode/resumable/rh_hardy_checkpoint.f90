! RH-owned restart controls for the pinned GPL-3.0 Hardy fastcode.
! Distributed with the derivative under GPL-3.0; see COPYING.GPL-3.0.txt.
module rh_hardy_checkpoint
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, parameter, public :: rh_dp = selected_real_kind(33)
  integer, parameter, public :: rh_i8 = selected_int_kind(16)
  integer, parameter, public :: rh_stage_zp = 0
  integer, parameter, public :: rh_stage_rs = 1
  integer, parameter, public :: rh_stage_complete = 2
  integer, parameter :: rh_serial_significand_digits = 41
  character(len=*), parameter :: rh_checkpoint_magic = 'RH_HARDY_CHECKPOINT_V1'

  type, public :: rh_checkpoint_config
     character(len=512) :: checkpoint_path = 'rh_hardy_checkpoint.dat'
     character(len=512) :: journal_path = 'rh_hardy_checkpoint.dat.jsonl'
     character(len=512) :: stop_path = 'rh_hardy_checkpoint.stop'
     character(len=64) :: run_id = ''
     logical :: resume = .false.
     real(kind=8) :: max_seconds = 0.0d0
     integer :: stop_after_stage = -1
     integer(rh_i8) :: stop_after_unit = -1_rh_i8
     integer(rh_i8) :: rs_chunk_size = 100000_rh_i8
     integer(kind=8) :: start_count = 0_8
     integer(kind=8) :: count_rate = 0_8
  end type rh_checkpoint_config

  type, public :: rh_checkpoint_state
     integer :: stage = rh_stage_zp
     integer(rh_i8) :: unit = 0_rh_i8
     integer :: comm_size = 1
     integer(rh_i8) :: numbercalc = 0_rh_i8
     integer(rh_i8) :: nchalf = 0_rh_i8
     integer :: totblock = 0
     integer(rh_i8) :: mt = 0_rh_i8
     integer(rh_i8) :: mto = 0_rh_i8
     integer(rh_i8) :: aenums = 0_rh_i8
     integer :: mmax = 0
     integer(rh_i8) :: nc = 0_rh_i8
     integer(rh_i8) :: nmax = 0_rh_i8
     integer(rh_i8) :: rs_chunk_size = 100000_rh_i8
     real(rh_dp) :: t = 0.0_rh_dp
     real(rh_dp) :: et = 0.0_rh_dp
     real(rh_dp) :: rn1 = 0.0_rh_dp
     real(rh_dp) :: zsum(15) = 0.0_rh_dp
     real(rh_dp) :: rszsumtot(15) = 0.0_rh_dp
     real(rh_dp) :: rszsum(15) = 0.0_rh_dp
  end type rh_checkpoint_state

  public :: rh_checkpoint_configure
  public :: rh_checkpoint_begin
  public :: rh_checkpoint_capture
  public :: rh_checkpoint_load
  public :: rh_checkpoint_commit

contains

  subroutine rh_checkpoint_configure(cfg, status, message)
    type(rh_checkpoint_config), intent(out) :: cfg
    integer, intent(out) :: status
    character(len=*), intent(out) :: message
    character(len=512) :: value
    integer :: env_status, ios

    cfg = rh_checkpoint_config()
    status = 0
    message = ''

    call get_environment_variable('RH_HARDY_CHECKPOINT', value, status=env_status)
    if (env_status == 0 .and. len_trim(value) > 0) cfg%checkpoint_path = trim(value)

    call get_environment_variable('RH_HARDY_JOURNAL', value, status=env_status)
    if (env_status == 0 .and. len_trim(value) > 0) then
       cfg%journal_path = trim(value)
    else
       cfg%journal_path = trim(cfg%checkpoint_path)//'.jsonl'
    endif

    call get_environment_variable('RH_HARDY_STOP_FILE', value, status=env_status)
    if (env_status == 0 .and. len_trim(value) > 0) then
       cfg%stop_path = trim(value)
    else
       cfg%stop_path = trim(cfg%checkpoint_path)//'.stop'
    endif

    call get_environment_variable('RH_HARDY_RUN_ID', value, status=env_status)
    if (env_status /= 0 .or. len_trim(value) /= 64 .or. &
        verify(trim(value), '0123456789abcdefABCDEF') /= 0) then
       status = 1
       message = 'RH_HARDY_RUN_ID must be a 64-character hexadecimal provenance hash'
       return
    endif
    cfg%run_id = trim(value)

    call get_environment_variable('RH_HARDY_RESUME', value, status=env_status)
    if (env_status == 0 .and. len_trim(value) > 0) then
       select case (adjustl(trim(value)))
       case ('1', 'true', 'TRUE', 'yes', 'YES')
          cfg%resume = .true.
       case ('0', 'false', 'FALSE', 'no', 'NO')
          cfg%resume = .false.
       case default
          status = 1
          message = 'RH_HARDY_RESUME must be 0/1, false/true, or no/yes'
          return
       end select
    endif

    call get_environment_variable('RH_HARDY_MAX_SECONDS', value, status=env_status)
    if (env_status == 0 .and. len_trim(value) > 0) then
       read(value, *, iostat=ios) cfg%max_seconds
       if (ios /= 0 .or. cfg%max_seconds < 0.0d0) then
          status = 1
          message = 'RH_HARDY_MAX_SECONDS must be a nonnegative real number'
          return
       endif
    endif

    call get_environment_variable('RH_HARDY_STOP_AFTER_STAGE', value, status=env_status)
    if (env_status == 0 .and. len_trim(value) > 0) then
       read(value, *, iostat=ios) cfg%stop_after_stage
       if (ios /= 0 .or. cfg%stop_after_stage < rh_stage_zp .or. &
           cfg%stop_after_stage > rh_stage_rs) then
          status = 1
          message = 'RH_HARDY_STOP_AFTER_STAGE must be 0 or 1'
          return
       endif
    endif

    call get_environment_variable('RH_HARDY_STOP_AFTER_UNIT', value, status=env_status)
    if (env_status == 0 .and. len_trim(value) > 0) then
       read(value, *, iostat=ios) cfg%stop_after_unit
       if (ios /= 0 .or. cfg%stop_after_unit < 0_rh_i8) then
          status = 1
          message = 'RH_HARDY_STOP_AFTER_UNIT must be a nonnegative integer'
          return
       endif
    endif

    call get_environment_variable('RH_HARDY_RS_CHUNK_SIZE', value, status=env_status)
    if (env_status == 0 .and. len_trim(value) > 0) then
       read(value, *, iostat=ios) cfg%rs_chunk_size
       if (ios /= 0 .or. cfg%rs_chunk_size < 1_rh_i8) then
          status = 1
          message = 'RH_HARDY_RS_CHUNK_SIZE must be a positive integer'
          return
       endif
    endif

    if ((cfg%stop_after_stage < 0) .neqv. (cfg%stop_after_unit < 0_rh_i8)) then
       status = 1
       message = 'RH_HARDY_STOP_AFTER_STAGE and RH_HARDY_STOP_AFTER_UNIT must be set together'
       return
    endif

    if (index(cfg%checkpoint_path, '"') > 0 .or. index(cfg%journal_path, '"') > 0) then
       status = 1
       message = 'checkpoint and journal paths may not contain a double quote'
       return
    endif

    call system_clock(count=cfg%start_count, count_rate=cfg%count_rate)
  end subroutine rh_checkpoint_configure


  subroutine rh_checkpoint_begin(cfg, status, message)
    type(rh_checkpoint_config), intent(in) :: cfg
    integer, intent(out) :: status
    character(len=*), intent(out) :: message
    logical :: exists

    status = 0
    message = ''
    if (cfg%resume) then
       inquire(file=trim(cfg%checkpoint_path), exist=exists)
       if (.not. exists) then
          status = 1
          message = 'resume requested but checkpoint snapshot does not exist'
       endif
       return
    endif

    call rh_delete_file(trim(cfg%checkpoint_path), status)
    if (status /= 0) then
       message = 'could not remove stale checkpoint snapshot'
       return
    endif
    call rh_delete_file(trim(cfg%checkpoint_path)//'.tmp', status)
    if (status /= 0) then
       message = 'could not remove stale checkpoint temporary file'
       return
    endif
    call rh_delete_file(trim(cfg%journal_path), status)
    if (status /= 0) message = 'could not remove stale checkpoint journal'
  end subroutine rh_checkpoint_begin


  subroutine rh_checkpoint_capture(state, stage, unit_number, comm_size, numbercalc, &
       nchalf, totblock, t, et, rn1, mt, mto, aenums, mmax, nc, nmax, zsum, &
       rszsumtot, rszsum, rs_chunk_size)
    type(rh_checkpoint_state), intent(out) :: state
    integer, intent(in) :: stage, comm_size, totblock, mmax
    integer(rh_i8), intent(in) :: unit_number, numbercalc, nchalf, mt, mto, aenums, nc, nmax
    integer(rh_i8), intent(in) :: rs_chunk_size
    real(rh_dp), intent(in) :: t, et, rn1
    real(rh_dp), intent(in) :: zsum(15), rszsumtot(15), rszsum(15)

    state%stage = stage
    state%unit = unit_number
    state%comm_size = comm_size
    state%numbercalc = numbercalc
    state%nchalf = nchalf
    state%totblock = totblock
    state%t = t
    state%et = et
    state%rn1 = rn1
    state%mt = mt
    state%mto = mto
    state%aenums = aenums
    state%mmax = mmax
    state%nc = nc
    state%nmax = nmax
    state%rs_chunk_size = rs_chunk_size
    state%zsum = zsum
    state%rszsumtot = rszsumtot
    state%rszsum = rszsum
  end subroutine rh_checkpoint_capture


  subroutine rh_checkpoint_load(cfg, expected_t, expected_et, expected_numbercalc, &
       expected_nchalf, expected_totblock, expected_comm_size, state, status, message)
    type(rh_checkpoint_config), intent(in) :: cfg
    real(rh_dp), intent(in) :: expected_t, expected_et
    integer(rh_i8), intent(in) :: expected_numbercalc, expected_nchalf
    integer, intent(in) :: expected_totblock, expected_comm_size
    type(rh_checkpoint_state), intent(out) :: state
    integer, intent(out) :: status
    character(len=*), intent(out) :: message
    character(len=64) :: magic
    character(len=64) :: stored_run_id
    integer :: unit_number, ios, i
    real(rh_dp) :: stored_checksum(3), actual_checksum(3)

    state = rh_checkpoint_state()
    status = 0
    message = ''
    open(newunit=unit_number, file=trim(cfg%checkpoint_path), status='old', &
         action='read', form='formatted', iostat=ios)
    if (ios /= 0) then
       status = 1
       message = 'could not open checkpoint snapshot'
       return
    endif

    read(unit_number, '(A)', iostat=ios) magic
    if (ios == 0) read(unit_number, '(A)', iostat=ios) stored_run_id
    if (ios == 0) read(unit_number, *, iostat=ios) state%stage, state%unit, &
         state%comm_size, state%numbercalc, state%nchalf, state%totblock
    if (ios == 0) read(unit_number, *, iostat=ios) state%t, state%et, state%rn1
    if (ios == 0) read(unit_number, *, iostat=ios) state%mt, state%mto, &
         state%aenums, state%mmax, state%nc, state%nmax, state%rs_chunk_size
    if (ios == 0) read(unit_number, *, iostat=ios) (state%zsum(i), i=1,15)
    if (ios == 0) read(unit_number, *, iostat=ios) (state%rszsumtot(i), i=1,15)
    if (ios == 0) read(unit_number, *, iostat=ios) (state%rszsum(i), i=1,15)
    if (ios == 0) read(unit_number, *, iostat=ios) stored_checksum
    close(unit_number)

    if (ios /= 0) then
       status = 1
       message = 'checkpoint snapshot is truncated or malformed'
       return
    endif
    if (trim(magic) /= rh_checkpoint_magic) then
       status = 1
       message = 'checkpoint source version does not match this executable'
       return
    endif
    if (trim(stored_run_id) /= trim(cfg%run_id)) then
       status = 1
       message = 'checkpoint provenance hash does not match this run'
       return
    endif
    if (state%comm_size /= expected_comm_size .or. &
        state%numbercalc /= expected_numbercalc .or. &
        state%nchalf /= expected_nchalf .or. state%totblock /= expected_totblock) then
       status = 1
       message = 'checkpoint integer metadata does not match this run'
       return
    endif
    if (state%t /= expected_t .or. state%et /= expected_et) then
       status = 1
       message = 'checkpoint height or tolerance does not match this run'
       return
    endif
    if (state%rs_chunk_size /= cfg%rs_chunk_size) then
       status = 1
       message = 'checkpoint RS chunk size does not match this run'
       return
    endif
    if (state%stage < rh_stage_zp .or. state%stage > rh_stage_complete) then
       status = 1
       message = 'checkpoint stage is outside the supported range'
       return
    endif
    if (state%stage == rh_stage_zp .and. &
        (state%unit < 0_rh_i8 .or. state%unit > int(state%totblock, rh_i8))) then
       status = 1
       message = 'checkpoint Gaussian-block index is outside the supported range'
       return
    endif
    if (.not. rh_state_is_finite(state)) then
       status = 1
       message = 'checkpoint contains a non-finite real value'
       return
    endif

    call rh_state_checksums(state, actual_checksum)
    if (any(actual_checksum /= stored_checksum)) then
       status = 1
       message = 'checkpoint checksum mismatch'
    endif
  end subroutine rh_checkpoint_load


  subroutine rh_checkpoint_commit(cfg, state, stop_requested, status, message)
    type(rh_checkpoint_config), intent(in) :: cfg
    type(rh_checkpoint_state), intent(in) :: state
    integer, intent(out) :: stop_requested, status
    character(len=*), intent(out) :: message

    stop_requested = 0
    status = 0
    message = ''
    if (.not. rh_state_is_finite(state)) then
       status = 1
       message = 'refusing to checkpoint non-finite state'
       return
    endif

    call rh_write_snapshot(cfg, state, status, message)
    if (status /= 0) return
    call rh_append_journal(cfg, state, status, message)
    if (status /= 0) return
    call rh_check_stop(cfg, state%stage, state%unit, stop_requested, message)
  end subroutine rh_checkpoint_commit


  subroutine rh_write_snapshot(cfg, state, status, message)
    type(rh_checkpoint_config), intent(in) :: cfg
    type(rh_checkpoint_state), intent(in) :: state
    integer, intent(out) :: status
    character(len=*), intent(out) :: message
    character(len=1024) :: temporary_path
    integer :: unit_number, ios, rename_status, i
    real(rh_dp) :: checksum(3)

    status = 0
    message = ''
    temporary_path = trim(cfg%checkpoint_path)//'.tmp'
    call rh_state_checksums(state, checksum)
    open(newunit=unit_number, file=trim(temporary_path), status='replace', &
         action='write', form='formatted', iostat=ios)
    if (ios /= 0) then
       status = 1
       message = 'could not open temporary checkpoint snapshot'
       return
    endif

    write(unit_number, '(A)', iostat=ios) rh_checkpoint_magic
    if (ios == 0) write(unit_number, '(A)', iostat=ios) trim(cfg%run_id)
    if (ios == 0) write(unit_number, *, iostat=ios) state%stage, state%unit, &
         state%comm_size, state%numbercalc, state%nchalf, state%totblock
    if (ios == 0) write(unit_number, '(3(ES49.40E4,1X))', iostat=ios) &
         state%t, state%et, state%rn1
    if (ios == 0) write(unit_number, *, iostat=ios) state%mt, state%mto, &
         state%aenums, state%mmax, state%nc, state%nmax, state%rs_chunk_size
    if (ios == 0) write(unit_number, '(15(ES49.40E4,1X))', iostat=ios) &
         (state%zsum(i), i=1,15)
    if (ios == 0) write(unit_number, '(15(ES49.40E4,1X))', iostat=ios) &
         (state%rszsumtot(i), i=1,15)
    if (ios == 0) write(unit_number, '(15(ES49.40E4,1X))', iostat=ios) &
         (state%rszsum(i), i=1,15)
    if (ios == 0) write(unit_number, '(3(ES49.40E4,1X))', iostat=ios) checksum
    flush(unit_number)
    close(unit_number, iostat=rename_status)
    if (ios /= 0 .or. rename_status /= 0) then
       status = 1
       message = 'could not flush temporary checkpoint snapshot'
       return
    endif

    call rh_sync_file(trim(temporary_path), status)
    if (status /= 0) then
       message = 'could not fsync temporary checkpoint snapshot'
       return
    endif
    call rename(trim(temporary_path), trim(cfg%checkpoint_path), rename_status)
    if (rename_status /= 0) then
       status = 1
       message = 'could not atomically replace checkpoint snapshot'
       return
    endif
    call rh_sync_file(trim(cfg%checkpoint_path), status)
    if (status /= 0) message = 'could not fsync checkpoint snapshot'
  end subroutine rh_write_snapshot


  subroutine rh_append_journal(cfg, state, status, message)
    type(rh_checkpoint_config), intent(in) :: cfg
    type(rh_checkpoint_state), intent(in) :: state
    integer, intent(out) :: status
    character(len=*), intent(out) :: message
    integer :: unit_number, ios, i

    status = 0
    message = ''
    open(newunit=unit_number, file=trim(cfg%journal_path), status='unknown', &
         position='append', action='write', form='formatted', iostat=ios)
    if (ios /= 0) then
       status = 1
       message = 'could not open append-only checkpoint journal'
       return
    endif

    write(unit_number, &
         '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)', &
         advance='no', iostat=ios) &
         '{"format":"RH_HARDY_CHECKPOINT_V1","run_id":"'//trim(cfg%run_id)// &
         '","real_radix":', radix(state%t), ',"real_digits":', digits(state%t), &
         ',"serialized_significand_digits":', rh_serial_significand_digits, &
         ',"stage":', state%stage, &
         ',"unit":', state%unit, ',"comm_size":', state%comm_size, &
         ',"numbercalc":', state%numbercalc, ',"nchalf":', state%nchalf, &
         ',"totblock":', state%totblock, ',"t":'
    if (ios == 0) write(unit_number, '(ES49.40E4,A)', advance='no', iostat=ios) &
         state%t, ',"et":'
    if (ios == 0) write(unit_number, '(ES49.40E4,A)', advance='no', iostat=ios) &
         state%et, ',"rn1":'
    if (ios == 0) write(unit_number, &
         '(ES49.40E4,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)', &
         advance='no', iostat=ios) state%rn1, ',"mt":', state%mt, &
         ',"mto":', state%mto, ',"aenums":', state%aenums, &
         ',"mmax":', state%mmax, ',"nc":', state%nc, ',"nmax":', &
         state%nmax, ',"rs_chunk_size":', state%rs_chunk_size, ',"zsum":['
    do i=1,15
       if (ios /= 0) exit
       if (i > 1) write(unit_number, '(A)', advance='no', iostat=ios) ','
       if (ios == 0) write(unit_number, '(ES49.40E4)', advance='no', iostat=ios) state%zsum(i)
    enddo
    if (ios == 0) write(unit_number, '(A)', advance='no', iostat=ios) '],"rszsumtot":['
    do i=1,15
       if (ios /= 0) exit
       if (i > 1) write(unit_number, '(A)', advance='no', iostat=ios) ','
       if (ios == 0) write(unit_number, '(ES49.40E4)', advance='no', iostat=ios) state%rszsumtot(i)
    enddo
    if (ios == 0) write(unit_number, '(A)', advance='no', iostat=ios) '],"rszsum":['
    do i=1,15
       if (ios /= 0) exit
       if (i > 1) write(unit_number, '(A)', advance='no', iostat=ios) ','
       if (ios == 0) write(unit_number, '(ES49.40E4)', advance='no', iostat=ios) state%rszsum(i)
    enddo
    if (ios == 0) write(unit_number, '(A)', iostat=ios) ']}'
    flush(unit_number)
    close(unit_number)
    if (ios /= 0) then
       status = 1
       message = 'could not append a complete checkpoint journal record'
       return
    endif
    call rh_sync_file(trim(cfg%journal_path), status)
    if (status /= 0) message = 'could not fsync checkpoint journal'
  end subroutine rh_append_journal


  subroutine rh_check_stop(cfg, stage, unit_number, stop_requested, message)
    type(rh_checkpoint_config), intent(in) :: cfg
    integer, intent(in) :: stage
    integer(rh_i8), intent(in) :: unit_number
    integer, intent(out) :: stop_requested
    character(len=*), intent(out) :: message
    logical :: stop_exists
    integer(kind=8) :: now_count
    real(kind=8) :: elapsed

    stop_requested = 0
    message = ''
    inquire(file=trim(cfg%stop_path), exist=stop_exists)
    if (stop_exists) then
       stop_requested = 1
       message = 'stop file observed after committed atomic unit'
       return
    endif
    if (cfg%stop_after_stage == stage .and. unit_number >= cfg%stop_after_unit) then
       stop_requested = 1
       message = 'configured stage/unit stop reached after committed atomic unit'
       return
    endif
    if (cfg%max_seconds > 0.0d0 .and. cfg%count_rate > 0_8) then
       call system_clock(count=now_count)
       elapsed = real(now_count-cfg%start_count, kind=8)/real(cfg%count_rate, kind=8)
       if (elapsed >= cfg%max_seconds) then
          stop_requested = 1
          message = 'runtime limit reached after committed atomic unit'
       endif
    endif
  end subroutine rh_check_stop


  subroutine rh_state_checksums(state, checksum)
    type(rh_checkpoint_state), intent(in) :: state
    real(rh_dp), intent(out) :: checksum(3)
    integer :: i

    checksum = 0.0_rh_dp
    do i=1,15
       checksum(1) = checksum(1) + real(i, rh_dp)*state%zsum(i)
       checksum(2) = checksum(2) + real(i+17, rh_dp)*state%rszsumtot(i)
       checksum(3) = checksum(3) + real(i+37, rh_dp)*state%rszsum(i)
    enddo
  end subroutine rh_state_checksums


  logical function rh_state_is_finite(state)
    type(rh_checkpoint_state), intent(in) :: state

    rh_state_is_finite = ieee_is_finite(state%t) .and. ieee_is_finite(state%et) .and. &
         ieee_is_finite(state%rn1) .and. all(ieee_is_finite(state%zsum)) .and. &
         all(ieee_is_finite(state%rszsumtot)) .and. all(ieee_is_finite(state%rszsum))
  end function rh_state_is_finite


  subroutine rh_sync_file(path, status)
    character(len=*), intent(in) :: path
    integer, intent(out) :: status
    character(len=1200) :: command

    command = 'sync -f -- "'//trim(path)//'"'
    call execute_command_line(trim(command), wait=.true., exitstat=status)
  end subroutine rh_sync_file


  subroutine rh_delete_file(path, status)
    character(len=*), intent(in) :: path
    integer, intent(out) :: status
    integer :: unit_number, ios
    logical :: exists

    status = 0
    inquire(file=trim(path), exist=exists)
    if (.not. exists) return
    open(newunit=unit_number, file=trim(path), status='old', action='readwrite', iostat=ios)
    if (ios /= 0) then
       status = ios
       return
    endif
    close(unit_number, status='delete', iostat=ios)
    if (ios /= 0) status = ios
  end subroutine rh_delete_file

end module rh_hardy_checkpoint
