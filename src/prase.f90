MODULE prase
!====================================================================================
! purpose: Character prasing
!
! 
!====================================================================================
! created by Menaka
! Menaka@MSU 2026
!====================================================================================
!$ use omp_lib

contains
subroutine parse_csv_string(in_string, max_items, out_count, out_names)
  implicit none
  
  ! --- Arguments ---
  character(len=256), intent(in)  :: in_string  ! The input comma-separated string
  integer,          intent(in)    :: max_items  ! Maximum size of the output array
  integer,          intent(out)   :: out_count  ! Total items successfully parsed
  character(len=256), intent(out) :: out_names(max_items) ! Target array for results
  
  ! --- Internal Private Variables ---
  character(len=len_trim(in_string)+1) :: buffer
  integer :: comma_pos, c_idx
  
  ! Initialize counters and wipe output array clean
  out_count = 0
  out_names = ""
  
  buffer = trim(adjustl(in_string))
  
  ! Loop until buffer is empty or we hit the maximum allowed array size
  do while (len_trim(buffer) > 0 .and. out_count < max_items)
    out_count = out_count + 1
    comma_pos = scan(buffer, ',')
    
    ! Extract string segment
    if (comma_pos > 0) then
      out_names(out_count) = buffer(1:comma_pos-1)
      buffer = adjustl(buffer(comma_pos+1:len(buffer)))
    else
      out_names(out_count) = buffer
      buffer = ""
    end if
    
    ! Force contents of the extracted token to lowercase
    do c_idx = 1, len_trim(out_names(out_count))
      if (out_names(out_count)(c_idx:c_idx) >= 'A' .and. &
          out_names(out_count)(c_idx:c_idx) <= 'Z') then
          out_names(out_count)(c_idx:c_idx) = char(iachar(out_names(out_count)(c_idx:c_idx)) + 32)
      end if
    end do
  end do
return
end subroutine parse_csv_string


END MODULE prase