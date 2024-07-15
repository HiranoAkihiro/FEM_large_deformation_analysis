module mod_debug
    use mod_util
    implicit none
    
contains
    
    subroutine std_error_string(string)
        implicit none
        !> [in] 出力ログ
        character(*), intent(in) :: string
        character :: esc*1 = char(27)

        write(*,"(a,a)")esc//"[31m"//"** ERROR: "//esc//"[0m", trim(string)
    end subroutine std_error_string

    subroutine std_error_stop()
        implicit none
        error stop 1
    end subroutine std_error_stop

    subroutine plot_time(header, time)
        implicit none
        real(8) :: time
        character(*) :: header

        write(*,"(a,1pe10.3,a)")"  - "//trim(header)//" elapse time: ", time, " [sec]"
    end subroutine plot_time
end module mod_debug