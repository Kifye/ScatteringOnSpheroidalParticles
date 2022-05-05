! Contains constants and functions for general logging

module logging
    use constants
    implicit none
    !  length of logged strings
    !  normal message
    integer, parameter :: MESSAGE_LENGTH = 1024

    !  LOGS_TO_WRITE(i) shows whether the log of level
    !  i should be written
    logical, parameter :: LOGS_TO_WRITE(ERROR:DETAIL) = (/ &
            .true., &  ! - ERROR
            .true., &  ! - WARNING
            .true., &  ! - DEBUG
            .false., &  ! - INFO
            .false. &   ! - DETAIL
        /)
    !  file descriptors for general messages of giver levels
    !  if fd /= 6 (stdout) the log file should be set
    integer, parameter :: FILE_DESCRIPTOR(ERROR:DETAIL) = (/ 6, 6, 6, 6, 6 /)
    integer, parameter :: FD_GENERAL = 15
    !  log files for general messages if the respecpective descriptor is not set to 6(stdout)
    character(128), parameter :: LOG_FILE(ERROR:DETAIL) = (/ '', '', '', '', '' /)
!            '/home/drakosha/Documents/pulkovo/programs/ScatteringOnSpheroids/logs/SpheroidalCalculation.log'
contains

    !  the base subroutine to write a string message data with a given level
    !  into a stream indicated by the file descriptor fd
    !  trims the data
    subroutine log_message(data, level, fd)
        character (*) :: data
        integer :: level, fd

        if (LOGS_TO_WRITE(level)) then
            write(fd, *) trim(data)
        end if
    end subroutine log_message

    subroutine log_new_line(level, fd)
        integer :: level, fd

        if (LOGS_TO_WRITE(level)) then
            write(fd, *)
        end if
    end subroutine log_new_line

    !  subroutine to open the files for general messeges
    !  should be called at the start of the program if
    !  there is an i in [1..5] such as LOGS_TO_WRITE(i) = .true.
    !  and FD_BY_LEVEL(i) /= 6
    subroutine log_init()
        integer :: level

        do level = ERROR, DETAIL
            if (FILE_DESCRIPTOR(level) /= 6) then
                open(FILE_DESCRIPTOR(level), file=LOG_FILE(level))
            end if
        end do
    end subroutine log_init

    !  closes the files for general messeges
    !  should be called at the end of the program if log_init() was
    !  called at the beginning
    subroutine log_close()
        integer :: level

        do level = ERROR, DETAIL
            if (FILE_DESCRIPTOR(level) /= 6) then
                close(FILE_DESCRIPTOR(level))
            end if
        end do
    end subroutine log_close

    subroutine log_array(fd, name, array, border)
        character(*), intent(in) :: name
        integer, intent(in) :: fd
        integer, optional, intent(in) :: border
        complex(knd), dimension(:), intent(in) :: array

        integer :: length

        length = size(array)

        if (present(border) .and. border >= 0) then
            length = border
        end if

        write(fd, *) name, ' = ', array(1:length)
    end subroutine log_array

    subroutine log_matrix(fd, name, matrix, by_columns, border)
        character(*), intent(in) :: name
        integer, intent(in) :: fd
        integer, optional, intent(in) :: border
        complex(knd), dimension(:,:), intent(in) :: matrix
        logical, intent(in) :: by_columns

        integer :: border_row, border_column, i

        border_row = size(matrix, 1)
        border_column = size(matrix, 2)
        if (present(border)) then
            if (border > 0) then
            border_row = border
            border_column = border

            end if
        end if

        write(fd, *) name, ':'
        if (by_columns) then
            do i = 1, border_column
                write(fd, *) matrix(1:border_row, i)
            end do
        else
            do i = 1, border_row
                write(fd, *) matrix(i, 1:border_column)
            end do
        end if
    end subroutine log_matrix

end module logging