module parameters_define
    implicit none
    character(LEN=40)::wrf_file


    contains
    subroutine read_nml(parameter_file)
        implicit none
        character(LEN=*)::parameter_file
        !set default value
        wrf_file="wrfinput_d01.2"

    end subroutine


end module parameters_define
