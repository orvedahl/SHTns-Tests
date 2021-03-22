!
! write out basic build information
!
! the idea was borrowed from the MAESTRO code and modified according to 
! my simple needs
!
! R. Orvedahl (2013-04-18)
!
module job_info

  contains

  subroutine write_job_info(dirname)

    ! write out some basic information about the way the job was run
    ! to a file called job_info in the directory dir_name.  Usually
    ! dir_name will be the name of the checkpoint or plotfile toplevel
    ! directory

    use data_types
    use input_params, only: job_name, inputs_file_used, runtime_print_params
    use build_info_module
  !, only: build_date, build_dir, build_machine, &
  !                               FCOMP, FCOMP_version, f90_compile_line, &
  !                               link_line, source_git_hash

    implicit none

    character (len=*), intent(in) :: dirname

    character (len=256) :: out_name, cwd
    character (len=16) :: date_in, time_in
    integer, dimension(8) :: values
    integer :: stat, io = 34

    call date_and_time(date_in, time_in, VALUES=values)
    stat = getcwd(cwd)
    if (stat /= 0) then
       write(*,*)
       write(*,*) "ERROR: could not get cwd in write_job_info.f90"
       write(*,*)
       stop
    endif
 
    out_name = trim(dirname) // "/job_info"

999  format(79('='))
1001 format(a,a)
1003 format(a,i4.4,'-',i2.2,'-',i2.2)
1004 format(a,i2.2,':',i2.2,':',i2.2)

    open(unit=io,file=out_name,form = "formatted", access = "sequential",&
          action="write")
     
    write (io,999)
    write (io,*) "Job Information"
    write (io,999)
    write (io,1001) "job name:    ", trim(job_name)
    write (io,1001) "inputs file: ", trim(inputs_file_used)
    write (io,*) " "     
    write (io,*) " "

    write (io,999)
    write (io,*) "Data/Time Information"
    write (io,999)
    write (io,1003) "output date:              ",values(1),values(2),values(3)
    write (io,1004) "output time:              ",values(5),values(6),values(7)
    write (io,1001) "output dir:               ", trim(cwd)

    write (io,*) " "
    write (io,*) " "

    write (io,999)
    write (io,*) "Build Information"
    write (io,999)
    write (io,1001) "build date:    ", trim(build_date)
    write (io,1001) "build machine: ", trim(build_machine)
    write (io,1001) "build dir:     ", trim(build_dir)
    write (io,*) " "
    write (io,1001) "Source     git hash: ", trim(source_git_hash)

    write (io,*) " "

    write (io,1001) "FCOMP:            ", trim(FCOMP)
    write (io,1001) "FCOMP version:    ", trim(FCOMP_version)
    write (io,*) " "
    write (io,1001) "F90 compile line: ", trim(f90_compile_line)
    write (io,*) " "
    write (io,1001) "linker line:      ", trim(link_line)

    write (io,*) " "
    write (io,*) " "

    write (io,999)
    write (io,*) "Runtime Parameter Information"
    write (io,999)

    call runtime_print_params(io)

    close(unit=io)

  end subroutine write_job_info

end module job_info

