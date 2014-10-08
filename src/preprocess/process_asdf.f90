module process_asdf_mod

contains
!> \Time-series Analysis for global adjoint tomography
!! Developer: Princeton Global Tomography Group(PGTG)
!! Group Member: James Smith jas11(@princeton.edu), Ebru Bozdag(bozdag@princeton.edu),
!! Wenjie Lei (lei@princeton.edu)
!! Bug Report: jas11@princeton.edu

subroutine process_asdf(observed_raw, synthetic_raw, observed_rotate, synthetic_rotate, rank, comm, nproc, adios_group, ierr)

  use mpi
  use asdf_data
  use asdf_read_subs
  use asdf_write_subs
  use process_var
  use process_subs
  use process_subs2
  use var_main

  implicit none

  integer                   :: comm, rank, ierr, nproc, adios_err
  integer(kind=8)           :: adios_handle, adios_group
  integer(kind=8)           :: adios_groupsize, adios_totalsize

  type (asdf_event),intent(inout) :: observed_raw, synthetic_raw
  type (asdf_event),intent(inout) :: observed_rotate, synthetic_rotate
  type (asdf_event)               :: observed_proc, synthetic_proc

  !local trace variables(In accordance with SAC lib, the data type is real)
  real(kind=4) :: b_obsd, b_synt, dt_obsd, dt_synt
  integer :: npts_obsd, npts_synt
  real(kind=4) :: obsd(MAX_TRACE_LENGTH), synt(MAX_TRACE_LENGTH)

  !after cut
  real(kind=4) :: b_obsd_cut, b_synt_cut, dt_obsd_cut, dt_synt_cut
  integer :: npts_obsd_cut, npts_synt_cut
  real(kind=4) :: obsd_cut(MAX_TRACE_LENGTH), synt_cut(MAX_TRACE_LENGTH)

  !remove trend and mean
  integer :: ipts
  real :: mean, yint, slope, siga, sigb, sig, cc, minnm, maxnm, increment, num

  !instrument response
  character(len=20) :: sta, net, locid, cha

  !after interpolation
  real(kind=4) :: b_obsd_interp, b_synt_interp, dt_obsd_interp, dt_synt_interp
  integer :: npts_obsd_interp, npts_synt_interp
  real(kind=4) :: obsd_interp(MAX_TRACE_LENGTH), synt_interp(MAX_TRACE_LENGTH)
  real :: end_time_syn, end_time_obs

  real :: epsi = 1e-8

  !after final cut
  real(kind=4) :: b_obsd_final, b_synt_final, dt_obsd_final, dt_synt_final
  integer :: npts_obsd_final, npts_synt_final
  real(kind=4) :: obsd_final(MAX_TRACE_LENGTH), synt_final(MAX_TRACE_LENGTH)

  integer :: nerr

!-----------------------------------------------------------------------------------
! READ CMT FILE AND COMPUTE JULIAN DAY AND ORIGIN OF EVENT
!-----------------------------------------------------------------------------------
print *, "ok"
  call read_CMT(observed_raw%event, rank, comm, ierr)
print *, "ok2"
  call event_origin(rank, comm, ierr)

!-----------------------------------------------------------------------------------
! ADJUST TIME OF OBSERVED DATA BASED ON CMT EVENT ORIGIN
!-----------------------------------------------------------------------------------
  print *, "observed origin time adjusted based on cmt event, rank:", rank
  call adjust_event_time(observed_raw, ierr)

  call init_asdf_data(observed_proc, observed_raw%nrecords, .false.)
  call init_asdf_data(synthetic_proc, synthetic_raw%nrecords, .false.)

!-----------------------------------------------------------------------------------
! LOOOP THROUGH PAIRS of OBSERVED AND SYNTHETIC TRACES
!-----------------------------------------------------------------------------------
  do irecord = 1, observed_raw%nrecords
    print *, "Processing record ", irecord, " of ", observed_raw%nrecords, " on processor " , rank, "."
    print *, "record_name: ", trim(observed_raw%receiver_name_array(irecord))&
      //"."//trim(observed_raw%network_array(irecord))&
      //"."//trim(observed_raw%receiver_id_array(irecord))&
      //"."//trim(observed_raw%component_array(irecord))

  !-----------------------------------------------------------------------------------
  ! CUT OBSERVED AND SYNTHETIC SEISMOGRAM
  !-----------------------------------------------------------------------------------
    print *, "Cut begin...rank:", rank
    !copy trace info to local var
    b_obsd=real(observed_raw%begin_value(irecord))
    dt_obsd=real(observed_raw%sample_rate(irecord))
    npts_obsd=observed_raw%npoints(irecord)
    obsd(1:npts_obsd)=real(observed_raw%records(irecord)%record(1:npts_obsd))

    b_synt=real(synthetic_raw%begin_value(irecord))
    dt_synt=real(synthetic_raw%sample_rate(irecord))
    npts_synt=synthetic_raw%npoints(irecord)
    if (npts_synt .eq. 1) cycle
    synt(1:npts_synt)=real(synthetic_raw%records(irecord)%record(1:npts_synt))

    sta=observed_raw%receiver_name_array(irecord)
    net=observed_raw%network_array(irecord)
    locid=observed_raw%receiver_id_array(irecord)
    cha=observed_raw%component_array(irecord)
    print *, locid, cha
    !call wsac1("adjust_event_time.obs", obsd, npts_obsd, &
    !      b_obsd, dt_obsd, nerr)
    !call wsac1("adjust_event_time.syn", synt, npts_synt, &
    !      b_synt, dt_synt, nerr)

    call cut_seis(obsd, npts_obsd, dt_obsd, b_obsd, synt, npts_synt, dt_synt, b_synt, &
    obsd_cut, npts_obsd_cut, dt_obsd_cut, b_obsd_cut, &
    synt_cut, npts_synt_cut, dt_synt_cut, b_synt_cut)
    print *, "Cut done...rank:", rank

    !call wsac1("cut.obs", obsd_cut, npts_obsd_cut, &
    !      b_obsd_cut, dt_obsd_cut, nerr)
    !call wsac1("cut.syn", synt_cut, npts_synt_cut, &
    !      b_synt_cut, dt_synt_cut, nerr)

  !-----------------------------------------------------------------------------------
  ! PREPARE FOR DECONVOLUTION: rmean/rtrend/taper
  !-----------------------------------------------------------------------------------
    mean=0
    slope=0
    yint=0
    !call extrma(obsd_cut, increment, num, minnm, maxnm, mean)
    !call rmean(obsd_cut, npts_obsd_cut, mean)
    call demean(obsd_cut, npts_obsd_cut, mean)
    !call lifite(b_obsd_cut, dt_obsd_cut, obsd_cut, npts_obsd_cut,&
    !      slope, yint, siga, sigb, sig, cc)
    call detrend(obsd_cut, npts_obsd_cut, yint, slope, &
            b_obsd_cut, dt_obsd_cut)
    call taper_width_to_points(.05, real(npts_obsd_cut), ipts)
    call taper(obsd_cut, npts_obsd_cut, 2, ipts)
    print *, "npoints:", npts_obsd_cut
    print *, "ipts:", ipts
    print *, "mean:", sum(obsd_cut(1:npts_obsd_cut))/npts_obsd_cut
    print *, "npts:", npts_obsd_cut
    print *, "mean, yint, slope:", mean, yint, slope

    mean=0
    slope=0
    yint=0
    !call extrma(synt_cut, increment, num, minnm, maxnm, mean)
    !call rmean(synt_cut, npts_synt_cut, mean)
    call demean(synt_cut, npts_synt_cut, mean)
    !call lifite(b_synt_cut, dt_synt_cut, synt_cut, npts_synt_cut,&
    !      slope, yint, siga, sigb, sig, cc)
    call detrend(synt_cut, npts_synt_cut, yint, slope, &
             b_synt_cut, dt_synt_cut)
    call taper_width_to_points(.05, real(npts_synt_cut), ipts)
    call taper(synt_cut, npts_synt_cut, 2, ipts)

    !print *, "npts:", npts_synt_cut
    !print *, "mean, yint, slope:", mean, yint, slope

    !call wsac1("rmean.obs", obsd_cut, npts_obsd_cut, &
    !      b_obsd_cut, dt_obsd_cut, nerr)
    !call wsac1("rmean.syn", synt_cut, npts_synt_cut, &
    !      b_synt_cut, dt_synt_cut, nerr)

    print *, "Taper done and ready to decon. Rank:", rank

    !call wsac1("before_filter.syn", synt_cut, npts_synt_cut, &
    !      b_synt_cut, dt_synt_cut, nerr)

    !-----------------------------------------------------------------------------------
    ! REMOVE INSTRUMENT RESPONSE
    !-----------------------------------------------------------------------------------
    call remove_response(obsd_cut, npts_obsd_cut, dt_obsd_cut, &
      sta, net, cha, locid, observed_raw%responses(irecord)%response_string, &
      observed_raw%responses(irecord)%response_length)
    call filter(synt_cut, npts_synt_cut, dt_synt_cut)

    !call wsac1("filter.obs", obsd_cut, npts_obsd_cut, &
    !      b_obsd_cut, dt_obsd_cut, nerr)
    !call wsac1("filter.syn", synt_cut, npts_synt_cut, &
    !      b_synt_cut, dt_synt_cut, nerr)
    print *, "Response removed. Rank:", rank

    !-----------------------------------------------------------------------------------
    ! rmean/rtrend/taper again in case instrument response introduces artifacts
    !-----------------------------------------------------------------------------------

    call demean(obsd_cut, npts_obsd_cut, mean)
    !call lifite(b_obsd_cut, dt_obsd_cut, obsd_cut, npts_obsd_cut,&
    !      slope, yint, siga, sigb, sig, cc)
    call detrend(obsd_cut, npts_obsd_cut, yint, slope, &
               b_obsd_cut, dt_obsd_cut)
    call taper_width_to_points(.05, real(npts_obsd_cut), ipts)
    call taper(obsd_cut, npts_obsd_cut, 2, ipts)

    !call demean(synt_cut, npts_synt_cut, mean)
    !call lifite(b_synt_cut, dt_synt_cut, synt_cut, npts_synt_cut,&
    !      slope, yint, siga, sigb, sig, cc)
    !call detrend(synt_cut, npts_synt_cut, yint, slope, &
    !            b_synt_cut, dt_synt_cut)
    !call taper_width_to_points(.05, sngl(npts_synt_cut), ipts)
    !call taper(synt_cut, npts_synt_cut, 2, ipts)
    !print *, "prepared for deconvolution"
    !call wsac1("taper.syn", synt_cut, npts_synt_cut,&
    !      b_synt_cut, dt_synt_cut, nerr)
    !call wsac1("taper.obs", obsd_cut, npts_obsd_cut,&
             !b_obsd_cut, dt_obsd_cut, nerr)

    print *, "Second Taper done. Rank:", rank

  !-----------------------------------------------------------------------------------
  ! RESAMPLE OBSERVED AND SYNTHETIC SEISMOGRAMS TO SAME SAMPLE RATE
  !-----------------------------------------------------------------------------------

    end_time_syn=b_synt_cut+dt_synt_cut*float(npts_synt_cut)
    end_time_obs=b_obsd_cut+dt_obsd_cut*float(npts_obsd_cut)
    b_synt_interp=b_synt_cut
    b_obsd_interp=b_obsd_cut
    dt_synt_interp = synthetic_raw%sample_rate(irecord)
    dt_obsd_interp = synthetic_raw%sample_rate(irecord)
    npts_synt_interp = int(npts_synt_cut*dt_synt_cut/dt_synt_interp)
    npts_obsd_interp = int(npts_obsd_cut*dt_obsd_cut/dt_obsd_interp)
    call interp(synt_cut, npts_synt_cut, synt_interp, npts_synt_interp,&
                b_synt_cut, end_time_syn, dt_synt_cut, &
                b_synt_interp, dt_synt_interp, epsi)
    call interp(obsd_cut, npts_obsd_cut, obsd_interp, npts_obsd_interp,&
                b_obsd_cut, end_time_obs, dt_obsd_cut, &
                b_obsd_interp, dt_obsd_interp, epsi)

    print *, "Resample complete. Rank:", rank
    !call wsac1("final.syn", synt_interp, npts_synt_interp,&
    !      b_synt_interp, dt_synt_interp, nerr)
    !call wsac1("final.obs", obsd_interp, npts_obsd_interp,&
    !         b_obsd_interp, dt_obsd_interp, nerr)

    !-----------------------------------------------------------------------------------
    ! FINAL CUT OF OBSERVED AND SYNTHETIC SEISMOGRAMS
    !-----------------------------------------------------------------------------------

    call cut_seis(obsd_interp, npts_obsd_interp, dt_obsd_interp, b_obsd_interp, &
      synt_cut, npts_synt_interp, dt_synt_interp, b_synt_interp, &
      obsd_final, npts_obsd_final, dt_obsd_final, b_obsd_final, &
      synt_final, npts_synt_final, dt_synt_final, b_synt_final)
    print *, "final cut finished"
!    
    observed_proc%begin_value(irecord)=b_obsd_interp
    observed_proc%sample_rate(irecord)=dt_obsd_interp
    observed_proc%npoints(irecord)=npts_obsd_interp
    allocate(observed_proc%records(irecord)%record(npts_obsd_interp))
    observed_proc%records(irecord)%record(1:npts_obsd_interp)=dble(obsd_interp(1:npts_obsd_interp))
!
    synthetic_proc%begin_value(irecord)=b_synt_interp
    synthetic_proc%sample_rate(irecord)=dt_synt_interp
    synthetic_proc%npoints(irecord)=npts_synt_interp
    allocate(synthetic_proc%records(irecord)%record(npts_synt_interp))
    synthetic_proc%records(irecord)%record(1:npts_synt_interp)=dble(synt_interp(1:npts_synt_interp))
    !endif
  enddo

  !-----------------------------------------------------------------------------------
  ! COPY METADATA TO ASDF CONTAINER
  !-----------------------------------------------------------------------------------
  call store_asdf_metadata(observed_raw, observed_proc)
  call store_asdf_metadata(synthetic_raw, synthetic_proc)

  !-----------------------------------------------------------------------------------
  ! Write out file before rotate
  !----------------------------------------------------------------------------------- 
  call write_asdf_file("obs.bp", observed_proc, adios_group, rank, nproc, comm, ierr)
  call write_asdf_file("syn.bp", synthetic_proc, adios_group, rank, nproc, comm, ierr)

  !-----------------------------------------------------------------------------------
  ! ROTATE
  !----------------------------------------------------------------------------------- 
  print *, "Rotation begin, rank:", rank
  call rotate_traces(observed_proc, observed_rotate)
  call rotate_traces(synthetic_proc, synthetic_rotate)
  print *, "Rotation done rank:", rank

end subroutine process_asdf

end module
