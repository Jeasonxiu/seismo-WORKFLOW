module process_subs2

  implicit none

  double precision :: FLT_MIN=1e-30

contains

!======================================================================
!> Removes the instrument response from the observed seismogram
!! \param observed The observed seismogram
!! \param synthetic The synthetichetic seismogram
subroutine remove_response(obsd, npts_obsd, dt_obsd, sta, net, cha, locid, &
  response, resp_length)

  use var_main
  use process_var
  use, intrinsic :: iso_c_binding
  implicit none

  real,dimension(*),intent(inout):: obsd
  integer, intent(in) :: npts_obsd
  real(kind=4),intent(in)   :: dt_obsd
  character(len=*) :: sta, net, cha, locid, response
  integer :: resp_length

  !local var
  integer :: LNPT, FFT_NPTS
  integer :: i, j, ierr

  !evresp var
  character*200        :: datime
  character(len=200)   :: temp_resp_fn
  character(len=10)    :: units, rtype
  integer              :: start_stage=-1,stop_stage=0,stdio_flag=0
  integer              :: lsta, lnet, lcha, llocid, ldatime, lunits, lfile, lrtype, lverbose
  integer              :: useTotalSensitivityFlag

  character*20         :: datime2,datime3

  integer              :: iflag, evresp
  character*10         :: verbose
  real(kind=4), allocatable    :: resp(:), freq(:)
  
  double precision :: denr

  !transfer
  real(kind=8), allocatable :: s_re(:), s_im(:)
  real(kind=8), allocatable :: x_re(:), x_im(:)
  real(kind=8)              :: f1, f2
  real(kind=8)              :: f3, f4
  real(kind=8), dimension(4):: F
  integer                   :: nfreq
  double precision          :: delfrq

  !filter
  real(kind=8)         :: unused=0.0, unused2=0.0
  real*8  :: dt_temp

  start_stage=-1
  stop_stage = 0
  stdio_flag=0
  if (npts_obsd .ge. 524288) then
    LNPT = 20
    FFT_NPTS = 2**LNPT
  elseif ( npts_obsd .ge. 262154) then
    LNPT = 19
    FFT_NPTS = 2**LNPT
  else
    LNPT = 18
    FFT_NPTS = 2**LNPT
  endif
  allocate(freq(FFT_NPTS))
  allocate(s_re(FFT_NPTS))
  allocate(s_im(FFT_NPTS))
  allocate(x_re(FFT_NPTS))
  allocate(x_im(FFT_NPTS))
  allocate(resp(FFT_NPTS*2))
  !allocate(observed_cmplx(FFT_NPTS))
  !allocate(synthetic_cmplx(FFT_NPTS))

  !observed_dt = observed_raw%sample_rate(irecord)

  !+++++++++++++++++++++++++++++++++++
  !evresp
  do i = 1, FFT_NPTS
    freq(i) = real(i-1)/FFT_NPTS/dt_obsd
  enddo
  freq(1) = eps

  !call create_temp_resp_file(sta, cha, net, locid, response, resp_length,&
  !        temp_resp_fn, ierr)
  !------------------------------------------------
  ! DECONVOLUTION: EVALRESP to compute instrument response
  ! FIX 32-64 bit conflict between SAC and EVALRESP: USE 64-bit SAC
  !-------------------------------------------------
  !length = observed_raw%responses(irecord)%response_length
  write(datime2,*) gmt_year
  write(datime3,*) gmt_day+1 ! need to make sure the event is afer this date
  datime = trim(adjustl(datime2))//","//trim(adjustl(datime3))
  units='DIS'
  verbose='-v'
  rtype = 'CS'
  lsta=len_trim(sta); lcha=len_trim(cha); lnet=len_trim(net)
  llocid=len_trim(locid)
  ldatime=len_trim(datime)
  lunits=len_trim(units)
  lfile=len_trim(temp_resp_fn)
  !lfile=16
  lrtype=len_trim(rtype)
  lverbose=len_trim(verbose)
  !print *,"response: ", trim(response) 
  !print *,"end response"
  !temp_resp_fn=adjustl(temp_resp_fn)
  !print *, "enter..."
  !print *, "sta:", trim(sta), trim(cha), trim(net), trim(locid), ' || ',trim(datime)
  !print *, "resp_file:", trim(temp_resp_fn), "|||"
  !print *, lsta, lcha, lnet, llocid, ldatime
  !print *, lunits, lfile, lrtype, lverbose
  !print *, "FFT_NPTS:", FFT_NPTS
  !print *, "stdio_flag:", stdio_flag

  nfreq=FFT_NPTS
  !********
  !test block
  !nfreq=91
  !do i=1,91
  !  freq(i)=0.01+0.001*(i-1)
  !enddo
print *, cha
  iflag = evresp(sta, cha, net, locid, datime, &
            units, response(1:resp_length), freq, nfreq, resp, &
            rtype, verbose, start_stage, stop_stage, &
            stdio_flag)!, lsta, lcha, lnet, llocid, &
  
  !iflag = evresp(sta, cha, net, locid, datime, &
  !          units, temp_resp_fn, freq, nfreq, resp, &
  !          rtype, verbose, start_stage, stop_stage, &
  !          stdio_flag)!, lsta, lcha, lnet, llocid, &
            !ldatime, lunits, lfile, lrtype, lverbose, useTotalSensitivityFlag)

  !call delete_temp_resp_file(temp_resp_fn)
  !write out spectrum
  !open(unit=111, file="spec.log")
  !do i=1,nfreq
  !  write(111,*) freq(i), resp(2*i-1), resp(2*i)
  !enddo
  !close(111)

  print *, "checkpoint: ",trim(sta), trim(cha),trim(net), trim(locid)
  !deallocate(response)

  if (iflag.ne.0) then
    write(*,*) 'ERROR computing instrument response:'
    do i = 1, FFT_NPTS
      resp(i)=cmplx(1.0,0.0)
    enddo
  endif
  j = 1
  do i = 1, FFT_NPTS
    s_re(i) = dble(resp(j))
    s_im(i) = dble(resp(j+1))
    j = j + 2
  enddo

  ! Compute transfer function
  ! 1.0/ 'FROM' FROM = Sre + iSim
  ! Set to zero if less than FLT_MIN
  do i=1, nfreq
    denr=(s_re(i)**2+s_im(i)**2)
    if(denr.lt.FLT_MIN)then
      s_re(i)=0.0
      s_im(i)=0.0
    else
      denr=1.0/denr
      s_re(i)=s_re(i)*denr
      s_im(i)=-s_im(i)*denr
    endif
  enddo

  !+++++++++++++++++++++
  !transfer
  f3 = 1.0/MIN_PERIOD
  f2 = 1.0/MAX_PERIOD
  f1 = f2 * 0.8
  f4 = f3 * 1.2
  F = (/ f1, f2, f3, f4 /)
  x_re(1:FFT_NPTS)=1.0
  x_im(:)=0.0
  delfrq = 1.0/(dble(FFT_NPTS)*dble(dt_obsd))
  !print *, npts_obsd
  !open(unit=111, file="xre1.log")
  !do i=1,FFT_NPTS
  !  write(111,*) x_re(i), x_im(i), s_re(i), s_im(i)
  !enddo
  !close(111)
  nfreq = FFT_NPTS/2+1
  !call ztransfer(obsd(1:npts_obsd), npts_obsd, dble(dt_obsd), s_re, s_im, &
  !        x_re, x_im, nfreq, FFT_NPTS, delfrq, F)
  !open(unit=111, file="xre2.log")
  !do i=1,FFT_NPTS
  !  write(111,*) x_re(i), x_im(i), s_re(i), s_im(i)
  !enddo
  !close(111)
  

  print *, "transfer done"
  !transfer from nanometer to meter
  obsd(1:npts_obsd) = s_re(1:npts_obsd)!*1e-9

  !++++++++++++++++++++++
  !filter data
  dt_temp=dt_obsd
  !call xapiir(obsd, npts_obsd, 'BU',unused, unused, 4, 'BP', &
  !        f2, f3, dt_temp, 2)

  !print *, "remove instrument response done!"

end subroutine remove_response

subroutine filter(trace, npts, dt)

  use var_main

  real(kind=4), dimension(*) :: trace
  integer :: npts
  real(kind=4) :: dt

  double precision :: dt_temp
  double precision :: f1, f4,f3, f2, unused=0.0, unused2=0.0
  real(kind=8), allocatable :: s_re(:), s_im(:)
  real(kind=8), allocatable :: x_re(:), x_im(:)
  real(kind=8), dimension(4):: F
  integer                   :: nfreq
  double precision          :: delfrq

  integer FFT_NPTS, LNPT
  
  if (npts .ge. 524288) then
    LNPT = 20
    FFT_NPTS = 2**LNPT
  elseif ( npts .ge. 262154) then
    LNPT = 19
    FFT_NPTS = 2**LNPT
  else
    LNPT = 18
    FFT_NPTS = 2**LNPT
  endif

  allocate(s_re(FFT_NPTS))
  allocate(s_im(FFT_NPTS))
  allocate(x_re(FFT_NPTS))
  allocate(x_im(FFT_NPTS))
  f3 = 1.0/MIN_PERIOD
  f2 = 1.0/MAX_PERIOD
  dt_temp=dble(dt)

  !+++++++++++++++++++++
  !transfer
  f3 = 1.0/MIN_PERIOD
  f2 = 1.0/MAX_PERIOD
  f1 = f2 * 0.8
  f4 = f3 * 1.2
  F = (/ f1, f2, f3, f4 /)
  x_re(1:FFT_NPTS)=1.0
  x_im(:)=0.0
  s_re(1:FFT_NPTS)=1.0
  s_im(:)=0.0
  delfrq = 1.0/(dble(FFT_NPTS)*dble(dt))
  !print *, npts
  nfreq = FFT_NPTS/2+1
  !call ztransfer(trace, npts, dble(dt), s_re, s_im, &
  !        x_re, x_im, nfreq, FFT_NPTS, delfrq, F)
  trace(1:npts) = s_re(1:npts)!*1e-9


  !print *, "sum synt:", sum(trace(1:npts))
  !print *, "npts:", npts
  !print *, "dt_synt:", dt
  !print *, "f2, f3:", f2, f3

  !call xapiir(trace, npts, "BU", &
  !        unused, unused2, 4, 'BP', f2, f3, dt_temp, 2)

end subroutine filter

!Creat temporary resp file
subroutine create_temp_resp_file(sta, cha, net, locid, response, resp_length, &
               temp_resp_fn, ierr)

  character(len=*) :: sta, cha, net, locid, response, temp_resp_fn
  integer :: ierr, resp_length

  character(len=60) :: temp_dir=''

  temp_resp_fn=''

  temp_resp_fn=trim(temp_dir)//'RESP.'//trim(net)//'.'//trim(sta)//&
    '.'//trim(locid)//'.'//trim(cha)

  !print *, "temp_file: ",trim(temp_resp_fn)
  !print *, "len:",len_trim(temp_resp_fn)
  open(unit=111, file=temp_resp_fn, iostat=ierr)
  write(111, '(a)') response(1:resp_length)
  close(111)

end subroutine create_temp_resp_file

!Delete temporary resp file
subroutine delete_temp_resp_file(resp_fn)

  character(len=*) :: resp_fn

  call system('rm '//trim(resp_fn)//'')

end subroutine delete_temp_resp_file

end module process_subs2
