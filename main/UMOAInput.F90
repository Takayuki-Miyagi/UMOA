module UMOAInput
  implicit none

  public :: InputParameters

  private :: InitInputParameters
  private :: PrintInputParameters
  type :: InputParameters
    integer :: umoa_rank
    integer :: emax
    integer :: lmax
    integer :: e2max
    integer :: e3max
    real(8) :: hw
    real(8) :: beta_cm
    real(8) :: alpha
    character(256), allocatable :: Ops(:)
    character(:), allocatable :: Nucl
    character(:), allocatable :: Core
    character(:), allocatable :: valence_list
    character(:), allocatable :: basis
    ! two-body file
    character(256) :: int_nn_file
    character(256), allocatable :: files_nn(:)
    integer :: emax_nn
    integer :: e2max_nn
    integer :: lmax_nn
    ! three-body file
    character(:), allocatable :: int_3n_file
    character(256), allocatable :: files_3n(:)
    integer :: emax_3n
    integer :: e2max_3n
    integer :: e3max_3n
    integer :: lmax_3n
    ! output files
    character(:), allocatable :: out_dir
    character(:), allocatable :: summary_file

    !
    logical :: is_Op_out
    character(256) :: Op_file_format
    logical :: X1_is_0
    logical :: X2_is_0
    logical :: X3_is_0
    logical :: S1_is_0
    logical :: S2_is_0
    logical :: S3_is_0
    logical :: is_NO2B

  contains
    procedure :: init => InitInputParameters
    procedure :: PrintInputParameters
  end type InputParameters
contains
  subroutine InitInputParameters(this, inputfile)
    use ClassSys, only: sys
    class(InputParameters), intent(inout) :: this
    character(*), intent(in) :: inputfile
    integer :: umoa_rank = 2
    integer :: emax=6
    integer :: e2max=-1
    integer :: e3max=-1
    integer :: lmax=-1
    real(8) :: hw=20.d0
    real(8) :: beta_cm = 0.d0 ! lowson's cm parameter, H => H + beta_cm * (hw/A) * Hcm
    real(8) :: alpha=0.7d0
    character(20) :: Nucl='O16'
    character(20) :: Core=""
    character(512) :: valence_list = ""
    character(20) :: basis = "HO"
    character(256) :: int_nn_file
    character(256) :: int_3n_file

    character(1024) :: optrs="none"
    character(1024) :: files_nn="none"
    integer :: emax_nn=6
    integer :: e2max_nn=12
    integer :: lmax_nn=-1
    ! three-body file
    character(1024) :: files_3n="none"
    integer :: emax_3n=6
    integer :: e2max_3n=6
    integer :: e3max_3n=6
    integer :: lmax_3n=-1

    character(256) :: out_dir = '.'
    character(256) :: summary_file = 'summary.out'
    logical :: is_Op_out = .false.
    character(256) :: Op_file_format = "snt"
    logical :: X1_is_0 = .false.
    logical :: X2_is_0 = .false.
    logical :: X3_is_0 = .false.
    logical :: S1_is_0 = .false.
    logical :: S2_is_0 = .false.
    logical :: S3_is_0 = .false.
    logical :: is_NO2B = .true.

    type(sys) :: s
    integer :: io
    namelist /input/ umoa_rank, emax, e2max, e3max, lmax, hw, &
        & Nucl, int_nn_file, files_nn, emax_nn, optrs, &
        & e2max_nn, lmax_nn, int_3n_file, files_3n, &
        & emax_3n, e2max_3n, e3max_3n, lmax_3n, alpha, &
        & summary_file, is_Op_out, &
        & beta_cm, out_dir, basis, &
        & Op_file_format, S1_is_0, S2_is_0, S3_is_0, &
        & X1_is_0, X2_is_0, X3_is_0, Core, valence_list, &
        & is_NO2B

    open(118, file=inputfile, action='read', iostat=io)
    if(io /= 0) then
      write(*,'(2a)') 'File opening error: ', trim(inputfile)
      stop
    end if
    read(118,nml=input)
    close(118)

    this%umoa_rank = umoa_rank
    this%emax = emax
    this%e2max = e2max
    if(e2max == -1) this%e2max = 2*emax
    this%e3max = e3max
    if(e3max == -1) this%e3max = 3*emax
    this%hw = hw
    this%alpha = alpha
    this%Nucl = Nucl
    this%Core = Core
    this%valence_list = valence_list
    this%int_nn_file = int_nn_file
    this%int_3n_file = int_3n_file
    this%basis = basis

    this%emax_nn = emax_nn
    this%e2max_nn = e2max_nn

    this%emax_3n = emax_3n
    this%e2max_3n = e2max_3n
    this%e3max_3n = e3max_3n

    this%lmax = lmax
    this%lmax_nn = lmax_nn
    this%lmax_3n = lmax_3n

    this%out_dir = out_dir
    this%S1_is_0 = S1_is_0
    this%S2_is_0 = S2_is_0
    this%S3_is_0 = S3_is_0
    this%X1_is_0 = X1_is_0
    this%X2_is_0 = X2_is_0
    this%X3_is_0 = X3_is_0
    this%is_NO2B = is_NO2B
    if(this%out_dir == '') this%out_dir = '.'
    call s%mkdir(this%out_dir)
    this%summary_file = trim(this%out_dir) // "/" // trim(summary_file)

    this%is_Op_out = is_Op_out
    this%Op_file_format = Op_file_format
    this%beta_cm = beta_cm

    if(lmax == -1) this%lmax = emax
    if(lmax_nn == -1) this%lmax_nn = emax_nn
    if(lmax_3n == -1) this%lmax_3n = emax_3n

    call s%split(optrs, ',', this%Ops)
    call s%split(files_nn, ',', this%files_nn)
    call s%split(files_3n, ',', this%files_3n)

  end subroutine InitInputParameters

  subroutine PrintInputParameters(this,iunit)
    class(InputParameters), intent(in) :: this
    integer, intent(in), optional :: iunit
    integer :: iut = 6
    integer :: n

    if(present(iunit)) iut = iunit

    write(iut,'(a)') "######  Input parameters  ####  "
    write(iut,'(2a)') "#  Target nuclide is ", trim(this%Nucl)
    write(iut,'(a,i3,a,i3,a,i3,a,i3)') "#  Model space: emax =", &
        & this%emax, ", e2max =", this%e2max, &
        & ", e3max =", this%e3max, ", lmax =", this%lmax
    write(iut,'(a,f8.3,a)') "#  HO basis parameter hw = ", this%hw, " MeV"

    write(iut,'(a)') "#"
    write(iut,'(a)') "#  NN files:"
    write(iut,'(2a)') "#  ", trim(this%int_nn_file)
    do n = 1, size(this%files_nn)
      if( this%files_nn(n) == 'none' .or. this%files_nn(n) == '') cycle
      write(iut,'(2a)') "#  ", trim(this%files_nn(n))
    end do
    write(iut,'(a,i3,a,i3,a,i3)') "#  File boundaries are emax =",this%emax_nn, &
        & ", e2max =",this%e2max_nn, ", lmax =",this%lmax_nn

    write(iut,'(a)') "#"
    write(iut,'(a)') "#  3N files:"
    write(iut,'(2a)') "#  ", trim(this%int_3n_file)
    do n = 1, size(this%files_3n)
      if( this%files_3n(n) == 'none' .or. this%files_3n(n) == '') cycle
      write(iut,'(2a)') "#  ", trim(this%files_3n(n))
    end do
    write(iut,'(a,i3,a,i3,a,i3,a,i3)') "#  File boundaries are emax =",this%emax_3n, &
        & ", e2max =",this%e2max_3n, ", e3max =", this%e3max_3n, &
        & ", lmax =",this%lmax_3n
    if(this%beta_cm > 1.d-4) write(iut,'(a,f6.3)') "#  Lawson's beta parameter is ", this%beta_cm
    !do n = 1, size(this%Ops)
    !  if( this%Ops(n) == 'none' .or. this%Ops(n) == '') cycle
    !  write(iut,'(a,a)') "#  Operator is ", trim(this%Ops(n))
    !end do
  end subroutine PrintInputParameters
end module UMOAInput
