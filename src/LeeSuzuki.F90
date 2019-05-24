module LeeSuzuki
  use omp_lib
  use LinAlgLib
  use Profiler, only: timer
  implicit none

  public :: LeeSuzukiSolver
  private :: Inils, finls, solve
  type :: LeeSuzukiSolver
    ! H: transformed Hamiltonian
    type(DMat) :: S, H, expS
    logical :: is_taylor = .true.
  contains
    procedure :: init => inils
    procedure :: fin => finls
    procedure :: solve
  end type LeeSuzukiSolver
contains
  subroutine inils(this, Hin, is_taylor)
    class(LeeSuzukiSolver), intent(inout) :: this
    type(DMat), intent(in) :: Hin
    logical, optional, intent(in) :: is_taylor
    integer :: n
    n = size(Hin%m, 1)
    if(present(is_taylor)) this%is_taylor = is_taylor
    call this%S%zeros(n,n)
    call this%expS%eye(n)
    this%H = Hin
  end subroutine inils

  subroutine finls(this)
    class(LeeSuzukiSolver), intent(inout) :: this
    call this%S%fin(); call this%H%fin(); call this%expS%fin()
  end subroutine finls

  subroutine solve(this, ndim_ex, npdim, nqdim)
    class(LeeSuzukiSolver), intent(inout) :: this
    integer, intent(in) :: ndim_ex, npdim, nqdim
    type(EigenSolSymD) :: sol
    type(DMat) :: hwork, omg
    integer(4) :: n, info
    real(8) :: ti

    ti = omp_get_wtime()
    n = ndim_ex + npdim + nqdim
    if(npdim < 1 .or. nqdim < 1) return
    n = npdim + nqdim
    call hwork%ini(n,n)
    hwork%m(:,:) = this%h%m(ndim_ex+1:, ndim_ex+1:)
    call sol%init(hwork)
    call sol%DiagSym(hwork)
    omg = Omega_LeeSuzuki(sol%vec, n, npdim, nqdim, info)

    call hwork%fin()
    if(info .ne. 0) return
    call sol%fin()
    this%S = Omega2S(ndim_ex, npdim, nqdim, omg)
    if(this%is_taylor) then
      this%expS = exp(this%S)
    end if

    if(.not. this%is_taylor) then
      this%expS = Omega2expS(ndim_ex, npdim, nqdim, omg)
    end if
    call omg%fin()
    this%H = this%expS%T() * this%H * this%expS

    call timer%add('Lee-Suzuki method', omp_get_wtime() - ti)
  contains
    type(DMat) function Omega_LeeSuzuki(wave, n, npdim, nqdim, info) result(omg)
      type(DMat), intent(in) :: wave
      integer(4), intent(in) :: n, npdim, nqdim
      integer(4), intent(out) :: info
      real(8), allocatable :: povlap(:)
      integer, allocatable :: lgt(:)
      real(8) :: ovlap
      type(DMat) :: a, ainv, ww_qp
      integer(4) :: i, j, k
      info = 0
      call omg%ini(nqdim, npdim)
      allocate(povlap(n), lgt(npdim))
      call a%Ini(npdim, npdim); call ww_qp%ini(nqdim, npdim)
      do i = 1, n
        ovlap = dot_product(wave%m(:npdim, i), wave%m(:npdim, i))
        povlap(i) = ovlap
      end do

      do i = 1, npdim
        k = 1
        do j = 1, n
          if(povlap(k) < povlap(j)) k = j
        end do
        lgt(i) = k
        povlap(k) = 0.d0
      end do

      do i = 1, npdim
        k = lgt(i)
        do j = 1, nqdim
          ww_qp%m(j,i) = wave%m(j+npdim, k)
        end do
      end do

      do i = 1, npdim
        do j = 1, npdim
          k = lgt(j)
          a%m(i,j) = wave%m(i,k)
        end do
      end do

      if(abs(a%Det()) < 1.d-8) then
        info = 1
        write(*,'(a, es18.6)') "in Omega_LeeSuzuki, linear dependence is not preserved: Det = ", a%Det()
        deallocate(povlap, lgt)
        call a%fin(); call ww_qp%fin()
        return
      end if
      ainv = a%Inv()
      omg = ww_qp * ainv
      deallocate(povlap, lgt)
    end function Omega_LeeSuzuki

    type(DMat) function Omega2S(nex, np, nq, omg) result(S)
      integer(4), intent(in) :: nex, np, nq
      type(DMat), intent(in) :: omg
      type(DMat) :: a, theta
      type(DVec) :: frac
      type(EigenSolSymD) :: sol
      integer(4) :: i, j

      a = omg%T() * omg
      call S%zeros(nex + np + nq, nex + np + nq)
      call sol%Init(a); call sol%DiagSym(a)
      call frac%ini(np)
      do i = 1, np
        if(abs(sol%eig%v(i)) < 1.d-8) then
          frac%v(i) = 1.d0
        else
          frac%v(i) = atan(sqrt(sol%eig%v(i))) / sqrt(sol%eig%v(i))
        end if
      end do
      call theta%DiagMat(frac)
      a = omg * sol%vec * theta * sol%vec%T()
      call sol%fin(); call frac%fin(); call theta%fin()
      do j = 1, np
        do i = 1, nq
          S%m(nex + np + i, nex + j) =   a%m(i, j)
          S%m(nex + j, nex + np + i) = - a%m(i, j)
        end do
      end do
    end function Omega2S

    function Omega2expS(nex, np, nq, omg) result(exS)
      integer(4), intent(in) :: nex, np, nq
      type(DMat), intent(in) :: omg
      type(DMat) :: a
      type(DMat) :: exS
      type(DVec) :: d
      type(DMat) :: pp, qp, qq, dd, qv
      integer :: i

      a = omg%T() * omg
      call exS%eye(nex+np+nq)
      call sol%Init(a); call sol%DiagSym(a)
      call d%ini(np)
      do i = 1, np
        d%v(i) = 1.d0 / sqrt(1.d0 + sol%eig%v(i))
      end do
      call dd%DiagMat(d)
      pp = sol%vec * dd * sol%vec%T()
      exS%m(nex+1:nex+np, nex+1:nex+np) = pp%m(:,:)

      pp = omg * sol%vec
      qp = pp * dd * sol%vec%T()
      exS%m(nex+np+1:, nex+1:nex+np) = qp%m(:,:)
      exS%m(nex+1:nex+np, nex+np+1:) = -transpose(qp%m)

      do i = 1, np
        if(abs(sol%eig%v(i)) < 1.d-8) then
          d%v(i) = - 0.5d0
        else
          d%v(i) = (1.d0/sqrt(1.d0 + sol%eig%v(i)) - 1.d0) / sol%eig%v(i)
        end if
      end do
      call dd%DiagMat(d)
      call qq%eye(nq)
      qq = qq + (pp * dd * pp%T())
      exS%m(nex+np+1:,nex+np+1:) = qq%m(:,:)

      call pp%fin()
      call qp%fin()
      call qq%fin()
      call qv%fin()
      call a%fin()
      call sol%fin()
      call dd%fin()
      call d%fin()
    end function Omega2expS
  end subroutine solve

end module LeeSuzuki

!program test
!  use Profiler, only: timer
!  use LinAlgLib
!  use LeeSuzuki
!  implicit none
!  type(DMat) :: H
!  type(LeeSuzukiSolver) :: LS
!  integer :: nex, np, nq, n
!  nex = 0
!  np = 4
!  nq = 10
!  n = nex + np + nq
!  call timer%init()
!  call H%Random(n,n)
!  H = H%T() + H
!  call H%prt("init")
!
!  call LS%init(H)
!  call LS%solve(nex,np,nq,H)
!  call LS%H%prt("LS")
!
!  call LS%fin()
!  call timer%fin()
!end program test
