module UMOA
  use omp_lib
  use LinAlgLib
  use Operators
  use Iteration, only: Optimizer
  implicit none

  public :: UMOASolver

  private :: InitUMOASolver
  private :: release_hamiltonian
  private :: FinUMOASolver
  private :: SolveUMOA
  private :: get_initial_hamil
  private :: one_particle_one_hole_decoupling
  private :: two_particle_two_hole_decoupling
  private :: three_particle_three_hole_decoupling
  private :: get_auxiliary_fields
  private :: fields_inverse_transformation
  private :: transformation_one_body
  private :: transformation_two_body
  private :: transformation_three_body
  private :: print_iteration_spe
  private :: initialize_optimizer
  private :: iteration_method
  private :: get_exp

  type :: UMOASolver
    type(Ops) :: Sop, Xop, Wop, h0, h
    type(Ops) :: expS, expX
    logical :: is_valence = .false.
    logical :: is_three_body_S = .true.
    integer :: n_iter_max = 1000
    integer :: rank = 2
    real(8) :: alpha = 1.d0
    real(8) :: tol = 1.d-10
    real(8) :: e0   = 0.d0 ! cm term
    real(8) :: e1   = 0.d0 ! one-body term
    real(8) :: e2   = 0.d0 ! two-body term
    real(8) :: e3   = 0.d0 ! three-body term
    real(8) :: etot = 0.d0 ! umoa energy
    logical :: S1_is_0 = .false.
    logical :: S2_is_0 = .false.
    logical :: S3_is_0 = .false.
    logical :: X1_is_0 = .false.
    logical :: X2_is_0 = .false.
    logical :: X3_is_0 = .false.
  contains
    procedure :: InitUMOASolver
    procedure :: FinUMOASolver
    procedure :: SolveUMOA
    procedure :: TransformUMOA
    procedure :: release_hamiltonian

    generic :: init => InitUMOASolver
    generic :: fin => FinUMOASolver
    generic :: solve => SolveUMOA
  end type UMOASolver

  type(Optimizer), private :: opt
  integer, allocatable, private :: n2ch(:)
  integer, allocatable, private :: n2i1(:)
  integer, allocatable, private :: n2i2(:)
  character(:), allocatable, private :: file_spe
contains

  subroutine release_hamiltonian(this)
    class(UMOASolver), intent(inout) :: this
    call this%h%fin()
    call this%h0%fin()
    call this%Wop%fin()
  end subroutine release_hamiltonian

  subroutine FinUMOASolver(this)
    class(UMOASolver), intent(inout) :: this
    call this%h%fin()
    call this%h0%fin()
    call this%Sop%fin()
    call this%Xop%fin()
    call this%expS%fin()
    call this%expX%fin()
    call this%Wop%fin()
  end subroutine FinUMOASolver

  subroutine InitUMOASolver(this, h, is_valence, &
        & S1_is_0, S2_is_0, S3_is_0, X1_is_0, X2_is_0, X3_is_0)
    class(UMOASolver), intent(inout) :: this
    type(Ops), intent(in) :: h
    logical, intent(in), optional :: is_valence
    logical, intent(in), optional :: S1_is_0
    logical, intent(in), optional :: S2_is_0
    logical, intent(in), optional :: S3_is_0
    logical, intent(in), optional :: X1_is_0
    logical, intent(in), optional :: X2_is_0
    logical, intent(in), optional :: X3_is_0
    this%h0 = h
    this%h = h
    if(present(is_valence)) this%is_valence = is_valence
    if(present(S1_is_0)) this%S1_is_0 = S1_is_0
    if(present(S2_is_0)) this%S2_is_0 = S2_is_0
    if(present(S3_is_0)) this%S3_is_0 = S3_is_0
    if(present(X1_is_0)) this%X1_is_0 = X1_is_0
    if(present(X2_is_0)) this%X2_is_0 = X2_is_0
    if(present(X3_is_0)) this%X3_is_0 = X3_is_0
    call this%Sop%init(0, 1, 0, "S", h%ms, h%rank)
    call this%Xop%init(0, 1, 0, "X", h%ms, h%rank)
    call this%Wop%init(0, 1, 0, "W", h%ms, h%rank)
    call this%expS%init(0, 1, 0, "expS", h%ms, h%rank)
    call this%expX%init(0, 1, 0, "expX", h%ms, h%rank)
    this%rank = h%rank
    this%h%is_normal_ordered = .True.
    file_spe = "spe_iteration.dat"
  end subroutine InitUMOASolver

  subroutine SolveUMOA(this)
    class(UMOASolver), intent(inout) :: this
    integer :: ite
    integer :: unit_spe=222

    if(this%h0%is_normal_ordered) then
      write(*,"(a)") "Warning: Inpur Hamiltonian of UMOA iteration"
      return
    end if

    call get_initial_hamil(this)
    call initialize_optimizer(this)

    write(*,"(a)") "## UMOA iteration starts "
    write(*,'(2x,a,9x,a,10x,a,9x,a,7x,a,12x,a,12x,a)') "iter", "zero-body", &
        & "one-body ", "two-body ", &
        & "three-body ", "total", " error"
    write(*,'(4x,a,5f18.6)') "  ", this%e0, this%e1, &
        &  this%e2, this%e3,this%etot

    open(unit_spe, file=file_spe, status="replace", action="write")
    do ite = 1, this%n_iter_max
      call print_iteration_spe(this, ite, unit_spe)
      call one_particle_one_hole_decoupling(this)
      call two_particle_two_hole_decoupling(this)
      if(this%rank == 3 .and. this%h%ms%is_three_body) then
        call three_particle_three_hole_decoupling(this)
      end if
      call get_auxiliary_fields(this)
      write(*,'(2x,i4,5f18.6,es18.6)') ite, this%e0, this%e1, &
          &  this%e2, this%e3, this%etot, opt%r

      if(this%tol > opt%r) exit
      if(ite == this%n_iter_max) then
        write(*,'(a,i5,a)') "UMOA iteration does not converge after ", &
            & this%n_iter_max, " iterations"
        write(*,'(es14.6)') opt%r
      end if
      call fields_inverse_transformation(this)
      call iteration_method(this, ite)
    end do
    this%h%zero = this%etot
    if(ite <  this%n_iter_max) write(*,"(a)") "## UMOA iteration converges."
    if(ite >= this%n_iter_max) write(*,"(a)") "## UMOA iteration does not converge!"
    call opt%fin()
    close(unit_spe)
  end subroutine SolveUMOA

  subroutine get_initial_hamil(this)
    class(UMOASolver), intent(inout) :: this
    type(OneBodyPart) :: w1from2, w1from3
    type(TwoBodyPart) :: w2from3

    if(this%rank == 3 .and. this%h%ms%is_three_body) then
      w2from3 = this%h%thr%NormalOrderingFrom3To2(this%h%ms%two)
      w1from3 = w2from3%NormalOrderingFrom2To1(this%h%ms%one)
    end if
    w1from2  = this%h%two%NormalOrderingFrom2To1(this%h%ms%one)
    this%Wop%one = w1from2
    if(this%rank == 3 .and. this%h%ms%is_three_body) then
      this%Wop%two = w2from3
      this%Wop%one = this%Wop%one - w1from3 * 0.5d0
    end if

    this%e0 = this%h0%zero
    this%e1 = this%h0%one%NormalOrderingFrom1To0() + &
        & this%Wop%one%NormalOrderingFrom1To0()
    this%e2 = - w1from2%NormalOrderingFrom1To0() * 0.5d0
    if(this%rank == 3 .and. this%h%ms%is_three_body) then
      this%e3 = w1from3%NormalOrderingFrom1To0() / 6.d0
    end if
    this%etot = this%e0 + this%e1 + this%e2 + this%e3
    this%h%zero = this%etot
  end subroutine get_initial_hamil

  subroutine one_particle_one_hole_decoupling(this)
    use LeeSuzuki, only: LeeSuzukiSolver
    class(UMOASolver), target, intent(inout) :: this
    type(Ops), pointer :: h0, h, Sop, Xop, Wop, expS, expX
    type(OneBodyChannel), pointer :: obs
    type(DMat) :: Hdec
    type(LeeSuzukiSolver) :: LS
    integer :: ch, nex, np, nq
    !real(8) :: ti
    !ti = omp_get_wtime()

    h0 => this%h0
    h  => this%h
    Sop=> this%Sop
    Xop=> this%Xop
    expS=> this%expS
    expX=> this%expX
    Wop=> this%Wop

    do ch = 1, h0%ms%one%NChan
      obs => h0%ms%one%jpz(ch)
      Hdec = h0%one%MatCh(ch,ch)%DMat + &
          & Wop%one%MatCh(ch,ch)%DMat

      if(.not. this%S1_is_0) then
        nex = 0
        np  = obs%n_h_state
        nq  = obs%n_v_state + obs%n_p_state
        call LS%init(Hdec, is_taylor=.false.)
        call LS%solve(nex, np, nq)
        Sop%one%MatCh(ch,ch)%DMat = LS%S
        expS%one%MatCh(ch,ch)%DMat = LS%expS
        HDec = LS%H
        call LS%fin()
      end if

      if(.not. this%X1_is_0 .and. this%is_valence) then
        nex = obs%n_h_state
        np  = obs%n_v_state
        nq  = obs%n_p_state
        call LS%init(HDec, is_taylor=.false.)
        call LS%solve(nex, np, nq)
        Xop%one%MatCh(ch,ch)%DMat = LS%S
        expX%one%MatCh(ch,ch)%DMat = LS%expS
        HDec = LS%H
        call LS%fin()
      end if
      h%one%MatCh(ch,ch)%DMat = HDec
    end do
    !call timer%add("1p1h decoupling", omp_get_wtime() - ti)
  end subroutine one_particle_one_hole_decoupling

  subroutine two_particle_two_hole_decoupling(this)
    use LeeSuzuki, only: LeeSuzukiSolver
    class(UMOASolver), target, intent(inout) :: this
    type(Ops), pointer :: h0, h, Sop, Xop, Wop, expS, expX
    type(LeeSuzukiSolver) :: LS
    type(TwoBodyPartChannel) :: h1, expsx
    type(DMat) :: Hdec
    type(TwoBodyChannel), pointer :: tbs
    integer :: ch, nex, np, nq
    real(8) :: ti
    ti = omp_get_wtime()

    h0 => this%h0
    h  => this%h
    Sop=> this%Sop
    Xop=> this%Xop
    expS=> this%expS
    expX=> this%expX
    Wop=> this%Wop

    do ch = 1, h%ms%two%NChan
      tbs => h0%ms%two%jpz(ch)

      call h1%init(tbs,tbs); call h1%set(h0%one+Wop%one)
      HDec = h0%two%MatCh(ch,ch)%DMat + Wop%two%MatCh(ch,ch)%DMat + h1%DMat

      call expsx%init(tbs,tbs)
      if(.not. this%S1_is_0) then
        call expsx%set(expS%one, "*")
        HDec = expsx%DMat%T() * HDec * expsx%DMat
      end if

      if(.not. this%X1_is_0 .and. this%is_valence) then
        call expsx%set(expX%one, "*")
        HDec = expsx%DMat%T() * HDec * expsx%DMat
      end if
      call expsx%fin()

      if(.not. this%S2_is_0) then
        nex = tbs%n_hp_state
        np  = tbs%n_hh_state + tbs%n_hv_state
        nq  = tbs%n_vv_state + tbs%n_vp_state + tbs%n_pp_state
        call LS%init(Hdec, is_taylor=.false.)
        call LS%solve(nex, np, nq)
        Sop%two%MatCh(ch,ch)%DMat = LS%S
        expS%two%MatCh(ch,ch)%DMat = LS%expS
        HDec = LS%H
        call LS%fin()
      end if

      if(.not. this%X2_is_0 .and. this%is_valence) then
        nex = tbs%n_hp_state + tbs%n_hv_state + tbs%n_hh_state
        np  = tbs%n_vv_state
        nq  = tbs%n_vp_state + tbs%n_pp_state
        call LS%init(HDec, is_taylor=.false.)
        call LS%solve(nex, np, nq)
        Xop%two%MatCh(ch,ch)%DMat = LS%S
        expX%two%MatCh(ch,ch)%DMat = LS%expS
        HDec = LS%H
        call LS%fin()
      end if
      call h1%set(h%one)
      h%two%MatCh(ch,ch)%DMat = HDec - h1%DMat
      call h1%fin()
    end do
    !call timer%add("2p2h decoupling", omp_get_wtime() - ti)
  end subroutine two_particle_two_hole_decoupling

  subroutine three_particle_three_hole_decoupling(this)
    use LeeSuzuki, only: LeeSuzukiSolver
    class(UMOASolver), target, intent(inout) :: this
    type(Ops), pointer :: h0, h, Sop, Xop, Wop
    type(LeeSuzukiSolver) :: LS
    type(ThreeBodyPartChannel) :: h1, h2, sx
    type(DMat) :: expsx, Hdec
    type(ThreeBodyChannel), pointer :: tbs
    integer :: ch, nex, np, nq
    !real(8) :: ti
    !ti = omp_get_wtime()
    h0 => this%h0
    h  => this%h
    Sop=> this%Sop
    Xop=> this%Xop
    Wop=> this%Wop

    do ch = 1, h%ms%thr%NChan
      tbs => h0%ms%thr%jpz(ch)

      call h1%init(tbs,tbs); call h1%set(h0%one+Wop%one)
      call h2%init(tbs,tbs); call h2%set(h0%two+Wop%two)
      Hdec = h1%DMat + h2%DMat + h%thr%MatCh(ch,ch)%DMat

      call sx%init(tbs,tbs)
      if(.not. this%S1_is_0) then
        call sx%set(Sop%one)
        expsx = get_exp(sx%DMat)
        Hdec = expsx%T() * Hdec * expsx
      end if

      if(.not. this%X1_is_0 .and. this%is_valence) then
        call sx%set(Xop%one)
        expsx = get_exp(sx%DMat)
        Hdec = expsx%T() * Hdec * expsx
      end if

      if(.not. this%S2_is_0) then
        call sx%set(Sop%two)
        expsx = get_exp(sx%DMat)
        Hdec = expsx%T() * Hdec * expsx
      end if

      if(.not. this%X2_is_0 .and. this%is_valence) then
        call sx%set(Xop%two)
        expsx = get_exp(sx%DMat)
        Hdec = expsx%T() * Hdec * expsx
      end if

      if(.not. this%S3_is_0) then
        nex = tbs%n_hhp_state + tbs%n_hpv_state + tbs%n_hpp_state
        np  = tbs%n_hhh_state + tbs%n_hhv_state + tbs%n_hvv_state
        nq  = tbs%n_vvv_state + tbs%n_pvv_state + tbs%n_ppv_state + tbs%n_ppp_state
        call LS%init(HDec, is_taylor=.false.)
        call LS%solve(nex, np, nq)
        Sop%thr%MatCh(ch,ch)%DMat = LS%S
        HDec = LS%H
        call LS%fin()
      end if

      if(.not. this%X3_is_0 .and. this%is_valence) then
        nex = tbs%n_hhp_state + tbs%n_hpv_state + tbs%n_hpp_state + &
            & tbs%n_hhh_state + tbs%n_hhv_state + tbs%n_hvv_state
        np  = tbs%n_vvv_state
        nq  = tbs%n_pvv_state + tbs%n_ppv_state + tbs%n_ppp_state
        call LS%init(HDec, is_taylor=.false.)
        call LS%solve(nex, np, nq)
        Xop%thr%MatCh(ch,ch)%DMat = LS%S
        HDec = LS%H
        call LS%fin()
      end if
      call h1%set(h%one); call h2%set(h%two)
      h%thr%MatCh(ch,ch)%DMat = HDec - h2%DMat - h1%DMat
      call h1%fin(); call h2%fin()
    end do
    !call timer%add("3p3h decoupling", omp_get_wtime() - ti)
  end subroutine three_particle_three_hole_decoupling

  subroutine get_auxiliary_fields(this)
    class(UMOASolver), target, intent(inout) :: this
    type(Ops), pointer :: Wop, h
    type(OneBodyPart) :: w1from3

    Wop=> this%Wop
    h  => this%h
    Wop%one = h%two%NormalOrderingFrom2To1(h%ms%one)
    this%e2 = - Wop%one%NormalOrderingFrom1To0() * 0.5d0
    call Wop%two%fin()
    call Wop%two%init(h%ms%two,.true.,"W",0,1,0)
    if(this%rank == 3 .and. this%h%ms%is_three_body) then
      Wop%two = h%thr%NormalOrderingFrom3To2(h%ms%two)
      w1from3 = Wop%two%NormalOrderingFrom2To1(h%ms%one)
      Wop%one = Wop%one - w1from3 * 0.5d0
      this%e3 = w1from3%NormalOrderingFrom1To0() / 6.d0
      call w1from3%fin()
    end if

    this%e1 = h%one%NormalOrderingFrom1To0()
    this%etot = this%e0 + this%e1 + this%e2 + this%e3
  end subroutine get_auxiliary_fields

  subroutine fields_inverse_transformation(this)
    class(UMOASolver), target, intent(inout) :: this
    type(Ops), pointer :: Sop, Xop, Wop, expX, expS
    type(OneBodyPart) :: w1_til
    type(OneBodyChannel), pointer :: obs
    type(TwoBodyChannel), pointer :: tbs
    type(DMat) :: expsx, wtr
    type(TwoBodyPartChannel) :: w1_, expsx1
    integer :: ch
    Sop=> this%Sop
    Xop=> this%Xop
    expS=> this%expS
    expX=> this%expX
    Wop=> this%Wop

    w1_til = Wop%one
    do ch = 1, Wop%ms%one%NChan
      obs => Wop%ms%one%jpz(ch)
      if(.not. this%X1_is_0 .and. this%is_valence) then
        expsx = expX%one%MatCh(ch,ch)%DMat
        Wop%one%MatCh(ch,ch)%DMat = expsx * Wop%one%MatCh(ch,ch)%DMat * expsx%T()
      end if

      if(.not. this%S1_is_0) then
        expsx = expS%one%MatCh(ch,ch)%DMat
        Wop%one%MatCh(ch,ch)%DMat = expsx * Wop%one%MatCh(ch,ch)%DMat * expsx%T()
      end if
    end do

    !do ch = 1, Wop%ms%two%NChan
    !  Wop%two%MatCh(ch,ch)%m = 0.d0
    !end do
    if(this%rank == 2) return
    if(.not. this%h%ms%is_three_body) return

    do ch = 1, Wop%ms%two%NChan
      tbs => Wop%ms%two%jpz(ch)

      call w1_%init(tbs,tbs); call w1_%set(w1_til)
      wtr = Wop%two%MatCh(ch,ch)%DMat + w1_%DMat

      if(.not. this%X2_is_0 .and. &
          & this%is_valence .and. tbs%n_vv_state * &
          & (tbs%n_vp_state + tbs%n_pp_state) > 0) then
        expsx = expX%two%MatCh(ch,ch)%DMat
        wtr = expsx * wtr * expsx%T()
      end if

      if(.not. this%S2_is_0 .and. &
          & (tbs%n_hh_state + tbs%n_hv_state) * &
          & (tbs%n_vv_state + tbs%n_vp_state + tbs%n_pp_state) > 0) then
        expsx = expS%two%MatCh(ch,ch)%DMat
        wtr = expsx * wtr * expsx%T()
      end if

      call expsx1%init(tbs,tbs)
      if(.not. this%X1_is_0 .and. this%is_valence) then
        call expsx1%set(expX%one, "*")
        wtr = expsx1%DMat * wtr * expsx1%DMat%T()
      end if

      if(.not. this%S1_is_0) then
        call expsx1%set(expS%one, "*")
        wtr = expsx1%DMat * wtr * expsx1%DMat%T()
      end if
      call expsx1%fin()

      call w1_%set(Wop%one)
      Wop%two%MatCh(ch,ch)%DMat = wtr - w1_%DMat
      call w1_%fin()
    end do
    call w1_til%fin()
  end subroutine fields_inverse_transformation

  subroutine TransformUMOA(this, op)
    class(UMOASolver), target, intent(inout) :: this
    type(Ops), intent(inout) :: op
    type(Ops) :: op_in

    if(op%is_normal_ordered) then
      write(*,"(a)") "Warning: Inpur operator of UMOA transformation"
      return
    end if
    op_in = op

    call transformation_one_body(this, op, op_in)
    call transformation_two_body(this, op, op_in)
    if(this%Sop%rank == 3) then
      call transformation_three_body(this, op, op_in)
    end if
  end subroutine TransformUMOA

  subroutine transformation_one_body(this, op, op_in)
    class(UMOASolver), target, intent(inout) :: this
    type(Ops), intent(inout) :: op
    type(Ops), intent(in) :: op_in
    type(Ops), pointer :: Sop, Xop
    type(OneBodySpace), pointer :: obs
    type(OneBodyChannel), pointer :: ch_one
    type(DMat) :: op1, s1
    integer :: ch

    Sop => this%Sop
    Xop => this%Xop
    obs => op_in%ms%one
    do ch = 1, obs%NChan
      ch_one => obs%jpz(ch)

      op1 = op_in%one%MatCh(ch,ch)%DMat
      s1 = Sop%one%MatCh(ch,ch)%DMat
      op%one%MatCh(ch,ch)%DMat = exp((-1.d0) * s1) * op1 * exp(s1)

      if(this%is_valence) then
        op1 = op_in%one%MatCh(ch,ch)%DMat
        s1 = Xop%one%MatCh(ch,ch)%DMat
        op%one%MatCh(ch,ch)%DMat = exp((-1.d0) * s1) * op1 * exp(s1)
      end if
    end do
  end subroutine transformation_one_body

  subroutine transformation_two_body(this, op, op_in)
    class(UMOASolver), target, intent(inout) :: this
    type(Ops), intent(inout) :: op
    type(Ops), intent(in) :: op_in
    type(Ops), pointer :: expS, expX
    type(TwoBodySpace), pointer :: tbs
    type(TwoBodyChannel), pointer :: ch_two
    type(TwoBodyPartChannel) :: op1, expsx
    type(DMat) :: op2, ex
    integer :: ch

    expS=> this%expS
    expX=> this%expX
    tbs => op_in%ms%two
    do ch = 1, tbs%NChan
      ch_two => tbs%jpz(ch)

      call op1%init(ch_two, ch_two); call op1%set(op_in%one)

      op2 = op_in%two%MatCh(ch,ch)%DMat + op1%DMat
      if(.not. this%S1_is_0) then
        call expsx%init(ch_two, ch_two)
        call expsx%set(expS%one, "*")
        op2 = expsx%DMat%T() * op2 * expsx%DMat
      end if

      if(.not. this%X1_is_0 .and. this%is_valence) then
        call expsx%set(expX%one, "*")
        op2 = expsx%DMat%T() * op2 * expsx%DMat
      end if
      call expsx%fin()

      if(.not. this%S2_is_0) then
        ex = expS%two%MatCh(ch,ch)%DMat
        op2 = ex%T() * op2 * ex
      end if

      if(this%is_valence) then
        ex = expX%two%MatCh(ch,ch)%DMat
        op2 = ex%T() * op2 * ex
      end if

      call op1%set(op%one)
      op%two%MatCh(ch,ch)%DMat = op2 - op1%DMat
    end do
  end subroutine transformation_two_body

  subroutine transformation_three_body(this, op, op_in)
    class(UMOASolver), target, intent(inout) :: this
    type(Ops), intent(inout) :: op
    type(Ops), intent(in) :: op_in
    type(Ops), pointer :: Sop, Xop
    type(ThreeBodySpace), pointer :: tbs
    type(ThreeBodyChannel), pointer :: ch_thr
    type(ThreeBodyPartChannel) :: s, op1, op2
    type(DMat) :: op3, exps
    integer :: ch

    Sop => this%Sop
    Xop => this%Xop
    tbs => op_in%ms%thr
    do ch = 1, tbs%NChan
      ch_thr => tbs%jpz(ch)

      call op1%init(ch_thr, ch_thr); call op1%set(op_in%one)
      call op2%init(ch_thr, ch_thr); call op2%set(op_in%two)
      op3 = op_in%thr%MatCh(ch,ch)%DMat + op1%DMat + op2%DMat

      call s%init(ch_thr, ch_thr); call s%set(Sop%one)
      exps = get_exp(s%DMat)
      op3 = exps%T() * op3 * exps

      if(this%is_valence) then
        call s%set(Xop%one)
        exps = get_exp(s%DMat)
        op3 = exps%T() * op3 * exps
      end if

      call s%set(Sop%two)
      exps = get_exp(s%DMat)
      op3 = exps%T() * op3 * exps

      if(this%is_valence) then
        call s%set(Xop%two)
        exps = get_exp(s%DMat)
        op3 = exps%T() * op3 * exps
      end if

      s%DMat = Sop%thr%MatCh(ch,ch)%DMat
      exps = get_exp(s%DMat)
      op3 = exps%T() * op3 * exps

      if(this%is_valence) then
        s%DMat = Xop%thr%MatCh(ch,ch)%DMat
        exps = get_exp(s%DMat)
        op3 = exps%T() * op3 * exps
      end if

      call op1%set(op%one)
      call op2%set(op%two)
      op%thr%MatCh(ch,ch)%DMat = op3 - op2%DMat - op1%DMat
    end do
  end subroutine transformation_three_body

  function get_exp(S, tol, show_process) result(R)
    type(DMat), intent(in) :: S
    real(8), optional, intent(in) :: tol
    logical, optional, intent(in) :: show_process
    type(DMat) :: R
    real(8) :: ti
    ti = omp_get_wtime()
    R = exp(S, ord=20, tol_in=tol, show_process=show_process)
    call timer%Add("Matrix Exponential", omp_get_wtime() - ti)
  end function get_exp

  subroutine initialize_optimizer(this)
    class(UMOASolver), intent(inout) :: this
    integer :: n = 0
    integer :: n_st, ch, bra, ket
    real(8), allocatable :: v(:)

    do ch = 1, this%h%ms%one%NChan
      n_st = this%h%ms%one%jpz(ch)%n_state
      n = n + n_st * (n_st + 1) / 2
    end do

    allocate(n2ch(n))
    allocate(n2i1(n))
    allocate(n2i2(n))
    allocate(v(n))

    n = 0
    do ch = 1, this%h%ms%one%NChan
      n_st = this%h%ms%one%jpz(ch)%n_state
      do bra = 1, n_st
        do ket = 1, bra
          n = n + 1
          n2ch(n) = ch
          n2i1(n) = bra
          n2i2(n) = ket
          v(n) = this%h%one%MatCh(ch,ch)%m(bra,ket)
        end do
      end do
    end do
    call opt%init(n=n, m=10, method="l-bfgs", a=0.7d0)
    !call opt%init(n=n, method="direct", a=0.7d0)
    call opt%set_init(v)
  end subroutine initialize_optimizer

  subroutine iteration_method(this, ite)
    class(UMOASolver), intent(inout) :: this
    integer, intent(in) :: ite
    integer :: n = 0, ch, bra, ket
    real(8), allocatable :: v(:)

    allocate(v(opt%n))
    do n = 1, opt%n
      ch = n2ch(n)
      bra= n2i1(n)
      ket= n2i2(n)
      v(n) = this%Wop%one%MatCh(ch,ch)%m(bra,ket)
    end do
    call opt%get_next(v, ite)
    do n = 1, opt%n
      ch = n2ch(n)
      bra= n2i1(n)
      ket= n2i2(n)
      this%Wop%one%MatCh(ch,ch)%m(bra,ket) = v(n)
      this%Wop%one%MatCh(ch,ch)%m(ket,bra) = v(n)
    end do
    deallocate(v)
  end subroutine iteration_method

  subroutine print_iteration_spe(this, ite, unt)
    class(UMOASolver), intent(inout) :: this
    integer, intent(in) :: ite, unt
    integer :: ch, n
    write(unt,"(i6)",advance="no") ite
    do ch = 1, this%h%ms%one%NChan
      do n = 1, this%h%ms%one%jpz(ch)%n_state
        write(unt,"(f16.6)",advance="no") this%h%one%MatCh(ch,ch)%m(n,n)
      end do
    end do
    write(unt,*)
  end subroutine print_iteration_spe

end module UMOA

!program test
!  use Profiler, only: timer
!  use UMOA
!  implicit none
!  type(MSpace) :: ms
!  type(Ops) :: h
!  type(UMOASolver) :: sol_umoa
!  character(:), allocatable :: fnn
!
!  call timer%init()
!  fnn = "none"
!  call ms%init("He4", 24.d0, 6, 12, is_three_body=.true.)
!  call h%init("hamil", ms, 2)
!  call h%set(fnn, "none", [15,15,15], [0,0,0,0])
!  call sol_umoa%init(h)
!  call sol_umoa%solve()
!  call sol_umoa%fin()
!  call timer%fin()
!end program test
