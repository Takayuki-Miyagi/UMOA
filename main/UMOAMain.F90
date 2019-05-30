program UMOAMain
  use omp_lib
  use Profiler, only: timer
  use UMOAInput
  use ModelSpace
  use Operators
  use HartreeFock
  use HFMBPT
  use UMOA
  use WriteOperator
  implicit none
  type(InputParameters) :: p
  type(MSpace) :: ms_init, ms
  type(Ops) :: h_init, h, op
  type(HFSolver) :: HF
  type(MBPTDMat) :: PTd
  type(UMOASolver) :: sol_umoa
  type(WriteFiles) :: w
  logical :: is_valence=.false.
  character(256) :: inputfile='none', conffile='none'
  integer :: istatus, n

  call timer%init()

  select case(command_argument_count())
  case(0)
    write(*,'(a)') "This code needs input file!"
    stop
  case(1)
    call get_command_argument(1,inputfile,status=istatus)
    write(*,'(2a)') "Input file: ", trim(inputfile)
  case(2)
    call get_command_argument(1,inputfile,status=istatus)
    call get_command_argument(2,conffile,status=istatus)
    write(*,'(4a)') "Input files: ", trim(inputfile), ", ", trim(conffile)
  case default
    write(*,'(a)') "Too many arguments!"
    stop
  end select

  call p%init(inputfile)
  call p%PrintInputParameters()

  call set_model_space(ms_init, p, conffile)
  ! Hamiltonian -----
  n = 2
  if(ms_init%is_three_body) n = 3
  if(ms_init%is_three_body_jt) n = 3
  call h_init%init('hamil',ms_init, n)
  call h_init%set(p%int_nn_file,p%int_3n_file,&
      & [p%emax_nn,p%e2max_nn,p%lmax_nn],&
      & [p%emax_3n,p%e2max_3n,p%e3max_3n,p%lmax_3n])

  select case(p%basis)
  case("HF", "hf")
    call HF%init(h_init,alpha=p%alpha)
    call HF%solve()
    h = HF%BasisTransform(h_init, p%is_NO2B)
    if(p%is_NO2B) then
      call h%UnNormalOrdering2B()
    end if

  case("NAT", "nat", "Nat")
    call HF%init(h_init,alpha=p%alpha)
    call HF%solve()
    h = HF%BasisTransform(h_init, p%is_NO2B)
    call PTd%init(HF, h)
    HF%C = PTd%C_HO2NAT
    call PTd%fin()
    h = HF%BasisTransform(h_init, p%is_NO2B)
    if(p%is_NO2B) then
      call h%UnNormalOrdering2B()
    end if

  case default
    h = h_init
    if(p%is_NO2B) then
      call h%NO2BApprox()
      call h%UnNormalOrdering2B()
    end if
  end select
  call h_init%fin()

  call set_model_space(ms, p, conffile, mode="UMOA")
  h = h%Truncate(ms)
  if(ms%Nucl /= ms%Core) is_valence=.true.
  ! Hamiltonian -----
  call sol_umoa%init(h, is_valence=is_valence, &
      & S1_is_0 = p%S1_is_0, S2_is_0 = p%S2_is_0, &
      & S3_is_0 = p%S3_is_0, X1_is_0 = p%X1_is_0, &
      & X2_is_0 = p%X2_is_0, X3_is_0 = p%X3_is_0)
  call sol_umoa%solve()

  call w%init(p%emax, p%e2max, ms%beta)
  if(is_valence) then
    call sol_umoa%h%ReNormalOrdering2B()
    call w%SetFileName(p%out_dir, p%Op_file_format, p%basis, sol_umoa%h)
    call w%WriteValenceFile(p, sol_umoa%h)
  end if

  call sol_umoa%release_hamiltonian()
  if(p%Ops(1) /= 'none' .and. p%Ops(1) /= '') then
    write(*,"(a)") "## Calculation for expectation values"
  end if

  do n = 1, size(p%Ops)
    if(p%Ops(n) == "none" .or. p%Ops(n) == "") cycle
    call op%init(p%Ops(n), ms, 2)
    call op%set()
    call sol_umoa%TransformUMOA(op)
    call op%NormalOrdering()
    write(*,"(2a,f18.8)") trim(p%Ops(n)), ": ", op%zero
    call op%fin()
  end do
  call sol_umoa%fin()

  select case(p%basis)
  case("HF", "hf")
    call HF%fin()
  end select
  call h%fin()
  call ms_init%fin()
  call ms%fin()

  call timer%fin()

contains
  subroutine set_model_space(ms, p_in, f_configuration, mode)
    type(MSpace), intent(out) :: ms
    type(InputParameters), intent(in) :: p_in
    type(InputParameters) :: p_set
    character(*), intent(in) :: f_configuration
    character(*), intent(in), optional :: mode

    p_set = p_in
    if(present(mode)) then

      select case(mode)
      case("umoa", "UMOA")
        p_set%emax =  p_in%e_umoa
        p_set%e2max = p_in%e2_umoa
        p_set%e3max = p_in%e3_umoa
        p_set%lmax =  p_in%l_umoa
        p_set%int_3n_file =  "none"
      case default
      end select

    end if

    select case(p_set%int_3n_file)
    case('none', 'None', 'NONE')

      if(p_set%umoa_rank==2) then

        if(f_configuration == 'none') then
          call ms%init(Nucl=p_set%Nucl, Core=p_set%Core, valence_orbits=p_set%valence_list, &
              & hw=p_set%hw, emax=p_set%emax, e2max=p_set%e2max, lmax=p_set%lmax, beta=p_set%beta_cm)
        end if

        if(f_configuration /= 'none') then
          call ms%init(filename=conffile, hw=p_set%hw, emax=p_set%emax, e2max=p_set%e2max, lmax=p_set%lmax, beta=p_set%beta_cm)
        end if

      end if

      if(p_set%umoa_rank==3) then

        if(f_configuration == 'none') then
          call ms%init(Nucl=p_set%Nucl, Core=p_set%Core, valence_orbits=p_set%valence_list, &
              & hw=p_set%hw, emax=p_set%emax, e2max=p_set%e2max, e3max=p_set%e3max, lmax=p_set%lmax, &
              & beta=p_set%beta_cm, is_three_body=.true.)
        end if

        if(f_configuration /= 'none') then
          call ms%init(filename=conffile, hw=p_set%hw, emax=p_set%emax, e2max=p_set%e2max, &
              & e3max=p_set%e3max, lmax=p_set%lmax, beta=p_set%beta_cm, is_three_body=.true.)
        end if

      end if

    case default

      if(p_set%umoa_rank==2) then
        if(f_configuration == 'none') then
          call ms%init(Nucl=p_set%Nucl, Core=p_set%Core, valence_orbits=p_set%valence_list, &
              & hw=p_set%hw, emax=p_set%emax, e2max=p_set%e2max, e3max=p_set%e3max, lmax=p_set%lmax, &
              & beta=p_set%beta_cm, is_three_body_jt=.true.)
        end if

        if(f_configuration /= 'none') then
          call ms%init(filename=conffile, hw=p_set%hw, emax=p_set%emax, e2max=p_set%e2max, &
              & e3max=p_set%e3max, lmax=p_set%lmax, beta=p_set%beta_cm, is_three_body_jt=.true.)
        end if

      end if

      if(p_set%umoa_rank==3) then
        if(f_configuration == 'none') then
          call ms%init(Nucl=p_set%Nucl, Core=p_set%Core, valence_orbits=p_set%valence_list, &
              & hw=p_set%hw, emax=p_set%emax, e2max=p_set%e2max, e3max=p_set%e3max, lmax=p_set%lmax, &
              & beta=p_set%beta_cm, is_three_body_jt=.true., is_three_body=.true.)
        end if

        if(f_configuration /= 'none') then
          call ms%init(filename=conffile, hw=p_set%hw, emax=p_set%emax, e2max=p_set%e2max, &
              & e3max=p_set%e3max, lmax=p_set%lmax, beta=p_set%beta_cm, is_three_body_jt=.true., is_three_body=.true.)
        end if
      end if

    end select
  end subroutine set_model_space

end program UMOAMain
