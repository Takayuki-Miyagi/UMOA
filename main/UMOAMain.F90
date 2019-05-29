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

  call set_model_space(ms_init, p)
  ! Hamiltonian -----
  n = 2
  if(ms%is_three_body) n = 3
  if(ms%is_three_body_jt) n = 3
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
    h = HF%BasisTransform(h_init, p%is_NO2B)
    if(p%is_NO2B) then
      call h%UnNormalOrdering2B()
    end if

  case default
    if(p%is_NO2B) then
      h = h_init
      call h%NO2BApprox()
      call h%UnNormalOrdering2B()
    end if
  end select
  call h_init%fin()

  call set_model_space(ms, p)
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
  subroutine set_model_space(ms, p)
    type(MSpace), intent(out) :: ms
    type(InputParameters), intent(in) :: p

    select case(p%int_3n_file)
    case('none', 'None', 'NONE')

      if(p%umoa_rank==2) then

        if(conffile == 'none') then
          call ms%init(Nucl=p%Nucl, Core=p%Core, valence_orbits=p%valence_list, &
              & hw=p%hw, emax=p%emax, e2max=p%e2max, lmax=p%lmax, beta=p%beta_cm)
        end if

        if(conffile /= 'none') then
          call ms%init(filename=conffile, hw=p%hw, emax=p%emax, e2max=p%e2max, lmax=p%lmax, beta=p%beta_cm)
        end if

      end if

      if(p%umoa_rank==3) then

        if(conffile == 'none') then
          call ms%init(Nucl=p%Nucl, Core=p%Core, valence_orbits=p%valence_list, &
              & hw=p%hw, emax=p%emax, e2max=p%e2max, e3max=p%e3max, lmax=p%lmax, &
              & beta=p%beta_cm, is_three_body=.true.)
        end if

        if(conffile /= 'none') then
          call ms%init(filename=conffile, hw=p%hw, emax=p%emax, e2max=p%e2max, &
              & e3max=p%e3max, lmax=p%lmax, beta=p%beta_cm, is_three_body=.true.)
        end if

      end if

    case default

      if(p%umoa_rank==2) then
        if(conffile == 'none') then
          call ms%init(Nucl=p%Nucl, Core=p%Core, valence_orbits=p%valence_list, &
              & hw=p%hw, emax=p%emax, e2max=p%e2max, e3max=p%e3max, lmax=p%lmax, &
              & beta=p%beta_cm, is_three_body_jt=.true.)
        end if

        if(conffile /= 'none') then
          call ms%init(filename=conffile, hw=p%hw, emax=p%emax, e2max=p%e2max, &
              & e3max=p%e3max, lmax=p%lmax, beta=p%beta_cm, is_three_body_jt=.true.)
        end if

      end if

      if(p%umoa_rank==3) then
        if(conffile == 'none') then
          call ms%init(Nucl=p%Nucl, Core=p%Core, valence_orbits=p%valence_list, &
              & hw=p%hw, emax=p%emax, e2max=p%e2max, e3max=p%e3max, lmax=p%lmax, &
              & beta=p%beta_cm, is_three_body_jt=.true., is_three_body=.true.)
        end if

        if(conffile /= 'none') then
          call ms%init(filename=conffile, hw=p%hw, emax=p%emax, e2max=p%e2max, &
              & e3max=p%e3max, lmax=p%lmax, beta=p%beta_cm, is_three_body_jt=.true., is_three_body=.true.)
        end if
      end if

    end select

  end subroutine set_model_space

end program UMOAMain
