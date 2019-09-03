module held_suarez_test_mod

  use mesh_mod
  use parallel_mod
  use params_mod
  use data_mod

  implicit none

  private

  public held_suarez_test_set_initial_condition


contains
  subroutine held_suarez_test_set_initial_condition()

    integer i, j

    write(6, *) '[Notice]: Use held suarez test initial condition.'

    static%ghs(:,:) = 0.0

    state(1)%u(:,:) = 0.0

    call parallel_fill_halo(state(1)%u, all_halo=.true.)

    state(1)%v(:,:) = 0.0


    call parallel_fill_halo(state(1)%v, all_halo=.true.)

    state(1)%gd(:,:) = 10**3 * g 

    call parallel_fill_halo(state(1)%gd, all_halo=.true.)


  end subroutine held_suarez_test_set_initial_condition

end module held_suarez_test_mod