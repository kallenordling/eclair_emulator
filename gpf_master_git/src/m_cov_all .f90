module m_cov_all
  use m_cov
  use m_cov_sqexp
  use m_cov_sqexp_param4
  use m_cov_lin
  use m_cov_linsqexp

  use m_util, only: max_name_len

contains

  function cov_fn_to_string(CovFunction) result (cov_fn_name)
    class(cov_fn), intent(in) :: CovFunction
    character(len=max_name_len) cov_fn_name

    select type(cf => CovFunction)
    type is (cov_sqexp) 
       cov_fn_name = 'SQEXP'
    type is (cov_sqexp_param4)
       cov_fn_name = 'SQEXP4PARAM'
    type is (cov_lin)
       cov_fn_name = 'LIN'
    type is (cov_linsqexp)
       cov_fn_name = 'LINSQEXP'
    class default
       cov_fn_name = 'UNKNOWN'
    end select
  end function cov_fn_to_string

  subroutine string_to_cov_fn(cov_fn_name, CovFunction)
    character(len=max_name_len), intent(in) :: cov_fn_name
    class(cov_fn), intent(out), allocatable :: CovFunction

    select case (cov_fn_name)
    case ('SQEXP')
       allocate(cov_sqexp :: CovFunction)
    case ('SQEXP4PARAM')
       allocate(cov_sqexp_param4 :: CovFunction)
    case ('LIN')
       allocate(cov_lin :: CovFunction)
    case ('LINSQEXP')
       allocate(cov_linsqexp :: CovFunction)
    case default
       print *, "unknown covariance function type, ", cov_fn_name
       stop 1
    end select
  end subroutine string_to_cov_fn

end module m_cov_all
