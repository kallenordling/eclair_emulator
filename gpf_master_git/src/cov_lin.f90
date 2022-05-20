module m_cov_lin
  use m_util
  use m_cov
  
  implicit none
  
  private
  public cov_lin
  
  type, extends(cov_fn) :: cov_lin
   contains
     procedure, nopass :: ntheta_required
     procedure, nopass :: cov_val
     procedure, nopass :: dcov_x1
     procedure, nopass :: dcov_x2
     procedure, nopass :: d2cov_xx
  end type cov_lin

contains
  pure function ntheta_required(dims)
    integer ntheta_required
    integer, intent(in) :: dims
 ! no parameters to be estimated
    ntheta_required = 0
  end function ntheta_required

  pure function cov_val(x,y,hypers)
    real(dp) :: cov_val
    real(dp), intent(in), dimension(:) :: x, y, hypers
    cov_val = 10 + 10 * (sum(x * y)) ! why don't we use the crossing parameter here? It would need to be MAP estimated
  end function cov_val

  pure function dcov_x1(n,x,y,hypers)
    real(dp) :: dcov_x1
    real(dp), intent(in), dimension(:) :: x, y, hypers
    integer, intent(in) :: n
    dcov_x1 = 10 * x(n) * y(n)
  end function dcov_x1

  pure function dcov_x2(n,x,y,hypers)
    real(dp) :: dcov_x2
    real(dp), intent(in), dimension(:) :: x, y, hypers
    integer, intent(in) :: n
    dcov_x2 = 10 * x(n) * y(n)
  end function dcov_x2

  pure function d2cov_xx(n,m,x,y,hypers)
    real(dp) :: d2cov_xx
    real(dp), intent(in), dimension(:) :: x, y, hypers
    integer, intent(in) :: n, m
    if (n .eq. m) then
    d2cov_xx = 10
    else
    d2cov_xx = 0
    end if
  end function d2cov_xx
end module m_cov_lin
