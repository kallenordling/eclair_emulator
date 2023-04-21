! train (and optionally optimize) a sparse GP from data files

program gp_in
  use m_util
  use m_gp
  use m_gp_sparse
  use m_gp_dense
  use m_gp_optim
  use m_cov_all
  use m_noise_all

  implicit none

!  type(SparseGP) :: gp
  class(BaseGP), allocatable :: gp
  
  ! number of training points, sparse points and dimension of the input
  integer :: n, nsparse, input_dimension
  ! unit number
  integer :: u
  ! loop counter
  integer :: i
  ! number of hyperparameters and noise parameters required
  integer :: ntheta, nnu
  
  !tmp variable
  real(dp), dimension(500) :: res
  real(dp), dimension(:), allocatable :: theta, nu, t, lbounds, ubounds
  real(dp), dimension(:,:), allocatable :: x
  integer, dimension(:), allocatable :: obs_type
  logical :: optimize

  integer :: optimize_max_iter
  real(dp) :: optimize_ftol

  character(len=max_name_len) :: covariance_function = 'LINSQEXP'
  character(len=max_name_len) :: noise_model_name    = 'VAL'

  class(cov_fn), allocatable :: cf
  class(noise_model), allocatable :: nm

  !namelist/DIMENSIONS/input_dimension, n, nsparse
  !namelist/MODEL/covariance_function, noise_model_name
  !namelist/HYPERPARAMETERS/nu, theta, lbounds, ubounds
  !namelist/CONTROL/optimize, optimize_max_iter, optimize_ftol
  
  input_dimension = 6
  n = 500
  !nsparse = 60
  
 

  call string_to_cov_fn(covariance_function, cf)
  call string_to_noise_model(noise_model_name, nm)

  nnu = nm%nparams_required(input_dimension)
  ntheta = cf%ntheta_required(input_dimension)

  allocate(real(dp) :: nu(nnu))
  allocate(real(dp) :: theta(ntheta))
  allocate(real(dp) :: lbounds(nnu + ntheta))
  allocate(real(dp) :: ubounds(nnu + ntheta))
  allocate(real(dp) :: t(n))
  allocate(real(dp) :: x(n,input_dimension))
  allocate(integer :: obs_type(n))
 
  nu = 0.001
  theta = (/ 0.9010,0.9650,0.6729,3.5576,4.7418,1.2722,4.0612 /)
  lbounds(1) = 0.001
  lbounds(2:) = 0.01
  ubounds(:) = 100.0
  ubounds(1) = 0.001

  optimize = .true.
  optimize_max_iter = 10000
  optimize_ftol = 1.0d-7


  open(newunit=u, file="./data_night.w/DATA")
  
  read (u,*) (x(i,:), obs_type(i), t(i), i=1,n)
  
  ! Transform design and response here

  t = standardize(t,n)
  res = standardize(x(:,1),n)
  x(:,1) = res
    res = standardize(x(:,2),n)
  x(:,2) = res
    res = standardize(x(:,3),n)
  x(:,3) = res
    res = standardize(x(:,4),n)
  x(:,4) = res
    res = standardize(x(:,5),n)
  x(:,5) = res
    res = standardize(x(:,6),n)
  x(:,6) = res


     allocate(gp, source=DenseGP(nu, theta, x, obs_type, t, cf, nm))

 
  if (optimize) then
     call log_lik_optim(nnu + ntheta, gp, lbounds, ubounds, optimize_max_iter, optimize_ftol)
  else
     gp%theta = (/ 0.9010,0.9650,0.6729,3.5576,4.7418,1.2722,4.0612  /)
     gp%nu = 0.001
     call gp%update_matrices
  end if
  
  print *, gp%nu,' and ', gp%theta

  call gp%write_out("emulators/out.gp.night.w_v2")

contains
	function mean(x,dmn) result(res)
		integer dmn
		real(dp) x(dmn)
		real(dp) :: res

		res = sum(x)/dmn
	end function mean

	function std(x,meanx,dmn) result(res)
		integer dmn
		real(dp) x(dmn)
		real(dp) :: meanx
		real(dp) :: res

		res = sqrt(sum((x - meanx)**2)/(dmn-1))
	end function std
	 
	function standardize(x,dmn) result(res)
		integer dmn
		real(dp) x(dmn)

		real(dp), dimension(dmn) :: res

		real(dp) :: meanx 
		real(dp) :: stdx


		meanx= mean(x,dmn)

		stdx = std(x,meanx,dmn)

		res = (x - meanx)/stdx
	end function standardize
	
	function logistic(x) result(res)
	real(dp), intent(in) :: x
	real(dp) :: res
	res = 0.9 * x + 0.05;
	res = log(res) - log(1-res)
	end function logistic
	
	function logistic_vector(x,dmn) result(res)
	integer dmn
	real(dp)  x(dmn)
	real(dp),dimension(dmn) :: res
	res = 0.9 * x + 0.05
	res = log(res) - log(1-res)
	end function logistic_vector
	
	function inv_logistic(x) result(res)
	real(dp), intent(in) :: x
	real(dp) :: res
	res =  ( (1.0 / (1.0 + exp(-x))) - 0.05) / 0.9
	end function inv_logistic

	function inv_logistic_vector(x,dmn) result(res)
	integer dmn
	real(dp), intent(in),dimension(dmn) :: x
	real(dp),dimension(dmn) :: res
	res = ( (1.0 / (1.0 + exp(-x))) - 0.05) / 0.9
	end function inv_logistic_vector

end program gp_in





