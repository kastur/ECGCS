function [x1,z1] = amp_iteration(A, x0, z0, sigma_t0, delta)
  z1 = y - A*x0 + (1/delta)*z0 * residual_function(soft_threshold_derivative(A*z0 + x0));
  x1 = soft_threshold(At*z0 + x, sigma_t0);
end