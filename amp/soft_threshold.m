function [x] = soft_threshold(u, lambda, sigma_t)
  [size_x, size_y] = size(u);
  thresh = lambda*sigma_t;
  thresh
  x = zeros(size(u));
  for i = 1:size_x,
    ui = u(i);
    if ui >= thresh,
      x(i) = ui - thresh;
    elseif ui <= -thresh,
      x(i) = ui + thresh;
    else
      x(i) = 0;
    end
  end
end