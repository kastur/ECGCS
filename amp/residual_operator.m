function r = residual_operator(u)
  [size_x, size_y] = size(u);
  r = sum(u)/size_x;
end