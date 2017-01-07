function [f, g, H] = my_objective(x, A, b, c)
  f = dot(c, x) - sum(log(x * A' + b));
  g = c - sum(A ./ (x * A' + b)');
  H = A' * (A ./ (x * A' + b).^2');
  %g = c;
  %H = zeros(size(A, 2));
  %for i = 1:size(A, 1)
    %g = g - A(i, :) / (dot(A(i, :), x) + b(i));
    %H = H + A(i, :)' * A(i, :) / (dot(A(i, :), x) + b(i))^2;
  %end
end

