function main(m, n)
  A = randn(m, n);
  b = rand(1, m);
  c = rand(1, n);
  x = zeros(1, n);

  alpha = 0.1;
  beta = 0.7;
  epsilon = 1e-10;

  all_est_ep = [];
  all_f = [];
  while true
    [f, g, H] = my_objective(x, A, b, c);
    step = - g * inv(H); % Newton's step
    dec = sqrt(-dot(g, step)); % Newton's decrement
    est_ep = dec^2 / 2; % estimate of f(x) - p*
    all_est_ep = [all_est_ep est_ep];
    all_f = [all_f f];
    if est_ep < epsilon
      break
    end
    
    % backtracking line search
    t = 1;
    while true
      if sum((x + t * step) * A' + b <= 0) == 0 % make sure x + t*step is in dom f
        [f_new, g_new, H_new] = my_objective(x + t * step, A, b, c);
        if f_new <= f - alpha * t * dec^2
          x = x + t * step; % update x
          break
        end
      end
      t = t * beta;
    end
    %disp(x);
    %fflush(stdout);
  end
  figure;
  %subplot(2, 1, 1);
  %plot(all_f);
  %title('f(x)');
  
  %subplot(2, 1, 2);
  plot(all_est_ep);
  title('0.5\lambda(x)^{2}');
end
