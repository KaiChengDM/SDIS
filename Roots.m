function [model_run, x_root, root, fval] = Roots(x_failure,g,sigma)

%% Find the root on specific direction

[m, n] = size(x_failure);  % failure sample size

for i = 1 : m

  e = diag(ones(1,n));

  a(i,:) = x_failure(i,:)./norm(x_failure(i,:));  % direction defined by failure samples

  w = a(i,:)'-e(:,1);

  v = e-2.*w.*w'./(w'*w);          % define the rotatation matrix

  g0 = @(x)g(sigma*x*v');          % rotated performance function

  g1 = @(u)g0([u zeros(1,n-1)]);   % one-dimensional performance function

  options = optimset('TolX',10^-2,'MaxIter',300,'TolFun',10^-2,'Display','off');

  [root(i,:), fval, exitflag, output] = fsolve(g1,norm(x_failure(i,:)),options);  % find root 

 % [root(i,:), fval, exitflag, output] = fsolve(g1,0,options);  % find root 

  x_root(i,:) = (v*[root(i,:) zeros(1,n-1)]')';  % compute the input parameter on the limit state surface corresponding the root 

  model_run(i) = output.funcCount;  % computational cost 

end

% outlier = find(abs(fval)>10^-3);
% root(outlier) = [];
% x_root(outlier,:) = [];

end

