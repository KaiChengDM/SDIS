
clear;  clc;

%% preparation

dim = [10 100 300 500 1000];          % number of dimensions
l = [ 3 5 7 9];

for kk = 1 : 4
for k = 1 : length(dim)

 d    = dim(k);
 Beta = 4;
 g    = @(x)Beta*sqrt(d)-sum(x');  % Performance function
 P    = normcdf(-Beta)             % True failure probability

%% Sequential directional importance sampling

 nf     = 100;      % Importance directions per level 
 len    = l(kk);    % Length of each Markov chain 
 sigma  = 3;        % Initial sigma
 tarCoV = 1.5;      % target coefficient of variation of important weight

for i = 1 : 300
   [pf(i), cov(i), n_cost(i), level(i)] = SDIS(g,nf,len,sigma,d,tarCoV);
end

  N(k,kk) = mean(n_cost');
  Pf(k,kk) = mean(pf');
  CV1(k,kk) = mean(cov');
  CV2(k,kk) = std(pf')./mean(pf');

end
end
