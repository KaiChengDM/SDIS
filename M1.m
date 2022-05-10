clear; clc; format long

%% preparation

d  = 2;          % number of dimensions
g1 = @(x)3+0.1.*(x(:,1)-x(:,2)).^2-(x(:,1)+x(:,2))./2^0.5;
g2 = @(x)3+0.1.*(x(:,1)-x(:,2)).^2+(x(:,1)+x(:,2))./2^0.5;
g3 = @(x)(x(:,1)-x(:,2))+6./2^0.5;
g4 = @(x)(x(:,2)-x(:,1))+6./2^0.5;
g  = @(x)min([g1(x),g2(x),g3(x),g4(x)]')'+2;  % limit state function  

%% Sequential directional importance sampling

nf     = 100;  % importance directions per level 
len    = 5;    % length of each Markov chain 
sigma  = 3;    % initial sigma
tarCoV = 1.5;  % target coefficient of variation of important weight
num    = 10;   % number of runs

for i = 1 : num                                                           % repeated runs
   [pf(i), cov(i), n_cost(i), level(i)] = SDIS(g,nf,len,sigma,d,tarCoV);  % run SDIS algorithm
end

n_m  = mean(n_cost')        % mean of computational costs
pf_m = mean(pf')            % mean of failure probability
cv_m = mean(cov')           % mean of coefficient of variation
cv   = std(pf')./mean(pf')  % coefficient of variation of multiple runs