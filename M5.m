
clear;  clc;

%% preparation

d    = 1000;                  % number of dimensions
Beta = 4;                     % reliability index
g1   = @(x)Beta*sqrt(d)-sum(x');
g2   = @(x)Beta*sqrt(d)+sum(x');
g    = @(x)min(g1(x),g2(x));  % limit state function  
Pf = 2*normcdf(-Beta)         % true failure probability

%% Sequential directional importance sampling

nf     = 100;  % importance directions per level 
len    = 5;    % length of each Markov chain 
sigma  = 3;    % initial sigma
tarCoV = 1.5;  % target coefficient of variation of important weight
num    = 1;    % number of runs

for i = 1 : num                                                           % repeated runs
   [pf(i), cov(i), n_cost(i), level(i)] = SDIS(g,nf,len,sigma,d,tarCoV);  % run SDIS algorithm
end

n_m  = mean(n_cost')        % mean of computational costs
pf_m = mean(pf')            % mean of failure probability
cv_m = mean(cov')           % mean of coefficient of variation
cv   = std(pf')./mean(pf')  % coefficient of variation of multiple runs