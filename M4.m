clear all ; clc; format long;

%% preparation

W0 = @(x)sqrt((x(:,2)+x(:,3))./x(:,1));
g  = @(x)3.*x(:,4)-abs(2.*x(:,5)./(x(:,1).*W0(x).^2).*sin(W0(x).*x(:,6)./2)); % limit state function 

mu    = [1 1 0.1 0.5 0.3 1];          % mean of input variable
sigma = [0.05 0.1 0.01 0.05 0.2 0.2]; % standard deviation of input variable

g = @(x)g(x.*sigma+mu);               % transformed limit state function in standard normal space
d = 6;                                % input dimension

%% Sequential directional importance sampling

nf     = 100;  % importance directions per level 
len    = 5;    % length of each Markov chain 
sigma  = 3;    % initial sigma
tarCoV = 1.5;  % target coefficient of variation of important weight
num    = 100;  % number of runs

for i = 1 : num                                                           % repeated runs
   [pf(i), cov(i), n_cost(i), level(i)] = SDIS(g,nf,len,sigma,d,tarCoV);  % run SDIS algorithm
end

n_m  = mean(n_cost')        % mean of computational costs
pf_m = mean(pf')            % mean of failure probability
cv_m = mean(cov')           % mean of coefficient of variation
cv   = std(pf')./mean(pf')  % coefficient of variation of multiple runs