function [pf, cov, n_cost, level] = SDIS(g,nf,len,sigma,d,tarCoV)

%% Sequential importance sampling algorithm for rare event estimation
%{
Created by: Kai Cheng (ckwow@mail.nwpu.edu.cn)
Based on: "Rare event estimation with sequential directional importance sampling" 
---------------------------------------------------------------------------
Input:
* nf     : important directions per level
* len    : length of each Markov chain
* g      : limit state function
* sigma  : initial simga
* d      : input dimension
* tarCoV : target coefficient of variation of the weights
---------------------------------------------------------------------------
Output:
* pf     : probability of failure
* cov    : coefficient of variation of SDIS estimator
* n_cost : total model evaluations
* level  : total number of intermediate levels
%}

%% Step 1 : Monte Carlo simulation 

% initilization
n                = 0;     % number of MCS population
n_failure        = 0;     % number of failure samples
model_evaluation = 0;     % number of model evaluations
sig              = sigma; % initial magnification factor of input standard deviation

mu = zeros(1,d); stdu = ones(1,d);  % mean and standard deviation of input variable

while n_failure < nf                % sequential enrichment of initial MCS population

   n      = n + 1;   
   x(n,:) = lhsnorm(mu,diag(stdu.^2),1);  % generate random samples with MCS
   y(n)   = g(sig.*x(n,:));               % evaluate the auxiliary limit state function 

   if y(n) < 0

     [model_run, x_f, root_f] = Roots(x(n,:),g,sig);   % find the roots of the failure samples
    
      model_evaluation = model_evaluation + model_run;  % update computational cost

      if (1-chi2cdf((root_f).^2,d)) ~= 0

         n_failure           = n_failure + 1;
         x_root(n_failure,:) = x_f;
         root(n_failure,:)   = root_f;

      end
   end
end

pf1       = n_failure/n;          % estimation of P_1 with MCS
cov1      = sqrt((1-pf1)/n/pf1);  % coefficient of variation of P_1
ind       = find( y < 0);         % find failure samples
x_failure = x(ind ,:);            % failure samples
model_evaluation = model_evaluation + n; % update computational cost

%% Sequential reduction of Sigma

tarWk = tarCoV;    % Target coefficient of variation of important weight
optimal_sig = sig;           
k = 1;

while optimal_sig > 1
     
  fun = @(x)abs(std((1-chi2cdf((sig(k)./x.*root).^2,d))./(1-chi2cdf((root).^2,d)))./...     % target function 
             mean((1-chi2cdf((sig(k)./x.*root).^2,d))./(1-chi2cdf((root).^2,d)))- tarWk);           

  options = optimoptions('fmincon','Display','off');            
  optimal_sig = fmincon(fun,1,[],[],[],[],1,sig(k),[],options);  % find optimal sigma


 if abs(optimal_sig-1) < 10^-3
    optimal_sig = 1;
 end

 k = k + 1;  sig(k) = optimal_sig;
  
 wk = (1-chi2cdf((sig(k-1)./optimal_sig.*root).^2,d))./(1-chi2cdf((root).^2,d));   % Compute the importance weight

 wk(find(isnan(wk) == 1)) = [];   
 sk(k-1) = mean(wk);                              % mean of importance weight 
 cv(k-1) = sqrt((std(wk)/sk(k-1))^2/length(wk));  % coefficient of variation of importance weight  

 if sig(k) == 1                                   % convergence criterion
     break;
 end

%%  Resampling 

  seeds = Resample(x_root, root, sig, k, d, nf);             % resampling 
 
  [x_mcmc, accepted_ratio] = MCMC(seeds, g, nf*len, sig(k)); % MCMC for sampling

  x_stable = x_mcmc(len:len:end,:);                          % select the desired stable samples
  
  model_evaluation = model_evaluation + nf*len;              % update computational cost 

  [model_run, x_new, root_new] = Roots(x_stable, g, sig(k)); % find the roots of the failure samples

  model_evaluation = model_evaluation + sum(model_run);      % update computational cost
 
  root      = root_new;
  x_failure = x_stable;
  x_root    = x_new;

 end

 pf     = pf1*prod(sk);               % failure probability 
 cov    = sqrt(cov1^2 + sum(cv.^2));  % coefficient of variation 
 n_cost = model_evaluation;           % total model evaluations
 level  = length(sk)+1;               % total number of intermediate levels

end




