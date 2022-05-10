
function  [x_mcmc, alpha_mu] = MCMC(seeds,g,ntotal,sigma)

%% MCMC sampling: conditional sampling M-H (CSM-H) algorithm

  [m, n]      = size(seeds);   % size of seeds
  seedslength = m;             % number of Markov chains 
  lenchain    = ntotal/m;      % length of each Markov chain
  beta        = 0.5;           % initial beta    
  alphak      = zeros(m,1);    % space for the standard deviation
  sigmaf      = 1;             % initial standard deviation
  counta      = 0; 
  lambda      = 0.6;
 
  ii = 0;

  for num = 1 : seedslength
  
      x_current = seeds(num,:);   % initial state
    
   for i = 1 : lenchain

     rhok = beta;
     x_candidate = normrnd(rhok*x_current', sqrt(1-rhok^2))';   % candidate state
     y_candidate = g(sigma.*x_candidate);    
   
     ii = ii + 1;
  
     if y_candidate <0                 % accept or reject
        x_current = x_candidate;
        alphak(num) = alphak(num)+1/lenchain;
     else
       x_current = x_current;
     end
     x_mcmc(ii,:) = x_current;
   end

   adapchains = 10;
   
   % check whether to adapt now
  if mod(num,adapchains) == 0 

     % mean acceptance rate of last adap_chains
     alpha_mu = mean(alphak(num-adapchains+1:num));
     counta   = counta + 1;
     gamma    = counta^(-0.5);
     lambda   = exp(log(lambda)+gamma*(alpha_mu-0.44));    
     % compute parameter rho
     sigmafk = min(lambda*sigmaf,1);
     beta   = sqrt(1-sigmafk.^2);

  end

end


