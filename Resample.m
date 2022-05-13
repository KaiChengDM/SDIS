
function seeds = Resample(x_root, root, sig, k, d, nchain)

%% Resample for important directions

 ind = find((1-chi2cdf((root).^2,d)) == 0);  root1 = root;  root1(ind) = [];
 wnork       = (1-chi2cdf((sig(k-1)./sig(k).*root1).^2,d))./(1-chi2cdf(root1.^2,d));   % weight
 nsamlev     = length(wnork);                          % number of samples
 ind         = randsample(nsamlev,nchain,true,wnork);  % resampling of directions 
 root_accept = root1(ind,:);                            % accepted roots along important directions

%% Resample of radius

 for i = 1 : nchain   
    
   r          = sqrt(chi2rnd(d,1,10^3));                       % random radius 
   r_accepted = r(find(r > sig(k-1)./sig(k)*root_accept(i)));  % accepted radius 

   while length(r_accepted) == 0
      r = sqrt(chi2rnd(d,1,10^4));
      r_accepted = r(find( r > sig(k-1)./sig(k)*root_accept(i)));
   end
   
   num  = size(r_accepted,2);     % number of accepted radius
   ind1 = randsample(num,1);      % resampling of radius

   seeds(i,:) = r_accepted(ind1).*x_root(ind(i),:)./norm(x_root(ind(i),:));  % transform radius and direction into Cartesian coordinates 

 end

end



