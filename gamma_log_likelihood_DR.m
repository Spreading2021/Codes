function y = gamma_log_likelihood(abc, t, dx)
   %% distribution parameters
   a = abc(1);
   b = abc(2);
   c = abc(3);
   % size of the data
   [nt,ns] = size(dx);
   % loop to calculate t^c for given value of c
   for j=1:ns
      for i = 2:nt
                  dtc(i,j) = t(i,j)^c - t(i-1,j)^c;
      end 
   end
   % evaluate pdf
   y = gampdf(dx(:),a*dtc(:),b); 
   % filter out NaN's, which could happen just for numerical reasons
   y = y(~isnan(y));
   % return -1 * log-likelihood (for minimization)
   y = -sum(log(y));
return
   
