

clear all;
diary('ITS_v2_results_log.txt')
pkg load statistics

%EXAMPLE DATA
data0=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
17	1.869	0.984	0.829	1.132	0.314	0.487	0.775	0.563	1.409	0.96	0.768	0.74	0.884	3.296
42	3.733	2.729	2.527	2.974	2.056	2.355	2.585	2.424	3.258	3.016	2.797	2.789	2.826	4.747
62	5.589	4.349	4.158	4.682	3.795	4.084	4.374	4.28	5.029	4.82	4.429	4.498	4.544	5.788
82	6.426	5.46	5.332	5.498	4.868	4.963	5.365	5.416	6.258	5.889	5.197	5.238	5.698	6.767
102	7.84	6.504	6.428	6.559	6.307	6.274	6.682	6.82	7.585	7.318	6.297	6.466	6.015	7.519
122	8.097	7.065	6.927	6.955	6.573	6.574	7.16	7.287	8.118	7.836	6.653	6.756	6.844	8.043
142	9.921	8.451	8.34	8.309	8.355	8.174	8.701	9.054	9.752	9.486	8.148	8.222	9.409	9.343
162	11.728	10.036	9.951	9.829	10.047	9.883	10.473	10.82	11.565	11.348	9.662	9.779	10.863	10.043
202	14.266	12.27	12.037	11.954	12.086	12.18	12.31	13.119	13.416	13.566	11.645	11.879	13.01	11.575
262	18.106	16.344	15.842	15.496	16.082	16.135	15.777	17.416	17.159	16.686	14.766	15.609	16.448	15.414
442	28.009	26.45	25.747	24.436	26.948	28.438	25.637	28.471	27.855	27.457	25.488	24.227	27.228	27.444
622	33.106	32.705	32.965	32.167	34.842	36.213	32.004	34.82	38.154	37.626	36.084	35.518	38.354	42.277
722	45.675	40.124	39.443	42.628	46.483	47.587	43.081	45.756	42.78	43.528	40.338	39.457	43.609	48.1
782	50.109	43.647	43.352	44.962	52.053	51.173	46.241	48.573	46.034	47.28	42.511	43.559	45.901	50.706
942	56.811	50.327	50.74	53.232	56.543	57.167	51.041	54.949	51.914	53.13	48.101	48.819	52.479	57.143
1102	64.112	59.607	59.577	62.632	64.703	66.741	58.864	64.144	59.956	61.475	57.284	55.946	59.771	65.58
];

% for convenience
treal = data0(:,1);
cf_real = data0(:,2:end);

%Pooled NHGP parameters of 14 cells (obtained by Daniela from dataO - original data)
abc_mle_estimates_real_data_pooled = [0.3537 0.4015 0.8758];
a = 0.3537;
b = 0.4015;
c = 0.8758;
%--------------------------------------------------


% number of samples (number of degradation trajectories to simulate)
ns = 100
% time sequence limits (measurements made up to 1200 days)
tlim = [min(treal) max(treal)];
% number of time points
nt = numel(treal)
% points in time at which the capacity was measured
% this could be just a regular sampling, which is what an actual experiment
% would do, e.g, measuring every xxx days
iRandom = false;
if false
   % generate measurement times randomly
   t = randi(tlim,nt,ns);
   t = [zeros(1,ns) ; t];
   t = sort(t);
   dt = diff(t);
else
   % generate measurement times deterministically, equal to the measurement
   % times in the real data (data0)
   t = repmat(treal,1,ns);
   dt = diff(t);
end

% to store degradation increments
dx = zeros(nt,ns);
% limits to degradation increment (is it a percentage? What units does the Gamma have?)
dxlim = [0 100];
# initial capacity
cinit = 100.0;
% initial degradation
xinit = 0.0;

fprintf('\n*** Simulating %d degradation trajectories from the non-homogeneous Gamma process (NHGP)...',ns)
% inverse transformation to sample from Gamma CDF
for j=1:ns
   for i = 2:nt
      % uniform random number used to generate the sample from the Gamma CDF
      Ui=unifrnd(0,1);
      
      % calculate parameter of the NHGP
      dtc(i,j) = t(i,j)^c - t(i-1,j)^c;
      % use built-in Gamma CDF
      F=@(x)( gamcdf(x,a*dtc(i,j),b) - Ui);
      % this one should have worked, but it does not and I do not know why
      %F=@(x)( ((1./gamma(a.*dtc(i,j))) .* gammainc(a.*dtc(i,j), b.*x) ) - Ui);

      % sample from the Gamma CDF.
      try
         dx(i,j)=fzero(F,dxlim);
      catch
         % because sometimes fzero fails with an error (no zeros found!)
         % we will filter the NaNs later
         dx(i,j) = NaN;
      end
      
      % For debugging. Do NOT activate, it will create ns*nt plots!
      %figure
      %plot(0:100,F(0:100)+Ui)
      %saveas(gcf,sprintf('cdf_%f_%d_%d.jpg',a*dtc(i,j),i,j))
   end
end
fprintf('Done!\n')


% Capacity vs time matrix: initial capacity minus accumulated degradation
capacities = cinit * ones(nt,ns) - cumsum(dx);

% because sometimes the inverse transform fails, some elements of dx are NaNs
% so we filter those from the output
% idx of columns to remove
idx = find(all(not(isnan(capacities))));
% remove the columns from the data
capacities = capacities(:,idx);
t = t(:,idx);
dtc = dtc(:,idx); % this dtc is valid ONLY for the given c of the pooled data
dx = dx(:,idx);
% and re-define the number of samples to match the remaining data
ns = size(capacities,2);

% visualize simulated degradation trajectories
figure
plot(t,capacities)
xlabel('time')
ylabel('capacity')
title('Simulated Data')
ylim([0 100])
xlim([tlim(1) tlim(2)])
saveas(gcf,'degradation_trajectories__simulated.jpg')
% visualize real degradation trajectories
figure
plot(data0(:,1),cinit-data0(:,2:end),'-o')
xlabel('time')
ylabel('capacity')
title('Real Data')
ylim([0 100])
xlim(tlim)
saveas(gcf,'degradation_trajectories__measured.jpg')

% overlay both in the same plot
% use only the bounds (max and min of the simulations)
figure
plot(data0(:,1),cinit-data0(:,2:end),'-o')
hold on
cmax = max(capacities,[],2);
cmin = min(capacities,[],2);
if iRandom
   for i = 1:numel(cmax)
      tmax(i) = t(i,find(capacities(i,:) == cmax(i),1));
      tmin(i) = t(i,find(capacities(i,:) == cmin(i),1));
   end
else
   tmax = t(:,1);
   tmin = t(:,1);
end

plot(tmax,cmax,'k-','LineWidth',2)
plot(tmin,cmin,'k-','LineWidth',2)
xlabel('time')
ylabel('capacity')
ylim([0 100])
xlim(tlim)
title('Real Data + Simulated Bounds')
legend('Real Data','Simulated - Max','Simulated - Min')
saveas(gcf,'degradation_trajectories__measured_plus_bounds.jpg')


%% PART 2 - MLE of the NHGP
% Determining NHGP Parameters via Maximum Likelihood
[n, m] = size(data0);
% Initial guess of the parameters
abc_guess = [1.0e-3 1.5 1.8]';

% optimization options
opt = optimset ( "TolX",1.e-6, "MaxFunEvals", 100000,"MaxIter", 10000,"display","iter");

% FIRST: REAL DATA - one parameter set per trajectory
warning off
fprintf('\n*** Estimating NHGP parameters via MLE for real data...')
tm = data0(:,1);
for i=2:m
    ym = data0(:,i);
   
    [abc_opt, fval, info, output] = fminunc(@(abc)gamma_log_likelihood_DG(abc,tm,ym,n),abc_guess,opt);
    %% PARAMETERs VALUES
    abc_mle_estimates_real_data(i-1,:)= abc_opt;
    % use a better initial guess (i.e. use last value found)
    abc_guess = abc_opt;
end
fprintf('Done!\n')
fprintf('\n************\n')
abc_mle_estimates_real_data
fprintf('\n************\n')

%% SECOND: SIMULATED DATA - one parameter set per trajectory
fprintf('\n*** Estimating NHGP parameters via MLE for simulated data...')
abc_guess = mean(abc_mle_estimates_real_data);
for i = 1:ns
[abc_opt, fval, info, output] = fminunc(@(abc)gamma_log_likelihood_DR(abc,t(:,i),dx(:,i)),abc_guess,opt);
abc_mle_estimates_sim_data(i,:) = abc_opt;
end
fprintf('Done!\n')
fprintf('\n************\n')
abc_mle_estimates_sim_data
fprintf('\n************\n')

%% THIRD: SIMULATED DATA - one parameter for all trajectories together
fprintf('\n*** Estimating NHGP parameters via MLE for pooled simulated data...')
abc_guess = mean(abc_mle_estimates_sim_data);
[abc_opt, fval, info, output] = fminunc(@(abc)gamma_log_likelihood_DR(abc,t,dx),abc_guess,opt);
fprintf('Done!\n')
fprintf('\n************\n')
abc_mle_estimates_sim_data_pooled = abc_opt
fprintf('\n************\n')
% added this here for comparison
fprintf('\n************\n')
abc_mle_estimates_real_data_pooled
fprintf('\n************\n')

save('ITS_v2_results.mat')

diary off


##a = parameter(:,1);
##b = parameter(:,2);
##c = parameter(:,3);
##w = 20.0;  %threshold --> Change according to pack 
##xo = 0.0;
##wb=(w-xo)./b;
##Epsilon= (wb./a+1./(2*a)).^(1./c)
