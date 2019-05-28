clear all
close all
clc

% SPGL1 for l1 minimization
addpath ./spgl1-1.9/ ./progress_bar
% Graph Signal Processing Toolbox
run('./gspbox/gsp_start.m')

opts = spgSetParms('verbosity',0,'iterations',10000);
warning('off','all')

vec = @(a) a(:);

%% Parameters
N = 100;   % signal length
k = 10;     % sparsity
supp = 1:k; % lowpass support
%supp = randperm(N); supp = supp(1:k);  % random sparsity support
supp_logi = false(N,1); supp_logi(supp) = true; % 1 in supp, 0 otherwise
notsupp_logi = ~supp_logi; % 0 in supp, 1 otherwise
graph_name = 'sensor';
supp_known = 1; % is the support known?

p = 0.5;
M = 30;
sigma_noise = 1e-5;

% create graph and graph signal
[x_bl,U,G,x_bl_noiseless] = create_graph_signal(graph_name, k, sigma_noise, N, notsupp_logi, p);
% compute dominating set
[D1, ~] = get_dominating_set(G);


%% Sampling with randomized local aggregations (two strategies)
Gcopy = G;
Gcopy2 = G;
[y_la1,hops1,sampling_nodes1,Phi1] = random_local_aggregations_fast(x_bl,M,Gcopy,D1,1);
[y_la2,hops2,sampling_nodes2,Phi2] = random_local_aggregations_fast(x_bl,M,Gcopy2,D1,2);


%% Reconstruction
if supp_known
    % pseudoinverse reconstruction
    
    X_hat_la1 = zeros(G.N,1);
    X_hat_la1(supp_logi) = ( Phi1*U(:,supp_logi) ) \ y_la1;
    x_rec_la1 = U*X_hat_la1;
    msedB_suppknown_la1 = 10*log10( mean( (x_bl_noiseless-x_rec_la1).^2 ) );
    
    X_hat_la2 = zeros(G.N,1);
    X_hat_la2(supp_logi) = ( Phi2*U(:,supp_logi) ) \ y_la2;
    x_rec_la2 = U*X_hat_la2;
    msedB_suppknown_la2 = 10*log10( mean( (x_bl_noiseless-x_rec_la2).^2 ) );
    
else
    % l1 minimization reconstruction
    
    X_hat_la1 = spg_bp(Phi1*U, y_la1, opts);
    x_rec_la1 = U*X_hat_la1;
    msedB_l1_la1 = 10*log10( mean( (x_bl_noiseless-x_rec_la1).^2 ) );
    
    X_hat_la2 = spg_bp(Phi2*U, y_la2, opts);
    x_rec_la2 = U*X_hat_la2;
    msedB_l1_la2 = 10*log10( mean( (x_bl_noiseless-x_rec_la2).^2 ) );
    
end
