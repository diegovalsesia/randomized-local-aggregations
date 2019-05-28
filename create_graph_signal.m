function [x_bl,U,G,x_bl_noiseless] = create_graph_signal(graph_name, k, sigma_noise, N, notsupp_logi, params)
%
% function [x_bl,U,G,x_bl_noiseless] = create_graph_signal(graph_name, k, sigma_noise, N, notsupp_logi, params)
%
% graph_name = {'erdos', 'minnesota', 'ring', 'bunny', 'full', 'sensor','scale_free'}
% k : sparsity
% sigma_noise : noise standard deviation
% N (optional) : no. of vertices
% notsupp_logi (optional) : logical vector with entries not in the support
% params : extra parameters for specific graphs
%
% If support is not specified, the first k compoennets are used
%

    if nargin < 3
        error('Too few parameters');
    end
    
    if nargin == 3
    	N = 256; 
    end
    
    %% Graph and GFT
    switch(graph_name)
        case 'erdos'
            G = gsp_erdos_renyi(N,params);
        case 'minnesota'
            G = gsp_minnesota(); N = G.N;
        case 'ring'
            G = gsp_ring(N);
        case 'bunny'
            G = gsp_bunny(); N = G.N;
        case 'full'
            G = gsp_full_connected(N);
        case 'sensor'
            G = gsp_sensor(N);
        case 'scale_free'
            G = gsp_barabasi_albert(N);
        case 'community'
            G = gsp_community(N);
        case '2dgrid'
            G = gsp_2dgrid(sqrt(N),sqrt(N));
    end
    
    G.A = double(G.A);
    
    G = gsp_create_laplacian(G,'normalized');
    G = gsp_compute_fourier_basis(G);
    U = G.U;

    %% Signal on the graph
    x = randn(N,1);
    X = U'*x;
    
    if nargin < 5
        notsupp_logi = true(N,1);
        notsupp_logi(1:k) = false;
    end
    
    X(notsupp_logi) = 0; % sparse
    x_bl_noiseless = U*X;
    x_bl = x_bl_noiseless + sigma_noise*randn(N,1);


end