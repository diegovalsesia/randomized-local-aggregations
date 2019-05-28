function [y,hops,sampling_nodes,Phi] = random_local_aggregations_fast(x,m,G,D1,strategy)
%
% function [y,hops,sampling_nodes,Phi] = random_local_aggregations(x,m,G,D1,strategy)
%
% x : signal
% m : no. of measurements
% G : graph structure as per GSP toolbox
% D1 : 1-hop dominating set (optional)
% strategy : 1 = repeat D nodes with greedy rank check, 2 = new nodes
%
% y : measurements
% hops: power of the adjacency matrix used
% sampling_nodes : indexes of the nodes with the random aggregations (
%    length(unique(sampling_nodes))) is the size of a dominating set of
%    A^hops )
% Phi : sensing matrix
%

if nargin<5
    strategy=1;
end

N = length(x);
G.A = double(full(G.A));

% determine how many communication hops are needed
hops = 1;
if nargin<4
    [D, size_D] = get_dominating_set(G);
else
    D=D1; size_D=length(D1);
end
while size_D > m
    G.A = G.A + G.A*G.A;
    G.A = G.A-diag(diag(G.A));
    G.A = double(G.A~=0);
    [D, size_D] = get_dominating_set(G);
    hops = hops+1;
end

%% OLD
% determine if extra measurements are needed, in case select a node in
% D neighbor of the node with lowest g (i.e. lowest no. of neighbors in D)
%     sampling_nodes = zeros(m,1);
%     sampling_nodes(1:size_D) = D;
%     cur_size = size_D+1;
%     A_dom = zeros(m,N);
%     A_tmp = G.A + eye(N);
%     A_dom(1:size_D,:) = A_tmp(D,:);
% %     while cur_size<=m
% %        % compute g
% %        g = sum(A_dom>0);
% %        [~,nn] = min(g);
% %        ii = find(A_dom(:,nn),1);
% %        A_dom(cur_size,:) = A_dom(ii,:);
% %        sampling_nodes(cur_size) = D(ii);
% %        cur_size = cur_size+1;
% %     end
%     while cur_size<=m
%        % compute g
%        g = sum(A_dom>0);
%        [~,nn] = min(g);
%        A_dom(cur_size,:) = A_tmp(nn,:);
%        sampling_nodes(cur_size) = nn;
%        cur_size = cur_size+1;
%     end

%% NEW
% determine if extra measurements are needed.
% Strategy 1: select a node in D that has the most neighbors with lowest g (i.e. number of times that node has been measured)
% Strategy 2: add any node (not necessarily in D) that has the most neighbors with lowest g
if strategy==1
    
    sampling_nodes = zeros(m,1);
    sampling_nodes(1:size_D) = D;
    eligible_nodes = D;
    nmeas_eligible_nodes = ones(size(eligible_nodes));
    cur_size = size_D+1;
    A_dom = zeros(m,N);
    A_tmp = G.A + eye(N);
    A_dom(1:size_D,:) = A_tmp(D,:);
    alert_flag=0;
    while cur_size<=m && ~isempty(eligible_nodes)
        % compute g and find gmin
        g = sum(A_dom>0);
        gmin = min(g);
        gmax = max(g);
        % select the node in D having most neighbors with g==gmin
        while(1)
            gmasked = g.*A_tmp(eligible_nodes,:);
            N_neig_at_gmin=sum(gmasked==gmin,2);
            if sum(N_neig_at_gmin)==0
                gmin=gmin+1;
                if gmin>gmax
                    alert_flag=1;
                    break
                end
            else
                break
            end
        end
        if alert_flag
            break
        end
        [~,nn] = max(N_neig_at_gmin);
        if rank([A_dom(1:cur_size-1,:);A_tmp(eligible_nodes(nn),:)].*randn(cur_size,N)) ~= min(cur_size,N)
            eligible_nodes(nn) = [];
            nmeas_eligible_nodes(nn) = [];
        else
            A_dom(cur_size,:) = A_tmp(eligible_nodes(nn),:);
            sampling_nodes(cur_size) = eligible_nodes(nn);
            nmeas_eligible_nodes(nn) = nmeas_eligible_nodes(nn)+1;
            if nmeas_eligible_nodes(nn) == sum(A_tmp(eligible_nodes(nn),:)>0)
                eligible_nodes(nn) = [];
                nmeas_eligible_nodes(nn) = [];
            end
            cur_size = cur_size+1;
        end
    end
    
else
    
    sampling_nodes = zeros(m,1);
    sampling_nodes(1:size_D) = D;
    eligible_nodes = setdiff(1:N,D);
    cur_size = size_D+1;
    A_dom = zeros(m,N);
    A_tmp = G.A + eye(N);
    A_dom(1:size_D,:) = A_tmp(D,:);
    alert_flag=0;
    while cur_size<=m && ~isempty(eligible_nodes)
        % compute g and find gmin
        g = sum(A_dom>0);
        gmin = min(g);
        gmax = max(g);
        % select the node in D having most neighbors with g==gmin
        while(1)
            gmasked = g.*A_tmp(eligible_nodes,:);
            N_neig_at_gmin=sum(gmasked==gmin,2);
            if sum(N_neig_at_gmin)==0
                gmin=gmin+1;
                if gmin>gmax
                    alert_flag=1;
                    break
                end
            else
                break
            end
        end
        if alert_flag
            break
        end
        [~,nn] = max(N_neig_at_gmin);
        
        A_dom(cur_size,:) = A_tmp(eligible_nodes(nn),:);
        sampling_nodes(cur_size) = eligible_nodes(nn);
        eligible_nodes(nn) = [];
        cur_size = cur_size+1;       
    end 
    
end

% random combinations
Phi = zeros(m,N);
g = sum(A_dom>0);
for ii=1:N
    Phi(:,ii) = (randn(m,1)./sqrt(g(ii))).*A_dom(:,ii);
end
y = Phi*x;
