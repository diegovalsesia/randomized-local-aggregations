%
% Estimated size of dominating set vs. number of hops
%

addpath ./spgl1-1.9/ ./progress_bar
run('./gspbox/gsp_start.m')

set(0,'DefaultAxesFontSize',13);

G = gsp_minnesota();
%G = gsp_sensor(1000);

G.A = double(full(G.A));
Gcopy = G;

Nhops = 10;
progBar = ProgressBar(Nhops, 'UpdateRate', 0.5);
size_D = zeros(Nhops,1);
Ntx = zeros(Nhops,1);
for hops=1:Nhops 
    [D, size_D(hops)] = get_dominating_set(G);
    Ntx(hops) = number_of_tx(Gcopy,hops,D);
    G.A = G.A + G.A*G.A;
    G.A = G.A - diag(diag(G.A));
    progBar([], [], []);
end

h=figure(1);
plot(1:Nhops,size_D(1:Nhops),'LineWidth',2), grid, xlabel('No. of hops'), ylabel('Size of dominating set'), axis square
saveas(h,'figures/dominating_size_vs_n_hops.fig')

% h2=figure(2);
% plot(1:Nhops,Ntx(1:Nhops)), grid, xlabel('No. of hops'), ylabel('Number of transmissions'), axis square
% saveas(h2,'figures/number_tx_vs_n_hops.fig')
% 
% save dominating_size_and_number_tx_vs_n_hops.mat