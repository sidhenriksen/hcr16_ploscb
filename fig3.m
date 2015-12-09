% This will recreate figure 3 from Henriksen, Cumming, & Read (2016).
% Note: This requires BEMToolbox - https://github.com/sidh0/BEMtoolbox

% create binocular energy model unit, and set parameters
bem = BEMunit();    
bem.Nx = 292; bem.Ny = 292;
bem.dx=-0.1;
bem.rescale(0.5);
bem = bem.update();
%bem.silent = 1;

% create pairedRDS object and set the stimulus parameters
rds = pairedRDS();                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  

rds.Nx = 292; rds.Ny = 292;

dx_pix = (bem.dx/bem.deg_per_pixel);
dxs = round((-10:10) - dx_pix); % disparities
corrs = [-1,0,1]; % anticorrelated, halfmatched and correlated

bootstrap_mode = 0;

N = 5e3;

tcs_C = zeros(length(corrs),length(dxs));
tcs_C2 = zeros(length(corrs),length(dxs));

for j = 1:length(corrs);
    rds.correlation = corrs(j);
    for k = 1:length(dxs);
        current_rds = rds;
        current_rds.dx = dxs(k);

        C = bem.simulate_spatial(current_rds,N,bootstrap_mode);
        
        tcs_C(j,k) = mean(C);
        tcs_C2(j,k) = mean(C.^2);
    end
    
end

tcs_C = tcs_C./max(tcs_C(:));
tcs_C2 = tcs_C2./max(tcs_C2(:));



figure();
subplot(1,2,1); hold on;
plot(dxs*bem.deg_per_pixel,tcs_C(1,:),'k - ^','linewidth',3,'markersize',10,'markerfacecolor','k');
plot(dxs*bem.deg_per_pixel,tcs_C(2,:),'b - ^','linewidth',3,'markersize',10,'markerfacecolor','b');
plot(dxs*bem.deg_per_pixel,tcs_C(3,:),'r - ^','linewidth',3,'markersize',10,'markerfacecolor','r');
xlabel('Disparity (deg)');
ylabel('Response');
set(gca,'ytick',0:.25:1,'xtick',dxs(1:4:21)*bem.deg_per_pixel)
xlim([0,0.6]); ylim([0,1]);
text(0.05,0.9,'a)','fontsize',16);

subplot(1,2,2); hold on;
plot(dxs*bem.deg_per_pixel,tcs_C2(1,:),'k - ^','linewidth',3,'markersize',10,'markerfacecolor','k');
plot(dxs*bem.deg_per_pixel,tcs_C2(2,:),'b - ^','linewidth',3,'markersize',10,'markerfacecolor','b');
plot(dxs*bem.deg_per_pixel,tcs_C2(3,:),'r - ^','linewidth',3,'markersize',10,'markerfacecolor','r');
xlabel('Disparity (deg)');
ylabel('Response');
set(gca,'ytick',0:.25:1,'xtick',dxs(1:4:21)*bem.deg_per_pixel)
xlim([0,0.6]); ylim([0,1]);
text(0.05,0.9,'b)','fontsize',16);

set(gcf,'color','white');