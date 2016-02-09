function fig3()
    % This will recreate Figure 3 from Henriksen, Cumming, & Read (2016).
    % Note: This requires BEMToolbox - https://github.com/sidh0/BEMtoolbox

    % create binocular energy model unit, and set parameters
    silent = 1;
    bem = BEMunit('silent',silent);
    bem.Nx = 292; bem.Ny = 292;
    bem.deg_per_pixel = 0.03;    
    bem.dx = 0.09;
    
    bem = bem.update();

    % create pairedRDS object and set the stimulus parameters
    rds = pairedRDS();                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
    rds.Nx = bem.Nx; rds.Ny = bem.Ny;

    dx_pix = (bem.dx/bem.deg_per_pixel);
    dxs = round((-12:1:12) + dx_pix); % disparities
    corrs = [-1,0,1]; % anticorrelated, halfmatched and correlated

    bootstrap_mode = 0;
    run_parallel = 1;

    N = 5e3;

    tcs_C = zeros(length(corrs),length(dxs));
    tcs_C2 = zeros(length(corrs),length(dxs));

    for j = 1:length(corrs);
        rds.correlation = corrs(j);
        for k = 1:length(dxs);
            current_rds = rds;
            current_rds.dx = dxs(k);

            C = bem.simulate_spatial(current_rds,N,bootstrap_mode,run_parallel);

            tcs_C(j,k) = mean(C);
            tcs_C2(j,k) = mean(C.^2);
        end

    end

    tcs_C = tcs_C./max(tcs_C(:));
    tcs_C2 = tcs_C2./max(tcs_C2(:));


    %%% Generate plots

    figure();
    % Figure 3a
    subplot(1,2,1); hold on;
    plot(dxs*bem.deg_per_pixel,tcs_C(1,:),'k - ^','linewidth',3,'markersize',10,'markerfacecolor','k');
    plot(dxs*bem.deg_per_pixel,tcs_C(2,:),'b - ^','linewidth',3,'markersize',10,'markerfacecolor','b');
    plot(dxs*bem.deg_per_pixel,tcs_C(3,:),'r - ^','linewidth',3,'markersize',10,'markerfacecolor','r');
    xlabel('Disparity (deg)');
    ylabel('Response');
    set(gca,'ytick',0:.25:1,'xtick',dxs(1:4:21)*bem.deg_per_pixel)
    xlim([min(dxs)*bem.deg_per_pixel-0.05,max(dxs)*bem.deg_per_pixel+0.05]); 
    ylim([0,1]);
    

    % Figure 3b
    subplot(1,2,2); hold on;
    plot(dxs*bem.deg_per_pixel,tcs_C2(1,:),'k - ^','linewidth',3,'markersize',10,'markerfacecolor','k');
    plot(dxs*bem.deg_per_pixel,tcs_C2(2,:),'b - ^','linewidth',3,'markersize',10,'markerfacecolor','b');
    plot(dxs*bem.deg_per_pixel,tcs_C2(3,:),'r - ^','linewidth',3,'markersize',10,'markerfacecolor','r');
    leg=legend('Anticorrelated','Half-matched','Correlated');
    set(leg,'fontsize',14);
    xlabel('Disparity (deg)');
    ylabel('Response');
    set(gca,'ytick',0:.25:1,'xtick',dxs(1:4:21)*bem.deg_per_pixel)
    xlim([min(dxs)*bem.deg_per_pixel-0.05,max(dxs)*bem.deg_per_pixel+0.05]); 
    ylim([0,1]);
       
    set_plot_params(gcf);
end