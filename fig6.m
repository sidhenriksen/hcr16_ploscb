
% Generating Figure 6 from Henriksen, Cumming, & Read (2016).
% Tuning curves for the different cells.
% Be warned: This takes a while to run.
function fig6()
    bem = BEMunit();
    bem.Nx = 292; bem.Ny = 292;
    bem.deg_per_pixel = 0.03;
    dxs = [-0.48,-0.03,0.03,0.48];
    dx_to_sx = @(x)(0.023 + abs(x)*0.4);
    mother_bem = bem;

    N = 20000; % Number of frames to average over.

    % We sample the tuning curves more finely in the range where the
    % fine disparity cells are active
    stimulus_dxs = [-50:2:-6,-5:5,6:2:50];

    rds = pairedRDS();
    rds.Nx = bem.Nx; rds.Ny = bem.Ny;

    tcs = zeros(length(dxs),length(stimulus_dxs));


    % Run this for each cell
    for j = 1:length(dxs);
        fprintf('Cell %i of %i.\n',j,length(dxs));
        sx = dx_to_sx(dxs(j));
        bem = mother_bem.rescale(sx/0.1);
        bem.dx = dxs(j);
        bem = bem.update();


        bem = bem.compute_normalization_constant(rds);


        for k = 1:length(stimulus_dxs);
            rds.dx = stimulus_dxs(k);

            C = bem.simulate_spatial(rds,N);

            tcs(j,k) = mean(C);
        end
    end

    figure(); hold on;
    red=[0.8,0.1,0.1];
    green=[0.1,0.8,0.1];
    cols = {red,green,green,red};
    for j = 1:length(dxs);
        plot(stimulus_dxs*bem.deg_per_pixel,tcs(j,:),'-','linewidth',2);
    end

    xlabel('Disparity (deg)','fontsize',16);
    ylabel('Normalized response','fontsize',16);
end