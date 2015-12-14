% This will recreate Figure 5 from Henriksen, Cumming, & Read (2016).

clear all; close all; clc;

bootstrap_mode = 1;
silent = 1;

bem = BEMunit('silent',silent);
bem.Nx = 292; bem.Ny = 292;
bem.deg_per_pixel = 0.03;
bem.dx = 0.1;
bem.temporal_kernel = 'gamma-cosine';
bem.outputNL = @(x)(x.^2);
bem = bem.update();


antibem = bem;
antibem.dx = bem.dx*-1;
antibem = antibem.update();

halfmatched_rds = pairedRDS(); halfmatched_rds.dx = bem.dx/bem.deg_per_pixel;
halfmatched_rds.Nx = 292; halfmatched_rds.Ny = 292;
halfmatched_rds.dotsize=3; % pixels
halfmatched_rds.correlation = 0;

sxs = linspace(0.025,0.30,21);
densities = logspace(log10(0.04),log10(4),11);
density = densities(5);
sx = sxs(4);

halfmatched_rds.density = density;
correlated_rds = halfmatched_rds; correlated_rds.correlation = 1;

mother_bem = bem;
mother_antibem = antibem;

cells_to_run = [3,7,15];


duration = 1; % 1 second trials
N = 1e5;


run_a = 0;
run_b = 1;


freqs = [1,10:10:100];
% Figure 5a
Ma = zeros(length(cells_to_run),length(freqs));

if run_a
    fprintf('Running simulations for Figure 5a.\n');
    for j = 1:length(cells_to_run);
        fprintf('Cell %i of %i\n',j,length(cells_to_run));
        k = cells_to_run(j);
        sx = sxs(k);

        bem = mother_bem.rescale(sx/0.1);
        antibem = mother_antibem.rescale(sx/0.1);


        if bootstrap_mode;
            cbem = bem.load_bootstrap(correlated_rds);
            bem = bem.load_bootstrap(halfmatched_rds);
            antibem = antibem.load_bootstrap(halfmatched_rds);

            % This will now automatically normalize the responses from
            % bem.simulate... .
            %bem = bem.compute_normalization_constant(halfmatched_rds,0,bootstrap_mode);
            %antibem = antibem.compute_normalization_constant(halfmatched_rds,0,bootstrap_mode);
        else

            bem = bem.compute_normalization_constant(halfmatched_rds);
            antibem = antibem.compute_normalization_constant(halfmatched_rds);
        end


        fprintf('Frequencies: ')
        for f = 1:length(freqs);
            fprintf('%i ',f);
            freq = freqs(f);
            HM_bem = zeros(1,N);
            HM_antibem = zeros(1,N);
            C = zeros(1,N);
            parfor i = 1:N;

                HM_bem(i) = mean(bem.simulate_spatiotemporal(halfmatched_rds,freq*duration+1,duration,bootstrap_mode));
                HM_antibem(i) = mean(antibem.simulate_spatiotemporal(halfmatched_rds,freq*duration+1,duration,bootstrap_mode));
                C(i) = mean(cbem.simulate_spatiotemporal(correlated_rds,freq*duration+1,duration,bootstrap_mode));
            end
            %C = mean(cbem.simulate_spatial(correlated_rds,0,1));

            Ma(j,f) = mean(HM_bem-HM_antibem)./mean(C(:));
        end
        fprintf('. Done.\n');    

    end
end

% Figure 5b
fs_to_run = [1,10,100];
Mb = zeros(length(sxs),length(fs_to_run));
if run_b
    fprintf('Running simulations for Figure 5b.\n')
    for j = 1:length(sxs);
        fprintf('Running cell %i of %i.\n',j,length(sxs));
        sx = bem.subunits(1).rf_params.left.sx;
        bem = bem.rescale(sxs(j)/sx);
        antibem = antibem.rescale(sxs(j)/sx);

        if bootstrap_mode
            cbem = bem.load_bootstrap(correlated_rds);
            bem = bem.load_bootstrap(halfmatched_rds);
            antibem = antibem.load_bootstrap(halfmatched_rds);

        end

        for f = 1:length(fs_to_run);
            HM_bem = zeros(1,N);
            HM_antibem = zeros(1,N);
            C = zeros(1,N);
            freq = fs_to_run(f);
            parfor k = 1:N;
                HM_bem(k) = mean(bem.simulate_spatiotemporal(halfmatched_rds,freq*duration+1,duration,bootstrap_mode));
                HM_antibem(k) = mean(antibem.simulate_spatiotemporal(halfmatched_rds,freq*duration+1,duration,bootstrap_mode));
                C(k) = mean(cbem.simulate_spatiotemporal(correlated_rds,freq*duration+1,duration,bootstrap_mode));
            end


            Mb(j,f) = mean(HM_bem-HM_antibem)./mean(C);

        end
    end
    fprintf('Done.\n\n');
end


sxs_pix = sxs/bem.deg_per_pixel;
dotsize_pix = round(correlated_rds.dotsize);

figure();
subplot(1,2,1); hold on;
leg_a = cell(1,length(cells_to_run));
for j = 1:length(cells_to_run);
    rrfsize = round(sxs_pix(cells_to_run(j))/dotsize_pix,1); %relative rf size
    plot(freqs,Ma(j,:),'-','linewidth',3);
    leg_a{j} = sprintf('\\sigma/r=%.1f',rrfsize);
end
legend_a=legend(leg_a);
xlabel('Pattern refresh rate (Hz)','fontsize',16)
ylabel('Normalized response','fontsize',16);

subplot(1,2,2); hold on;
leg_b = cell(1,length(fs_to_run));
for j = 1:length(fs_to_run);    
    plot(sxs_pix(3:end)/dotsize_pix,Mb(3:end,j),'-','linewidth',3);
    leg_b{j} = sprintf('%i Hz',fs_to_run(j));
end
legend_b = legend(leg_b);
xlabel('RF size/dot size','fontsize',16);
ylabel('Normalized response','fontsize',16);

