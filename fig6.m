clear all;

%% Define run-time parameters
% If you want to save the responses to disk (note: this takes up some
% space)
run_bootstrap = 0;
bootstrap_mode = 1;
build_bootstrap_arrays = 1;

scaling_factor=1;

%% Define the parameters of the mother BEMunit; all the other energy model
%% units are derived from this  mother unit.
silent = 1;
bem = BEMunit('silent',silent,'x0',0,'y0',0,'bootstrap_dir','/sid/Modelling/hcr16_ploscb/fig6_data/');
bem.Nx = 292*scaling_factor; bem.Ny=292*scaling_factor;
bem.deg_per_pixel = 0.03/scaling_factor;
bem.temporal_kernel = 'gamma-cosine';
bem.outputNL = @(x)(x.^2); % squaring output nonlinearity

% Set temporal properties of bem unit
bem.tk.tau = 0.035;
bem.tk.omega = 4;
bem.dt=1/(85*2);
 

f = 0.3125; % cycles per SD
%f=0.3;
% Set frequency of Gabor
for j = 1:length(bem.subunits);
    bem.subunits(j).rf_params.left.f=f/bem.subunits(j).rf_params.left.sx;
    bem.subunits(j).rf_params.right.f=f/bem.subunits(j).rf_params.right.sx;
end
bem = bem.update();


%% Create the daughter energy model units.
% There are four types of cells: 
% near_fine, far_fine, coarse_fine, coarse_far
% near_fine has negative disparity, small dx/RF.

% function to map dx to rf size
dx_to_sx = @(x)(0.023 + abs(x)*0.41);
sx = bem.subunits(1).rf_params.left.sx;

fine_near_bem = bem; fine_near_bem.dx = -0.03;
fine_sx = dx_to_sx(fine_near_bem.dx);
fine_near_bem = fine_near_bem.rescale(fine_sx/sx);

coarse_near_bem = bem; coarse_near_bem.dx=-0.48;
coarse_sx = dx_to_sx(coarse_near_bem.dx);
coarse_near_bem = coarse_near_bem.rescale(coarse_sx/sx);

% These are exactly the same except that their disparities are the opposite
% sign.
fine_far_bem = fine_near_bem; 
fine_far_bem.dx = fine_near_bem.dx*-1;
fine_far_bem = fine_far_bem.update();

coarse_far_bem = coarse_near_bem; 
coarse_far_bem.dx = coarse_near_bem.dx*-1;
coarse_far_bem = coarse_far_bem.update();



%% Create stimulus generator object
rds = pairedRDS(); rds.dotsize=3*scaling_factor;
rds.Nx = bem.Nx; rds.Ny = bem.Ny;


%% Stimulus parameters
dxs = [-0.48,-0.03,0.03,0.48]./bem.deg_per_pixel;
correlation_levels = -1:0.2:1; 
freqs = [85/16,85/4,85/2]; % rds pattern refresh rates
duration = 1.5; % stimulus duration in seconds

nt = round(duration/bem.dt); % number of time points
n_bootstrap = 4e4; % this many samples to save
run_parallel = 1;
%% Precompute the monocular responses; this lets us compute a large
%% number of dynamic RDS responses very quickly
if run_bootstrap          
    for j = 1:length(correlation_levels);
        rds.correlation = correlation_levels(j);
        fprintf('Correlation %i of %i\n',j,length(correlation_levels));
        for k = 1:length(dxs);
            rds.dx = dxs(k);
            
            temp_bem = coarse_near_bem.load_bootstrap(rds);
            if isfield(temp_bem.subunits(1),'V_L');
                n_temp = length(temp_bem.subunits(1).V_L);
            else
                n_temp = 0;
            end
            if n_temp < n_bootstrap
                Nb = n_bootstrap-n_temp;
                
                % This guarantees that they see the same image sequence
                seed = randi(1e9,1);
                fine_far_bem.simulate_spatial(rds,Nb,bootstrap_mode,run_parallel,seed);
                fine_near_bem.simulate_spatial(rds,Nb,bootstrap_mode,run_parallel,seed);
                coarse_far_bem.simulate_spatial(rds,Nb,bootstrap_mode,run_parallel,seed);
                coarse_near_bem.simulate_spatial(rds,Nb,bootstrap_mode,run_parallel,seed);
            end
        end
        
    end
end

% Run normalisation if normailsation constant is already computed
run_norm = ~fine_far_bem.check_normalization_constant(rds);
if run_norm;
    fprintf('Calculating normalisation constants... ');
    fine_far_bem = fine_far_bem.compute_normalization_constant(rds,0,bootstrap_mode);
    fine_near_bem = fine_near_bem.compute_normalization_constant(rds,0,bootstrap_mode);
    coarse_near_bem = coarse_near_bem.compute_normalization_constant(rds,0,bootstrap_mode);
    coarse_far_bem = coarse_far_bem.compute_normalization_constant(rds,0,bootstrap_mode);
    fprintf('Done.\n');
end



%% Build bootstrap arrays
% First we build the bootstrap arrays. We make the decision model by
% resampling responses of cells to a sequence of RDSs. The rationale for
% doing this is that generating the stimuli, computing the cells'
% spatiotemporal response and generating the decision would take forever
% with 200 cells. Instead, we get the distribution of responses for each
% cell per trial (given some noise model), assume that their responses 
% are independent (as they see
% independent parts of the image), and take random picks from this
% distribution. Once we've built the bootstrap arrays, we can quickly find 
% the optimal number of cells + noise combination that fits the observed
% data.
n_bootstrap_trials = 1e4; % number of samples in bootstrap distribution
conditions = CombVec(freqs,correlation_levels,dxs);
noise_levels = [0,linspace(5,30,41)];

n_bootstrap=2e4;
if build_bootstrap_arrays
    ff_resp = zeros(size(conditions,2),n_bootstrap_trials,length(noise_levels));
    fn_resp = zeros(size(conditions,2),n_bootstrap_trials,length(noise_levels));
    cn_resp = zeros(size(conditions,2),n_bootstrap_trials,length(noise_levels));
    cf_resp = zeros(size(conditions,2),n_bootstrap_trials,length(noise_levels));
    
    norm_constants = zeros(size(conditions,2),4);
    completed_freqs = zeros(3,4);
    for j = 1:size(conditions,2);
        fprintf('Condition %i of %i\n',j,size(conditions,2));
        freq = conditions(1,j);
        c = conditions(2,j);
        dx = conditions(3,j);
        n_frames = round(duration*freq);
        rds.correlation = c;
        rds.dx = dx;
        
        
                      
        fine_indices = randi(n_bootstrap,1,n_frames*n_bootstrap_trials);
        fine_indices = reshape(fine_indices,[n_frames,n_bootstrap_trials]);
        
        coarse_indices = randi(n_bootstrap,1,n_frames*n_bootstrap_trials);
        coarse_indices = reshape(coarse_indices,[n_frames,n_bootstrap_trials]);

        % Load the data for the subunits
        bems = {fine_far_bem,fine_near_bem,coarse_far_bem,coarse_near_bem};                
        for i = 1:length(bems);
            bem = bems{i};
            
            if c == -1 && j <= 3;                                
                fprintf('freq: %.2f\n', freq);
                crds = rds;
                crds.correlation = 1;
                
                crds.dx = bem.dx/bem.deg_per_pixel;
                
                cbem = bem.load_bootstrap(crds);
                
                this_f = conditions(1,:) == freq;
                
                resps = zeros(1,n_bootstrap_trials);
                for k = 1:n_bootstrap_trials;
                    resps(k) = mean(cbem.simulate_spatiotemporal(crds,n_frames,duration,bootstrap_mode));
                end
                norm_constants(this_f,i) = mean(resps);
            end
                    
                                    
            bem = bem.load_bootstrap(rds);

            current_norm = norm_constants(j,i);
            parfor k = 1:n_bootstrap_trials;
                
                if i < 3
                    idx=fine_indices(:,k);
                else
                    idx=coarse_indices(:,k);
                end
                
                C = bem.simulate_spatiotemporal(rds,n_frames,duration,bootstrap_mode,idx);
                C = C/current_norm;
                
                Cn = zeros(1,length(noise_levels));
                
                for n_i = 1:length(noise_levels);                    
                    Cn(n_i) = sum(C + randn(size(C)).*sqrt(C).*noise_levels(n_i));                    
                end
                
                switch i 
                    case 1
                        ff_resp(j,k,:) = Cn;
                    case 2
                        fn_resp(j,k,:) = Cn;
                    case 3
                        cf_resp(j,k,:) = Cn;
                    case 4
                        cn_resp(j,k,:) = Cn;
                end
            end

        end
        
    end

    %save('fig7_v5_bootstrap_data.mat','ff_resp','fn_resp','cf_resp','cn_resp')
else
    %load('fig7_v5_bootstrap_data.mat');
end

all_n_cells = 50; all_ps = 0.4;

N = 1e3;


%% Run the decision model
fprintf('Running conditions...')
sub_dims = [length(freqs),length(correlation_levels),length(dxs)];
Psi = zeros([length(all_ps),length(all_n_cells),length(noise_levels),sub_dims]);
for j = 1:size(conditions,2);   
    [f,g,k] = ind2sub(sub_dims,j);


    for n_i = 1:length(noise_levels);                
        % Get the bootstrap distribution for the current condition
        current_ff = squeeze(ff_resp(j,:,n_i));
        current_fn = squeeze(fn_resp(j,:,n_i));
        current_cf = squeeze(cf_resp(j,:,n_i));
        current_cn = squeeze(cn_resp(j,:,n_i));

        for nk = 1:length(all_n_cells);
            for p_i = 1:length(all_ps);
                p = all_ps(p_i);
                n_cells = all_n_cells(nk);
                K = round(2*n_cells*p);


                idx = randi(length(current_ff),[n_cells,N]);
                ff = current_ff(idx) - current_fn(idx);                        
                fn = -ff;

                idx = randi(length(current_ff),[n_cells,N]);            
                cf = current_cf(idx) - current_cn(idx);            
                cn = -cf;

                near = [fn;cn];
                far = [ff;cf];

                near_sorted = sort(near,'descend');            
                far_sorted = sort(far,'descend');

                psi = sum(near_sorted(1:K,:)) > sum(far_sorted(1:K,:));


                Psi(p_i,nk,n_i,f,g,k) = mean(psi);
            end
        end
    end
end
fprintf('Done.\n');








%% Find best-fitting model and plot
data = load('doi_data.mat');
x = linspace(0,1,9);
x2 = linspace(0,1,11);
doi_dx(1,:) = interp1(x,data.doi_dx(1,:),x2);
doi_dx(2,:) = interp1(x,data.doi_dx(2,:),x2);
doi_hz(1,:) = interp1(x,data.doi_hz(1,:),x2);
doi_hz(2,:) = interp1(x,data.doi_hz(2,:),x2);

dims = size(Psi);
n_conds = length(noise_levels)*length(all_n_cells);
allSS = zeros(1,n_conds);
for n_i = 1:length(noise_levels);
    for nk = 1:length(all_n_cells);
        for p_i = 1:length(all_ps);

            j = sub2ind(dims(1:3),p_i,nk,n_i);

            temp_Psi = squeeze(Psi(p_i,nk,n_i,:,:,:));
            current_Psi(:,1,:) = (temp_Psi(:,:,2)+(1-temp_Psi(:,:,3)))'/2;
            current_Psi(:,2,:) = (temp_Psi(:,:,1)+(1-temp_Psi(:,:,4)))'/2;
            obs_dx(2,:) = current_Psi(:,1,2);
            obs_dx(1,:) = current_Psi(:,2,2);
            obs_hz(1,:) = current_Psi(:,1,1);
            obs_hz(2,:) = current_Psi(:,1,3);
            SS = sum(abs((obs_dx(:)-doi_dx(:))).^1.0) + ...
                 sum(abs((obs_hz(:)-doi_hz(:))).^1.0);
             
             diff_ac = 15*(obs_hz(1,1)-obs_hz(2,1)).^2;
             diff_dx = 15*(obs_dx(1,1)-0.2).^2;

            allSS(j) = SS+diff_ac;
        end
        
    end
end

[~,idx] = min(allSS);
[~,idxs] = sort(allSS);




%    [p_i,i,j] = ind2sub([length(all_ps),length(all_n_cells),length(noise_levels)],idxs(k));
i = 1; j = 20;
p_i = 1;

        
    

clear PsiM
PsiM(:,:,:,2) = squeeze(Psi(p_i,i,:,:,:,1)+(1-Psi(p_i,i,:,:,:,4)))/2;
PsiM(:,:,:,1) = squeeze(Psi(p_i,i,:,:,:,2)+(1-Psi(p_i,i,:,:,:,3)))/2;
brown = [0.6,0.3,0.1];
%

psi_fine = squeeze(PsiM(j,2,:,1)); psi_coarse = squeeze(PsiM(j,2,:,2));
psi_slow = squeeze(PsiM(j,1,:,1)); psi_fast = squeeze(PsiM(j,3,:,1));

psi_func  = @(x,params)(1-(1-params(1))*exp(-(x/params(2)).^params(3)));
x = linspace(0,1,501); x2 = (x-0.5)*2;
params_fine = fit_psych(psi_fine); params_coarse = fit_psych(psi_coarse);
params_slow = fit_psych(psi_slow); params_fast = fit_psych(psi_fast);
fit_fine = psi_func(x,params_fine);
fit_coarse = psi_func(x,params_coarse);
fit_slow = psi_func(x,params_slow);
fit_fast = psi_func(x,params_fast);

fig=figure(); 
sub1=subplot(1,2,1); hold on;
plot(correlation_levels,psi_fine,'g o','markerfacecolor','g','linewidth',3,'markersize',8);
plot(correlation_levels,psi_coarse,'o','color',brown,'markerfacecolor',brown,'linewidth',3,'markersize',8);

L1=legend('\Deltax=\pm 0.03^\circ','\Deltax=\pm 0.48^\circ');
set(L1,'location','southeast');
plot(x2,fit_fine,'g -','linewidth',2);
plot(x2,fit_coarse,'-','linewidth',2,'color',brown);

xlabel('Correlation');
ylabel('Proportion correct');
xlim([-1,1]); ylim([0,1]);
plot([-1,1],[0.5,0.5],'k --')
set(gca,'ytick',0:0.25:1);

red = [0.8,0.1,0.1];
sub2=subplot(1,2,2); hold on;
plot(correlation_levels,psi_slow,'o' ,'color',red,'markerfacecolor',red,'linewidth',3,'markersize',8);
plot(correlation_levels,psi_fast,'b o','markerfacecolor','b','linewidth',3,'markersize',8);
L2=legend('5.3 Hz', '42.5 Hz');
set(L2,'location','southeast')
plot(x2,fit_slow,'-','linewidth',2,'color',red);
plot(x2,fit_fast,'b -','linewidth',2)


xlabel('Correlation');
ylabel('Proportion correct');
xlim([-1,1]); ylim([0,1]);
set(gca,'ytick',0:0.25:1);
plot([-1,1],[0.5,0.5],'k --')


set_plot_params(gcf);




% 
%     subplot(1,2,1);
%     title(sprintf('ncells: %i (%i); noise %.2f (%i)',all_n_cells(i),i,noise_levels(j),j))
%     subplot(1,2,2);
%     title(sprintf('p=%.2f',p));

A_dx = [fractional_area(fit_fine),fractional_area(fit_coarse)];
a = figure(); hold on;    
plot([0.03,0.48],A_dx,'k o -','markerfacecolor','k','linewidth',2,...
    'markersize',8);

xlabel('\Deltax');
ylabel(sprintf('Fractional\n area'));
ylim([0.0,1.0]); xlim([0,0.55]);
set(gca,'ytick',0:0.5:1,'xtick',[0.05,0.5]);

%set_plot_params(a);

inset_size = 0.15;
inset_handle = gca;
new_handle = copyobj(inset_handle,fig);
close(a);
ax=get(sub1,'Position');
set(new_handle,'Position',[1.8*ax(1)+ax(3)-inset_size 0.1*ax(2)+ax(4)-inset_size inset_size*0.5 inset_size])



subplot(1,2,2);
% Make inset for subplot 2
a = figure(); hold on;    
A_hz = [fractional_area(fit_slow),fractional_area(fit_fast)];
plot([85/16,85/2],A_hz,'k o -','markerfacecolor','k','linewidth',2,...
    'markersize',8);

xlabel('Hz');
ylabel(sprintf('Fractional\n area'));
ylim([0.15,0.45]); 
set(gca,'ytick',[0.2,0.3,0.4],'xtick',round([85/16,85/2],2),'xscale','linear');
xlim([0,85/2+85/16])
%set_plot_params(a);

inset_size = 0.15;
inset_handle = gca;
new_handle = copyobj(inset_handle,fig);
close(a);
ax=get(sub2,'Position');
set(new_handle,'Position',[1.2*ax(1)+ax(3)-inset_size 0.1*ax(2)+ax(4)-inset_size inset_size*0.5 inset_size])


