clear all;

%% Define run-time parameters
% If you want to save the responses to disk (note: this takes up some
% space)
run_bootstrap = 0;
bootstrap_mode = 1;
build_bootstrap_arrays = 1;


%% Define the parameters of the mother BEMunit; all the other energy model
%% units are derived from this  mother unit.
silent = 1;
bem = BEMunit('silent',silent,'x0',0,'y0',0);
bem.Nx = 292; bem.Ny=292;
bem.deg_per_pixel = 0.03;
bem.temporal_kernel = 'gamma-cosine';
bem.outputNL = @(x)(x.^2); % squaring output nonlinearity

% Set temporal properties of bem unit
bem.tk.tau = 0.035;
bem.tk.omega = 4;
 

f = 0.3125; % cycles per SD
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
rds = pairedRDS(); rds.dotsize=3;
rds.density = 2;
rds.Nx = bem.Nx; rds.Ny = bem.Ny;


%% Stimulus parameters
dxs = [-0.03,0.03]./bem.deg_per_pixel;
correlation_levels = [-1,1];
freq = 120;
alternation_rates = 120./(2.^(5:-1:1));

duration = 0.5; % stimulus duration in seconds
N = 500; % number of frames
nt = round(duration/bem.dt); % number of time points

%% Precompute the monocular responses; this lets us compute a large
%% number of dynamic RDS responses very quickly
if run_bootstrap
    run_parallel = 1;
    n_bootstrap = 2e4; % this many samples to save        
    
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
                fine_far_bem.simulate_spatial(rds,Nb,bootstrap_mode,run_parallel);
                fine_near_bem.simulate_spatial(rds,Nb,bootstrap_mode,run_parallel);
                coarse_far_bem.simulate_spatial(rds,Nb,bootstrap_mode,run_parallel);
                coarse_near_bem.simulate_spatial(rds,Nb,bootstrap_mode,run_parallel);
            end
        end
        
    end
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
M = 500; % number of trials per cell
n_cells = 50; % number of cells
conditions = CombVec(alternation_rates,dxs);

noise_levels = [0,linspace(5,30,41)];
norm_consts = zeros(1,4);

if build_bootstrap_arrays
    ff_resp = zeros(size(conditions,2),N*n_cells,length(noise_levels));
    fn_resp = zeros(size(conditions,2),N*n_cells,length(noise_levels));
    cn_resp = zeros(size(conditions,2),N*n_cells,length(noise_levels));
    cf_resp = zeros(size(conditions,2),N*n_cells,length(noise_levels));
    tic
    for j = 1:size(conditions,2);
        
        fprintf('Condition %i of %i\n',j,size(conditions,2));
        alternation_rate = conditions(1,j);        
        dx = conditions(2,j);
        
        n_frames = round(duration*freq);
        
        
        rds.dx = dx;


        % Load the data for the subunits
        bems = {fine_far_bem,fine_near_bem,coarse_far_bem,coarse_near_bem};                
        for i = 1:length(bems);
            bem = bems{i};

            rds.correlation = 1;
            cbem = bem.load_bootstrap(rds);
            
            rds.correlation=-1;
            acbem = bem.load_bootstrap(rds);
            
            %% Norm constants
            if j == 1;
                norm_rds = rds;
                norm_rds.correlation = 1;
                norm_rds.dx = bem.dx/bem.deg_per_pixel;
                if ~bem.check_bootstrap(norm_rds);
                    fprintf('Running bootstrap...')
                    bem.simulate_spatial(norm_rds,5e3,1);                
                    fprintf(' Done.\n');
                end
                bem = bem.load_bootstrap(norm_rds);
                Cnorms = zeros(1,M);
                parfor k = 1:M;
                    Cnorms(k) = mean(compute_alternating(bem,bem,duration,freq,freq));
                end
                norm_consts(i) = mean(Cnorms);
                
            end
            
            parfor k = 1:(M*n_cells)
                C = compute_alternating(cbem,acbem,duration,freq,alternation_rate)./norm_consts(i);
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

    save('fig10_bootstrap_data.mat','ff_resp','fn_resp','cf_resp','cn_resp')
    
else
    load('fig10_bootstrap_data.mat');
end



n_cells = [50]; p = 0.4;

N = 1000;
dims = size(ff_resp);


%% Run the decision model
fprintf('Running conditions...')
sub_dims = [length(alternation_rates),length(dxs)];
Psi = zeros([length(noise_levels),sub_dims]);
for j = 1:size(conditions,2);            
    for n_i = 1:length(noise_levels);
        % Get the bootstrap distribution for the current condition
        current_ff = squeeze(ff_resp(j,:,n_i));
        current_fn = squeeze(fn_resp(j,:,n_i));
        current_cf = squeeze(cf_resp(j,:,n_i));
        current_cn = squeeze(cn_resp(j,:,n_i));
        
        
        [f,k] = ind2sub(sub_dims,j);
        
        K = round(2*n_cells*p);

        replace=1;
        ff = randsample(current_ff,n_cells*N,replace)-randsample(current_fn,n_cells*N,replace);
        ff=reshape(ff,[n_cells,N]);
        fn = -ff;

        cf = randsample(current_cf,n_cells*N,replace)-randsample(current_cn,n_cells*N,replace);
        cf = reshape(cf,[n_cells,N]);
        cn = -cf;

        near = [fn;cn];
        far = [ff;cf];

        near_sorted = sort(near,'descend');            
        far_sorted = sort(far,'descend');

        psi = sum(near_sorted(1:K,:)) > sum(far_sorted(1:K,:));


        Psi(n_i,f,k) = mean(psi);

    end
end
fprintf('Done.\n');


PsiM = (Psi(:,:,1)+(1-Psi(:,:,2)))/2;
figure();
imagesc(PsiM);

fig = alternating_psych_analysis();
subplot(1,2,1); hold on;
plot(log(alternation_rates),PsiM(16,:),'o m -','linewidth',4,'markersize',10);

subplot(1,2,2); hold on;
plot(log(alternation_rates),PsiM(16,:),'o m -','linewidth',4,'markersize',10);