clear
warning('off','all')

%% clear figures
for j = 1:11
    figure(j)
    clf
end

%% settings/parameters
num_competitors = 250;% number of total competitors

num_traits = 1; % we want to generalize to more dimensions
f_mode = 2;% Which of the performance functions is to be used
fit_mode = 'cubic'; % Whether to use interpolation, cubic, or quintic fit

genetic_drift = 1*10^(-3);% Amount that the traits can change with reproduction
softmax_power = 5;%Exponent we raise each probability to in the softmax function
num_epochs = 100;% Number of evolutionary steps in each trial

num_experiments = 500;% Number of total trials
games_per_competitor = 100;% Number of games each competitor plays to determine fitness


num_conv_rate_boundary = 10;%Number of times we run convergence rate analysis for boundary cases
num_conv_rate_interior = 10;%Number of times we run convergence rate analysis for interior cases
% num_conv_rate = 10^6;% Number of competitors in convergence rate analysis
std_bounds = [-2.5,-1]; % bounds on range of stds to test for convergence (log base 10)
num_stds = 40;% Number of standard deviations that we run convergence rate analysis over
network_sample_size = 100; % size of network to sample per trial in convergence test
epoch_bounds = [10,10^3]; % max and min number of trials to run in convergence test
tol = 10^(-6); % desired tolerance on numerical differentiation of gradient/Hessian in chicken interpolation

tic;

%% Save parameters to results struct
results.parameters.num_competitors = num_competitors;
results.parameters.num_traits = num_traits;

results.parameters.f_mode = f_mode;
results.parameters.fit_mode = fit_mode;

results.parameters.genetic_drift = genetic_drift;
results.parameters.num_epochs = num_epochs;
results.parameters.num_experiments = num_experiments;
results.parameters.softmax_power = softmax_power;

results.parameters.games_per_competitor = games_per_competitor;
results.parameters.num_conv_rate_boundary = num_conv_rate_boundary;
results.parameters.num_conv_rate_interior = num_conv_rate_interior;
% results.parameters.num_conv_rate = num_conv_rate;
results.parameters.std_bounds = std_bounds;
results.parameters.network_sample_size = network_sample_size;
results.epoch_bounds = [10,10^3];
results.tol = tol;

%% generate edge to endpoint mapping
edge_to_endpoints = NaN(num_competitors*(num_competitors - 1)/2,2);
k = 0;
for i = 1:num_competitors-1
    for j = i+1:num_competitors
        k = k+1;
        edge_to_endpoints(k,:) = [i,j];
    end
end


%% define performance function
if f_mode == 1
     %f = @(x,y) Moran_performance_function_interp(x,y,'log_odds_moran_pd');
     f = @(x,y) performance_pd_moran_24_individuals(x,y,fit_mode);
elseif f_mode == 2
    f = @(x,y) performance_stag_moran_24_individuals(x,y,fit_mode);
elseif f_mode == 3
    f = @(x,y) Moran_performance_function_interp(x,y,'log_odds_moran_chicken_10000_midway');
end


%% preallocate (generate data array that tracks data over multiple experiments)
%experiment_array = zeros(num_experiments,12); % will store experimental results at end of evolution
step_by_step = NaN(num_epochs,7,num_experiments); % will store experimental results over each step
%competitors_step_by_step = NaN(num_competitors,num_epochs,num_experiments); %will save competitor information by step
mu_step_by_step = NaN(num_epochs,round(num_competitors/10),num_experiments);
sigma_step_by_step = NaN(num_epochs,round(num_competitors/10),num_experiments);


grad = NaN(num_traits,num_experiments); % stores gradient in performance at end of evolution
Hessian_xx = NaN(num_traits,num_traits,num_experiments); % stores on diagonal block of the Hessian at end of evolution
Hessian_xy = NaN(num_traits,num_traits,num_experiments); % stores off diagonal block of the Hessian at end of evolution
mus = NaN(round(num_competitors/10),num_epochs,num_experiments); % stores mus of clustering algorithm at end of evolution
sigmas = NaN(round(num_competitors/10),num_epochs,num_experiments); % stores sigmas of clustering algorithm at end of evolution

epsilon_prediction = NaN(num_experiments,1); % stores predicted value of epsilon at end of each experiment
rhos.prediction = NaN(num_experiments,1); % stores predicted value of rho (correlation coefficient) at end of each experiment
rhos.empirical = NaN(num_experiments,1); % stores the empirical value of rho at the end of each experiment
rel_intransitivity_prediction = NaN(num_experiments,1); % stores the predicted relative intransitivity at the end of each experiment
rhos.convergence.boundary.analytic = NaN(num_experiments,num_stds); % stores predicted correlation coefficient as a function of standard deviation around final centroid
rhos.convergence.interior.analytic = NaN(num_experiments,num_stds); % stores predicted correlation coefficient as a function of standard deviation around final centroid
rhos.convergence.boundary.mean = NaN(num_experiments,num_stds); % stores empirical correlation coefficient as a function of standard deviation around final centroid
rhos.convergence.interior.mean = NaN(num_experiments,num_stds); % stores empirical correlation coefficient as a function of standard deviation around final centroid
rhos.convergence.boundary.std = NaN(num_experiments,num_stds); % stores std in empirical correlation coefficient as a function of standard deviation around final centroid
rhos.convergence.interior.std = NaN(num_experiments,num_stds); % stores std in empirical correlation coefficient as a function of standard deviation around final centroid
boundary_correlation_coefficient = NaN(num_experiments,num_stds); % stores empirical correlation coefficient as a function of standard deviation around final centroid if on boundary
interior_correlation_coefficient = NaN(num_experiments,num_stds); % stores empirical correlation coefficient as a function of standard deviation around final centroid if in interior


%% loop over experiments
parfor experiment = 1:num_experiments
    warning('off','all')
    
    %% randomly generate traits
    %[x,y,z] = rand_pick_sphere(num_competitors,0,1); %Call rand_pick_sphere function
    competitor_traits = rand(num_competitors,num_traits);
    
%     %% display scatter of strategies
%     figure(1)
%     clf;
%     hold on
%     histogram(competitor_traits,20,'normalization','pdf')
%     grid on
%     axis square
%     xlim([-0.1,1.1])
%     %     ylim([-0.1,0.1])
%     drawnow
    
    %% generate data array that tracks transitivity/intransitivity over time (epochs)
    epoch_array = zeros(num_epochs,7);
    mu_array = zeros(num_epochs,round(num_competitors/10));
    %traits_array = zeros(num_competitors,num_epochs);
    sigma_array = zeros(num_epochs,round(num_competitors/10));
    cluster_centroid = zeros(1,num_traits);
    evolution_over = 0;%Dummy variable to check if we can stop the evolution process
    on_boundary = 0;%Dummy variable to state if the final location is on the boundary
    final_epoch = num_epochs;%final_epoch will output the right row of the epoch_array to experiment_array, this starts at the maximal epoch and changes if evolution_over changes
    epoch = 1;
    
    %% loop over epochs
    while epoch < final_epoch && evolution_over == 0
        
        %traits_array(:,epoch) = competitor_traits;
        %% Calculate performance using the performance function
        X = competitor_traits(edge_to_endpoints(:,1),:);
        Y = competitor_traits(edge_to_endpoints(:,2),:);
        competition = sparse(edge_to_endpoints(:,1), edge_to_endpoints(:,2),...
            f(X,Y),num_competitors,num_competitors) - ...
            sparse(edge_to_endpoints(:,2), edge_to_endpoints(:,1),...
            f(X,Y),num_competitors,num_competitors);
        competition = full(competition);
        
        %% Sample event outcomes
        n_events = games_per_competitor*num_competitors;
        
        competitors = NaN([n_events,2]);
        
        Z = rand([n_events,1]); % this is all the random numbers we need
        p_win = NaN([n_events,1]);
        
        num_games = 0;
        for first_competitor = 1:num_competitors
            for events = 1:games_per_competitor
                %% count games
                num_games = num_games+1;
                
                %% pick opponent
                stop = 0;
                while stop == 0
                    second_competitor = randperm(num_competitors,1);
                    if second_competitor ~= first_competitor
                        stop = 1;
                    end
                end
                
                %% store competitors
                competitors(num_games,1) = first_competitor;
                competitors(num_games,2) = second_competitor;
                
                %% get win probability
                p_win(num_games) = (1 + exp(-competition(first_competitor,second_competitor)))^(-1); % p_win = logistic(performance)
                
            end
            
        end
        
        Outcomes = (Z <= p_win); % = 1 if i wins, 0 if j wins
        
        %% Compute win frequencies (only left competitor)
        win_freq = sum(reshape(Outcomes, [games_per_competitor,num_competitors])',2)/games_per_competitor;
        indices = (1:num_competitors)';
        winfreq_with_indices = [indices,win_freq];
        
        %% Compute Kendall intransitivity
        ratings_1 = competition > 0;
        ratings_2 = competition == 0;
        kendall_ratings = ratings_1 + 1/2*ratings_2 - 1/2 * eye(num_competitors);
        row_totals = sum(kendall_ratings,2);
        triads = 1/12*num_competitors*(num_competitors - 1)*(2*num_competitors - 1) - 1/2 * sumsqr(row_totals);
        degree_linearity = 1 - 24*triads/(num_competitors^3-4*num_competitors);
        
        %% performing the HHD for a complete graph ***only works for complete graphs
        ratings = (1/num_competitors)*sum(competition,2);
        
        
        %% decompose log odds
        F_t = ratings - ratings';
        F_c = competition - F_t;
        
        %% sizes of components
        Transitivity = norm(F_t,'fro')/sqrt(2);
        Intransitivity = norm(F_c,'fro')/sqrt(2);
        
        
        %% covariance calculations
        cov_at_step = trace(cov(competitor_traits));
        
        %% choose competitors to reproduce based on softmax fitness
        winfreq_with_indices = sortrows(winfreq_with_indices,2, 'descend');
        selectprobs = (exp(winfreq_with_indices(:,2)*softmax_power))/(sum(exp(winfreq_with_indices(:,2)*softmax_power)));
        % selectprobs = selectprobs/norm(selectprobs);
        selectedrowindices = randsample(height(winfreq_with_indices), num_competitors, true, selectprobs);
        reprorows = winfreq_with_indices(selectedrowindices,:);
        highest_fitness = winfreq_with_indices(1:num_competitors/10,:);
        
        %% analyze for clusters      
        AIC = zeros(1,round(num_competitors/10));
        GMModels = cell(1,round(num_competitors/10));
        options = statset('MaxIter',100);
        for num_clusters = 1:(round(num_competitors/10))%Getting an ill-conditioned covariance error
            GMModels{num_clusters} = fitgmdist(competitor_traits,num_clusters,'Options',options,'CovarianceType','full', 'RegularizationValue',genetic_drift/100);
            AIC(num_clusters)= GMModels{num_clusters}.AIC;
        end
        
        [minAIC,numComponents] = min(AIC);
        n_classes = numComponents;
        
        BestModel = GMModels{numComponents};
        mu = BestModel.mu;
        mu(end+1:round(num_competitors/10)) = 0;
        sigma = BestModel.Sigma;
        sigma(end+1:round(num_competitors/10)) = 0;
        
        mu_array(epoch,:) = mu;
        sigma_array(epoch,:) = sigma;
        
        %mus(:,epoch,experiment) = mu;
        %sigmas(:,epoch,experiment) = sigma;
        
        
%         %% display scatter of strategies
%         figure(1)
%         clf;
%         hold on
%         histogram(competitor_traits,20,'normalization','pdf')
%         grid on
%         axis square
%         xlim([-0.1,1.1])
%         
%         plot(linspace(min(competitor_traits),max(competitor_traits),500)',...
%                 pdf(BestModel,linspace(min(competitor_traits),max(competitor_traits),500)'),...
%                 'k-','Linewidth',1.5)
%         for k = 1:n_classes
%             cluster_mu = BestModel.mu(k);
%             cluster_var = BestModel.Sigma(:,:,k);
%             plot(linspace(min(competitor_traits),max(competitor_traits),500),...
%                 BestModel.ComponentProportion(k)*normpdf(linspace(min(competitor_traits),max(competitor_traits),500),cluster_mu,sqrt(cluster_var)), ...
%                 'Linewidth',1.5);
%         end
%         drawnow
        
        
        
        %% cluster centroid (I moved the analyze for clusters function earlier to put this)
        %if cov_at_step < 4*num_traits*genetic_drift^2 %covariance is within twice range we'd expect if it was a single point that expanded outwards
        new_cluster_centroid = mean(competitor_traits);
        if (norm(cluster_centroid - new_cluster_centroid)) < genetic_drift/5 && n_classes == 1%Cluster centroid has stopped moving
            evolution_over = 1;
            final_epoch = epoch;
        end
        cluster_centroid = new_cluster_centroid;
        %end
        
        %% reproduction process
        new_competitor_traits = NaN(num_competitors,num_traits);
        parent_traits = competitor_traits(reprorows(:,1),:);
        for i=1:num_competitors
            drift_vector = genetic_drift*(randn(1,num_traits));
            child_vector = competitor_traits(reprorows(i,1),:) + drift_vector;
            %if norm(child_vector) > 1 (code for sphere only)
            %    child_vector = 1/norm(child_vector)*child_vector;
            %end
            if child_vector > 1
                child_vector = 1;
            elseif child_vector < 0
                child_vector = 0;
            end
            new_competitor_traits(i,:) = child_vector;
        end
        competitor_traits = new_competitor_traits;
       
         
        
        
        %% add to arrays
        for epoch_fill = epoch:num_epochs%Fixed to now fill subsequent epochs up to max
            epoch_array(epoch_fill,1) = Transitivity;
            epoch_array(epoch_fill,2) = Intransitivity;
            epoch_array(epoch_fill,3) = Intransitivity/sqrt(Transitivity^2 + Intransitivity^2);
            epoch_array(epoch_fill,4) = cov_at_step;
            epoch_array(epoch_fill,5) = n_classes;
            epoch_array(epoch_fill,6) = triads;
            epoch_array(epoch_fill,7) = degree_linearity;
            mu_array(epoch_fill,:) = mu;
            sigma_array(epoch_fill,:) = sigma;
        end
        epoch = epoch + 1;
        
        
    end
    step_by_step(:,:,experiment) = epoch_array;
    mu_step_by_step(:,:,experiment) = mu_array;
    sigma_step_by_step(:,:,experiment) = sigma_array;
    %competitors_step_by_step(:,:,experiment) = traits_array;
    
    %% Flag if final centroid is on the boundary
    final_gen_traits = max(parent_traits,[],2);
    max_coordinate = quantile(final_gen_traits,0.25,1);
    if max_coordinate > 1-1.96*genetic_drift
        on_boundary = 1;
    end
    
    % Compute gradient and Hessian at centroid
    if f_mode == 1
        [g,H] = Hessian_for_performance_pd_moran_24_individuals(cluster_centroid,fit_mode);
    elseif f_mode == 2
        [g,H] = Hessian_for_performance_stag_moran_24_individuals(cluster_centroid,fit_mode);
    elseif f_mode == 3
         [g,H] = Hessian_for_performance_chicken_moran_24_individuals(f,cluster_centroid,tol);
    end
    
    grad(:,experiment) = g;
    Hessian_xx(:,:,experiment) = H.xx;
    Hessian_xy(:,:,experiment) = H.xy;
    hessian_xx_norm = norm(H.xx,'fro');
    hessian_xy_norm = norm(H.xy,'fro');
    
    
    %% predict relative size of intransitivity, use to compare to empirical result for the final sample
    Sigma = cov(competitor_traits);
    final_std(experiment) = sqrt(trace(Sigma)/num_traits);
    epsilon_prediction(experiment) = trace(-H.xy*Sigma*H.xy*Sigma)/(2*g'*Sigma*g + trace(H.xx*Sigma*H.xx*Sigma));
    rhos_prediction(experiment) = 1/(2*(1 + epsilon_prediction(experiment)));
    
    E = num_competitors*(num_competitors - 1)/2;
    L = E - (num_competitors - 1);
    rel_intransitivity_prediction(experiment) = sqrt((1 - 2*rhos_prediction(experiment))*(L/E));
    
    
    
    %% empirical correlation
     rhos_empirical(experiment) = (1 - (E/L)*(Intransitivity^2/(Transitivity^2 + Intransitivity^2)))/2;
%     
%        %% predict rho as a function of standard deviation 
%     epsilon = @(std_traits) norm(H.xy,'fro')^2./(2*(norm(g)./std_traits).^2 + norm(H.xx,'fro')^2);
%     rho = @(std_traits) 1./(2*(1 + epsilon(std_traits))); %we will want to compare this to the results from the empirical approach defined below
    
    
%     %% Calculate convergence rate given deviations of competitors from final cluster centroid (boundary)
%     if  num_conv_rate_boundary > 0 && on_boundary == 1 
%         boundary_std_convergence_list = 10.^(linspace(std_bounds(2),std_bounds(1),num_stds));
%         for repopulate_step = 1:length(boundary_std_convergence_list)
%             %% compute theoretical rho and epsilon
%             epsilons(experiment,repopulate_step) = epsilon(boundary_std_convergence_list(repopulate_step));
%             rhos.convergence.boundary.analytic(experiment,repopulate_step) = rho(boundary_std_convergence_list(repopulate_step)); % compare to rhos.mean from below
%             
%             if num_conv_rate_boundary > 0 %run 10 total convergence rate tests empirically
%                 %% estimate rho empirically
%                 [rhos.convergence.boundary.mean(experiment,repopulate_step),rhos.convergence.boundary.std(experiment,repopulate_step)] = estimate_rho_Gauss(f,cluster_centroid,boundary_std_convergence_list(repopulate_step)^2,tol,epoch_bounds,network_sample_size);
%                 boundary_correlation_coefficient(experiment,repopulate_step) = rhos.convergence.boundary.mean(experiment,repopulate_step);
%             end
%             
%         end
%         num_conv_rate_boundary = num_conv_rate_boundary-1;
%     end
%     
%     %% Calculate convergence rate given deviations of competitors from final cluster centroid (interior)
%     if num_conv_rate_interior > 0 && on_boundary == 0 %run 10 total convergence rate tests
%         interior_std_convergence_list = 10.^(linspace(std_bounds(2),std_bounds(1),num_stds));
%         for repopulate_step = 1:length(interior_std_convergence_list)
%             %% compute theoretical rho and epsilon
%             epsilons(experiment,repopulate_step) = epsilon(interior_std_convergence_list(repopulate_step));
%             rhos.convergence.interior.analytic(experiment,repopulate_step) = rho(interior_std_convergence_list(repopulate_step)); % compare to rhos.mean from below
%             
%             %% estimate rho empirically
%             [rhos.convergence.interior.mean(experiment,repopulate_step),rhos.convergence.interior.std(experiment,repopulate_step)] = estimate_rho_Gauss(f,cluster_centroid,interior_std_convergence_list(repopulate_step)^2,tol,epoch_bounds,network_sample_size);
%             interior_correlation_coefficient(experiment,repopulate_step) = rhos.convergence.interior.mean(experiment,repopulate_step);
%         end
%         num_conv_rate_interior = num_conv_rate_interior - 1;
%     end
%     
    %% Add parameters from this to experiments array and structs
%     experiment_array(experiment,1) = epoch_array(1,3);
%     experiment_array(experiment,2) = epoch_array(1,4);
%     experiment_array(experiment,3) = epoch_array(2,5);
%     experiment_array(experiment,4) = epoch_array(final_epoch,3);
%     experiment_array(experiment,5) = epoch_array(final_epoch,4);
%     experiment_array(experiment,6) = epoch_array(final_epoch,5);
%     experiment_array(experiment,7) = final_epoch;
%     experiment_array(experiment,8) = max_coordinate;
%     experiment_array(experiment,9) = on_boundary;
%     experiment_array(experiment,11) = hessian_xx_norm;
%     experiment_array(experiment,12) = hessian_xy_norm;
    
%    experiment_array(experiment,:) = [epoch_array(1,3), epoch_array(1,4), epoch_array(2,5),...
%        epoch_array(final_epoch,3), epoch_array(final_epoch,4), epoch_array(final_epoch,5),...
%        final_epoch, max_coordinate, on_boundary, NaN, hessian_xx_norm, hessian_xy_norm];
    
    %% print
    fprintf('\n Trial %d of %d complete \n',experiment,num_experiments)
end

%% Add parameters to structs
grad_norm = vecnorm(grad);
grad_norm = grad_norm';

results.stepbysteparray = step_by_step;
%results.competitorarray = competitors_step_by_step;
%results.intransitivity.initial = experiment_array(:,1);
%results.covariance.initial = experiment_array(:,2);
%results.clusters.initial = experiment_array(:,3);
%results.intransitivity.final = experiment_array(:,4);
%results.covariance.final = experiment_array(:,5);
%results.clusters.final = experiment_array(:,6);
%results.num.steps = experiment_array(:,7);%(Struct)
%results.maxcoordinate = experiment_array(:,8);
results.norms.gradient = grad_norm;
%results.norms.xxhessian = experiment_array(:,11);
%results.norms.xyhessian = experiment_array(:,12);

results.analysis.grad = grad;
Hess.xx = Hessian_xx;
Hess.xy = Hessian_xy;
results.analysis.Hessian = Hess;
results.analysis.mus = mu_step_by_step;
results.analysis.sigmas = sigma_step_by_step;
results.analysis.epsilon = epsilon_prediction;
results.analysis.rho = rhos_prediction;
results.analysis.rho_empirical = rhos_empirical;
results.analysis.rel_intransitivity = rel_intransitivity_prediction; 
% results.analysis.epsilon_function = epsilon;


% results.convergence_test.boundary.stds = boundary_std_convergence_list;
% results.convergence_test.boundary.rhos = boundary_correlation_coefficient;
% results.convergence_test.boundary.rho_mean = nanmean(boundary_correlation_coefficient);
% results.convergence_test.boundary.rho_predictions = rhos.convergence.boundary.analytic;
% results.convergence_test.interior.stds = interior_std_convergence_list;
% results.convergence_test.interior.rhos = interior_correlation_coefficient;
% results.convergence_test.interior.rho_mean = nanmean(interior_correlation_coefficient);
% results.convergence_test.interior.rho_predictions = rhos.convergence.interior.analytic;

%% Save step-by-step results
step_by_step_array = NaN(num_epochs,6);
step_by_step_array = nansum(step_by_step,3);
step_by_step_array(:,(1:5)) = step_by_step_array(:,(1:5))./step_by_step_array(:,6);

step_by_step_std = NaN(num_epochs,6);
step_by_step_std = std(step_by_step,[],3,'omitnan');

%% Quantiles for step-by-step
step_by_step_first_quantile = quantile(step_by_step, 0.25, 3);
step_by_step_second_quantile = quantile(step_by_step, 0.75, 3);

%% Save step-by-step array to structs
results.bystep.intransitivity.mean = step_by_step_array(:,3);
results.bystep.covariance.mean = step_by_step_array(:,4);
results.bystep.clusters.mean = step_by_step_array(:,5);
results.bystep.remaining.mean = step_by_step_array(:,6);
results.bystep.intransitivity.std = step_by_step_std(:,3);
results.bystep.covariance.std = step_by_step_std(:,4);
results.bystep.clusters.std = step_by_step_std(:,5);
results.bystep.remaining.std = step_by_step_std(:,6);
results.bystep.intransitivity.firstquant = step_by_step_first_quantile(:,3);
results.bystep.covariance.firstquant = step_by_step_first_quantile(:,4);
results.bystep.clusters.firstquant = step_by_step_first_quantile(:,5);
results.bystep.remaining.firstquant = step_by_step_first_quantile(:,6);
results.bystep.intransitivity.secondquant = step_by_step_second_quantile(:,3);
results.bystep.covariance.secondquant = step_by_step_second_quantile(:,4);
results.bystep.clusters.secondquant = step_by_step_second_quantile(:,5);
results.bystep.remaining.secondquant = step_by_step_second_quantile(:,6);

%% Find overall information about our experiment's runtime and save to struct
runtime = toc;
results.parameters.runtime = runtime;

%% save structs
save('evolution_test_moran_results.mat', 'results');

%% visualize data
PlotResults_Moran('evolution_test_moran_results.mat');

