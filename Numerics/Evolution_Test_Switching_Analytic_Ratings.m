%% Test gaussian-quadratic-replicator_with_diffusion theory using agent based switching model
% Goals: 
% 1. confirm the equations for epsilon, rho, and intransitivity when
% quadratic. DONE
% 2. confirm the Gaussian approximation to ^^ when quadratic. DONE
% 3. confirm that, if we start ~ Gaussian, we stay approximately Gaussian.
% DONE
% 4. confirm steady state covariance eqn. DONE (up to scaling)
% 5. confirm moment dynamics. DONE

% Coding tasks: update so that all parameters are saved to results, save
% video to results, save results, revise so that it can handle repeated
% experiments

%% Clear

clear
warning('off','all')

%% settings/parameters
num_competitors = 2*10^4;% number of total competitors: DEFAULT: 10^4

num_traits = 2; % we want to generalize to more dimensions: DEFAULT: 2

a = 1;%Tuning parameter for matrix xx-Hessian (quadratic performance function): DEFAULT: 1
b = 1;%Tuning parameter for matrix xy-Hessian (quadratic performance function): DEFAULT: 1

f_mode = 4;% Which of the performance functions is to be used: DEFAULT: 4

genetic_drift = 2*10^(-1);% Amount that the traits can change with reproduction, : DEFAULT: 0.2
switching_w = 1/10; % controls the biasing of the switch probability given the rating difference of the agents, want small enough so that w*max(rating diff) < 1/2: DEFAULT: 1/10

tau = 1/num_competitors; % time constant to match ODE

num_steps = round(200*num_competitors/2);% Number of evolutionary steps in each trial
expensive_calculate_rate = round(num_competitors/2); % sampling rate per epoch

num_experiments = 1;% Number of total trials

initial_shift = 5; % distance of initial centroid from the origin: DEFAULT: 5
initial_trait_std = 0.1; % standard deviation (per component) in initial trait distribution: DEFAULT: 1
    
compute_complete_tournament = 0; % 1 to complete full tournaments, 0 to avoid. Suggest 0 for num_competitors > 2*10^3


%% Save parameters to results struct
results.parameters.num_competitors = num_competitors;
results.parameters.num_traits = num_traits;

results.parameters.a = a;
results.parameters.b = b;

results.parameters.f_mode = f_mode;

results.parameters.genetic_drift = genetic_drift;
results.parameters.switching_likelihood = switching_w;
results.parameters.tau = tau;

results.parameters.num_steps = num_steps;
results.parameters.num_experiments = num_experiments;
results.parameters.expensivecalc = expensive_calculate_rate;
results.parameters.initialshift = initial_shift;
results.parameters.initialsd = initial_trait_std;
results.parameters.completetournament = compute_complete_tournament;


%% graph theory basics: generate edge to endpoint mapping (for complete graph evaluation)
E = num_competitors*(num_competitors - 1)/2;
L = E - (num_competitors - 1);

edge_to_endpoints = NaN(num_competitors*(num_competitors - 1)/2,2);
k = 0;
for i = 1:num_competitors-1
    for j = i+1:num_competitors
        k = k+1;
        edge_to_endpoints(k,:) = [i,j];
    end
end


%% define parameters for performance function
% Generate xx-Hessian block A
U = orth(randn(num_traits));
v = exprnd(1,1,num_traits);
S = diag(v);
A = -U*S*U'; % Make matrix symmetric & negative-definite
A = a*A/norm(A, 'fro'); % Normalize

% Generate xy-Hessian block B
V = orth(randn(num_traits));
v2 = exprnd(1,1,num_traits);
for i = 1:floorDiv(num_traits,2)
    v2(2*i) = v2(2*i-1); % Set even-indexed singular values to be equal to preceding odd-indexed s.v.
end
if mod(num_traits,2) == 1
    v2(num_traits) = 0; % If odd, set last singular value to 0
end
S2 = diag(v2);
w = zeros(num_traits);
for i = 1:floorDiv(num_traits,2) % Create rotation matrix w
    w(2*i-1,2*i) = 1;
    w(2*i,2*i-1) = -1;
end
S2 = w*S2;
B = V*S2*V'; %Generate random matrix B with decaying singular values
B = b*B/norm(B, 'fro'); % Normalize

quad_rating_function = @(x) (1/2)*x*A*(x'); % assumes the input is a row vector

%% define performance function
f = @(x,y) example_performance_quad(x,y,A,B);

%% solve for steady state covariance (note: we only need this up to scaling)
[U,S,~] = svd(-A);
A_sqrt = U*sqrt(S)*U';
D = (1/2)*genetic_drift*eye([2,2]); % covariance for genetic drift
M = inv(A_sqrt)*D*inv(A_sqrt);
[U,S,~] = svd(M);
M_sqrt = U*sqrt(S)*U';
steady_state_cov = inv(A_sqrt)*M_sqrt*A_sqrt;


%% predict steady state intransitivity
numerator = trace(-B*steady_state_cov*B*steady_state_cov);
denominator = trace(A*steady_state_cov*A*steady_state_cov);

epsilon.steady_state.gauss = numerator/denominator;
rho.steady_state.gauss = (1/2)*(1/(1 + epsilon.steady_state.gauss));

intransitivity.relative.steady_state.gauss = sqrt((1 - 2*rho.steady_state.gauss)*(L/E));



%% preallocate (generate data array that tracks data over multiple experiments)
expensive_calculation_count = 0;

%% loop over evolution experiments
for experiment = 1:num_experiments
    warning('off','all')
    
    %% randomly generate traits
    % start from a multivariate Gaussian, isovariant, shifted off the
    % origin
    initial_trait_centroid = zeros([1,num_traits]);
    initial_trait_centroid(1) = initial_shift;
        
    competitor_traits = initial_trait_std*randn([num_competitors,num_traits]) + initial_trait_centroid;    

    %% generate data array that tracks transitivity/intransitivity over time (epochs)
    centroid_array = zeros(num_steps,num_traits); % track means
    covariance_array = zeros(num_steps,num_traits,num_traits); % track covariance matrices
    skewness_array = zeros(round(num_steps/expensive_calculate_rate), num_traits); % track skewness
    kurtosis_array = zeros(round(num_steps/expensive_calculate_rate), num_traits); % track kurtosis

    centroid = mean(competitor_traits); 

    evolution_over = 0;%Dummy variable to check if we can stop the evolution process
    on_boundary = 0;%Dummy variable to state if the final location is on the boundary
    final_epoch = num_steps;%final_epoch will output the right row of the epoch_array to experiment_array, this starts at the maximal epoch and changes if evolution_over changes
    step = 1;

    %% compute initial moments
    centroid = mean(competitor_traits);
    covariance = cov(competitor_traits);

    centroid_array(step,:) = centroid;
    covariance_array(step,:,:) = covariance;
    
    loops = final_epoch;

    %% Preallocate for movie
    v = VideoWriter("evolutiontest.mp4",'MPEG-4');
    open(v)

    %% loop over epochs
    while step < final_epoch && evolution_over == 0
        

        %% update competitors
        new_competitor_traits = competitor_traits;

        % pick competitors who could swap
        swap_one = randi(num_competitors); %Pick competitor for switching strategy
        swap_two = randi(num_competitors); %Pick second competitor for switching strategy

        while swap_one == swap_two
            swap_two = randi(num_competitors); % Ensure competitors have different indices
        end

        % determine rating difference (note: these are true, current
        % population ratings, not based on sample outcomes. Models a
        % scenario where the agents are constantly playing games and
        % have a good estimate of their current quality, but swap
        % rarely).
        x = competitor_traits(swap_one,:);
        y = competitor_traits(swap_two,:);
        z = B*centroid';
        delta_rating = (1/2)*(x*A*(x') - y*A*(y')) + (x- y)*z;

        % determine whether to swap
        randintone = rand;
        randinttwo = rand;
        if randintone < (1/2 + switching_w *delta_rating) % swap 2 to 1 (2 swaps to 1 if increase rating)
            new_competitor_traits(swap_two,:) = competitor_traits(swap_one,:);
        end
        if randinttwo < (1/2 - switching_w *delta_rating) % swap 1 to 2 (1 swaps to 2 if increase rating)
            new_competitor_traits(swap_one,:) = competitor_traits(swap_two,:);
        end

        % add drift
        new_competitor_traits = new_competitor_traits + genetic_drift/sqrt(num_competitors)*randn(num_competitors,num_traits); %Add genetic drift

        % update traits
        competitor_traits = new_competitor_traits;
        

        %% Update step count
        step = step + 1; 

        %% update the moments of the competitor distribution
        centroid = mean(competitor_traits);
        centroid_array(step,:) = centroid;
        

        %% compute gradient of f at the centroid
        J = A+ B;
        gradient = J*(centroid');
        norm_of_grad(step) = norm(gradient);

        %% compute covariance
        covariance = cov(competitor_traits);

        norm_of_covariance(step) = norm(covariance,'fro');
        cov_at_step = trace(covariance); 

        covariance_array(step,:,:) = covariance;


        %% compute the Gaussian approximation to epsilon
        epsilon.gauss_approx(step) = trace(-B*covariance*B*covariance)/(2*gradient'*covariance*gradient + trace(A*covariance*A*covariance));

        %% compute rho
        % rho.exact = (1/2)*(1/(1 + epsilon.exact));
        rho.gauss_approx(step) = (1/2)*(1/(1 + epsilon.gauss_approx(step)));

        %% estimate transitivity and intransitivity
        intransitivity.relative.gauss_approx(step) = sqrt((1 - 2*rho.gauss_approx(step))*(L/E));

        %% expensive calculations
        if mod(step,expensive_calculate_rate) == 1
            expensive_calculation_count = expensive_calculation_count + 1;

            expensive_calculation_steps(expensive_calculation_count) = step;

            %% higher order moments for Gaussianity check
            skewness_array(expensive_calculation_count,:) = skewness(competitor_traits);
            kurtosis_array(expensive_calculation_count,:) = kurtosis(competitor_traits);

            %% Calculate performance using the performance function (Complete tournament)
            % ***WARNING: EXPENSIVE, SCALES QUADRATICALLY IN COMPETITOR
            % COUNT
            if compute_complete_tournament == 1
                X = competitor_traits(edge_to_endpoints(:,1),:);
                Y = competitor_traits(edge_to_endpoints(:,2),:);
                performance_values = f(X,Y);
                competition = sparse(edge_to_endpoints(:,1), edge_to_endpoints(:,2),...
                    performance_values,num_competitors,num_competitors) - ...
                    sparse(edge_to_endpoints(:,2), edge_to_endpoints(:,1),...
                    performance_values,num_competitors,num_competitors);
                competition = full(competition);

                %% Compute components
                ratings.pop = mean(competition,2);
                F_t = ratings.pop - ratings.pop';
                F_c = competition - F_t;

                %% sizes of components and intransitivity
                Transitivity = norm(F_t,'fro')/sqrt(2);
                Intransitivity = norm(F_c,'fro')/sqrt(2);

                intransitivity.relative.exact(expensive_calculation_count) = Intransitivity/sqrt(Transitivity^2 + Intransitivity^2);
            end

            %% Compute third and fourth moments
            [Third,Fourth] = third_and_fourth_moments(competitor_traits);

            %% compute exact value for epsilon 
            numerator = trace(-B*covariance*B*covariance);
            
            denominator = gradient'*covariance*gradient +...
                tensorprod(gradient,tensorprod(A,Third,[1 2], [2 3]),1,1) + ...
                (1/4)*(tensorprod(A,tensorprod(A,Fourth,[1 2], [3 4]), [1 2], [1 2]) -...
                tensorprod(A,covariance,[1 2],[1 2])^2);

            epsilon.exact(expensive_calculation_count) = numerator/(2*denominator);

            %% compute exact value for rho and for intransitivity
            rho.exact(expensive_calculation_count) = (1/2)*(1/(1 + epsilon.exact(expensive_calculation_count)));
            intransitivity.relative.exact_moments(expensive_calculation_count) = sqrt((1 - 2*rho.exact(expensive_calculation_count))*(L/E));

            
            %% predict moment changes
            d_mean_dt(expensive_calculation_count,:) = 2*(2*switching_w)*covariance*gradient;
            d_cov_dt(expensive_calculation_count,:,:) = 2*(2*switching_w)*(D + covariance*A*covariance);

            %% compute observed moment changes
            % should update to a central moment
            stride = round(expensive_calculate_rate/2);
            del_mean(expensive_calculation_count,:) = (centroid - centroid_array(step-stride,:))/stride;
            del_cov(expensive_calculation_count,:,:) = (covariance - squeeze(covariance_array(step-stride,:,:)))/stride;

        end

        %% Save parameters to structs

        %% display average quality of agents and agent distribution
        if mod(step,round(num_competitors/4)) == 2
            if num_traits == 2

                %% Compute agent ratings given the current population for the whole population
                X = competitor_traits;
                Y = competitor_traits*A;
                ratings.quad = nan([num_competitors,1]);
                for j = 1:num_competitors
                    x = competitor_traits(j,:);
                    ratings.quad(j) = (1/2)*X(j,:)*(Y(j,:)');
                end
                mean_rating = mean(ratings.quad);

                ratings.pop = ratings.quad - mean_rating + competitor_traits*(B*(centroid_array(step,:))');

                %% plot
                h=figure(1);
                clf

                subplot(4,3,[1, 2, 4, 5, 7, 8])
                hold on
                scatter3(0,0,0,20,'k','filled')
                agents_to_plot = randperm(num_competitors,min(1000,num_competitors));
                scatter3(competitor_traits(agents_to_plot,1),competitor_traits(agents_to_plot,2),ratings.pop(agents_to_plot),5,ratings.pop(agents_to_plot),'filled')
                scatter3(centroid(1),centroid(2),0,12,'k','filled')
                if expensive_calculation_count > 0
                    quiver(centroid_array(expensive_calculation_steps(1:3:expensive_calculation_count),1),centroid_array(expensive_calculation_steps(1:3:expensive_calculation_count),2),...
                        d_mean_dt((1:3:expensive_calculation_count),1),d_mean_dt((1:3:expensive_calculation_count),2),'k','filled','Linewidth',0.5)
                end
                plot(centroid_array(1:step-1,1),centroid_array(1:step-1,2),'k-','Linewidth',1.5)
                grid on
                axis equal
                axis square
                colorbar
                xlim([-2,6])
                ylim([-4,4])
                xlabel('Trait 1','FontSize',16,'interpreter','latex')
                ylabel('Trait 2','FontSize',16,'interpreter','latex')
                zlabel('Population Rating','FontSize',16,'interpreter','latex')
                title_string = strcat('Epoch:', {' '}, num2str(2*step/num_competitors,3));
                title( title_string,'FontSize',16,'interpreter','latex')

                subplot(4,3,10)
                hold on
                histogram(competitor_traits(:,1),30,'Normalization','pdf');
                histogram(competitor_traits(:,2),30,'Normalization','pdf');
                grid on
                xlabel('Traits','FontSize',16,'interpreter','latex')
                l = legend('Trait 1','Trait 2');
                set(l,'FontSize',12,'interpreter','latex','Location','Best')

                if expensive_calculation_count > 0
                    subplot(4,3,11)
                    hold on
                    plot(2*(expensive_calculation_steps(1:expensive_calculation_count)-1)/num_competitors,skewness_array((1:expensive_calculation_count),:),'Linewidth',1.5)
                    plot(2*(expensive_calculation_steps(1:expensive_calculation_count)-1)/num_competitors,kurtosis_array((1:expensive_calculation_count),:),'Linewidth',1.5)
                    grid on
                    %xlabel('Epoch','FontSize',16,'interpreter','latex')
                    ylabel('Moments','FontSize',16,'interpreter','latex')
                    l = legend('Skew 1','Skew 2','Kurt 1','Kurt 2','NumColumns',2);
                    set(l,'FontSize',12,'interpreter','latex','Location','Best')
                    axis tight
                end

                subplot(4,3,3)
                hold on
                plot(2*(2:step)/num_competitors,norm_of_grad(2:step),'b','Linewidth',1.5)
                grid on
                axis tight
                %xlabel('Epochs','FontSize',16,'interpreter','latex')
                ylabel('$\|$ Grad $\|$','FontSize',16,'interpreter','latex')


                subplot(4,3,6)
                hold on
                plot([0,2*step/num_competitors],norm(steady_state_cov,'fro')*[1,1],'k--','Linewidth',1)
                plot(2*(2:step)/num_competitors,norm_of_covariance(2:step),'b','Linewidth',1.5)
                grid on
                axis tight
                %xlabel('Epochs','FontSize',16,'interpreter','latex')
                ylabel('$\|$ Covariance $\|$','FontSize',16,'interpreter','latex')
                l = legend('Steady State','Covariance');
                set(l,'FontSize',12,'interpreter','latex','location','best')


                subplot(4,3,9)
                hold on
                plot([0,2*step/num_competitors],rho.steady_state.gauss*[1,1],'k--','Linewidth',1)
                if expensive_calculation_count > 0
                    plot(2*(expensive_calculation_steps(1:expensive_calculation_count)-1)/num_competitors,rho.exact((1:expensive_calculation_count)),'k','Linewidth',1.5)
                end
                plot(2*(2:step)/num_competitors,rho.gauss_approx(2:step),'b','Linewidth',1.5)
                grid on
                axis tight
                ylim([0,0.5])
                %xlabel('Epochs','FontSize',16,'interpreter','latex')
                ylabel('Correlation $\rho$','FontSize',16,'interpreter','latex')
                l = legend('Steady State','Exact (Moments)','Gauss Approx');
                set(l,'FontSize',12,'interpreter','latex','location','best')


                subplot(4,3,12)
                hold on
                plot([0,2*step/num_competitors],intransitivity.relative.steady_state.gauss*[1,1],'k--','Linewidth',1)
                if expensive_calculation_count > 0
                    if compute_complete_tournament == 1
                        plot(2*(expensive_calculation_steps(1:expensive_calculation_count)-1)/num_competitors,intransitivity.relative.exact((1:expensive_calculation_count)),'r','Linewidth',1.5)
                    end
                    plot(2*(expensive_calculation_steps(1:expensive_calculation_count)-1)/num_competitors,intransitivity.relative.exact_moments((1:expensive_calculation_count)),'k','Linewidth',1.5)
                end
                plot(2*(2:step)/num_competitors,intransitivity.relative.gauss_approx(2:step),'b','Linewidth',1.5)
                grid on
                axis tight
                ylim([0,1])
                xlabel('Epochs','FontSize',16,'interpreter','latex')
                %ylabel('Relative Intransitivity','FontSize',16,'interpreter','latex')
                ylabel('Intransitivity','FontSize',16,'interpreter','latex')
                if expensive_calculation_count > 0
                    if compute_complete_tournament == 1
                        l = legend('Steady State','Exact (Tourney)','Exact (Moments)','Gauss Approx');
                    else
                        l = legend('Steady State','Exact (Moments)','Gauss Approx');
                    end
                    set(l,'FontSize',12,'interpreter','latex','location','best')
                end
                frame = getframe(gcf);
                writeVideo(v,frame)
                % saveas(h,sprintf('FIRSTFIG%d.png',step));


                j = figure(2);
                clf

                if expensive_calculation_count > 0
                    subplot(1,2,1)
                    hold on
                    plot(1.1*[min(min([d_mean_dt*tau,del_mean])),max(max([d_mean_dt*tau,del_mean]))],1.1*[min(min([d_mean_dt*tau,del_mean])),max(max([d_mean_dt*tau,del_mean]))],'k-','Linewidth',1)
                    scatter(d_mean_dt((1:expensive_calculation_count),1)*tau,del_mean((1:expensive_calculation_count),1),10,'r','filled');
                    scatter(d_mean_dt((1:expensive_calculation_count),2)*tau,del_mean((1:expensive_calculation_count),2),10,'b','filled');
                    grid on
                    axis square
                    axis tight
                    xlabel('Predicted','FontSize',16,'interpreter','latex')
                    ylabel('Observed','FontSize',16,'interpreter','latex')
                    title('Rate of Change in Centroid','FontSize',16,'interpreter','latex')
                    l = legend('Parity','Trait 1', 'Trait 2');
                    set(l,'FontSize',16,'interpreter','latex','Location','NorthWest')

                    subplot(1,2,2)
                    hold on
                    plot(1.1*[min(min(min([d_cov_dt*tau,del_cov]))),max(max(max([d_cov_dt*tau,del_cov])))],1.1*[min(min(min([d_cov_dt*tau,del_cov]))),max(max(max([d_cov_dt*tau,del_cov])))],'k-','Linewidth',1)
                    scatter(d_cov_dt((1:expensive_calculation_count),1,1)*tau,del_cov((1:expensive_calculation_count),1,1),10,'r','filled');
                    scatter(d_cov_dt((1:expensive_calculation_count),2,1)*tau,del_cov((1:expensive_calculation_count),2,1),10,'m','filled');
                    scatter(d_cov_dt((1:expensive_calculation_count),2,2)*tau,del_cov((1:expensive_calculation_count),2,2),10,'b','filled');
                    grid on
                    axis square
                    axis tight
                    xlabel('Predicted','FontSize',16,'interpreter','latex')
                    ylabel('Observed','FontSize',16,'interpreter','latex')
                    title('Rate of Change in Covariance','FontSize',16,'interpreter','latex')
                    l = legend('Parity','Var in 1', 'Cov 1,2', 'Var in 2');
                    set(l,'FontSize',16,'interpreter','latex','Location','NorthWest')
                end



                drawnow
                % saveas(j,sprintf('SECONDFIG%d.png',step));
            end
        end
        
        
    end

    %% Save parameters to structs
    results.bystep.centroid = centroid_array;
    results.bystep.covariance = covariance_array;
    results.bystep.skewness = skewness_array;
    results.bystep.kurtosis = kurtosis_array;
    results.intransitivity = intransitivity;

    % step_by_step(:,:,experiment) = epoch_array;
   
    
    %% print
    fprintf('\n Trial %d of %d complete \n',experiment,num_experiments)

end

%% save structs
close(v)
save('evolution_test_sw_ar_results.mat', 'results');
