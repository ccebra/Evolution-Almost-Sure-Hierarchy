%% Clear figures
for j = 1:6
    figure(j);
    clf;
end
clear;

%% Load data from staghunt trials
load('evolution_test_moran_results_stag_lowGED.mat')
final_covariance_1 = sqrt(squeeze(results.stepbysteparray(results.parameters.num_epochs,4,:)));
empirical_rho_1 = (0.5-results.analysis.rho_empirical)';
grad_1 = abs(results.analysis.grad)';
drift_1 = results.parameters.genetic_drift*ones(results.parameters.num_experiments,1);
load('evolution_test_moran_results_stag_highGD.mat')
final_covariance_2 = sqrt(squeeze(results.stepbysteparray(results.parameters.num_epochs,4,:)));
empirical_rho_2 = (0.5-results.analysis.rho_empirical)';
grad_2 = abs(results.analysis.grad)';
drift_2 = results.parameters.genetic_drift*ones(results.parameters.num_experiments,1);
load('evolution_test_moran_results_stag_1000.mat')
final_covariance_3 = sqrt(squeeze(results.stepbysteparray(results.parameters.num_epochs,4,:)));
empirical_rho_3 = (0.5-results.analysis.rho_empirical)';
grad_3 = abs(results.analysis.grad)';
drift_3 = results.parameters.genetic_drift*ones(results.parameters.num_experiments,1);
load('evolution_test_moran_results_stag_5E-3.mat')
final_covariance_4 = sqrt(squeeze(results.stepbysteparray(results.parameters.num_epochs,4,:)));
empirical_rho_4 = (0.5-results.analysis.rho_empirical)';
grad_4 = abs(results.analysis.grad)';
drift_4 = results.parameters.genetic_drift*ones(results.parameters.num_experiments,1);
load('evolution_test_moran_results_stag_vvhighGD.mat')
final_covariance_5 = sqrt(squeeze(results.stepbysteparray(results.parameters.num_epochs,4,:)));
empirical_rho_5 = (0.5-results.analysis.rho_empirical)';
grad_5 = abs(results.analysis.grad)';
drift_5 = results.parameters.genetic_drift*ones(results.parameters.num_experiments,1);

%% Combine vectors
final_cov_analysis = cat(1,final_covariance_1,final_covariance_2,final_covariance_3,final_covariance_4,final_covariance_5);
empirical_rho_analysis = cat(1,empirical_rho_1,empirical_rho_2,empirical_rho_3,empirical_rho_4,empirical_rho_5);
grad_analysis = cat(1,grad_1,grad_2,grad_3,grad_4,grad_5);
drift_analysis = cat(1,drift_1,drift_2,drift_3,drift_4,drift_5);
data_frame = cat(2,final_cov_analysis,empirical_rho_analysis,grad_analysis,drift_analysis);
data_frame_1 = cat(2,final_covariance_1,empirical_rho_1,grad_1,drift_1);
data_frame_2 = cat(2,final_covariance_2,empirical_rho_2,grad_2,drift_2);
data_frame_3 = cat(2,final_covariance_3,empirical_rho_3,grad_3,drift_3);
data_frame_4 = cat(2,final_covariance_4,empirical_rho_4,grad_4,drift_4);
data_frame_5 = cat(2,final_covariance_5,empirical_rho_5,grad_5,drift_5);

%% Generate convex hulls for each trial
% idx = data_frame_1(:,3) < 0.02;
% data_frame_1 = data_frame_1(idx,:);
% idx = data_frame_2(:,3) < 0.02;
% data_frame_2 = data_frame_2(idx,:);
% idx = data_frame_3(:,3) < 0.02;
% data_frame_3 = data_frame_3(idx,:);
% idx = data_frame_4(:,3) < 0.02;
% data_frame_4 = data_frame_4(idx,:);
% idx = data_frame_5(:,3) < 0.02;
% data_frame_5 = data_frame_5(idx,:);

full_matrix = cat(1,data_frame_1,data_frame_2,data_frame_3,data_frame_4,data_frame_5);
full_matrix = [log(full_matrix(:,1)),log(full_matrix(:,2))];
matrix_1 = [data_frame_1(:,1),data_frame_1(:,2)];
matrix_2 = [data_frame_2(:,1),data_frame_2(:,2)];
matrix_3 = [data_frame_3(:,1),data_frame_3(:,2)];
matrix_4 = [data_frame_4(:,1),data_frame_4(:,2)];
matrix_5 = [data_frame_5(:,1),data_frame_5(:,2)];
% hull_matrix = cat(1,data_frame_1,data_frame_2,data_frame_3,data_frame_4,data_frame_5);
% hull_matrix = [log(hull_matrix(:,1)),log(hull_matrix(:,2))];
% [k,av] = convhull(hull_matrix);
[k,av] = convhull(full_matrix);
[k2,av2] = convhull(matrix_1);
[k3,av3] = convhull(matrix_2);
[k4,av4] = convhull(matrix_3);
[k5,av5] = convhull(matrix_4);
[k6,av6] = convhull(matrix_5);


%% 
figure(8)
clf
plot(full_matrix(:,1),full_matrix(:,2),'*');
hold on
plot(full_matrix(k,1),full_matrix(k,2));

k7 = [4829,4165,555,336]';
figure(7)
clf
plot(full_matrix(:,1),full_matrix(:,2),'*');
hold on
plot(full_matrix(k7,1),full_matrix(k7,2));
drawnow

f = fit(full_matrix(k7,1),full_matrix(k7,2),'poly1');
figure(9)
clf
plot(full_matrix(:,1),full_matrix(:,2),'*');
hold on
plot(f,full_matrix(k7,1),full_matrix(k7,2));
drawnow

%% Draw lines
%% Plot results
figure(1)
clf
hold on
scatter(data_frame_1(:,1),data_frame_1(:,2),8,data_frame_1(:,3));
hold on
fill(matrix_1(k2,1),matrix_1(k2,2),'r','facealpha',0.05);
hold on
scatter(data_frame_2(:,1),data_frame_2(:,2),8,data_frame_2(:,3));
hold on
fill(matrix_2(k3,1),matrix_2(k3,2),'m','facealpha',0.05);
hold on
scatter(data_frame_3(:,1),data_frame_3(:,2),8,data_frame_3(:,3));
hold on
fill(matrix_3(k4,1),matrix_3(k4,2),'b','facealpha',0.05);
hold on
scatter(data_frame_4(:,1),data_frame_4(:,2),8,data_frame_4(:,3));
hold on
fill(matrix_4(k5,1),matrix_4(k5,2),'g','facealpha',0.05);
hold on
scatter(data_frame_5(:,1),data_frame_5(:,2),8,data_frame_5(:,3));
hold on
fill(matrix_5(k6,1),matrix_5(k6,2),'p','facealpha',0.15);
colorbar;
grid on
set(gca,'xscale','log','yscale','log')
title('Concentration vs. $\rho$ Predicted','FontSize',16,'interpreter','latex')
xlabel('$\sqrt{\frac{1}{T}E[||X - \bar{x}||^2]}$','FontSize',16,'interpreter','latex')
ylabel('Empirical $0.5 - \rho$','FontSize',16,'interpreter','latex')
axis square
drawnow

stag_data_frame = cat(1,data_frame_1,data_frame_2,data_frame_3,data_frame_4,data_frame_5);
[stag_y_max,stag_max_index] = max(stag_data_frame(:,2));
stag_x_max = (stag_data_frame(stag_max_index,1));
[stag_x_min,stag_min_index] = min(stag_data_frame(:,1));
stag_y_sqr = stag_y_max*(stag_x_min/stag_x_max)^2;
stag_y_quart = stag_y_max*(stag_x_min/stag_x_max)^4;

figure(2)
clf
hold on
scatter(data_frame_1(:,1),data_frame_1(:,2),8,data_frame_1(:,3));
hold on
scatter(data_frame_2(:,1),data_frame_2(:,2),8,data_frame_2(:,3));
hold on
scatter(data_frame_3(:,1),data_frame_3(:,2),8,data_frame_3(:,3));
hold on
scatter(data_frame_4(:,1),data_frame_4(:,2),8,data_frame_4(:,3));
hold on
scatter(data_frame_5(:,1),data_frame_5(:,2),8,data_frame_5(:,3));
plot([stag_x_min,stag_x_max],[stag_y_sqr,stag_y_max],'k','LineWidth',1);
plot([stag_x_min,stag_x_max],[stag_y_quart,stag_y_max],'k','LineWidth',1);
colorbar;
grid on
set(gca,'xscale','log','yscale','log','FontSize',24)
title('Stag','FontSize',32,'interpreter','latex')
xlabel('Concentration','FontSize',32,'interpreter','latex')
axis square
drawnow

figure(3)
clf
hold on
scatter(data_frame(:,1),data_frame(:,2),10,data_frame(:,3),'filled');
colorbar;
grid on
set(gca,'xscale','log','yscale','log')
title('Concentration vs. $\rho$ Predicted (Gradient coloring)','FontSize',16,'interpreter','latex')
xlabel('$\sqrt{\frac{1}{T}E[||X - \bar{x}||^2]}$','FontSize',16,'interpreter','latex')
ylabel('Empirical $0.5 - \rho$','FontSize',16,'interpreter','latex')
axis square
drawnow

%% Load data from chicken trials
load('evolution_test_moran_results_chicken_5E-5.mat')
final_covariance_1 = sqrt(squeeze(results.stepbysteparray(results.parameters.num_epochs,4,:)));
empirical_rho_1 = (0.5-results.analysis.rho_empirical)';
grad_1 = abs(results.analysis.grad)';
hessian_1 = squeeze(abs(results.analysis.Hessian.xx));
drift_1 = results.parameters.genetic_drift*ones(results.parameters.num_experiments,1);
mus_1 = squeeze(results.analysis.mus(50,1,:));
clusters_1 = squeeze(results.stepbysteparray(50,5,:));
load('evolution_test_moran_results_chicken_5E-4.mat')
final_covariance_2 = sqrt(squeeze(results.stepbysteparray(results.parameters.num_epochs,4,:)));
empirical_rho_2 = (0.5-results.analysis.rho_empirical)';
grad_2 = abs(results.analysis.grad)';
drift_2 = results.parameters.genetic_drift*ones(results.parameters.num_experiments,1);
mus_2 = squeeze(results.analysis.mus(50,1,:));
clusters_2 = squeeze(results.stepbysteparray(50,5,:));
hessian_2 = squeeze(abs(results.analysis.Hessian.xx));
load('evolution_test_moran_results_chicken_1E-4.mat')
final_covariance_3 = sqrt(squeeze(results.stepbysteparray(results.parameters.num_epochs,4,:)));
empirical_rho_3 = (0.5-results.analysis.rho_empirical)';
grad_3 = abs(results.analysis.grad)';
drift_3 = results.parameters.genetic_drift*ones(results.parameters.num_experiments,1);
mus_3 = squeeze(results.analysis.mus(50,1,:));
clusters_3 = squeeze(results.stepbysteparray(50,5,:));
hessian_3 = squeeze(abs(results.analysis.Hessian.xx));
load('evolution_test_moran_results_chicken_1E-3.mat')
final_covariance_4 = sqrt(squeeze(results.stepbysteparray(results.parameters.num_epochs,4,:)));
empirical_rho_4 = (0.5-results.analysis.rho_empirical)';
grad_4 = abs(results.analysis.grad)';
drift_4 = results.parameters.genetic_drift*ones(results.parameters.num_experiments,1);
mus_4 = squeeze(results.analysis.mus(50,1,:));
clusters_4 = squeeze(results.stepbysteparray(50,5,:));
hessian_4 = squeeze(abs(results.analysis.Hessian.xx));
load('evolution_test_moran_results_chicken_1E-2.mat')
final_covariance_5 = sqrt(squeeze(results.stepbysteparray(results.parameters.num_epochs,4,:)));
empirical_rho_5 = (0.5-results.analysis.rho_empirical)';
grad_5 = abs(results.analysis.grad)';
drift_5 = results.parameters.genetic_drift*ones(results.parameters.num_experiments,1);
mus_5 = squeeze(results.analysis.mus(50,1,:));
clusters_5 = squeeze(results.stepbysteparray(50,5,:));
hessian_5 = squeeze(abs(results.analysis.Hessian.xx));
load('evolution_test_moran_results_chicken_5E-3.mat')
final_covariance_6 = sqrt(squeeze(results.stepbysteparray(results.parameters.num_epochs,4,:)));
empirical_rho_6 = (0.5-results.analysis.rho_empirical)';
grad_6 = abs(results.analysis.grad)';
drift_6 = results.parameters.genetic_drift*ones(results.parameters.num_experiments,1);
mus_6 = squeeze(results.analysis.mus(50,1,:));
clusters_6 = squeeze(results.stepbysteparray(50,5,:));
hessian_6 = squeeze(abs(results.analysis.Hessian.xx));

%% Combine vectors
cmatrix_1 = cat(2,final_covariance_1,empirical_rho_1,grad_1,drift_1,mus_1,clusters_1,hessian_1,0*hessian_1);
cmatrix_2 = cat(2,final_covariance_2,empirical_rho_2,grad_2,drift_2,mus_2,clusters_2,hessian_2,0*hessian_2);
cmatrix_3 = cat(2,final_covariance_3,empirical_rho_3,grad_3,drift_3,mus_3,clusters_3,hessian_3,0*hessian_3);
cmatrix_4 = cat(2,final_covariance_4,empirical_rho_4,grad_4,drift_4,mus_4,clusters_4,hessian_4,0*hessian_4);
cmatrix_5 = cat(2,final_covariance_5,empirical_rho_5,grad_5,drift_5,mus_5,clusters_5,hessian_5,0*hessian_5);
cmatrix_6 = cat(2,final_covariance_6,empirical_rho_6,grad_6,drift_6,mus_6,clusters_6,hessian_6,0*hessian_6);
[k8,av8] = convhull(cmatrix_1(:,1),cmatrix_1(:,2));
[k9,av9] = convhull(cmatrix_2(:,1),cmatrix_2(:,2));
[k10,av10] = convhull(cmatrix_3(:,1),cmatrix_3(:,2));
[k11,av11] = convhull(cmatrix_4(:,1),cmatrix_4(:,2));
[k12,av12] = convhull(cmatrix_5(:,1),cmatrix_5(:,2));
[k13,av13] = convhull(cmatrix_6(:,1),cmatrix_6(:,2));
final_cov_analysis = cat(1,final_covariance_1,final_covariance_2,final_covariance_3,final_covariance_4,final_covariance_5,final_covariance_6);
empirical_rho_analysis = cat(1,empirical_rho_1,empirical_rho_2,empirical_rho_3,empirical_rho_4,empirical_rho_5,empirical_rho_6);
grad_analysis = cat(1,grad_1,grad_2,grad_3,grad_4,grad_5,grad_6);
drift_analysis = cat(1,drift_1,drift_2,drift_3,drift_4,drift_5,drift_6);
mus_analysis = cat(1,mus_1,mus_2,mus_3,mus_4,mus_5,mus_6);
clusters_analysis = cat(1,clusters_1,clusters_2,clusters_3,clusters_4,clusters_5,clusters_6);
hessian_analysis = cat(1,hessian_1,hessian_2,hessian_3,hessian_4,hessian_5,hessian_6);
zerovector = 0*hessian_analysis;
chicken_data_frame = cat(2,final_cov_analysis,empirical_rho_analysis,grad_analysis,drift_analysis,mus_analysis,clusters_analysis,hessian_analysis,zerovector);
chicken_data_frame(:,8) = chicken_data_frame(:,7).^2/chicken_data_frame(:,3).^2*chicken_data_frame(:,1).^2;

%% Clean results for one final cluster only
idx = chicken_data_frame(:,6)==1;
cleaned_data = chicken_data_frame(idx,:);

[chicken_y_max,chicken_max_index] = max(chicken_data_frame(:,2));
chicken_x_max = (chicken_data_frame(chicken_max_index,1));
[chicken_x_min,chicken_min_index] = min(chicken_data_frame(:,1));
chicken_y_sqr = chicken_y_max*(chicken_x_min/chicken_x_max)^2;
chicken_y_quart = chicken_y_max*(chicken_x_min/chicken_x_max)^4;

%% Plot results
figure(4)
clf
hold on
scatter(chicken_data_frame(:,1),chicken_data_frame(:,2),10,chicken_data_frame(:,4));
colorbar;
grid on
set(gca,'xscale','log','yscale','log','ColorScale','log')
title('Concentration vs. $\rho$ Predicted (Drift coloring, Chicken)','FontSize',16,'interpreter','latex')
xlabel('$\sqrt{\frac{1}{T}E[||X - \bar{x}||^2]}$','FontSize',16,'interpreter','latex')
ylabel('Empirical $0.5 - \rho$','FontSize',16,'interpreter','latex')
axis square
drawnow

figure(5)
clf
hold on
%scatter(chicken_data_frame(:,1),chicken_data_frame(:,2),10,chicken_data_frame(:,8));
scatter(cleaned_data(:,1),cleaned_data(:,2),10,cleaned_data(:,8),'filled'); % keeps one cluster only
plot([chicken_x_min,chicken_x_max],[chicken_y_sqr,chicken_y_max],'k');
plot([chicken_x_min,chicken_x_max],[chicken_y_quart,chicken_y_max],'k');
colorbar;
grid on
set(gca,'xscale','log','yscale','log','Colorscale','log','FontSize',24)
title('Chicken','FontSize',32,'interpreter','latex')
xlabel('Concentration','FontSize',32,'interpreter','latex')
axis square
drawnow

figure(6)
clf
hold on
scatter(cleaned_data(:,1),cleaned_data(:,2),10,cleaned_data(:,8),'filled');
colorbar;
grid on
set(gca,'xscale','log','yscale','log','Colorscale','log')
title('Concentration vs. $\rho$ Predicted ($\mu$ coloring, Chicken)','FontSize',16,'interpreter','latex')
xlabel('$\sqrt{\frac{1}{T}E[||X - \bar{x}||^2]}$','FontSize',16,'interpreter','latex')
ylabel('Empirical $0.5 - \rho$','FontSize',16,'interpreter','latex')
axis square
drawnow

%% Load data from PD trials
load('evolution_test_moran_results_PD_comptraits.mat')
final_covariance_1 = sqrt(squeeze(results.stepbysteparray(results.parameters.num_epochs,4,:)));
empirical_rho_1 = (0.5-results.analysis.rho_empirical)';
grad_1 = abs(results.analysis.grad)';
hessian_1 = squeeze(abs(results.analysis.Hessian.xx));
drift_1 = results.parameters.genetic_drift*ones(results.parameters.num_experiments,1);
mus_1 = squeeze(results.analysis.mus(50,1,:));
clusters_1 = squeeze(results.stepbysteparray(50,5,:));
load('evolution_test_moran_results_PD_highGD.mat')
final_covariance_2 = sqrt(squeeze(results.stepbysteparray(results.parameters.num_epochs,4,:)));
empirical_rho_2 = (0.5-results.analysis.rho_empirical)';
grad_2 = abs(results.analysis.grad)';
drift_2 = results.parameters.genetic_drift*ones(results.parameters.num_experiments,1);
mus_2 = squeeze(results.analysis.mus(50,1,:));
clusters_2 = squeeze(results.stepbysteparray(50,5,:));
hessian_2 = squeeze(abs(results.analysis.Hessian.xx));
load('evolution_test_moran_results_PD_midGD.mat')
final_covariance_3 = sqrt(squeeze(results.stepbysteparray(results.parameters.num_epochs,4,:)));
empirical_rho_3 = (0.5-results.analysis.rho_empirical)';
grad_3 = abs(results.analysis.grad)';
drift_3 = results.parameters.genetic_drift*ones(results.parameters.num_experiments,1);
mus_3 = squeeze(results.analysis.mus(50,1,:));
clusters_3 = squeeze(results.stepbysteparray(50,5,:));
hessian_3 = squeeze(abs(results.analysis.Hessian.xx));
load('evolution_test_moran_results_PD_mhighGD.mat')
final_covariance_4 = sqrt(squeeze(results.stepbysteparray(results.parameters.num_epochs,4,:)));
empirical_rho_4 = (0.5-results.analysis.rho_empirical)';
grad_4 = abs(results.analysis.grad)';
drift_4 = results.parameters.genetic_drift*ones(results.parameters.num_experiments,1);
mus_4 = squeeze(results.analysis.mus(50,1,:));
clusters_4 = squeeze(results.stepbysteparray(50,5,:));
hessian_4 = squeeze(abs(results.analysis.Hessian.xx));
load('evolution_test_moran_results_PD_mmhighGD.mat')
final_covariance_5 = sqrt(squeeze(results.stepbysteparray(results.parameters.num_epochs,4,:)));
empirical_rho_5 = (0.5-results.analysis.rho_empirical)';
grad_5 = abs(results.analysis.grad)';
drift_5 = results.parameters.genetic_drift*ones(results.parameters.num_experiments,1);
mus_5 = squeeze(results.analysis.mus(50,1,:));
clusters_5 = squeeze(results.stepbysteparray(50,5,:));
hessian_5 = squeeze(abs(results.analysis.Hessian.xx));
load('evolution_test_moran_results_PD_lowGD.mat')
final_covariance_6 = sqrt(squeeze(results.stepbysteparray(results.parameters.num_epochs,4,:)));
empirical_rho_6 = (0.5-results.analysis.rho_empirical)';
grad_6 = abs(results.analysis.grad)';
drift_6 = results.parameters.genetic_drift*ones(results.parameters.num_experiments,1);
mus_6 = squeeze(results.analysis.mus(50,1,:));
clusters_6 = squeeze(results.stepbysteparray(50,5,:));
hessian_6 = squeeze(abs(results.analysis.Hessian.xx));

final_cov_analysis = cat(1,final_covariance_1,final_covariance_2,final_covariance_3,final_covariance_4,final_covariance_5,final_covariance_6);
empirical_rho_analysis = cat(1,empirical_rho_1,empirical_rho_2,empirical_rho_3,empirical_rho_4,empirical_rho_5,empirical_rho_6);
grad_analysis = cat(1,grad_1,grad_2,grad_3,grad_4,grad_5,grad_6);
drift_analysis = cat(1,drift_1,drift_2,drift_3,drift_4,drift_5,drift_6);
mus_analysis = cat(1,mus_1,mus_2,mus_3,mus_4,mus_5,mus_6);
clusters_analysis = cat(1,clusters_1,clusters_2,clusters_3,clusters_4,clusters_5,clusters_6);
hessian_analysis = cat(1,hessian_1,hessian_2,hessian_3,hessian_4,hessian_5,hessian_6);
zerovector = 0*hessian_analysis;
pd_data_frame = cat(2,final_cov_analysis,empirical_rho_analysis,grad_analysis,drift_analysis,mus_analysis,clusters_analysis,hessian_analysis,zerovector);
pd_data_frame(:,8) = pd_data_frame(:,7).^2/pd_data_frame(:,3).^2*pd_data_frame(:,1).^2;
%[k14,av14] = convhull(pd_data_frame(:,1),pd_data_frame(:,2));
%f = fit(pd_data_frame(k14,1),pd_data_frame(k14,2),'poly1');

[pd_y_max,pd_max_index] = max(pd_data_frame(:,2));
pd_x_max = (pd_data_frame(pd_max_index,1));
[pd_x_min,pd_min_index] = min(pd_data_frame(:,1));
pd_y_sqr = pd_y_max*(pd_x_min/pd_x_max)^2;
pd_y_quart = pd_y_max*(pd_x_min/pd_x_max)^4;
figure(7)
clf
hold on
scatter(pd_data_frame(:,1),pd_data_frame(:,2),8,pd_data_frame(:,3));
plot([pd_x_min,pd_x_max],[pd_y_sqr,pd_y_max],'k');
plot([pd_x_min,pd_x_max],[pd_y_quart,pd_y_max],'k');
%plot(f,pd_data_frame(k14,1),pd_data_frame(k14,2));
colormap parula
colorbar;
grid on
set(gca,'xscale','log','yscale','log','FontSize',24)
title('PD','FontSize',32,'interpreter','latex')
xlabel('Concentration','FontSize',32,'interpreter','latex')
ylabel('Empirical $0.5 - \rho$','FontSize',32,'interpreter','latex')
axis square
drawnow

%% Load data for random performance functions
load('Test_1_Results.mat')
rho_1_rpf = results.analysis.rho;
rho_1_empirical_rpf = 0.5 - results.analysis.rho_empirical;
rho_1_covariance_rpf = sqrt(results.covariance.final);
load('Test_2_Results_small_drift.mat')
rho_2_rpf = results.analysis.rho;
rho_2_empirical_rpf = 0.5 - results.analysis.rho_empirical;
rho_2_covariance_rpf = sqrt(results.covariance.final);
load('Test_3_Results_big_drift.mat')
rho_3_rpf = results.analysis.rho;
rho_3_empirical_rpf = 0.5 - results.analysis.rho_empirical;
rho_3_covariance_rpf = sqrt(results.covariance.final);
rho_rpf = cat(1,rho_1_rpf,rho_2_rpf,rho_3_rpf);
rho_empirical_rpf = cat(1,rho_1_empirical_rpf,rho_2_empirical_rpf,rho_3_empirical_rpf);
rho_covariance_rpf = cat(1,rho_1_covariance_rpf,rho_2_covariance_rpf,rho_3_covariance_rpf);
rpf_data_frame = cat(2,rho_rpf,rho_empirical_rpf,rho_covariance_rpf);

%% Calculate max
[rpf_y_max,rpf_max_index] = max(rpf_data_frame(:,2));
rpf_x_max = (rpf_data_frame(rpf_max_index,1));
[rpf_x_min,rpf_min_index] = min(rpf_data_frame(:,3));
rpf_y_sqr = rpf_y_max*(rpf_x_min/rpf_x_max)^2;
rpf_y_quart = rpf_y_max*(rpf_x_min/rpf_x_max)^4;

%% Graph
figure(8)
clf
hold on
colormap parula
%cbh = colorbar ; %Create Colorbar
colorbar
cbh.Ticks = logspace(10^-2, 1, 8) ; %Create 8 ticks from zero to 1
cbh.TickLabels = [10^-2 5*10^-1 10^-1] ;    %Replace the labels of these 8 ticks with the numbers 1 to 8
scatter3(0.5 - rho_rpf,rho_empirical_rpf,sqrt(rho_covariance_rpf/results.parameters.num_traits),...
   15,sqrt(rho_covariance_rpf),'fill')
plot3([0.9*min(rho_empirical_rpf),1],[0.9*min(rho_empirical_rpf),1],[10^-2,10^(-2)],'k-','Linewidth',1)
scatter((rho_covariance_rpf/results.parameters.num_traits),rho_empirical_rpf,15,rho_covariance_rpf);
plot([rpf_x_min,rpf_x_max],[rpf_y_sqr,rpf_y_max],'k');
plot([rpf_x_min,rpf_x_max],[rpf_y_quart,rpf_y_max],'k');
plot([0.9*min(0.5-rho_empirical_rpf),1],[0.9*min(0.5-rho_empirical_rpf),1],[10^-2,10^(-2)],'k-','Linewidth',1)
grid on
set(gca,'xscale','log','yscale','log','zscale','log','ColorScale','log','FontSize',24)
xlim([10^-7 0.5]);
xticks([10^-6 10^-4 10^-2])
ylim([10^-7 0.5]);
title('$\rho$ Predicted vs. Actual','FontSize',32,'interpreter','latex')
xlabel('Concentration','FontSize',32,'interpreter','latex')
ylabel('Empirical $0.5 - \rho$','FontSize',32,'interpreter','latex')
zlabel('$\sqrt{\frac{1}{T}E[||X - \bar{x}||^2]}$','FontSize',32,'interpreter','latex')
axis square
drawnow