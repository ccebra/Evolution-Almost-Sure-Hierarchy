function PlotResults_RPFSM(filename)

%% Load input data
load(filename)

column_three = squeeze(results.stepbysteparray(:,3,:))';
column_four = squeeze(results.stepbysteparray(:,4,:))';
column_five = squeeze(results.stepbysteparray(:,5,:))';
grad_norm = abs(results.analysis.grad);
grad_norm = grad_norm';
%hessian_xx_norm = abs(squeeze(results.analysis.Hessian.xx)');
%hessian_xy_norm = abs(squeeze(results.analysis.Hessian.xy)');
final_covariance = squeeze(results.stepbysteparray(results.parameters.num_epochs,4,:));
% hessian_xx_norm = norm(results.analysis.Hessian.xx,'fro')';
% hessian_xy_norm = norm(results.analysis.Hessian.xy,'fro')';


%% Manipulate data for norms of gradient and Hessian
% gradient_norm = [results.norms.gradient,results.maxcoordinate];
% hessian_xx_norm = [results.norms.xxhessian,results.maxcoordinate];
% hessian_xy_norm = [results.norms.xyhessian,results.maxcoordinate];
% gradient_interior = gradient_norm(gradient_norm(:,2)<1-1.96*results.parameters.genetic_drift, :);
% gradient_boundary = gradient_norm(gradient_norm(:,2)>=1-1.96*results.parameters.genetic_drift, :);
% hessian_xx_interior = hessian_xx_norm(hessian_xx_norm(:,2)<1-1.96*results.parameters.genetic_drift, :);
% hessian_xx_boundary = hessian_xx_norm(hessian_xx_norm(:,2)>=1-1.96*results.parameters.genetic_drift, :);
% hessian_xy_interior = hessian_xy_norm(hessian_xy_norm(:,2)<1-1.96*results.parameters.genetic_drift, :);
% hessian_xy_boundary = hessian_xy_norm(hessian_xy_norm(:,2)>=1-1.96*results.parameters.genetic_drift, :);


%% Plot intransitivity on interior and boundary

figure(2)
clf
hold on
boxplot(column_three);%Grid on, set(gca,'yscale','log'),increase font size of axes
med = median(column_three(:,[1:results.parameters.num_epochs]));
xcoords = (1:results.parameters.num_epochs);
pbaspect([1 1 1])
set(gca,'FontSize',20)
xticks([0 20 40 60 80 100 120 140 160 180 200])
xticklabels({'0', '20', '40', '60', '80', '100', '120', '140' ...
    , '160', '180', '200', '70', '75', '80', '85', '90', '95', '100'})
xtickangle(45)
plot(xcoords, med,'b-','LineWidth',1.75)
xlabel('Evolutionary Steps', 'FontSize', 36, 'interpreter', 'latex')
%ylabel('Proportion Intransitivity', 'FontSize', 36, 'interpreter', 'latex')
axis tight
title('Intransitivity Step by Step', 'FontSize', 36, 'interpreter', 'latex')
grid on
drawnow

%% Plot covariance on interior and boundary
figure(3)
clf
hold on
boxplot(column_four);%Grid on, set(gca,'yscale','log'),increase font size of axes
med = median(column_four(:,[1:results.parameters.num_epochs]));
xcoords = (1:results.parameters.num_epochs);
pbaspect([1 1 1])
set(gca,'FontSize',20)
xticks([0 20 40 60 80 100 120 140 160 180 200])
xticklabels({'0', '20', '40', '60', '80', '100', '120', '140' ...
    , '160', '180', '200', '70', '75', '80', '85', '90', '95', '100'})
xtickangle(45)
plot(xcoords, med,'b-','LineWidth',1.75)
xlabel('Evolutionary Steps', 'FontSize', 36, 'interpreter', 'latex')
ylabel('Covariance', 'FontSize', 36, 'interpreter', 'latex')
axis tight
title('Covariance Step by Step', 'FontSize', 36, 'interpreter', 'latex')
grid on
drawnow

%% Plot number of clusters on interior and boundary
figure(4)
clf
hold on
boxplot(column_five);%Grid on, set(gca,'yscale','log'),increase font size of axes
med = median(column_five(:,[1:results.parameters.num_epochs]));
xcoords = (1:results.parameters.num_epochs);
pbaspect([1 1 1])
set(gca,'FontSize',20)
xticks([0 20 40 60 80 100 120 140 160 180 200])
xticklabels({'0', '20', '40', '60', '80', '100', '120', '140' ...
    , '160', '180', '200', '70', '75', '80', '85', '90', '95', '100'})
xtickangle(45)
plot(xcoords, med,'b-','LineWidth',1.75)
xlabel('Evolutionary Steps', 'FontSize', 36, 'interpreter', 'latex')
ylabel('Number of Clusters', 'FontSize', 36, 'interpreter', 'latex')
axis tight
title('Clusters Step by Step', 'FontSize', 36, 'interpreter', 'latex')
grid on
drawnow

%% plot histograms for gradient and Hessian norms
figure(5)
clf
hist(grad_norm);
xlabel('Norm', 'FontSize', 18, 'interpreter', 'latex')
ylabel('', 'FontSize', 18, 'interpreter', 'latex')
title('Norm of the Gradient', 'FontSize', 18, 'interpreter', 'latex')
axis tight
grid on

% figure(6)
% clf
% hist(hessian_xx_norm);
% xlabel('Norm', 'FontSize', 18, 'interpreter', 'latex')
% ylabel('', 'FontSize', 18, 'interpreter', 'latex')
% title('Norm of the xx-Hessian', 'FontSize', 18, 'interpreter', 'latex')
% axis tight
% grid on
% 
% figure(7)
% clf
% hist(hessian_xy_norm);
% xlabel('Norm', 'FontSize', 18, 'interpreter', 'latex')
% ylabel('', 'FontSize', 18, 'interpreter', 'latex')
% title('Norm of the xy-Hessian', 'FontSize', 18, 'interpreter', 'latex')
% axis tight
% grid on

%% display rhos (experimental versus predicted)
% figure(8)
% clf
% hold on
% scatter(sqrt(final_covariance),0.5-results.analysis.rho_empirical,20,abs(results.analysis.grad),'filled');
% colorbar;
% grid on
% set(gca,'xscale','log','yscale','log')
% title('Concentration vs. $\rho$ Predicted','FontSize',16,'interpreter','latex')
% xlabel('$\sqrt{\frac{1}{T}E[||X - \bar{x}||^2]}$','FontSize',16,'interpreter','latex')
% ylabel('Empirical $0.5 - \rho$','FontSize',16,'interpreter','latex')
% axis square
% drawnow

%% mus analysis for chicken tests
% figure(9)
% clf
% hold on
% histogram(results.analysis.mus(50,1,:));
% title('Final $\mu$ of the competitor clusters','FontSize',16,'interpreter','latex')
% axis square
% drawnow
% 
% figure(10)
% clf
% hold on
% scatter(squeeze(results.analysis.mus(50,1,:)),grad_norm,20,squeeze(results.stepbysteparray(50,5,:)),'filled')
% colorbar;
% grid on
% title('Gradient norm against final $\mu$ of the data cluster','FontSize',16,'Interpreter','latex')
% xlabel('Final $\mu$','Interpreter','latex')
% ylabel('Norm')
% axis square
% drawnow
% 
% figure(11)
% clf
% hold on
% scatter(squeeze(results.analysis.mus(50,1,:)),hessian_xx_norm,20,squeeze(results.stepbysteparray(50,5,:)),'filled')
% colorbar;
% grid on
% title('Hessian norm against final $\mu$ of the data cluster','FontSize',16,'Interpreter','latex')
% xlabel('Final $\mu$','Interpreter','latex')
% ylabel('Norm')
% axis square
% drawnow

% % 
% % %% display to test
% % figure(12)
% % clf
% % hold on
% % plot(results.convergence_test.boundary.stds,0.5 - mean(results.convergence_test.boundary.rho_predictions,'omitnan'),'k-','Linewidth',1.5)
% % plot(results.convergence_test.boundary.stds,0.5 - results.convergence_test.boundary.rho_predictions,'Linewidth',1)
% % grid on
% % set(gca,'yscale','log','xscale','log')
% % xlabel('Standard Deviation','FontSize',16,'interpreter','latex')
% % ylabel('$0.5 - \rho$','FontSize',16,'interpreter','latex')
% % title('Predicted $\rho$ on Boundary','FontSize',16,'interpreter','latex')
% % l = legend('Average over Trials','Per Trial');
% % set(l,'FontSize',14,'interpreter','latex','location','best');
% % drawnow
% % 
% % figure(14)
% % clf
% % hold on
% % plot(results.convergence_test.boundary.stds,mean(abs(results.convergence_test.boundary.rho_mean - results.convergence_test.boundary.rho_predictions)./results.convergence_test.boundary.rho_predictions,'omitnan'),...
% %     'k-','Linewidth',2)
% % plot(results.convergence_test.boundary.stds,abs(results.convergence_test.boundary.rho_mean - results.convergence_test.boundary.rho_predictions),...
% %     'Linewidth',1)
% % grid on
% % set(gca,'yscale','log','xscale','log')
% % xlabel('Standard Deviation','FontSize',16,'interpreter','latex')
% % ylabel('Relative Error','FontSize',16,'interpreter','latex')
% % title('Error in Estimated $\rho$ on Boundary','FontSize',16,'interpreter','latex')
% % l = legend('Average over Trials','Per Trial');
% % set(l,'FontSize',14,'interpreter','latex','location','best');
% % drawnow
% % 
% % %% display to test
% % figure(13)
% % clf
% % hold on
% % plot(results.convergence_test.interior.stds,0.5 - mean(results.convergence_test.interior.rho_predictions,'omitnan'),'k-','Linewidth',1.5)
% % plot(results.convergence_test.interior.stds,0.5 - results.convergence_test.interior.rho_predictions,'Linewidth',1)
% % grid on
% % set(gca,'yscale','log','xscale','log')
% % xlabel('Standard Deviation','FontSize',16,'interpreter','latex')
% % ylabel('$0.5 - \rho$','FontSize',16,'interpreter','latex')
% % title('Predicted $\rho$ on Interior','FontSize',16,'interpreter','latex')
% % l = legend('Average over Trials','Per Trial');
% % set(l,'FontSize',14,'interpreter','latex','location','best');
% % drawnow
% % 
% % figure(15)
% % clf
% % hold on
% % plot(results.convergence_test.interior.stds,mean(abs(results.convergence_test.interior.rho_mean - results.convergence_test.interior.rho_predictions)./results.convergence_test.interior.rho_predictions,'omitnan'),...
% %     'k-','Linewidth',2)
% % plot(results.convergence_test.interior.stds,abs(results.convergence_test.interior.rho_mean - results.convergence_test.interior.rho_predictions),...
% %     'Linewidth',1)
% % grid on
% % set(gca,'yscale','log','xscale','log')
% % xlabel('Standard Deviation','FontSize',16,'interpreter','latex')
% % ylabel('Relative Error','FontSize',16,'interpreter','latex')
% % title('Error in Estimated $\rho$ on Interior','FontSize',16,'interpreter','latex')
% % l = legend('Average over Trials','Per Trial');
% % set(l,'FontSize',14,'interpreter','latex','location','best');
% % drawnow