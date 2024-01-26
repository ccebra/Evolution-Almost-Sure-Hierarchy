function [f] = performance_pd_moran_24_individuals(X,Y,mode)
%% Inputs
% 1. x and y are column vectors storing traits (probability of defecting)
% 2. mode: either 'cubic' or 'quintic'

%% outputs
% 1. f: the performance
% 2. g: gradient
% 3. H: Hessian

%% get dimensions and preallocate
[pairs,~] = size(X);
f = nan([pairs,1]);

%% Loop over input pairs
for k = 1:pairs
    x = X(k,:);
    y = Y(k,:);
    
    %% cubic
    if strcmp(mode,'cubic')
        % Coefficients (with 95% confidence bounds):
        p00 =   0;
        p10 =   -5.72;  %(-5.841, -5.599)
        p01 =    5.72;  %(5.599, 5.841)
        p20 =    2.623;  %(2.383, 2.862)
        p11 =    0;
        p02 =    -2.623;  %(-2.862, -2.383)
        p30 =    -1.509;  %(-1.66, -1.358)
        p21 =     2.277;  %(2.146, 2.409)
        p12 =    -2.277;  %(-2.409, -2.146)
        p03 =     1.509;  %(1.358, 1.66)
        
        % Linear model Poly33:
        f(k) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y +...
            p02*y^2 + p30*x^3 + p21*x^2*y + p12*x*y^2 + p03*y^3;
        
        % Goodness of fit:
        %  SSE: 0.5212
        %  R-square: 0.9996
        %  Adjusted R-square: 0.9996
        %  RMSE: 0.03477
    elseif strcmp(mode,'quintic')
        %% quintic
        
        % Coefficients (with 95% confidence bounds):
        p00 =  0;
        p10 =      -5.732;  %(-6.155, -5.309)
        p01 =       5.732;  %(5.309, 6.155)
        p20 =       3.369;  %(1.231, 5.508)
        p11 =  0;  
        p02 =      -3.369;  %(-5.508, -1.231)
        p30 =      -3.993;  %(-8.908, 0.9211)
        p21 =       3.794;  %(0.2581, 7.33)
        p12 =      -3.794;  %(-7.33, -0.2581)
        p03 =       3.993;  %(-0.9211, 8.908)
        p40 =       2.376;  %(-2.804, 7.556)
        p31 =      -1.272;  %(-5.116, 2.572)
        p22 =  0;
        p13 =       1.272;  %(-2.572, 5.116)
        p04 =      -2.376;  %(-7.556, 2.804)
        p50 =      -0.661;  %(-2.694, 1.372)
        p41 =      0.7926;  %(-0.9276, 2.513)
        p32 =      -1.678;  %(-3.334, -0.02282)
        p23 =       1.678;  %(0.02282, 3.334)
        p14 =     -0.7926;  %(-2.513, 0.9276)
        p05 =       0.661;  %(-1.372, 2.694)
        
        
        f(k) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y + ...
        p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y + p22*x^2*y^2 + ...
        p13*x*y^3 + p04*y^4 + p50*x^5 + p41*x^4*y + p32*x^3*y^2 + ...
        p23*x^2*y^3 + p14*x*y^4 + p05*y^5;
        
        %Goodness of fit:
        %  SSE: 0.4066
        %  R-square: 0.9997
        %  Adjusted R-square: 0.9997
        %  RMSE: 0.03112
        
    end
end