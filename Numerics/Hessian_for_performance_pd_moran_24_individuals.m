function [g,H] = Hessian_for_performance_pd_moran_24_individuals(x,mode)
%% Inputs
% 1. x is the trait (probability of defecting)
% 2. mode: either 'cubic' or 'quintic'

%% outputs
% 1. g: gradient
% 2. H: Hessian

%% set y = x
y = x;

    
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
    
    % gradient w.r.t x:
    g = p10 + 2*p20*x + p11*y + 3*p30*x^2 + 2*p21*x*y + p12*y^2;
    
    % Hessian w.r.t xx and xy:
     H.xx = 2*p20 + 6*p30*x + 2*p21*y;
     H.xy = p11 + 2*p21*x + 2*p12*y;
    
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
    
    
    % gradient w.r.t x:
    g = p10 + ...
        2*p20*x + p11*y + ... 
        3*p30*x^2 + 2*p21*x*y + p12*y^2 + ...
        4*p40*x^3 + 3*p31*x^2*y + 2*p22*x*y^2 + p13*y^3 + ...
        5*p50*x^4 + 4*p41*x^3*y + 3*p32*x^2*y^2 + 2*p23*x*y^3 + p14*y^4;
    
    % Hessian w.r.t xx and xy:
     H.xx = 2*p20 + ... 
        6*p30*x + 2*p21*y  + ...
        12*p40*x^2 + 6*p31*x*y + 2*p22*y^2 + ...
        20*p50*x^3 + 12*p41*x^2*y + 6*p32*x*y^2 + 2*p23*y^3;
     H.xy = p11 + ... 
        2*p21*x + 2*p12*y + ...
        3*p31*x^2 + 4*p22*x*y + 3*p13*y^2 + ...
        4*p41*x^3 + 6*p32*x^2*y + 6*p23*x*y^2 + 4*p14*y^3;
    
    %Goodness of fit:
    %  SSE: 0.4066
    %  R-square: 0.9997
    %  Adjusted R-square: 0.9997
    %  RMSE: 0.03112
    
end
