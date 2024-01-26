function [f] = performance_stag_moran_24_individuals(X,Y,mode)
%% Inputs
% 1. x and y are column vectors storing traits (probability of hunting the stag)
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
        %Coefficients (with 95% confidence bounds):
        p00 =   0;
        p10 =       1.475;  %(1.441, 1.509)
        p01 =      -1.475;  %(-1.509, -1.441)
        p20 =      -3.206;  %(-3.273, -3.139)
        p11 =  0;
        p02 =       3.206;  %(3.139, 3.273)
        p30 =      0.9129;  %(0.8705, 0.9553)
        p21 =       1.096;  %(1.059, 1.133)
        p12 =      -1.096;  %(-1.133, -1.059)
        p03 =     -0.9129;  %(-0.9553, -0.8705)
        
        %Linear model Poly33:
        f(k) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y ...
        + p12*x*y^2 + p03*y^3;
        
        % Goodness of fit:
        %   SSE: 0.04099
        %   R-square: 0.9992
        %   Adjusted R-square: 0.9992
        %   RMSE: 0.009752
        
    elseif strcmp(mode,'quintic')
        %% quintic
        %Coefficients (with 95% confidence bounds):
       p00 =  0;
       p10 =       1.443;  %(1.322, 1.564)
       p01 =      -1.443;  %(-1.564, -1.322)
       p20 =      -2.898;  %(-3.511, -2.284)
       p11 =  0;
       p02 =       2.898;  %(2.284, 3.511)
       p30 =       0.214;  %(-1.196, 1.624)
       p21 =     -0.3664;  %(-1.381, 0.6484)
       p12 =      0.3664;  %(-0.6484, 1.381)
       p03 =      -0.214;  %(-1.624, 1.196)
       p40 =      0.6237;  %(-0.8628, 2.11)
       p31 =       2.228;  %(1.124, 3.331)
       p22 =  0;
       p13 =      -2.228;  %(-3.331, -1.124)
       p04 =     -0.6237;  %(-2.11, 0.8628)
       p50 =     -0.2054;  %(-0.7887, 0.3779)
       p41 =     -0.8083;  %(-1.302, -0.3146)
       p32 =     -0.8159;  %(-1.291, -0.3408)
       p23 =      0.8159;  %(0.3408, 1.291)
       p14 =      0.8083;  %(0.3146, 1.302)
       p05 =      0.2054;  %(-0.3779, 0.7887)

        %Linear model Poly55:
     f(k) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y ... 
                    + p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y + p22*x^2*y^2 ... 
                    + p13*x*y^3 + p04*y^4 + p50*x^5 + p41*x^4*y + p32*x^3*y^2 ... 
                    + p23*x^2*y^3 + p14*x*y^4 + p05*y^5;

% Goodness of fit:
%   SSE: 0.03349
%   R-square: 0.9994
%   Adjusted R-square: 0.9993
%   RMSE: 0.008929


        
        
    end
end