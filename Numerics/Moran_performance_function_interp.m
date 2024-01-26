function performance = Moran_performance_function_interp(x,y,dataset)
%% Load dataset
load(dataset);

%% Get size
[n,~] = size(log_odds_moran);

%% Interpolate grid
performance = interp2(log_odds_moran,y*(n-1)+1,x*(n-1)+1,'spline');%Can change method to cubic, makima, spline, default is linear