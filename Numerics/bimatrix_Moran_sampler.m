%% Parameters
num_moran = 24;% number of entrants in Moran model
num_interpolation = 50000;% number of Moran games to play for grid interpolation
grid_width = 20;% grid width for the Moran interpolation
start_row = 6;% if computer crashes part of the way through, can start at a later row
current_cell = 0;% tracks current cell

payout = cell(2,2);
payout{1,1} = [2,2];
payout{1,2} = [3,0];
payout{2,1} = [0,3];
payout{2,2} = [1,1];


%% grid trait space for performance function interpolation
log_odds_moran = NaN(grid_width,grid_width);
for i = start_row:grid_width
    for j = i:grid_width
        current_cell = current_cell + 1;
        num_i_wins = 0;
        parfor k = 1:num_interpolation
            if i == j
                num_i_wins = num_i_wins+1/2;
            else
                moran_result = moran(payout,i/grid_width,j/grid_width,num_moran);
                if moran_result == i/grid_width
                    num_i_wins = num_i_wins+1;
                end
            end
        end
        log_odds_moran(i+1,j+1) = num_i_wins;
        log_odds_moran(j+1,i+1) = num_interpolation - num_i_wins;
        save('log_odds_moran_temp.mat','log_odds_moran');
        fprintf('\n Trial %d of %d complete \n',i*grid_width+j,1/2*grid_width*(grid_width+1));
    end
end
log_odds_moran = log((log_odds_moran + 1)./(log_odds_moran' + 1));%Log odds

%% Perform interpolation

save('log_odds_moran.mat','log_odds_moran');