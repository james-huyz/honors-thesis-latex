%% ******************************* STEP 02 ********************************
% Identify the coordinates of second analyst's traced craters; generate
% list of craters in original traced set that are close to those in second;
% if there are multiple craters in this list choose the one with the
% highest resemblance; check how different the craters are in terms of
% azimuth (add this to an array); make an array with crater azimuth
% difference.
%                                                            James Hu, 2022
% *************************************************************************

e_crit = 1.05;

for mode = 1:2
    if mode == 1
        load(sprintf('craters_obs/RetraceKite%.2f.mat',e_crit),'a')
    elseif mode == 2
        load(sprintf('craters_obs/TraceHolo%.2f.mat',e_crit),'a')
    end
    
    b = a;
    
    load(sprintf('craters_obs/lAv%.2f.mat',e_crit))
    lAv = a;
    load(sprintf('craters_obs/AHv%.2f.mat',e_crit))
    AHv = a;
    a = cat(1,lAv,AHv);
    
    x = 0.1;
    y = 0.1;
    
    a1 = zeros(size(b));
    
    for i = 1:length(b)
        a0 = a(abs(a(:,1) - b(i,1)) <= x & abs(a(:,5) - b(i,5)) <= x &...
               abs(a(:,3) - b(i,3)) <= y, ...
             :);
        
        if isempty(a0)
            disp('Analog not found for crater.');
        elseif height(a0) > 1
            disp('Multiple analogs found.');
        else
            disp('Unique analog found.');
            a1(i,:) = a0;
        end
    end
    
    % Find diffs
    azimuth_diff = a1(:,2) - b(:,2);
    
    % Filter for where unique analog found
    azimuth_diff = azimuth_diff(find(a1(:,2)),:);
    
    % Duplicate and get opposite diffs
    azimuth_diff = cat(1,azimuth_diff,azimuth_diff.*-1);
    
    % Export
    if mode == 1
        fname = sprintf('iaerr_perts/iaerr_Kite%.2f.mat', e_crit);
        save(fname,'azimuth_diff')
    elseif mode == 2
        fname = sprintf('iaerr_perts/iaerr_Holo%.2f.mat', e_crit);
        save(fname,'azimuth_diff')
    end
end