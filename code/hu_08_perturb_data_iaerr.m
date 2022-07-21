%% ******************************* STEP 08 ********************************
% Create inter-analyst-perturbed crater azimuth arrays.
%                                                            James Hu, 2022
% *************************************************************************

e_crit = 1.05;
N = 1000;  % Number of instances of inter-analyst perturbation

% Load crater data and inter-analyst perturbation options
for mode = 1:2  % 1 = lAv; 2 = AHv
    
    if mode == 1
        load(sprintf('craters_obs/lAv%.2f.mat',e_crit))
    elseif mode == 2
        load(sprintf('craters_obs/AHv%.2f.mat',e_crit))
    end
    
    load(sprintf("iaerr_perts/iaerr_Kite%.2f.mat",e_crit))
    kite = azimuth_diff;
    load(sprintf("iaerr_perts/iaerr_Holo%.2f.mat",e_crit))
    holo = azimuth_diff;
    azimuth_diff = cat(1,kite,holo);
    
    cratdata = a(:,2);
    clear('a')
    
    % Set random seed
    rng('default');
    rng(1);
    
    % Create perturbation array
    l = length(cratdata);
    perts = reshape(datasample(azimuth_diff,N*l),[l,N]);
    
    % Perturb
    cratdata_p = zeros(l,N);
    for i = 1:N
        cratdata_p(:,i) = cratdata + perts(:,i);
        
        % Refold azimuths above 90° and below 0°
        for j = 1:l
            if (cratdata_p(j,i) > 90 || cratdata_p(j,i) < 0)
                cratdata_p(j,i) = 90 - mod(cratdata_p(j,i),90);
            end
        end
    end
    
    if mode == 1
        fname = sprintf('iaerr_perts/lAv%.2f_pert%d.mat',e_crit,N);
        save(fname,'cratdata_p')
    elseif mode == 2
        fname = sprintf('iaerr_perts/AHv%.2f_pert%d.mat',e_crit,N);
        save(fname,'cratdata_p')
    end
end