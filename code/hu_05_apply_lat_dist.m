%% ******************************* STEP 05 ********************************
% Modify the output of Sam Holo's forward model to conform to the
% latitudinal distribution of actually observed craters.
%                                                            James Hu, 2022
% *************************************************************************

for e_crit = 1.00:0.05:1.15
    for mode = 1:2  % 1 = lAv; 2 = AHv
        load('forward_model/obpreds4.mat','orientationdata','latdata',...
             'diamdata');
        
        if mode == 1
            load(sprintf('craters_obs/lAv%.2f.mat',e_crit));
        elseif mode == 2
            load(sprintf('craters_obs/AHv%.2f.mat',e_crit));
        end
        
        % Fold observed latitudes by assumption of N/S hemispheric symmetry
        folded_lats = abs(a(:,1));
        clear a;
        
        % Determine the empirical and simulated latitudinal  distribution
        emp_lat_dist = hu_05a_find_lat_dist(folded_lats);
        sim_lat_dists = zeros(17,91);
        
        for i = 0:90
            sim_lat_dists(:,i+1) = hu_05a_find_lat_dist(latdata(:,i+1));
        end
        
        % Normalize to bin with greatest simulation/empirical disparity
        scale_factors = emp_lat_dist ./ sim_lat_dists;
        scale_factors = scale_factors / max(scale_factors(:));
        
        % Get indicies by proportionally resampling sim craters wrt 5°
        % latitude bin; note: method obtains indicies before applying those
        % to lats, ors, diams
        size = 1e6;
        idx = zeros(size,91);

        for i = 0:90
            resampled = [];

            for bin = 1:17
                latdata_at_lat = latdata(:,i+1);
                id = find(latdata_at_lat > 5*(bin-1)&...
                          latdata_at_lat <= 5*(bin-1) + 5);
                sample = datasample(id,round(height(id)*...
                                    scale_factors(bin,i+1)),'Replace',...
                                    false);
                resampled = cat(1,resampled,sample);
            end
            idx(:,i+1) = datasample(resampled,size);
        end
        
        % Apply indicies to latitudes, orientations, and diameters
        resampled_lats = zeros(size,91);
        resampled_orientations = zeros(size,91);
        resampled_diams = zeros(size,91);
        
        % Loop through each obliquity
        for i = 0:90
            lat = latdata(:,i+1);
            resampled_lats(:,i+1) = lat(idx(:,i+1));
            orientation = orientationdata(:,i+1);
            resampled_orientations(:,i+1) = orientation(idx(:,i+1));
            diam = diamdata(:,i+1);
            resampled_diams(:,i+1) = diam(idx(:,i+1));
        end
        
        % Export
        latdata = resampled_lats;
        orientationdata = resampled_orientations;
        diamdata = resampled_diams;
        
        if mode == 1
            fname = sprintf('forward_model/obpreds4_norm_lAv%.2f.mat',...
                            e_crit);
        elseif mode == 2
            fname = sprintf('forward_model/obpreds4_norm_AHv%.2f.mat',...
                            e_crit);
        end    
        
        save(fname,'latdata','orientationdata','diamdata','-v7.3');
    end
end