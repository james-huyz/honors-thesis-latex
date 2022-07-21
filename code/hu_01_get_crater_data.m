%% ******************************* STEP 01 ********************************
% Read crater coordinate data from shapefiles and trace ellipses using
% Gal (2003)'s script. Export data as in five columns (latitude, azimuth,
% diameter, ellipticity, longitude).
%                                                            James Hu, 2022
% *************************************************************************

for option = 1:3
    if option == 1
        shapefile = 'traces/lAv.shp';
    elseif option == 2
        shapefile = 'traces/AHv.shp';
    elseif option == 3
        shapefile = 'traces/RetraceKite.shp';
    end
    
    S = shaperead(shapefile);
    
    a = zeros(length(S),4);
    R = 3389.5;
    digits(32);
    
    for i = 1:length(S)
        long = S(i).X;
        lat = S(i).Y;
        
        long = long(~isnan(long));
        lat = lat(~isnan(lat));
        
        cen_long = fit_ellipse(long,lat).X0_in;
        
        if isempty(cen_long)
            a(i,:) = [NaN,NaN,NaN,NaN,NaN];
            continue;
        else
            ellipse = fit_ellipse((long-cen_long).*cosd(lat)+cen_long,lat);
            if isempty(ellipse.status)
                % Latitude
                a(i,1) = ellipse.Y0_in;
                % Azimuth
                a(i,2) = 90-abs(rad2deg(ellipse.phi));
                % Diameter
                a(i,3) = 2*sqrt(2*pi*R/360*0.5*ellipse.long_axis...
                         *2*pi*R/360*0.5*ellipse.short_axis);
                % Ellipticity
                a(i,4) = ellipse.long_axis/ellipse.short_axis;
                % Longitude
                a(i,5) = ellipse.X0_in;
            else
                a(i,:) = [NaN,NaN,NaN,NaN,NaN];
            end
        end
    end
    
    % Filter for craters with ellipticity > e_crit
    for e_crit = 1:0.05:1.15
        if e_crit > 1
            a = a(a(:,4) >= e_crit,:);
        end

        if option == 1
            fname = sprintf('craters_obs/lAv%.2f.mat', e_crit);
            save(fname,'a')
        elseif option == 2
            fname = sprintf('craters_obs/AHv%.2f.mat', e_crit);
            save(fname,'a')
        elseif option == 3
            fname = sprintf('craters_obs/RetraceKite%.2f.mat', e_crit);
            save(fname,'a')
        end
    end
end

% Filter for craters with ellipticity > e_crit
load('craters_obs/TraceHolo1.00.mat')

for e_crit = 1.05:0.05:1.15
    a = a(a(:,4) >= e_crit,:);
    fname = sprintf('craters_obs/TraceHolo%.2f.mat', e_crit);
    save(fname,'a')
end