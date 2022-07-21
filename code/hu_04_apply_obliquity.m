%% ******************************* STEP 04 ********************************
% This is the second script in Sam Holo's forward model for the effect of
% obliquity on elliptic crater orientations. It takes the output ensemble
% from get_pois.m (get "pre-obliquity-impacts"), applies each integer
% degree obliquity to the ensemble, and saves the key parameters: latitude,
% diameter, and orientation of the elliptic craters.
%                                                            Sam Holo, 2018
% *************************************************************************

% Load the impact data from get_pois.m
load('forward_model/14mag_pois.mat','impactinfo');

% Set the diameter (in km) below which you throw out model craters
diam_cutoff = 4;

% Filter out craters that are too small for the analysis
R = 3.3895e6;
impactinfo = impactinfo(impactinfo(:,8)>= diam_cutoff,:);
len = length(impactinfo);

% Pre-allocate space for the centers/counts of the model output
latdata = zeros(len,91);
orientationdata = zeros(len,91);
diamdata = zeros(len,91);

% Loop through the possible values for obliquity
for i = 0:90
    % Set the obliquity NOTE: here we are working in degrees
    obl = i;

    % Loop through each elliptic-crater producing impact
    for j = 1:len
        % Define the obliquity rotation matrix
        cobl = cosd(obl);
        sobl = sind(obl);
        Mobl = [1,0,0;0,cobl,sobl;0,-1*sobl,cobl];

        % Apply the obliquity rotation to the position and velocity vectors
        pos = Mobl*impactinfo(j,1:3)';
        vimp = Mobl*impactinfo(j,4:6)';

        % Calculate latitude
        latdata(j,i+1) = abs(asind(pos(3)/norm(pos)));

        % Calculate the impact orientation by projecting northbound and
        % velocity vectors onto the planet's tangent plane
        poshat = pos/norm(pos);
        north = [0;0;R] - pos;
        northproj = north - dot(north,poshat)*poshat;
        trajectoryproj = vimp - dot(vimp,poshat)*poshat;

        if norm(northproj) ~= 0
            northproj = northproj/norm(northproj); % Normalize 
        end
        
        if norm(trajectoryproj) ~= 0
            trajectoryproj = trajectoryproj/norm(trajectoryproj); % N'lize
        end
        
        azimuth = acosd(dot(northproj,trajectoryproj));
        azimuth(azimuth > 90) = 180-azimuth;
        orientationdata(j,i+1) = azimuth;
        diamdata(j,i+1) = impactinfo(j,8);
    end
    % Display progress
    disp(i);
end

% Save the data
newfile = 'forward_model/obpreds4.mat';
save(newfile,'latdata','orientationdata','diamdata');