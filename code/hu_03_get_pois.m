%% ******************************* STEP 03 ********************************
% This is the first script in Sam Holo's forward model for the effect of
% obliquity on elliptic crater orientations. The input file is name
% '14mag_imp_info' and contains an ensemble of encounter inclinations and
% speeds from n-body output. The output of this script is an ensemble of
% elliptic crater locations and diameters, as well as impactor velocity
% vectors at the time of impact and the impact angle. This output can be
% ingested and processed by apply_obliquity.m.
%                                                            Sam Holo, 2018
% *************************************************************************

% Set variables for physical constants
G = 6.67408e-11;
R = 3.3895e6;
M = 6.4171e23;
vesc = 5027;
g = 3.71;

% Load impact info file - This needs two values: encounter
% inclinations and speeds from n-body output
load('forward_model/14mag_imp_info.mat','inclinations','speeds');

% Set parameters for the size frequency distribution
% NOTE: the SFD parameters have been tuned for a particular impactor
% population, this will NOT immediately translate to different populations
lmin = 21;
alpha = 0.65;
n = 1e6;

% Pre-assign random numbers and place for impact info
[speeds,id] = datasample(speeds,n);
inclinations = inclinations(id);
impactinfo = zeros(n,9);
randos = rand(n,4);

% Loop through each of the seeded inclination-speed pairs 
for i = 1:n
    % set the impact parameters 
    speed = speeds(i);
    tau = R*sqrt(1 + (2*G*M)/(R*speed*speed));
    inc = inclinations(i);

    % Sample from our assumed SFD and apply crater diameter scaling prior
    % to correction for impact angle
    prob = randos(i,4);
    imp_diam = lmin*(prob)^(1/(-1*alpha)); impdiamstore(i) = imp_diam;
    d90 = 1.161*(g^(-0.22))*(speed^(0.44))*(imp_diam^(0.78));
    

    % Calculate the critical angle from cratering efficiency
    % Later we will "flag" if the generated crater is elliptical or not
    theta_c = 45*(d90/imp_diam)^(-0.52) + 77*(d90/imp_diam)^(-1.85);
    
    % Here we create impacts with angles < 45 degrees
    bsq = (tau^2)*(0.5 + 0.5*randos(i,1)); 
    b = sqrt(bsq);
    theta = acos(b/tau);
    del = randos(i,2)*2*pi;
    phi = randos(i,3)*2*pi;

    % Semi major axis of hyperbolic orbit
    a = G*M/(speed^2);

    % Eccentricity
    ecc = sqrt(1 + (b^2)*(speed^4)/ ((G^2)*(M^2)) );

    % Cosine and sine of true anomaly at launch point
    cinf = -1/ecc;
    sinf = -1*sqrt(1 - (cinf^2));

    % Cosine and sine of true anomaly at impact
    cimp = (((a/R)*(ecc^2 - 1))-1)/ecc;
    simp = -1*sqrt(1 - (cimp^2));

    % Calculate impact location in projectile coordinates
    yi = -R*(cimp*cinf + simp*sinf);
    xi = R*(simp*cinf - cimp*sinf);
    zi = 0;
    pos = [xi;yi;zi];
    rotangle = pi/2 + theta;

    % Calculate impact velocities from momentum conservation
    speed = sqrt((speed)^2 + vesc^2);
    vimp = speed*[cos(rotangle),-1*sin(rotangle),0; sin(rotangle),...
           cos(rotangle),0; 0,0,1]*pos/norm(pos);

    % Rotatation matrix for the effect of impact argument
    cd = cos(del);
    sd = sin(del);
    Mdel = [cd,0,-1*sd; 0,1,0; sd,0,cd];

    % Rotatation matrix for the effect of impactor inclination
    ci = cos(inc);
    si = sin(inc);
    Minc = [1,0,0;0,ci,si;0,-1*si,ci];

    % Rotatation matrix for the effect of randomized precessional season
    sphi = sin(phi);
    cphi = cos(phi);
    Mprec = [cphi,-1*sphi,0;sphi,cphi,0;0,0,1];

    % Update position and velocity by applying rotations
    pos = (Mprec*(Minc*(Mdel*pos)));
    vimp = (Mprec*(Minc*(Mdel*vimp)));

    % Store position and velocity vectors, along with impact angles
    impactinfo(i,1:6) = [pos;vimp]';
    impactinfo(i,7) = theta;
    impactinfo(i,9) = deg2rad(theta_c);

    % Correct for the effect of impact angle in crater diameter and store
    crater_diam = d90*(sin(theta)^(1/3))/1000;  % Remember scaling is in m!
    impactinfo(i,8) = crater_diam;
end

impactinfo = impactinfo(impactinfo(:,7) < impactinfo(:,9),:);

% Save the data
newfile = 'forward_model/14mag_pois.mat';
save(newfile,'impactinfo','-v7.3');