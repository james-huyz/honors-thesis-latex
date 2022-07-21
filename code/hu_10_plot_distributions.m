%% ******************************* STEP 10 ********************************
% Plot distributions of crater latitudes, azimuths, and ellipticities.
%                                                            James Hu, 2022
% *************************************************************************

e_crit = 1.05;

% Latitude distribution
load('forward_model/obpreds4.mat','orientationdata','latdata','diamdata')

fig1 = figure;
t = tiledlayout(1,2,'TileSpacing','compact');

nexttile
load(sprintf('craters_obs/lAv%.2f.mat',e_crit))
folded_lats = abs(a(:,1));
lAv = a;
hist(folded_lats,2.5:5:87.5)
h = findobj(gca,'Type','patch');
h.FaceColor = '#74C1E4';
xlim([0 90])
ylim([0 300])
title('(a)')

nexttile
load(sprintf('craters_obs/AHv%.2f.mat',e_crit))
folded_lats = abs(a(:,1));
AHv = a;
hist(folded_lats,2.5:5:87.5)
h = findobj(gca,'Type','patch');
h.FaceColor = '#74C1E4';
xlim([0 90])
ylim([0 300])
title('(b)')

ylabel(t, 'Frequency')
xlabel(t, 'Latitude (°)')

fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0 0 24 16];
print(fig1,'output/fig4_1','-dpng','-r500')

% Azimuth distribution
a = cat(1,lAv,AHv);
azimuths = a(:,2);

fig2 = figure;

hist(azimuths,2.5:5:87.5)
h = findobj(gca,'Type','patch');
h.FaceColor = '#74C1E4';

ylabel('Frequency')
xlabel('Azimuth (°)')
xlim([0 90])

fig2.PaperUnits = 'centimeters';
fig2.PaperPosition = [0 0 24 16];
print(fig2,'output/fig4_4','-dpng','-r500')

% Ellipticity distribution

load('craters_obs/lAv1.00.mat')
lAv100 = a;
load('craters_obs/AHv1.00.mat')
AHv100 = a;
a = cat(1,lAv100,AHv100);
ellips = a(:,4);
clear('a')

fig2 = figure;

hist(ellips,1.0125:0.025:1.9875)
h = findobj(gca,'Type','patch');
h.FaceColor = '#74C1E4';

ylabel('Frequency')
xlabel('Ellipticity')
xlim([1 2])

fig2.PaperUnits = 'centimeters';
fig2.PaperPosition = [0 0 24 16];
print(fig2,'output/fig2_4','-dpng','-r500')