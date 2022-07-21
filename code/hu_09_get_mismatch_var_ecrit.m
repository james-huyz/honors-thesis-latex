%% ******************************* STEP 09 ********************************
% Redo step 6 (with simplifications) for various e_crit values, i.e., 1.00,
% 1.05, 1.10, and 1.15.
%                                                            James Hu, 2022
% *************************************************************************

% 1. Load the appropriate crater data and latitude-normalized obliquity
%    predictions.
scenario = 2;  % 1 = constant obliquity;
               % 2 = Gaussian variation around mean obliquity
if scenario == 2
    sigma = 5.0326;  % Standard deviation from Laskar 2004, past 2.5 Ma:
                     % load('laskar/laskar_nonchaotic.mat')
                     % rad2deg(std(laskar(1:2501,2)))
end

for mode = 1:2  % 1 = lAv; 2 = AHv
    if mode == 1
        load('craters_obs/lAv1.00.mat')
        a100 = a;
        load('craters_obs/lAv1.05.mat')
        a105 = a;
        load('craters_obs/lAv1.10.mat')
        a110 = a;
        load('craters_obs/lAv1.15.mat')
        a115 = a;
        clear('a')
        load('forward_model/obpreds4_norm_lAv1.00.mat','orientationdata')
        orientationdata100 = orientationdata;
        load('forward_model/obpreds4_norm_lAv1.05.mat','orientationdata')
        orientationdata105 = orientationdata;
        load('forward_model/obpreds4_norm_lAv1.10.mat','orientationdata')
        orientationdata110 = orientationdata;
        load('forward_model/obpreds4_norm_lAv1.15.mat','orientationdata')
        orientationdata115 = orientationdata;
        clear('orientationdata')
    elseif mode == 2
        load('craters_obs/AHv1.00.mat')
        a100 = a;
        load('craters_obs/AHv1.05.mat')
        a105 = a;
        load('craters_obs/AHv1.10.mat')
        a110 = a;
        load('craters_obs/AHv1.15.mat')
        a115 = a;
        clear('a')
        load('forward_model/obpreds4_norm_AHv1.00.mat','orientationdata')
        orientationdata100 = orientationdata;
        load('forward_model/obpreds4_norm_AHv1.05.mat','orientationdata')
        orientationdata105 = orientationdata;
        load('forward_model/obpreds4_norm_AHv1.10.mat','orientationdata')
        orientationdata110 = orientationdata;
        load('forward_model/obpreds4_norm_AHv1.15.mat','orientationdata')
        orientationdata115 = orientationdata;
        clear('orientationdata')
    end
    
    cratdata100 = a100(:,2);
    cratdata105 = a105(:,2);
    cratdata110 = a110(:,2);
    cratdata115 = a115(:,2);
    
% 2. Take the actual crater orientations and construct crude
%    probability distribution function (histogram)
    bin_centers = 5:10:85;
    num_bins = width(bin_centers);
    [N100,X] = hist(cratdata100,bin_centers);  % Raw frequencies
    [N105,X] = hist(cratdata105,bin_centers);
    [N110,X] = hist(cratdata110,bin_centers);
    [N115,X] = hist(cratdata115,bin_centers);
    num_crats100 = height(cratdata100);
    num_crats105 = height(cratdata105);
    num_crats110 = height(cratdata110);
    num_crats115 = height(cratdata115);
    N_norm100 = N100 / num_crats100 / num_bins;  % Normalized distribution
    N_norm105 = N105 / num_crats105 / num_bins;
    N_norm110 = N110 / num_crats110 / num_bins;
    N_norm115 = N115 / num_crats115 / num_bins;
    
%% Compare the PDFs of the actual and simulated crater orientations
    
% 3. Calculate the standard error on the actual crater probability
%    distribution function; plot with error bars
    fig1 = figure('color',[1 1 1]);
    subplot(2,1,mode)
    hold on
    plot(X,N_norm100,'Color','#74C1E4','LineWidth',1.5,'LineStyle',':');
    E105 = sqrt(N105) / num_crats105 / num_bins;
    errorbar(X,N_norm105,E105,'Color','#8D6135','LineWidth',1.5)
    plot(X,N_norm110,'Color','#74C1E4','LineWidth',1.5,'LineStyle','-.');
    plot(X,N_norm115,'Color','#74C1E4','LineWidth',1.5,'LineStyle','--');
    
% 4. Take the simulated crater orientation data and construct crude
%    probability distribution functions; plot
    sims100 = zeros(91,num_bins);
    sims105 = zeros(91,num_bins);
    sims110 = zeros(91,num_bins);
    sims115 = zeros(91,num_bins);

    num_simcrats = height(orientationdata100);  % This is the same for all e_crits
    
    for i = 0:1:90  % Loop through all obliquity scenarios
        if scenario == 1
            odata100 = orientationdata100(:,i+1);
            odata105 = orientationdata105(:,i+1);
            odata110 = orientationdata110(:,i+1);
            odata115 = orientationdata115(:,i+1);
        elseif scenario == 2
            odata100 = hu_06b_gaussamp(orientationdata100,i,sigma);
            odata105 = hu_06b_gaussamp(orientationdata105,i,sigma);
            odata110 = hu_06b_gaussamp(orientationdata110,i,sigma);
            odata115 = hu_06b_gaussamp(orientationdata115,i,sigma);
        end
        
        f100 = hist(odata100,bin_centers) / num_simcrats / num_bins;
        sims100(i+1,:) = f100;
        f105 = hist(odata100,bin_centers) / num_simcrats / num_bins;
        sims105(i+1,:) = f105;
        f110 = hist(odata100,bin_centers) / num_simcrats / num_bins;
        sims110(i+1,:) = f110;
        f115 = hist(odata100,bin_centers) / num_simcrats / num_bins;
        sims115(i+1,:) = f115;
    end
    
    
%% Conduct 91 chi-squared tests to determine best-match constant obliquity

% 5. Calculate reduced chi-squared values for each obliquity
    chi2_sims100 = hu_06a_chi2test(sims100'*num_crats100*num_bins,N100');
    chi2_sims105 = hu_06a_chi2test(sims105'*num_crats105*num_bins,N105');
    chi2_sims110 = hu_06a_chi2test(sims110'*num_crats110*num_bins,N110');
    chi2_sims115 = hu_06a_chi2test(sims115'*num_crats115*num_bins,N115');

    chi2_sims_rdc100 = chi2_sims100 / (num_bins-1);
    chi2_sims_rdc105 = chi2_sims105 / (num_bins-1);
    chi2_sims_rdc110 = chi2_sims110 / (num_bins-1);
    chi2_sims_rdc115 = chi2_sims115 / (num_bins-1);
    
    %% 6. Plot
    fig2 = figure('color',[1 1 1]);
    hold on
    ylim([0 25])
    plot(0:90,chi2_sims_rdc100,'LineWidth',1.5,'Color','#74C1E4',...
         'LineStyle',':');
    plot(0:90,chi2_sims_rdc105,'LineWidth',1.5,'Color','#8D6135');
    plot(0:90,chi2_sims_rdc110,'LineWidth',1.5,'Color','#74C1E4',...
         'LineStyle','-.');    
    plot(0:90,chi2_sims_rdc115,'LineWidth',1.5,'Color','#74C1E4',...
         'LineStyle','--');
        
    % 7. Draw lines on plot - commented out; use to find uncertainty
    [yMin_rdc,xIndex_rdc] = min(chi2_sims_rdc100);  % Observed
    xline(xIndex_rdc-1,'Color','#8D6135')
    y_uncert_rdc = yMin_rdc + 1/(num_bins-1);
    yline(y_uncert_rdc,'Color','#8D6135')
    y_uncert_rdc_arr = y_uncert_rdc + zeros(1,91);
    [xIndex_uncert_rdc,~] = intersections(0:90,chi2_sims_rdc100,0:90,...
                                          y_uncert_rdc_arr);
    xline([min(xIndex_uncert_rdc) max(xIndex_uncert_rdc)],'Color',...
         '#8D6135')  % Uncertainty (x)

% 8. Post-processing
    clearvars -except mode scenario sigma fig*
    
    if mode == 1
        fig11 = fig1;
        fig21 = fig2;
    elseif mode == 2
        fig12 = fig1;
        fig22 = fig2;
    end
    
    clear fig1 fig2
end

%% Tile plots together

% 9. Figure 3.1
fig1 = figure;
t1 = tiledlayout(1,2,'TileSpacing','compact');

fig1a = nexttile;
copyobj(allchild(get(fig11,'CurrentAxes')),fig1a)
title('(a)')
xlim([0 90])
ylim([0.006 0.03])

fig1b = nexttile;
copyobj(allchild(get(fig12,'CurrentAxes')),fig1b)
title('(b)')
xlim([0 90])
ylim([0.006 0.03])

ylabel(t1, 'Probability density')
xlabel(t1, 'Azimuth (° from due north)')
l1 = legend('1.00','1.05','1.10','1.15');
l1.Orientation = 'vertical';

fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0 0 24 16];
print(fig1,'output/fig4_2','-dpng','-r500')

% 10. Figure 3.2
fig2 = figure;
t2 = tiledlayout(1,2,'TileSpacing','compact');

fig2a = nexttile;
copyobj(allchild(get(fig21,'CurrentAxes')),fig2a)
title('(a)')
xlim([0 90])
ylim([0 30])
set(gca,'Ydir','reverse')

fig2b = nexttile;
copyobj(allchild(get(fig22,'CurrentAxes')),fig2b)
title('(b)')
xlim([0 90])
ylim([0 30])
set(gca,'Ydir','reverse')

ylabel(t2, 'Mismatch (reduced \chi^2 statistic)')
xlabel(t2, 'Obliquity (°)')

l2 = legend('1.00','1.05','1.10','1.15');
l2.Orientation = 'vertical';

fig2.PaperUnits = 'centimeters';
fig2.PaperPosition = [0 0 24 16];
print(fig2,'output/fig4_3','-dpng','-r500')