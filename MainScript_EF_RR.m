% Code for Bayesian estimation of advanced warning time of precipitation emergence 
% Lickley and Fletcher, submission 1, August 24, 2023 to Earth's Future
% Revision 1, December 7th 2023
% The 


clear all
close all

nModels = 16;
nLat = 90;
nLon = 180; 

HomeDir = pwd;

yrs = 1850.08:1/12:2300;
windowSize = 30;

load_cru = 0; 
create_pr_6x6 = 0; 


%% loading GCM data, spatially averaged at 6d x 6d resolution. 
if create_pr_6x6
    preprocess_EFSubmission
else
    load('pr_6x6_grid.mat'); 
end

%%  loading CRU data
% Harris, I., Osborn, T.J., Jones, P. et al. Version 4 of the CRU TS monthly 
% high-resolution gridded multivariate climate dataset. Sci Data 7, 109 (2020). 
% https://doi.org/10.1038/s41597-020-0453-3

if load_cru
    preprocess_EFSubmission
else
    load('CRU_regridded.mat');
end

%% Calibrating model hyperparameters and validating model 
calibrating_model = 0; 
validating_noisy_dist = 0; 
validating_future_mean = 1; 
model_calibration_validation(calibrating_model, validating_noisy_dist, validating_future_mean)

%% Calculating precip moving averages

pr_tmp =  log(pr_6x6); % log transform of precip
ind = find(isinf(pr_tmp));
pr_tmp(ind) = NaN;

histP = nanmean(pr_tmp(:,:,1:100,:),3);  % historical precip
histStDevP = std(pr_tmp(:,:,1:100,:), 0, 3, 'omitnan');  % historical noise

%Five Year anomalies relative to historical mean
for yy = 1:90
    FiveYr(:,:,yy,:) = (nanmean(pr_tmp(:,:,5*(yy-1)+1:5*yy,:),3) - histP)./histStDevP; 
end

% 50 year moving average
for yy = 1:451-50
    smoothedMA(:,:,yy,:) = (nanmean(pr_tmp(:,:,yy:yy+50,:),3)- histP)./histStDevP;  
end

% Five year anomaly relative to 23rd century precip
FiveYr_anoms = FiveYr(:,:,5:end-5,:) - smoothedMA(:,:,1:5:end,:); 

%% Figure 1: Illustration of methods for at Mombasa for one OOS GCM
% a) prior distribution and individual 50-year smoothed GCM 
% b) posterior mean and sigma given t_obs = 2000,
% c) posterior mean and sigma given t_obs = ToC
% d) posterior mean and sigma given t_obs = ToE


% years corresponding to Five Yr averages and anomalies
yrs_5yr = [1875:5:2275];  

% Location: Mombasa 4°S and 40°E, model 7
lon_ii = 20; 
lat_ii = 45 - 2; 


% Load in hyperparameters from model calibration
load('LOOCV_test_hyperparameter_selectionfull_Dec2023.mat'); 

sigma1 = Max_Likelihood_sigma1(lon_ii,lat_ii); 
sigma2 = Max_Likelihood_sigma2(lon_ii,lat_ii); 
sigma3 = Max_Likelihood_sigma3(lon_ii,lat_ii);

% smoothed values for this illustration
pr_mombasa_smoothed = squeeze(smoothedMA(lon_ii, lat_ii,:,:));  % for plotting smoothed data

mod_ii = 12; %  OOS_GCM for g = 7.  model 12 also provides a useful illustration for positive AWT
modinds = [1:16]; % model indices
mod_indtmp = modinds; 
mod_indtmp(mod_ii) = []; % removing oos_gcm from ensemble for developing prior

oos_gcm = squeeze(FiveYr(lon_ii, lat_ii, 5:end-5, mod_ii)); 
sigma_n = nanvar(squeeze(FiveYr_anoms(lon_ii,lat_ii, 1:20, mod_ii))); % Noise in observations
GCM_data = squeeze(smoothedMA(lon_ii,lat_ii, 1:5:end, mod_indtmp));

% Estimating ToC and ToE for our OOS_GCM  
[ToC ToC_sign ToE ToE_sign] = ToC_ToE_function(oos_gcm, sigma1, sigma2, sigma3, sigma_n, GCM_data, yrs_5yr, 1)


fighandle = figure(1)

% Plotting specifications
Year_end = 2150; % xlim max for plotting purposes
n = 6; % color spacing for blue shades
cmap = [linspace(.9,0,n)', linspace(.9447,.447,n)', linspace(.9741,.741,n)'];

subplot(4, 1, 1)
                
nobs = find(yrs_5yr == 2000);   
[prior_mu, prior_cov, post_mu, post_cov] = GPR_fun(oos_gcm, nobs, GCM_data, sigma1, sigma2, sigma3, sigma_n);
   
% Plotting priors on each subplot
% prior_cov = K.  To plot the prior of noisy data add sigma_n*I

prior_noisy_cov = prior_cov + sigma_n*eye(size(prior_cov));
bounds_opts = diag(prior_noisy_cov).^(0.5);

x = yrs_5yr;
curve1 = prior_mu - bounds_opts';
curve2 = prior_mu + bounds_opts';

for ii = 1:4
    subplot(4, 1, ii)
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween, [0.8, 0.8, 0.8],'EdgeColor','none');
end

subplot(4, 1, 1)
% Plotting individual smoothed GCMs
hold on; plot([1875:2275], pr_mombasa_smoothed(:,mod_indtmp)','r');  
hold on; plot(x, prior_mu,'b', 'LineWidth', 2)

xlim([1950,Year_end]); 
ylim([-1.5,2.5]); 
ylabel('\Delta P [s.d. log(P)]');


% Plotting posterior for nobs up to 2000
% for post covariance of noisy prediction add sigma_n*I to diagonal
post_noisy_cov = post_cov + sigma_n*eye(size(post_cov));
bounds_opts = diag(post_noisy_cov ).^(0.5);

curve1 = post_mu(nobs:end)' - bounds_opts(nobs:end)';
curve2 = post_mu(nobs:end)' + bounds_opts(nobs:end)';

x_tmp = x(nobs:end);
x2 = [x_tmp, fliplr(x_tmp)];

for ii = 2:4
    subplot(4, 1, ii)
    inBetween = [curve1, fliplr(curve2)];
    hold on;    fill(x2, inBetween, cmap(2,:),'EdgeColor','none')
end

% On second panel, including oos_gcm up to nobs along with posterior mean
subplot(4, 1, 2)
hold on; plot(x(1:nobs), oos_gcm(1:nobs),'k', 'LineWidth',2); 
hold on; plot(x(nobs:end), post_mu(nobs:end),'Color',[.5 0 .5], 'LineWidth',2); 
hold on; plot([1900,2150],[0,0],'k'); 
hold on; plot([x(nobs), x(nobs)], [-1.5, 2.5], '--k');

xlim([1950,Year_end]); 
ylim([-1.5,2.5]); 
ylabel('\Delta P [s.d. log(P)]');

% Plotting posterior at time of confidence
nobs = find(yrs_5yr == ToC);  

[prior_mu, prior_cov, post_mu, post_cov] = GPR_fun(oos_gcm, nobs, GCM_data, sigma1, sigma2, sigma3, sigma_n);

post_noisy_cov = post_cov + sigma_n*eye(size(post_cov));
bounds_opts = diag(post_noisy_cov).^(0.5);

curve1 = post_mu(nobs:end)' - bounds_opts(nobs:end)';
curve2 = post_mu(nobs:end)' + bounds_opts(nobs:end)';  

x_tmp = x(nobs:end);
x2 = [x_tmp, fliplr(x_tmp)];

for ii = 3:4
    subplot(4, 1, ii)
    inBetween = [curve1, fliplr(curve2)];
    hold on; fill(x2, inBetween, cmap(4,:),'EdgeColor','none');

end

% Adding ToC line and posterior time series from ToC onwards.
subplot(4, 1, 3)

hold on; plot(x(1:nobs), oos_gcm(1:nobs),'k', 'LineWidth',2)
hold on; plot([1900,2150],[0,0],'k'); 
hold on; plot(x(nobs:end), post_mu(nobs:end),'Color',[.5 0 .5], 'LineWidth',2)

hold on; plot([ToC, ToC], [-1.5, 2.5], '--k'); 
xlim([1950,Year_end]); 
ylim([-1.5, 2.5]);
ylabel('\Delta P [s.d. log(P)]');

nobs = find(yrs_5yr == 2150);   

[prior_mu, prior_cov, post_mu, post_cov] = GPR_fun(oos_gcm, nobs, GCM_data, sigma1, sigma2, sigma3, sigma_n);
post_noisy_cov = post_cov + sigma_n*eye(size(post_cov));
 
% Adding ToE line and posterior time series conditioned on all
% observations.

subplot(4, 1, 4)

bounds_opts = diag(post_noisy_cov).^(0.5);
curve1 = post_mu' - bounds_opts';
curve2 = post_mu' + bounds_opts';
 
x_tmp = x(1:end);   
x2 = [x_tmp, fliplr(x_tmp)];

inBetween = [curve1, fliplr(curve2)];
hold on; fill(x2, inBetween, cmap(6,:),'EdgeColor','none');
hold on; plot(x, oos_gcm,'k', 'LineWidth',2)
hold on; plot(x, post_mu,'Color',[.5 0 .5], 'LineWidth',2)

hold on; plot([1900,2150],[0,0],'k'); 
ind = find(x == ToE);
ToE_line = interp1([curve1(ind-1), curve1(ind)],[ToE - 5, ToE], 0);
hold on; plot([ToE_line, ToE_line], [-1.5, 2.5], '--k')
hold on; plot([ToC, ToC], [-1.5, 2.5], '--k'); 
xlim([1950,Year_end]); 
ylim([-1.5,2.5]);
ylabel('\Delta P [s.d. log(P)]');

Fig_str = strcat('Figures/Figure1.pdf');
figure_width = 6; 
figure_height = 14; 
savefig(fighandle, Fig_str, figure_width,figure_height)

%% Making Fig 2 of paper: % Delta P, stippling indicates uncertainty across GCMs
%  Red stippling indicates the uncertainty will be resolved by 2020, 
%  Yellow stippling indicates the uncertainty will be resolved by 2060
%  Grey stippling indicates not resolved by 2060

% This section calculates the prior (ind_CI = 1) and posterior distribution at 2100 for each
% oos_gcm given t_obs in 2020, 2060, and 2100 (ind_CI = 2, 2, and 4 respectively) 

load('LOOCV_test_hyperparameter_selectionfull_Dec2023.mat'); 

modinds = [1:16];
ChangeUpdatedCI = NaN(180,88,16,4);
ChangeUpdatedCI_noise_free = NaN(180,88,16,4);
for lon_ii = 1:180
    for lat_ii = 1:88
        sigma1 = Max_Likelihood_sigma1(lon_ii,lat_ii); 
        sigma2 = Max_Likelihood_sigma2(lon_ii,lat_ii); 
        sigma3 = Max_Likelihood_sigma3(lon_ii,lat_ii);

        for mod_ii = 1:16

            mod_indtmp = modinds; 
            mod_indtmp(mod_ii) = [];

            oos_gcm = squeeze(FiveYr(lon_ii, lat_ii, 5:end-5, mod_ii)); 
            nobs = find(yrs_5yr == 2020); 
            sigma_n = nanvar(squeeze(FiveYr_anoms(lon_ii,lat_ii, 1:20, mod_ii))); 
            GCM_data = squeeze(smoothedMA(lon_ii,lat_ii, 1:5:end, mod_indtmp));

            [prior_mu, prior_cov, post_mu, post_cov] = GPR_fun(oos_gcm, nobs, GCM_data, sigma1, sigma2, sigma3, sigma_n);
            prior_noisy_cov = prior_cov + sigma_n*eye(size(prior_cov));
            posterior_noisy_cov = post_cov + sigma_n*eye(size(prior_cov));

            % Year of prediction target: 2100
            ind_2100 = find(yrs_5yr == 2100); 
            tmp_mu = prior_mu(ind_2100);
            tmp_std = prior_noisy_cov(ind_2100,ind_2100)^0.5;
            tmp_std_noise_free = prior_cov(ind_2100,ind_2100)^0.5;
                        
            ind_CI = 1; % Noise and Signal for prior

            if tmp_mu > 0

                ChangeUpdatedCI(lon_ii, lat_ii, mod_ii, ind_CI) = 1 - normcdf(0, tmp_mu, tmp_std);
                ChangeUpdatedCI_noise_free(lon_ii, lat_ii, mod_ii, ind_CI) = 1 - normcdf(0, tmp_mu, tmp_std_noise_free);

            elseif tmp_mu < 0 

                ChangeUpdatedCI(lon_ii, lat_ii, mod_ii, ind_CI) = normcdf(0, tmp_mu, tmp_std);
                ChangeUpdatedCI_noise_free(lon_ii, lat_ii, mod_ii, ind_CI) = normcdf(0, tmp_mu, tmp_std_noise_free);

            end
   
            tmp_mu = post_mu(ind_2100);
            tmp_std = posterior_noisy_cov(ind_2100,ind_2100)^0.5;
            tmp_std_noise_free = post_cov(ind_2100,ind_2100)^0.5;
            
            ind_CI = 2; % Noise and Signal for posterior updated to 2020
            
            if tmp_mu > 0

                ChangeUpdatedCI(lon_ii, lat_ii, mod_ii, ind_CI) = 1 - normcdf(0, tmp_mu, tmp_std);
                ChangeUpdatedCI_noise_free(lon_ii, lat_ii, mod_ii, ind_CI) = 1 - normcdf(0, tmp_mu, tmp_std_noise_free);

            elseif tmp_mu < 0

                ChangeUpdatedCI(lon_ii, lat_ii, mod_ii, ind_CI) = normcdf(0, tmp_mu, tmp_std);
                ChangeUpdatedCI_noise_free(lon_ii, lat_ii, mod_ii, ind_CI) = normcdf(0, tmp_mu, tmp_std_noise_free);

            end

            nobs = find(yrs_5yr == 2060); % Uncertainty at 2050
            [prior_mu, prior_cov, post_mu, post_cov] = GPR_fun(oos_gcm, nobs, GCM_data, sigma1, sigma2, sigma3, sigma_n);
            posterior_noisy_cov = post_cov + sigma_n*eye(size(prior_cov));

            tmp_mu = post_mu(ind_2100);
            tmp_std = posterior_noisy_cov(ind_2100,ind_2100)^0.5;
            tmp_std_noise_free = post_cov(ind_2100,ind_2100)^0.5;
            
            ind_CI = 3; % Noise and Signal for posterior updated to 2050        

            if tmp_mu > 0

                ChangeUpdatedCI(lon_ii, lat_ii, mod_ii, ind_CI) = 1 - normcdf(0, tmp_mu, tmp_std);
                ChangeUpdatedCI_noise_free(lon_ii, lat_ii, mod_ii, ind_CI) = 1 - normcdf(0, tmp_mu, tmp_std_noise_free);

            elseif tmp_mu < 0

                ChangeUpdatedCI(lon_ii, lat_ii, mod_ii, ind_CI) = normcdf(0, tmp_mu, tmp_std);
                ChangeUpdatedCI_noise_free(lon_ii, lat_ii, mod_ii, ind_CI) = normcdf(0, tmp_mu, tmp_std_noise_free);

            end


            nobs = find(yrs_5yr == 2100); % Emergence by 2100
            [prior_mu, prior_cov, post_mu, post_cov] = GPR_fun(oos_gcm, nobs, GCM_data, sigma1, sigma2, sigma3, sigma_n);
            posterior_noisy_cov = post_cov + sigma_n*eye(size(prior_cov));

            tmp_mu = post_mu(ind_2100);
            tmp_std = posterior_noisy_cov(ind_2100,ind_2100)^0.5;
            tmp_std_noise_free = post_cov(ind_2100,ind_2100)^0.5;
            
            ind_CI = 4; % Noise and Signal for posterior updated to 2050        

            if tmp_mu > 0

                ChangeUpdatedCI(lon_ii, lat_ii, mod_ii, ind_CI) = 1 - normcdf(0, tmp_mu, tmp_std);
                ChangeUpdatedCI_noise_free(lon_ii, lat_ii, mod_ii, ind_CI) = 1 - normcdf(0, tmp_mu, tmp_std_noise_free);

            elseif tmp_mu < 0

                ChangeUpdatedCI(lon_ii, lat_ii, mod_ii, ind_CI) = normcdf(0, tmp_mu, tmp_std);
                ChangeUpdatedCI_noise_free(lon_ii, lat_ii, mod_ii, ind_CI) = normcdf(0, tmp_mu, tmp_std_noise_free);

            end

        end

    end
    lon_ii
end

%%  Figure 2 map for when at least half of GCMs are confident about 2100 precip being significantly different

for fig_ii = 2 %:3
    if fig_ii == 2
        ConfidentUpdate = squeeze(sum(ChangeUpdatedCI > 0.84 , 3));  %% Number of models where the posterior 1-sigma CI does not contain zero.        
        fig_str = 'Figures/Figure2map_uncertainty.pdf';
    else
        ConfidentUpdate = squeeze(sum(ChangeUpdatedCI_noise_free > 0.84 , 3));  %% Number of models where the posterior 1-sigma CI does not contain zero. 
        fig_str = 'Figures/FigureS2map_uncertainty_noise_free.pdf'
    end
    
    % calculating the fraction of land that is not confident about 2100
    % change in 2100 for prior and tobs = 2020, 2060 2100
    latgrid = flipud(linspace(-87,87, 88)'); 
    Lat_grid = repmat(-latgrid',180,1);
    Lat_weights = cosd(Lat_grid');
    UncertainPrior_84CI = 1 - sum(sum((squeeze(ConfidentUpdate(:,:,1))'> 7).*Lat_weights))/sum(sum(Lat_weights))
    UncertainPostrior2020_84CI = 1 - sum(sum((squeeze(ConfidentUpdate(:,:,2))'> 7).*Lat_weights))/sum(sum(Lat_weights))
    UncertainPostrior2060_84CI = 1 - sum(sum((squeeze(ConfidentUpdate(:,:,3))'> 7).*Lat_weights))/sum(sum(Lat_weights))
    UncertainPostrior2100_84CI = 1 - sum(sum((squeeze(ConfidentUpdate(:,:,4))'> 7).*Lat_weights))/sum(sum(Lat_weights))


    fighandle = figure(fig_ii)
    
    yrs = 1850:2300;
    
    ind0 = find((yrs > 1900) & (yrs < 1950));
    HistYrs = squeeze(nanmean(pr_6x6(:, :, ind0, :),3));
    
    ind1 = find((yrs > 2070) &  (yrs < 2100));
    BaseYrs = squeeze(nanmean(pr_6x6(:, :, ind1, :),3));
    
    Delta1 = 100*(BaseYrs - HistYrs)./HistYrs;  % percentage change in precip. 
    
    latgrid = flipud(linspace(-87,87,size(pr_6x6,2))'); 
    longrid = linspace(1,359,size(pr_6x6,1));
    
    m_proj('robinson','long',[-180 180],'lat',[-80 80]); 
    [lat, long] = meshgrid(linspace(-180, 180,360),linspace(-90, 90,360));
    Lon_grid = repmat(longrid',1,88);
    Lat_grid = repmat(-latgrid',180,1);
    
    % Mapping MMM percentage change in precip
    TempGrid = griddata([Lon_grid-360, Lon_grid], [Lat_grid, Lat_grid], [mean(Delta1,3), mean(Delta1,3)], lat, long); 
    m_pcolor(lat, long, TempGrid); shading flat;  
    
    m_coast('color',[0.2 .2 0.2]); 
    m_grid('linest','none','xticklabels',[],'yticklabels',[]);  
    
    clrmp = cbrewer('div','RdBu',50);
    clrmp(clrmp<0)=0;
    clrmp = (1/max(max(clrmp)))*clrmp;
    colormap(clrmp);  colorbar; caxis([-30, 30]); 
    
    
    for pp = 1:3
    
        LAT = [];
        LON = [];
        % For plotting purposes, only calculate every fourth grid cell
        for ii = 1:4:180
            for jj = 2:4:89
                if  ConfidentUpdate(ii, jj-1, pp) < 0.5*16 
    
                    LAT = [LAT; Lat_grid(ii,jj)];
                    LON = [LON; Lon_grid(ii,jj)];
                end
            end
        end
    
        hold on;
    
        LON(LON>180) = LON(LON>180)-360;
    
        if pp == 1
            m_plot(LON, LAT, '.','Color', [1.0000    0    0],'MarkerSize',12);  
            hold on; 
        elseif pp == 2
    
            m_plot(LON, LAT, '.','Color', [1    0.35    0.7],'MarkerSize',12);  
            hold on;
        elseif pp == 3
    
            m_plot(LON, LAT, '.','Color', [0.7 0.7 0.7], 'MarkerSize',12);  
            hold on;
        end
    
    end
    
    title('2070-2100 %\Delta P')
    
    figure_width = 6; 
    figure_height = 6; 
    savefig(fighandle, fig_str, figure_width,figure_height)
end


Will_emerge_by2100 = 1 - UncertainPostrior2100_84CI;
UncertainPrior_84CI
Certain_by_2020 = 1 - UncertainPostrior2020_84CI;
Certain_by_2060 = 1 - UncertainPostrior2060_84CI;
(Certain_by_2020-UncertainPostrior2100_84CI)/Will_emerge_by2100
(Certain_by_2060-UncertainPostrior2100_84CI)/Will_emerge_by2100

(Certain_by_2020-UncertainPostrior2100_84CI)
(Certain_by_2060-UncertainPostrior2100_84CI)
%%  Subplots for Figure 2: Submitted manuscript includes first three. 

titleopts{1}='Mombasa, Kenya'; % Mombasa, 4S, 40E
titleopts{2}='Damascus, Syria'; % Damascus, Syria, 34N, 46E
titleopts{3}='Nuuk, Greenland'; % Greenland 64N, 51W
titleopts{4}='Uvira, DR Congo'; % Uvira, Dr Congo, 3S, 29E
titleopts{5}='Islamabad Pakistan'; % Islamabad, Pakistan, 34N, 73E
titleopts{6}='Brasilia, Brazil'; % Brasilia, 15.8S, 48W

Lat_opts = 46 + floor(0.5*[-4, 34, 64, -3, 34, -16]);
Lon_opts = floor(0.5*[40, 46, 360 - 51, 29, 73, 360 - 48]);

n = 3;
cmap = [linspace(.9,0,n)', linspace(.9447,.447,n)', linspace(.9741,.741,n)'];

fighandle = figure(4)

mod_ii = 1;
mod_indtmp = modinds; 
mod_indtmp(mod_ii) = [];

for tt = 1:6
    lon_ii = Lon_opts(tt);
    lat_ii = Lat_opts(tt);

    sigma1 = Max_Likelihood_sigma1(lon_ii,lat_ii); 
    sigma2 = Max_Likelihood_sigma2(lon_ii,lat_ii); 
    sigma3 = Max_Likelihood_sigma3(lon_ii,lat_ii);
    
    Year_opts = [2020, 2050, 2080];
    
    oos_gcm = squeeze(FiveYr(lon_ii, lat_ii, 5:end-5, mod_ii)); 
    sigma_n = nanvar(squeeze(FiveYr_anoms(lon_ii,lat_ii, 1:20, mod_ii)));   
    GCM_data = squeeze(smoothedMA(lon_ii,lat_ii, 1:5:end, mod_indtmp));

    for yy = 1:3

        nobs = find(yrs_5yr == Year_opts(yy));
        [prior_mu, prior_cov, post_mu, post_cov] = GPR_fun(oos_gcm, nobs, GCM_data, sigma1, sigma2, sigma3, sigma_n);
        prior_cov = prior_cov + sigma_n*eye(size(prior_cov));
        post_cov = post_cov + sigma_n*eye(size(post_cov));
        if yy == 1
            bounds_opts = diag(prior_cov).^(0.5);
            subplot(2, 3, tt); 
            x = yrs_5yr;
            curve1 = prior_mu - bounds_opts';
            curve2 = prior_mu + bounds_opts';
            plot(x, curve1, 'k', 'LineWidth', 2);
            hold on;  
            plot(x, curve2, 'k', 'LineWidth', 2);
            x2 = [x, fliplr(x)];
            inBetween = [curve1, fliplr(curve2)];
            fill(x2, inBetween, [0.8, 0.8, 0.8]);
        end

        hold on; 
        bounds_opts = diag(post_cov).^(0.5);
        x = yrs_5yr(nobs:end);
        curve1 = [oos_gcm(nobs), post_mu(nobs+1:end)' - bounds_opts(nobs+1:end)'];
        curve2 = [oos_gcm(nobs), post_mu(nobs+1:end)' + bounds_opts(nobs+1:end)'];
        plot(x(2:end), curve1(2:end), 'k', 'LineWidth', 2);
        hold on;  plot(x(2:end), curve2(2:end), 'k', 'LineWidth', 2);
        hold on;  plot(x(1:2), curve1(1:2), '--k', 'LineWidth', 2);
        hold on;  plot(x(1:2), curve2(1:2), '--k', 'LineWidth', 2);
        x2 = [x, fliplr(x)];
        inBetween = [curve1, fliplr(curve2)];
        hold on;  fill(x2, inBetween, cmap(yy,:),'LineStyle','none','FaceAlpha',0.7);
        hold on;  plot(x(1), oos_gcm(nobs),'.k','MarkerSize',20)
    end
        hold on;  plot(yrs_5yr, oos_gcm, 'k', 'LineWidth', 2);
        xlim([1900,2150]);
        box on
        hold on;  plot(yrs_5yr, 0*yrs_5yr,'k', 'LineWidth', 1);
        title(titleopts{tt})
end
Fig_str = strcat('Figures/Fig2subplots_mod1.pdf');
figure_width = 12; 
figure_height = 6; 
savefig(fighandle, Fig_str, figure_width,figure_height)


%% Figure S1 c) and d)  Correlations between near-term and long-term change  % Bottom: fraction of world that never emerges

fighandle = figure(5)

yrs = 1850:2300;
ind0 = find((yrs > 1900) & (yrs < 1950));
HistYrs = squeeze(nanmean(pr_6x6(:, :, ind0, :),3));
ind1 = find((yrs > 2200) & (yrs < 2300));
FutYrs = squeeze(nanmean(pr_6x6(:, :, ind1, :),3));
DeltaFut = (FutYrs - HistYrs)./HistYrs;

cmin = -0.2;
cmax = 1; 

ind2 = find((yrs > 2020) &  (yrs < 2050));
IntermYrs = squeeze(nanmean(pr_6x6(:, :, ind2, :),3));
DeltaInterm = (IntermYrs - HistYrs)./HistYrs;

clrmp = cbrewer('seq','PuRd',60);
clrmp = clrmp(1:end-10,:); 

for lat_ii = 1:88
    for lon_ii = 1:180
        CorrMap(lon_ii, lat_ii) = corr(squeeze(DeltaInterm(lon_ii,lat_ii,:)),squeeze(DeltaFut (lon_ii,lat_ii,:)));
    end
end
str = '2020-2050';
subplot(2,1,1)
globalmapping_withcoords(CorrMap', -latgrid, longrid, 0, cmin, cmax, str, 0, 360, -80, 80, clrmp)

for yy = 1:250
    tmp_weights = reshape(Lat_weights, 180*88, 1); 
    
    ind2 = find((yrs > 1950 + yy) &  (yrs < 1980 + yy));
    
    IntermYrs = squeeze(nanmean(pr_6x6(:, :, ind2, :),3));
    DeltaInterm = (IntermYrs - HistYrs)./HistYrs;

    for lat_ii = 1:88
        for lon_ii = 1:180
            CorrMap(lon_ii, lat_ii) = corr(squeeze(DeltaInterm(lon_ii,lat_ii,:)),squeeze(DeltaFut (lon_ii,lat_ii,:)));
        end
    end
    
    tmp_corr = reshape(CorrMap, 180*88, 1); 
    
    [corr_out,idx] = sort(tmp_corr); 
    
    indx = ~isnan(tmp_weights(idx)); 
    tmp_weights2 = cumsum(tmp_weights(idx(indx)))/sum(tmp_weights(idx(indx))); 
    corr_out2 = corr_out(indx); 
    
    for ii = 1:100
        indtmp = find(tmp_weights2 > 0.005+0.01*(ii-1),1);
        Corr(ii, yy) = corr_out2(indtmp);
    end
    yy
end

subplot(2,1,2)
contourf([1966:1965+250],[0.5:99.5],Corr,50,'LineStyle','none'); 
colormap(clrmp); xlim([1980, 2200]); caxis([-0.2, 1]); 
ylabel('% Surface Area'); colorbar

Fig_str = strcat('Figures/FigS1c_d.pdf');
figure_width = 6; 
figure_height = 8; 
savefig(fighandle, Fig_str, figure_width,figure_height)


%% Estimating ToE and ToC at each grid cell and for each model


run_ToE_ToC = 1; 
ci_level = 1; 
if run_ToE_ToC
    if ci_level == 1
        % in ToC_ToE_function.m change fraction_past_emergence = 0.5
        str_name = 'ToE_ToC_p_1sigma_Dec2023.mat'

        % in ToC_ToE_function.m change fraction_past_emergence to 0.25
        % str_name = 'ToE_ToC_p_1sigma_25percentfuture.mat'; 

        % in ToC_ToE_function.m change fraction_past_emergence to 0.25
        % str_name = 'ToE_ToC_p_1sigma_75percentfuture.mat'; 
        sig = 1; 
    else 
        str_name = 'ToE_ToC_p_2sigma_Dec2023.mat'
        sig = 2;

    end

    modinds = [1:16];

    ToE_wet = NaN(180, 88, 16); 
    ToC_wet = NaN(180, 88, 16); 

    ToE_dry = NaN(180, 88, 16); 
    ToC_dry = NaN(180, 88, 16); 

    for lon_ii = 1:180
        for lat_ii = 1:88

        sigma1 = Max_Likelihood_sigma1(lon_ii,lat_ii); 
        sigma2 = Max_Likelihood_sigma2(lon_ii,lat_ii); 
        sigma3 = Max_Likelihood_sigma3(lon_ii,lat_ii);


        for mod_ii = 1:16
            

            mod_indtmp = modinds; 
            mod_indtmp(mod_ii) = [];
            
            oos_gcm = squeeze(FiveYr(lon_ii, lat_ii, 5:end-5, mod_ii)); 
            sigma_n = nanvar(squeeze(FiveYr_anoms(lon_ii,lat_ii, 1:20, mod_ii)));    
            GCM_data = squeeze(smoothedMA(lon_ii,lat_ii, 1:5:end, mod_indtmp));

            yrs_5yr = [1875:5:2275];

            [ToC ToC_sign ToE ToE_sign] = ToC_ToE_function(oos_gcm, sigma1, sigma2, sigma3, sigma_n, GCM_data, yrs_5yr, sig);

            if ToE_sign == -1
                ToE_dry(lon_ii, lat_ii, mod_ii) = ToE; 
            elseif ToE_sign == 1
                ToE_wet(lon_ii, lat_ii, mod_ii) = ToE; 
            end
            
            if ToC_sign == -1
                ToC_dry(lon_ii, lat_ii, mod_ii) = ToC; 
            elseif ToE_sign == 1
                ToC_wet(lon_ii, lat_ii, mod_ii) = ToC; 
            end
        end
        end
        lon_ii
    end
    save(str_name , 'ToE_dry','ToE_wet','ToC_dry','ToC_wet')

end

%% Testing sensitivity of ToE and ToC definition emerged for at least 25%, 50% and 75% of years beyond initial emergence
yrs = 1975:2150;
% For an estimate of the sensitivity, the above section of code was run on 
% every 3rd lat/lon value.  Weights here are modified accordingingly.
% ToC_ToE_function was modified such that the initial window of emergence
% was identified for three scenarios when at least 25%, 50%, and 75% of 
% future years were significant.  

LatWeights = repmat(cos((pi/180)*latgrid),1,180);
LatWeights_tmp = NaN(size(LatWeights));
LatWeights_tmp(1:3:end,1:3:end) = LatWeights(1:3:end,1:3:end);
LatWeights_tmp = LatWeights_tmp';

for ii = 1:3
    if ii == 1
        load('ToE_ToC_p_1sigma_25percentfuture.mat')
        ToC_25wet = ToC_wet; 
        ToE_25wet = ToE_wet; 
    elseif ii == 2
        load('ToE_ToC_p_1sigma_Dec2023.mat')
        ToC_50wet = ToC_wet; 
        ToE_50wet = ToE_wet; 
    elseif ii == 3
        load('ToE_ToC_p_1sigma_75percentfuture.mat')
        ToC_75wet = ToC_wet; 
        ToE_75wet = ToE_wet; 
    end
    

    for yy = 1:length(yrs)
        tmp = nanmean(ToE_dry,3);
        ind = find(tmp < yrs(yy)); 
        ToE_dry_sensitivity(yy,ii) = nansum(LatWeights_tmp(ind))./nansum(nansum(LatWeights_tmp)); 
        
        if yy == length(yrs)
            store_mean(1, ii)  = nansum(nansum(LatWeights_tmp(ind).*tmp(ind)))./nansum(nansum(LatWeights_tmp(ind)));  
        end

        tmp = nanmean(ToE_wet,3);
        ind = find(tmp < yrs(yy)); 
        ToE_wet_sensitivity(yy,ii) = nansum(LatWeights_tmp(ind))./nansum(nansum(LatWeights_tmp)); 
        
        if yy == length(yrs)
            store_mean(2, ii)  = nansum(nansum(LatWeights_tmp(ind).*tmp(ind)))./nansum(nansum(LatWeights_tmp(ind)));  
        end

        tmp = nanmean(ToC_dry,3);
        ind = find(tmp < yrs(yy)); 
        ToC_dry_sensitivity(yy,ii) = nansum(LatWeights_tmp(ind))./nansum(nansum(LatWeights_tmp)); 
        
        if yy == length(yrs)
            store_mean(3, ii)  = nansum(nansum(LatWeights_tmp(ind).*tmp(ind)))./nansum(nansum(LatWeights_tmp(ind)));  
        end

        tmp = nanmean(ToC_wet,3);
        ind = find(tmp < yrs(yy)); 
        ToC_wet_sensitivity(yy,ii) = nansum(LatWeights_tmp(ind))./nansum(nansum(LatWeights_tmp)); 

        if yy == length(yrs)
            store_mean(4, ii)  = nansum(nansum(LatWeights_tmp(ind).*tmp(ind)))./nansum(nansum(LatWeights_tmp(ind)));  
        end

    end


end

fighandle = figure; 
subplot(2,2,1)
plot(yrs, 100*ToE_dry_sensitivity, 'LineWidth', 2);
xlim([1976, 2150]);
ylim([0,100])
ylabel('Fraction global surface area')
title('ToE^C dry models')

subplot(2,2,2)
plot(yrs, 100*ToE_wet_sensitivity, 'LineWidth', 2)
xlim([1976, 2150])
ylim([0,100])
ylabel('Fraction global surface area')
title('ToE^C wet models')

subplot(2,2,3)
plot(yrs, 100*ToC_dry_sensitivity, 'LineWidth', 2)
xlim([1976, 2150])
ylim([0,100])
ylabel('Fraction global surface area')
title('ToC dry models')

subplot(2,2,4)
plot(yrs, 100*ToC_wet_sensitivity, 'LineWidth', 2)
xlim([1976, 2150])
ylim([0,100])
ylabel('Fraction global surface area')
title('ToC wet models')
legend('25% sig. years > t_{obs}', '50% sig. years > t_{obs}', '75% sig. years > t_{obs}')

Fig_str = strcat('Figures/Sensitivity_futureyears.pdf'); 
figure_width = 8; 
figure_height = 6; 
savefig(fighandle, Fig_str, figure_width,figure_height)

% average ToE/ToC values for different thresholds of future years
        
VarNames = {'25% future years sig', '50% future years sig', '75% future years sig'};
RowNames = {'ToE dry,' 'ToE wet','ToC dry,' 'ToC wet'};
T = table(store_mean(:, 1), store_mean(:, 2), store_mean(:, 3), 'VariableNames',VarNames,'RowNames',RowNames)



%%  Figure S3: Number of models that are wet vs dry. 

load('ToE_ToC_p_1sigma_Dec2023.mat')
clrmp = cbrewer('seq','PuRd',60);
clrmp = clrmp(1:end-10,:); 
fighandle = figure(6); 

num_dry = sum(~isnan(ToE_dry),3); 
num_wet = sum(~isnan(ToE_wet),3); 
num_change = num_dry + num_wet; 

num_dry_ToC = sum(~isnan(ToC_dry),3); 
num_wet_ToC = sum(~isnan(ToC_wet),3); 
num_change_ToC = num_dry_ToC + num_wet_ToC; 

subplot(1,3,1)
Data_plot = [num_wet(91:180,:); num_wet(1:90,:)];
globalmapping_withcoords(Data_plot', -latgrid,[longrid(91:180)-360, longrid(1:90)], 0, 0, 16, 'Number Wet models' , -180, 180, -80, 80, clrmp)

subplot(1,3,2)
Data_plot = [num_dry(91:180,:); num_dry(1:90,:)];
globalmapping_withcoords(Data_plot', -latgrid,[longrid(91:180)-360, longrid(1:90)], 0, 0, 16, 'Number Dry models' , -180, 180, -80, 80, clrmp)

subplot(1,3,3)
Data_plot = [num_change(91:180,:); num_change(1:90,:)];
globalmapping_withcoords(Data_plot', -latgrid,[longrid(91:180)-360, longrid(1:90)], 0, 0, 16, 'Number Wet + Dry models' , -180, 180, -80, 80, clrmp)

Fig_str = strcat('Figures/NumberOfModels.pdf'); 
figure_width = 20; 
figure_height = 4; 
savefig(fighandle, Fig_str, figure_width,figure_height)

%% Figure 3, ToE and LoC figures and lag figures

load('ToE_ToC_p_1sigma_Dec2023.mat');

nLon = 180; 
longrid = linspace(1,359,nLon);

for ii = 1:6
    if ii == 1
        tmp_Metric = ToC_wet;
        str = 'E(ToC_{G_W})';
        str1 = 'Fig3a_ToC_wet.pdf';

        clrmp = cbrewer('div','RdBu',50);
        clrmp = flipud(clrmp(26:47,:));
    elseif ii == 2
        tmp_Metric = ToC_dry;
        str = 'E(ToC_{G_D})';
        str1 = 'Fig3b_ToC_dry.pdf';

        clrmp = cbrewer('div','RdBu',50);
        clrmp = clrmp(4:25,:);
    elseif ii == 3
        tmp_Metric = nan(size(ToC_dry)); 
        tmp_Metric(~isnan(ToC_dry)) = ToC_dry(~isnan(ToC_dry));
        tmp_Metric(~isnan(ToC_wet)) = ToC_wet(~isnan(ToC_wet));
        ToC_change = nanmean(tmp_Metric,3);
        str = 'E(ToC_{G_W\cup G_D})';   
        str1 = 'Fig3c_ToC_both.pdf';
        
        clrmp = cbrewer('seq','RdPu',50);
        clrmp = flipud(clrmp(1:45,:));
    elseif ii == 4
        tmp_Metric = ToE_wet;
        str = 'E(ToE|\DeltaP_{G_W, T})';
        str1 = 'Fig3d_ToE_wet.pdf';

        clrmp = cbrewer('div','RdBu',50);
        clrmp = flipud(clrmp(26:47,:));

    elseif ii == 5
        tmp_Metric = ToE_dry;
        str = 'E(ToE|\DeltaP_{G_D, T})';
        str1 = 'Fig3e_ToE_dry.pdf';
        
        clrmp = cbrewer('div','RdBu',50);
        clrmp = clrmp(4:25,:);
    elseif ii == 6
        tmp_Metric = nan(size(ToE_dry)); 
        tmp_Metric(~isnan(ToE_dry)) = ToE_dry(~isnan(ToE_dry));
        tmp_Metric(~isnan(ToE_wet)) = ToE_wet(~isnan(ToE_wet));
        ToE_change = nanmean(tmp_Metric, 3);
        str = 'E(ToE|\DeltaP_{G_W\cup G_D, T})';
        str1 = 'Fig3f_ToE_both.pdf';
        clrmp = cbrewer('seq','RdPu',50);
        clrmp = flipud(clrmp(1:45,:));
    end
    
    % averaging over GCM dimension
    tmp_mean = nanmean(tmp_Metric,3);
      
    fighandle = figure(ii+6); 

    Data_plot = [tmp_mean(91:180,:); tmp_mean(1:90,:)]; % Reshaping for plotting purposes
    globalmapping_withcoords(Data_plot', -latgrid,[longrid(91:180)-360, longrid(1:90)], 0,2020,2100,str, -180, 180, -80, 80, clrmp)

    Fig_str = strcat('Figures/',str1); 
    figure_width = 4; 
    figure_height = 4; 
    savefig(fighandle, Fig_str, figure_width,figure_height)

end


%%  Advanced Warning Time Maps
fighandle = figure; 

AWT_wet = ToE_wet - ToC_wet; 
AWT_dry = ToE_dry - ToC_dry; 
AWT_change = ToE_change - ToC_change;

AWT_wet_mean = nanmean(AWT_wet,3);
AWT_dry_mean = nanmean(AWT_dry,3);
AWT_change_mean = nanmean(AWT_change,3);


subplot(1,3,1); 
clrmp = cbrewer('div','BrBG',50);
clrmp = clrmp(5:45,:); 
Data_plot = [AWT_wet_mean(91:180,:); AWT_wet_mean(1:90,:)];
globalmapping_withcoords(Data_plot', -latgrid,[longrid(91:180)-360, longrid(1:90)],0,0,50,'E(AWT_{G_W})', -180, 180, -80, 80, clrmp(20:end,:))

subplot(1,3,2); 
Data_plot = [AWT_dry_mean(91:180,:); AWT_dry_mean(1:90,:)];
globalmapping_withcoords(Data_plot', -latgrid,[longrid(91:180)-360, longrid(1:90)],0,0,50,'E(AWT_{G_D})', -180, 180, -80, 80, clrmp(20:end,:))

subplot(1,3,3); 
Data_plot = [AWT_change_mean(91:180,:); AWT_change_mean(1:90,:)];
globalmapping_withcoords(Data_plot', -latgrid,[longrid(91:180)-360, longrid(1:90)],0,0,50,'E(AWT_{G_W\cup G_D})', -180, 180, -80, 80, clrmp(20:end,:))


Fig_str = strcat('Figures/Figure3_bottom.pdf'); 
figure_width = 15; 
figure_height = 4; 
savefig(fighandle, Fig_str, figure_width,figure_height)

%% Figure, summary of Signal and Noise for Figure 3.  Noise needed for Giorgi and Bi
% For noise, take the multi-model mean of the variance of decadal average
% (following G&B)
for lon_ii = 1:180
    for lat_ii = 1:88
        for dec = 1:442
           decPLG(dec,:) = nanmean(pr_6x6(lon_ii, lat_ii, dec:dec+9,:),3);
        end

        tmp = decPLG(1:10:100,:); 
        Noise(lon_ii, lat_ii) = sqrt(mean(nanvar(tmp)));
        
    end
    lon_ii
end


yrs = 1850:2300; 
yr1 = find(yrs == 2070); 
yr2 = find(yrs == 2100);
clear Signal
Signal = nanmean(nanmean(pr_6x6(:,:,yr1:yr2,:),3) - nanmean(pr_6x6(:,:,50:100,:),3),4); 


Signal = [Signal(91:180,:); Signal(1:90,:)];
Noise = [Noise(91:180,:); Noise(1:90,:)];
Signal_to_Noise = Signal./Noise;

fighandle = figure; 
globalmapping_withcoords(Signal', -latgrid,[longrid(91:180)-360, longrid(1:90)],0,-15,15,'', -180, 180, -80, 80, clrmp)
Fig_str = strcat('Figures/Signal.pdf'); 
figure_width = 4; 
figure_height = 4; 
savefig(fighandle, Fig_str, figure_width,figure_height)


fighandle = figure
globalmapping_withcoords(Noise', -latgrid,[longrid(91:180)-360, longrid(1:90)],0,0,15,'', -180, 180, -80, 80, clrmp(20:end,:))
Fig_str = strcat('Figures/Noise.pdf'); 
figure_width = 4; 
figure_height = 4; 
savefig(fighandle, Fig_str, figure_width,figure_height)

fighandle = figure
globalmapping_withcoords(Signal_to_Noise', -latgrid,[longrid(91:180)-360, longrid(1:90)],0,-5,5,'', -180, 180, -80, 80, clrmp)
Fig_str = strcat('Figures/Signal_to_Noise.pdf'); 
figure_width = 4; 
figure_height = 4; 
savefig(fighandle, Fig_str, figure_width,figure_height)


% Spatial correlation between each map and MMM Signal_to_Noise 
tmp = [AWT_change_mean(91:180,:); AWT_change_mean(1:90,:)];
[rho, pval] = corr(reshape(abs(Signal_to_Noise),180*88,1),reshape(tmp,180*88,1),'Rows','Complete','Type','Spearman')

tmp = [ToC_change(91:180,:); ToC_change(1:90,:)];
[rho, pval] = corr(reshape(abs(Signal_to_Noise),180*88,1),reshape(tmp,180*88,1),'Rows','Complete','Type','Spearman')

tmp = [ToE_change(91:180,:); ToE_change(1:90,:)];
[rho, pval] = corr(reshape(abs(Signal_to_Noise),180*88,1),reshape(tmp, 180*88,1),'Rows','Complete','Type','Spearman')


%% Figure 5b: ToE using Giorgi and Bi approach
calculate_GiorgiBi = 0; 

yrs = [1850+windowSize:2300];

% Calculating moving signal
clear Signal
for yy = 1:451 - windowSize
    Signal(:,:,yy) = nanmean(nanmean(pr_6x6(:,:,yy:yy+windowSize,:),3) - nanmean(pr_6x6(:,:,50:100,:),3),4); 
end

save('SignalandNoise.mat','Signal','Noise'); 

if calculate_GiorgiBi
    for lat_ii = 1:88
        for lon_ii = 1:180
             ind = find(abs(squeeze(Signal(lon_ii,lat_ii,:))) > Noise(lon_ii,lat_ii), 1);
            if isempty(ind)
                ToE_GiorgiBi(lon_ii, lat_ii) = NaN;
            else
                ToE_GiorgiBi(lon_ii, lat_ii) = yrs(ind); 
            end
        end
    end 

    for lat_ii = 1:88
        for lon_ii = 1:180
             ind = find(abs(squeeze(Signal(lon_ii,lat_ii,:))) > 2*Noise(lon_ii,lat_ii),1);
            if isempty(ind)
                ToE_GiorgiBi2sigma(lon_ii, lat_ii) = NaN;
            else
                ToE_GiorgiBi2sigma(lon_ii, lat_ii) = yrs(ind); 
            end
        end
    end 

    save('GiorgiBi_ToE.mat','ToE_GiorgiBi2sigma','ToE_GiorgiBi','Signal');
else
    load('GiorgiBi_ToE.mat')
end


fighandle = figure

str = 'ToE^{MMM}';
clrmp = cbrewer('seq','RdPu',50);
clrmp = flipud(clrmp(1:45,:));

tmp  = [ToE_GiorgiBi(91:180,:); ToE_GiorgiBi(1:90,:)];
tmp(ind_ocean) = nan;

globalmapping_withcoords(tmp', -latgrid, [longrid(1,91:180)-360,longrid(1,1:90)],0,2020,2100,str, -180, 180, -80, 80, clrmp)

Fig_str = strcat('Figures/Fig5b_ToE_GiorgiBi.pdf');
figure_width = 6; 
figure_height = 6; 
savefig(fighandle, Fig_str, figure_width,figure_height)


%%

LatWeights_ToE = repmat(cosd(latgrid),1,180)';

load('ToE_ToC_p_1sigma_Dec2023.mat')
ToE_wetMMM_1sigma = nanmean(ToE_wet,3); 
ToE_dryMMM_1sigma = nanmean(ToE_dry,3); 

ind = isnan(ToE_wet); 
ToE_Conditional = ToE_wet;
ToE_Conditional(ind) = ToE_dry(ind); 
ToE_Conditional_1sigma = nanmean(ToE_Conditional, 3); 


ToC_wetMMM_1sigma = nanmean(ToC_wet, 3); %Changing name from learn of change to time of confidence
ToC_dryMMM_1sigma = nanmean(ToC_dry, 3); 

ind = isnan(ToC_wet); 
ToC_Conditional = ToC_wet;
ToC_Conditional(ind) = ToC_dry(ind); 
ToC_Conditional_1sigma = nanmean(ToC_Conditional, 3); 


AWT_wetMMM_1sigma = nanmean(ToE_wet - ToC_wet,3); %Changing name from learn of change to time of confidence
AWT_dryMMM_1sigma = nanmean(ToE_dry - ToC_dry,3); 

tmp_wet = ToE_wet - ToC_wet;
tmp_dry = ToE_dry - ToC_dry;
ind = isnan(tmp_wet); 
AWT_Conditional = tmp_wet; 
AWT_Conditional(ind) = tmp_dry(ind); 
AWT_Conditional_1sigma = nanmean(AWT_Conditional, 3); 

load('SignalandNoise.mat'); 
StoN = mean(Signal(:,:,220:250),3)./Noise;

load('ToE_ToC_p_2sigma_Dec2023.mat')
ToE_wetMMM_2sigma = nanmean(ToE_wet,3); 
ToE_dryMMM_2sigma = nanmean(ToE_dry,3); 

ind = isnan(ToE_wet); 
ToE_Conditional = ToE_wet; % conditional means both wet and dry AWT
ToE_Conditional(ind) = ToE_dry(ind); % combining ToE_wet and ToE_dry values for WUD union
ToE_Conditional_2sigma = nanmean(ToE_Conditional, 3); 

ToC_wetMMM_2sigma = nanmean(ToC_wet,3); %Changing name from learn of change to time of confidence
ToC_dryMMM_2sigma = nanmean(ToC_dry,3); 

ind = isnan(ToC_wet); 
ToC_Conditional = ToC_wet; % conditional means both wet and dry AWT
ToC_Conditional(ind) = ToC_dry(ind); % combining ToC_wet and ToC_dry values for WUD union
ToC_Conditional_2sigma = nanmean(ToC_Conditional, 3); 

AWT_wetMMM_2sigma = nanmean(ToE_wet - ToC_wet,3); %Changing name from learn of change to time of confidence
AWT_dryMMM_2sigma = nanmean(ToE_dry - ToC_dry,3); 

tmp_wet = ToE_wet - ToC_wet;
tmp_dry = ToE_dry - ToC_dry;
ind = isnan(tmp_wet); 
AWT_Conditional = tmp_wet;  % conditional means both wet and dry AWT
AWT_Conditional(ind) = tmp_dry(ind); % combining AWT_wet and AWT_dry values for WUD union
AWT_Conditional_2sigma = nanmean(AWT_Conditional, 3); 

yrs = [1850+30:2300];
for yy = 1:421
    ind = find(ToE_GiorgiBi <= yrs(yy)); 
    ToE_fracchange1sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE)); 
    
    ind = find(ToE_GiorgiBi2sigma <= yrs(yy)); 
    ToE_fracchange2sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE)); 
    
    ind = find(logical(ToE_GiorgiBi <= yrs(yy).*(squeeze(Signal(:,:,yy))>0))); 
    ToE_fracwet1sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE)); 
    
    ind = find(logical(ToE_GiorgiBi2sigma <= yrs(yy).*(squeeze(Signal(:,:,yy))>0))); 
    ToE_fracwet2sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE)); 
    
    ind = find(logical(ToE_GiorgiBi <= yrs(yy).*(squeeze(Signal(:,:,yy))<0))); 
    ToE_fracdry1sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE)); 
    
    ind = find(logical(ToE_GiorgiBi2sigma <= yrs(yy).*(squeeze(Signal(:,:,yy))<0))); 
    ToE_fracdry2sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE)); 
    
    ind = find(ToE_wetMMM_1sigma <= yrs(yy)); 
    ToE_wet1sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE)); 
    
    ind = find(ToE_dryMMM_1sigma <= yrs(yy)); 
    ToE_dry1sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE)); 
    
    ind = find(ToE_wetMMM_2sigma <= yrs(yy)); 
    ToE_wet2sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE)); 
    
    ind = find(ToE_dryMMM_2sigma <= yrs(yy)); 
    ToE_dry2sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE)); 
    
    % conditional means both wet and dry AWT
    ind = find(ToE_Conditional_1sigma <= yrs(yy)); 
    ToE_TotalCond_1sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE));
    
    ind = find(ToE_Conditional_2sigma <= yrs(yy)); 
    ToE_TotalCond_2sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE)); 
    
    ind = find(ToC_wetMMM_1sigma <= yrs(yy)); 
    ToC_wet1sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE)); 
    
    ind = find(ToC_dryMMM_1sigma <= yrs(yy)); 
    ToC_dry1sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE)); 
    
    ind = find(ToC_wetMMM_2sigma <= yrs(yy)); 
    ToC_wet2sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE)); 
    
    ind = find(ToC_dryMMM_2sigma <= yrs(yy)); 
    ToC_dry2sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE)); 
    
    % conditional means both wet and dry AWT
    ind = find(ToC_Conditional_1sigma <= yrs(yy)); 
    ToC_TotalCond_1sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE));
    
    ind = find(ToC_Conditional_2sigma <= yrs(yy)); 
    ToC_TotalCond_2sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE)); 

end


for yy = 1:842
    ind = find(AWT_wetMMM_1sigma <= -421+yy); 
    AWT_wet1sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE)); 
    
    ind = find(AWT_dryMMM_1sigma <= -421+yy); 
    AWT_dry1sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE)); 
    
    ind = find(AWT_wetMMM_2sigma <= -421+yy); 
    AWT_wet2sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE)); 
    
    ind = find(AWT_dryMMM_2sigma <= -421+yy); 
    AWT_dry2sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE)); 
    
    ind = find(AWT_Conditional_1sigma <= -421+yy); % conditional means both wet and dry AWT
    AWT_TotalCond_1sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE));
    
    ind = find(AWT_Conditional_2sigma <= -421+yy); 
    AWT_TotalCond_2sigma(yy) = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE)); 
end

%% Figure 4 ToE panels: Plotting cumulative percentage of global surface area that has emerged by year y. 

fighandle = figure; 
subplot(1,3,1); 
p1 = plot(yrs, 100*ToE_fracwet1sigma, '--b','LineWidth', 2); hold on; 
p2 = plot(yrs, 100*ToE_wet1sigma, 'b','LineWidth', 2); hold on; 
p3 = plot(yrs, 100*ToE_fracwet2sigma, '--k','LineWidth', 2); hold on; 
p4 = plot(yrs, 100*ToE_wet2sigma, 'k','LineWidth', 2); hold on; 
xlim([1990, 2150]); ylim([0,100]); ylabel('% Global Surface Area'); 

legend([p1 p2 p3 p4], 'ToE^{MMM} (1\sigma)','E(ToE^C_{G_W}) (1\sigma)', ...
    'ToE^{MMM} (2\sigma)','E(ToE^C_{G_W}) (2\sigma)','Location','northwest')


subplot(1,3,2); 
p1 = plot(yrs, 100*ToE_fracdry1sigma, '--r','LineWidth', 2); hold on; 
p2 = plot(yrs, 100*ToE_dry1sigma, 'r','LineWidth', 2); hold on; 
p3 = plot(yrs, 100*ToE_fracdry2sigma, '--k','LineWidth', 2); hold on; 
p4 = plot(yrs, 100*ToE_dry2sigma, 'k','LineWidth', 2); hold on; 
xlim([1990, 2150]); ylim([0,100]); ylabel('% Global Surface Area'); 

legend([p1 p2 p3 p4], 'ToE^{MMM} (1\sigma)','E(ToE^C_{G_D}) (1\sigma)', ...
    'ToE^{MMM} (2\sigma)','E(ToE^C_{G_D}) (2\sigma)','Location','northwest')


subplot(1,3,3); 
p1 = plot(yrs, 100*ToE_fracchange1sigma, '--','Color',[0.9290, 0.6940, 0.1250],'LineWidth', 2); hold on; 
p2 = plot(yrs, 100*ToE_TotalCond_1sigma,'Color',[0.9290, 0.6940, 0.1250],'LineWidth', 2); hold on; 
p3 = plot(yrs, 100*ToE_fracchange2sigma, '--k','LineWidth', 2); hold on; 
p4 = plot(yrs, 100*ToE_TotalCond_2sigma,'k','LineWidth', 2); hold on; 
xlim([1990, 2150]); ylim([0,100]); ylabel('% Global Surface Area'); 

legend([p1 p2 p3 p4], 'ToE^{MMM} (1\sigma)','E(ToE^C_{G_{W \cup D}}) (1\sigma)', ...
    'ToE^{MMM} (2\sigma)','E(ToE^C_{G_{W \cup D}}) (2\sigma)','Location','southeast')

Fig_str = strcat('Figures/Figure4_ToEpanels.pdf'); 
figure_width = 15; 
figure_height = 3; 
savefig(fighandle, Fig_str, figure_width,figure_height)


% Estimating how much earlier ToEC is than ToE, 
% must fill in nans for ToEMMM so we can take the difference 
% For discussion comparison of ToE^MMM vs ToE^C

ToE_GB_tmp = ToE_GiorgiBi; % If ToE_MMM is nan, set to 2300.
ToE_GB_tmp(isnan(ToE_GB_tmp)) = 2300; 

ToE_C_tmp = ToE_Conditional_1sigma;
ToE_C_tmp(isnan(ToE_C_tmp)) = 2300;

ind = ToE_C_tmp < ToE_GB_tmp;
Global_fraction_ToEC_lessthan_ToEMMM = nansum(LatWeights_ToE(ind))./nansum(nansum(LatWeights_ToE))

ind = isnan(ToE_GiorgiBi);
ToE_noemergenceToEMMM = nansum(LatWeights_ToE(ind).*ToE_Conditional_1sigma(ind))./nansum(nansum(LatWeights_ToE(ind)))

%mean ToEc
ind = ~isnan(ToE_Conditional_1sigma);
ToE_C_mean = nansum(LatWeights_ToE(ind).*ToE_Conditional_1sigma(ind))./nansum(nansum(LatWeights_ToE(ind)))


ind = find(yrs == 2050);
ToE_fracchange1sigma(ind) % fraction of Earth emerged at 1-sigma level for G&B 
ToE_TotalCond_1sigma(ind) % fraction of Earth emerged at 1-sigma level for ToE^C

ind = find(yrs == 2100);
ToE_fracchange1sigma(ind) % fraction of Earth emerged at 1-sigma level for G&B 
ToE_TotalCond_1sigma(ind) % fraction of Earth emerged at 1-sigma level for ToE^C

ind = find(yrs == 2050);
ToE_fracdry1sigma(ind) % fractin of G&B expected to have dried by 2050
ToE_dry1sigma(ind) % fractin of ToEC expected to have dried by 2050

% Fraction could be wetter or drier
ind = find(yrs == 2150);
ToE_wet1sigma(ind)
ToE_dry1sigma(ind)

%% Figure 4 ToC panels: Plotting cumulative percentage of global surface area that is confident by year y. 

fighandle = figure; 
subplot(1,3,1); 
p2 = plot(yrs, 100*ToC_wet1sigma, 'b','LineWidth', 2); hold on; 
p4 = plot(yrs, 100*ToC_wet2sigma, 'k','LineWidth', 2); hold on; 

xlim([1990, 2150]); ylim([0,100]); ylabel('% Global Surface Area'); 
legend([p2 p4], 'E(ToC_{G_W}) (1\sigma)', ...
   'E(ToC_{G_W})  (2\sigma)','Location','southeast')
%title('ToC for \DeltaP > 0 ')

subplot(1,3,2); 
p2 = plot(yrs, 100*ToC_dry1sigma, 'r','LineWidth', 2); hold on; 
p4 = plot(yrs, 100*ToC_dry2sigma, 'k','LineWidth', 2); hold on; 
xlim([1990, 2150]); ylim([0,100]); ylabel('% Global Surface Area'); 
legend([p2 p4], 'E(ToC_{G_D})  (1\sigma)', ...
    'E(ToC_{G_D})  (2\sigma)','Location','southeast')
%title('ToC for \DeltaP < 0 ')

subplot(1,3,3); 
p2 = plot(yrs, 100*ToC_TotalCond_1sigma,'Color',[0.9290, 0.6940, 0.1250],'LineWidth', 2); hold on; 
p4 = plot(yrs, 100*ToC_TotalCond_2sigma,'k','LineWidth', 2); hold on; 
xlim([1990, 2150]); ylim([0,100]); ylabel('% Global Surface Area'); 
legend([p2 p4], 'E(ToC_{G_{W \cup D}})  (1\sigma)', ...
   'E(ToC_{G_{W \cup D}})  (2\sigma)','Location','southeast')
%title('ToC for \DeltaP \neq 0 ')

Fig_str = strcat('Figures/Figure4_ToCpanels.pdf'); 
figure_width = 15; 
figure_height = 3; 
savefig(fighandle, Fig_str, figure_width,figure_height)

ind = find(yrs == 2020);
ToC_TotalCond_1sigma(ind) % fraction of Earth where ToC defined by 2020 at 1-sigma level for G&B 

ind = find(yrs == 2050);
ToC_TotalCond_1sigma(ind) % fraction of Earth where ToC defined by 2050 at 1-sigma level for G&B 

%% Figure 4: AWT panels: PLotting cumulative % global surface area where AWT < y years 
fighandle = figure; 
subplot(1,3,1); 
p2 = plot(-420:421, 100*AWT_wet1sigma, 'b','LineWidth', 2); hold on; 
p4 = plot(-420:421, 100*AWT_wet2sigma, 'k','LineWidth', 2); hold on; 
xlim([0, 100]); ylim([0,100]); ylabel('% Global Surface Area'); xlabel('Advanced Warning Time (years)');
legend([p2 p4], 'E(AWT_{G_W}) (1\sigma)', ...
    'E(AWT_{G_W}) (2\sigma)','Location','southeast')
%title('Expected AWT for \DeltaP > 0 ')

subplot(1,3,2); 
p2 = plot(-420:421, 100*AWT_dry1sigma, 'r','LineWidth', 2); hold on; 
p4 = plot(-420:421, 100*AWT_dry2sigma, 'k','LineWidth', 2); hold on; 
xlim([0, 100]); ylim([0,100]); ylabel('% Global Surface Area'); xlabel('Advanced Warning Time (years)');
legend([p2 p4], 'E(AWT_{G_D}) (1\sigma)', ...
    'E(AWT_{G_D}) (2\sigma)','Location','southeast')
%title('AWT for \DeltaP < 0 ')

subplot(1,3,3); 
p2 = plot(-420:421, 100*AWT_TotalCond_1sigma,'Color',[0.9290, 0.6940, 0.1250],'LineWidth', 2); hold on; 
p4 = plot(-420:421, 100*AWT_TotalCond_2sigma,'k','LineWidth', 2); hold on; 
xlim([0, 100]); ylim([0,100]); ylabel('% Global Surface Area'); xlabel('Advanced Warning Time (years)');
legend([p2 p4], 'E(AWT_{G_{W\cup D}})  (1\sigma)', ...
    'E(AWT_{G_{W\cup D}}) (2\sigma)','Location','southeast')
%title('AWT for \DeltaP \neq 0 ')


Fig_str = strcat('Figures/Figure4_AWTpanels.pdf'); 
figure_width = 15; 
figure_height = 3; 
savefig(fighandle, Fig_str, figure_width,figure_height)


% For last paragraph of paper: 
AWT_wet1sigma(421) % fraction of Earth where AWT is <= 0 
AWT_wet1sigma(end) - AWT_wet1sigma(421) % fraction of Earth with positive AWT for wetter future
AWT_TotalCond_1sigma(421)

AWT_dry1sigma(421) 
% mean advanced warning including negative values
Mean_AWT_1sigma = sum([zeros(1, 421), 1:420].*(AWT_TotalCond_1sigma(2:end)-AWT_TotalCond_1sigma(1:end-1)))
Mean_AWT_2sigma = sum([zeros(1, 421), 1:420].*(AWT_TotalCond_2sigma(2:end)-AWT_TotalCond_2sigma(1:end-1)))

%% Adding CRU data to analysis

% shifting E and W hemispheres to be the same as processed GCMs
Observational_data = cat(1,pr_cru(91:180,:,:), pr_cru(1:90,:,:)); 
yrs_5yr = [1900:5:2275];
modinds = [1:16];

ExpToE_data = nan(180,88); 
Expected_change = nan(180,88); 
Expected_changeCI = nan(180,88); 
for lon_ii = 1:180
    for lat_ii = 1:88
        
        Expected_ToE = nan;
        if lon_ii == 1
            tmp = cat(1,Observational_data(180, lat_ii:lat_ii+2, :),Observational_data(1:2, lat_ii:lat_ii+2, :));
            Obs_tmp = log(squeeze(nanmean(nanmean(tmp,1),2)));
        elseif lon_ii == 180
            tmp = cat(1,Observational_data(179:180, lat_ii:lat_ii+2, :),Observational_data(1, lat_ii:lat_ii+2, :));
            Obs_tmp = log(squeeze(nanmean(nanmean(tmp,1),2)));
        else
            Obs_tmp = log(squeeze(nanmean(nanmean(Observational_data(lon_ii-1:lon_ii+1, lat_ii:lat_ii+2, :),1),2)));  
        end
        
        if isnan(Observational_data(lon_ii, lat_ii, 1))
            ExpToE_data(lon_ii,lat_ii) = nan;
            Expected_change(lon_ii, lat_ii) = nan; 
        else
            
            sigma1 = Max_Likelihood_sigma1(lon_ii,lat_ii); 
            sigma2 = Max_Likelihood_sigma2(lon_ii,lat_ii); 
            sigma3 = Max_Likelihood_sigma3(lon_ii,lat_ii);

            HistP_obs = nanmean(Obs_tmp(1:50)); 
            HistStdP_obs = nanstd(Obs_tmp(1:50)); 

            oos_cru = (1/HistStdP_obs)*(mean(reshape(Obs_tmp,5,24),1) - HistP_obs); 
            sigma_n = nanvar(oos_cru(1:20)); 
    
            GCM_data = squeeze(smoothedMA(lon_ii,lat_ii, 6*5:5:end, :));  % since CRU starts in 1900, we start the prior then as well. 

            yrs_5yr = [1905:5:2275];

            nobs = find(yrs_5yr == 2020);
            t_end = find(yrs_5yr == 2150);
            [prior_mu, prior_cov, post_mu, post_cov] = GPR_fun(oos_cru', nobs, GCM_data, sigma1, sigma2, sigma3, sigma_n);
            
            post_cov = post_cov + sigma_n*eye(size(post_cov));
            upperLimit = post_mu(1:t_end) + diag(post_cov(1:t_end, 1:t_end)).^(0.5); 
            lowerLimit = post_mu(1:t_end) - diag(post_cov(1:t_end, 1:t_end)).^(0.5); 

            First_anom_indx = find(lowerLimit > 0, 1); 
            % If there is a future lowerlimit greater than 0, then
            % check time in which at least 50% of future lowerlimits
            % are greater than zero.  
            if ~isempty(First_anom_indx)
                tmp = cumsum(flipud(lowerLimit) > 0); 
                tmp = tmp'./[1:length(tmp)]; 
                tmp = flipud(tmp');  % Fraction of future values greater than 0. 

                % Finding initial window when at least 50% of values from then on are greater than 0 
                ToC_metric = find(tmp(First_anom_indx : end) > 0.5, 1);  
                if ~isempty(ToC_metric)

                    ToC_metric = ToC_metric + First_anom_indx - 1;
                    Expected_ToE = yrs_5yr(ToC_metric);

                end

            else
                %Same as above but for upper limit less than 0
                First_anom_indx = find(upperLimit < 0, 1); 
                if ~isempty(First_anom_indx)

                    tmp = cumsum(flipud(upperLimit) < 0); 
                    tmp = tmp'./[1:length(tmp)]; 
                    tmp = flipud(tmp');

                    ToC_metric = find(tmp(First_anom_indx :end) > 0.5, 1);

                    if ~isempty(ToC_metric)

                        ToC_metric = ToC_metric + First_anom_indx - 1;
  
                        Expected_ToE = yrs_5yr(ToC_metric);

                    end
                else
                    Expected_ToE = nan; 
                end

            end
            
            ind1 = find(yrs_5yr == 2085);
            
            ExpToE_data(lon_ii,lat_ii) = Expected_ToE;
            
            tmp = exp(post_mu(ind1)*HistStdP_obs+HistP_obs);
            hist_P = exp(Obs_tmp(1:50)); 
            Expected_change(lon_ii, lat_ii) = (tmp - mean(hist_P))/mean(hist_P); 
            
            if post_mu(ind1)>0
                
                sig_tmp = diag(post_cov(1:end,1:end).^(0.5));
                Expected_changeCI(lon_ii, lat_ii) = 1-normcdf(0, post_mu(ind1), sig_tmp(ind1)); 
            
            elseif post_mu(ind1)<0
                
                sig_tmp = diag(post_cov(1:end,1:end).^(0.5));
                Expected_changeCI(lon_ii, lat_ii) = normcdf(0, post_mu(ind1), sig_tmp(ind1)); 
            
            end
       end
        
    end
    lon_ii
end

%% Figure 5 
% c) E(%Delta P|CRU observations)
% d) E(ToE|CRU) 

fighandle = figure

clrmp = cbrewer('seq','RdPu',50);
clrmp = flipud(clrmp(1:45,:));

globalmapping_withcoords([ExpToE_data(91:180,:); ExpToE_data(1:90,:)]', [-87:2:87]',[longrid(1,91:180)-360, longrid(1,1:90)] ,0,2020,2100,'E(ToE|obs)', -180, 180, -80, 80, clrmp)

Fig_str = strcat('Figures/Figure5d_CRU_ToE.pdf'); 
figure_width = 6; 
figure_height = 6; 
savefig(fighandle, Fig_str, figure_width,figure_height)

  
% Figure c, expected delta P with stippling indicating agreement.

Lon_grid_tmp = [Lon_grid(91:180,:) - 360; Lon_grid(1:90,:)];
LAT = [];
LON = [];
for ii = 1:4:180
    for jj = 1:4:88
        if ii>90
        if  Expected_changeCI(ii-90, jj) < 0.84 
            %mapField(ii,jj) = NaN;
            LAT = [LAT; Lat_grid(ii,jj)];
            LON = [LON; Lon_grid_tmp(ii,jj)];
        end
        elseif ii<90
            if  Expected_changeCI(ii+90, jj) < 0.84 
            %mapField(ii,jj) = NaN;
            LAT = [LAT; Lat_grid(ii,jj)];
            LON = [LON; Lon_grid_tmp(ii,jj)];
            end
        end
    end
end

fighandle = figure

clrmp = cbrewer('div','RdBu',50);
clrmp = (1/max(max(clrmp)))*clrmp;
clrmp(clrmp<0) = 0;  

m_proj('robinson','long',[-180 180],'lat',[-80 80]); 
[lat, long] = meshgrid(linspace(-180, 180,360),linspace(-90, 90,360));
TempGrid = griddata(Lon_grid_tmp, Lat_grid, 100*[Expected_change(91:180,:); Expected_change(1:90,:)], lat, long); 

ind_ocean = isnan(100*[Expected_change(91:180,:); Expected_change(1:90,:)]);
Bayesian_change = 100*[Expected_change(91:180,:); Expected_change(1:90,:)];
m_pcolor(lat, long, TempGrid); shading flat;  
hold on; 
m_plot(LON, LAT, '.k');  
hold on; 

m_coast('color',[0.2 .2 0.2]); 
m_grid('linest','none','xticklabels',[],'yticklabels',[]);  
colormap(clrmp);  
colorbar; 
caxis([-30, 30]); 

Fig_str = strcat('Figures/Figure5c_CRU_GPR.pdf'); 
figure_width = 6; 
figure_height = 6; 
savefig(fighandle, Fig_str, figure_width,figure_height)


%% Figure 5. a Level of uncertainty in %Delta P using IPCC method
% Not including agreement in "no change" here, since we are providing a
% parallel with ToC which is only defined if confident about a change.

fighandle = figure
cmin = -30;
cmax = 30; 
cbarsign = 1;
yrs = 1850:1:2300;

ind0 = find((yrs > 1900) & (yrs < 1950));
HistYrs = squeeze(nanmean(pr_6x6(:, :, ind0, :),3));

ind1 = find((yrs > 2070) &  (yrs < 2100));
BaseYrs = squeeze(nanmean(pr_6x6(:, :, ind1, :),3));

Delta1 = 100*(BaseYrs - HistYrs)./HistYrs;
midrange = cat(1, Delta1(179, :, :), Delta1(1, :, :));
tmp = cat(1, Delta1(91:179, :, :), mean(midrange,1));
tmp = cat(1, tmp, Delta1(1:90, :, :));


for mod_ii = 1:16
    tmp2 = tmp(:,:,mod_ii);
    tmp2(ind_ocean) = nan;
    tmp(:,:,mod_ii) = tmp2;
end

Field = tmp;
median_field = nanmedian(Field,3);
nochange = 0; 
Rval = 0.84;

GT = Field > nochange;
LT = Field < nochange;  
%NC = (Field < nochange & Field > -nochange);

LAT = [];
LON = [];
for ii = 1:4:180
    for jj = 1:4:88
        N = sum(~isnan(Field(ii,jj,:)));
        if (median_field(ii,jj)> nochange && sum(GT(ii,jj,:))>(Rval*N)) || (median_field(ii,jj)<-nochange && sum(LT(ii,jj,:))>(Rval*N)) 
            mapField(ii,jj) = median_field(ii,jj);
        %elseif median_field(ii,jj)>-nochange && median_field(ii,jj)< nochange && sum(NC(ii,jj,:))>(Rval*N)
        %    mapField(ii,jj) = median_field(ii,jj);
        elseif  ~isnan(median_field(ii,jj))
            
            LAT = [LAT; Lat_grid(1,jj)];
            LON = [LON; Lon_grid_tmp(ii,1)];
        end
    end
end

m_proj('robinson','long',[-180 180],'lat',[-80 80]); 
[lat, long] = meshgrid(linspace(-180, 180,360),linspace(-90, 90,360));
TempGrid = griddata(Lon_grid_tmp, Lat_grid, nanmean(Field,3), lat, long); 

MMM_change = nanmean(tmp,3);
m_pcolor(lat, long, TempGrid); shading flat;  
hold on; 
m_plot(LON, LAT, '.k');  
hold on; 

clrmp = cbrewer('div','RdBu',50);
clrmp = (1/max(max(clrmp)))*clrmp;
clrmp(clrmp<0) = 0;

m_coast('color',[0.2 .2 0.2]); 
m_grid('linest','none','xticklabels',[],'yticklabels',[]);  
colormap(clrmp);  colorbar ; caxis([-30, 30]); 


Fig_str = strcat('Figures/Fig5a_IPCC.pdf');
figure_width = 6; 
figure_height = 6; 
savefig(fighandle, Fig_str, figure_width,figure_height)


%% Difference between Bayesian and MMM expected change

diff_methods = MMM_change - Bayesian_change;
TempGrid = griddata(Lon_grid_tmp, Lat_grid, diff_methods, lat, long); 

fighandle = figure
m_pcolor(lat, long, TempGrid); shading flat;  
hold on; 

clrmp = cbrewer('div','RdBu',50);
clrmp = (1/max(max(clrmp)))*clrmp;
clrmp(clrmp<0) = 0;

m_coast('color',[0.2 .2 0.2]); 
m_grid('linest','none','xticklabels',[],'yticklabels',[]);  
colormap(clrmp);  colorbar ; caxis([-20, 20]); 


Fig_str = strcat('Figures/Fig5e.pdf');
figure_width = 6; 
figure_height = 6; 
savefig(fighandle, Fig_str, figure_width,figure_height)

%% Uncertainty with MMM direction change  
% repeating above for uncertainty estimate between methods conditioned on CRU
% Reporting value on global model agreement.

Delta1 = 100*(BaseYrs - HistYrs)./HistYrs;
midrange = cat(1, Delta1(179, :, :), Delta1(1, :, :));
tmp = cat(1, Delta1(91:179, :, :), mean(midrange,1));
Delta1  = cat(1, tmp, Delta1(1:90, :, :));

Delta1_noocean = Delta1;

for mod_ii = 1:16
    tmp2 = Delta1_noocean(:,:,mod_ii);
    tmp2(ind_ocean) = nan;
    Delta1_noocean(:,:,mod_ii) = tmp2;
end


% Number of models that agree on direction of change globally
Field = Delta1;
 [l1,l2,N] = size(Field);

median_field = nanmedian(Field,3);
nochange = 0;
Rval = 0.84;

GT = Field > 0;
LT = Field < -0;  
NC = (Field < nochange & Field > -nochange);

uncertain_mmm = zeros(l1, l2);
for ii = 1:l1
    for jj = 1:l2
        N = sum(~isnan(Field(ii,jj,:)));
        if (median_field(ii,jj)> nochange && sum(GT(ii,jj,:))>(Rval*N)) || (median_field(ii,jj)<-nochange && sum(LT(ii,jj,:))>(Rval*N)) 
            mapField(ii,jj) = median_field(ii,jj);
        elseif median_field(ii,jj)>-nochange && median_field(ii,jj)< nochange && sum(NC(ii,jj,:))>(Rval*N)
            mapField(ii,jj) = median_field(ii,jj);
        else 
            uncertain_mmm(ii,jj) = 1; 
        end
    end
end

% Number of models that agree on direction of change over land

Field = Delta1_noocean;
 [l1,l2,N] = size(Field);

median_field = nanmedian(Field,3);
nochange = 0;
Rval = 0.84;

GT = Field > nochange;
LT = Field < -nochange;  
NC = (Field < nochange & Field > -nochange); %this term is redundant


uncertain_mmm_noocean = zeros(l1, l2);
for ii = 1:l1
    for jj = 1:l2
        
        N = sum(~isnan(Field(ii,jj,:)));
        if ii <= 90
            nochange = Noise(ii+90,jj);
        elseif ii > 90
            nochange = Noise(ii-90,jj);
        end

        if (median_field(ii,jj)> nochange && sum(GT(ii,jj,:))>(Rval*N)) || (median_field(ii,jj)<-nochange && sum(LT(ii,jj,:))>(Rval*N)) 
            mapField(ii,jj) = median_field(ii,jj);
        elseif median_field(ii,jj)>-nochange && median_field(ii,jj)< nochange && sum(NC(ii,jj,:))>(Rval*N)
            mapField(ii,jj) = median_field(ii,jj);
        else 
            uncertain_mmm_noocean(ii,jj) = 1; 
        end
    end
end



LatWeights = cosd(repmat([-87:2:87],180,1));


% Multi-model mean fraction of surface area where <84% of GCM agree on
% direction of change
uncertain_globalextent = sum(sum(LatWeights.*uncertain_mmm))/sum(sum(LatWeights));

Lat_Weights_land = LatWeights;
Lat_Weights_land(ind_ocean) = nan;
% Same as above but over land
uncertain_landextent = nansum(nansum(Lat_Weights_land.*uncertain_mmm))/nansum(nansum(Lat_Weights_land));

% fraction of land area where predictive model has less than 84% confidence
% in direction of change conditioned on CRU.
GPR_certainty_give_CRU = [Expected_changeCI(91:180,:); Expected_changeCI(1:90,:)] < 0.84;
uncertain_landextent_GPRupdated = nansum(nansum(Lat_Weights_land.*GPR_certainty_give_CRU))/nansum(nansum(Lat_Weights_land));




