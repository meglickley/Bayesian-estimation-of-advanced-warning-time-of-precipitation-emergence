function [ToC ToC_sign ToE ToE_sign] = ToC_ToE_function(oos_gcm, sigma1, sigma2, sigma3, sigma_n, GCM_data, yrs, sig)
    
    % This function reads in the oos_gcm and GCM_data and estimates the GPR
    % prior.  It then iteratively updates the posterior for each additional
    % year in the oos_gcm time series.  It identifies the ToC, the first year in
    % which the posterior is at some point different than zero for at least
    % half of future time periods.  It then identifies the ToE, the first
    % year in which the posterior is significantly different from zero
    % during the final year of updating
    %
    % Input: 
    %   oos_gcm: time series of 5-yr anomalies of out of sample gcm 
    %   sigma_n: noise in observations
    %   GCM_data: GCM ensemble (without the oos_gcm member) used to construct prior
    %   yrs: corresponding years for oos_gcm and GCM_data time series
    %
    % Output:
    %   ToC: Time of Confidence, the year in which the posterior (at some point in the future) is statistically different from zero
    %   ToC_sign: +1 = wetter future, -1 = drier future
    %   ToE:  The time of emergence,the year in which the posterior is statistically different from zero
    %   ToE_sign: +1= wetter future, -1 = drier future
    %
    %       

   ii = 0; 
   confidence = 0; 
   ToC = NaN;
   ToC_sign = 0; 
   ToE_sign = 0;          
   ToE = NaN;
   fraction_past_emergence = 0.75; 

   
   last_yr = find(yrs == 2150);

    while confidence == 0
         
        % Starting in 1975
        ii = ii + 1;
        nobs = find(yrs == 1975 + 5*(ii-1));

        [prior_mu, prior_cov, post_mu, post_cov] = GPR_fun(oos_gcm, nobs, GCM_data, sigma1, sigma2, sigma3, sigma_n);
        
        post_noisy_cov = post_cov + sigma_n*eye(size(post_cov));

        upperLimit = post_mu(1:last_yr) + sig*diag(post_noisy_cov(1:last_yr, 1:last_yr).^(0.5)); 
        lowerLimit = post_mu(1:last_yr) - sig*diag(post_noisy_cov(1:last_yr, 1:last_yr).^(0.5)); 


        % if LoC has not yet been defined, then identify first year in
        % which lower limit is greater than zero (for wetter future)
        % Or upperLimit is less than zero
        yr_emerge = find(lowerLimit> 0, 1); 
           
        if ~isempty(yr_emerge)

            % time-dependent fraction of future values greater than 0. 
            % each value, x, in the tmp vector indicates fraction of values
            % following x that is greater than zero. 
            tmp = cumsum(flipud(lowerLimit) > 0); 
            tmp = tmp'./[1:length(tmp)]; % flipping the 
            tmp = flipud(tmp');  % Fraction of future values greater than 0. 

            % Finding initial window when at least 50% of values from then on are greater than 0 
            ToC_metric = find(tmp(yr_emerge : end) > fraction_past_emergence, 1);  
            
            if ~isempty(ToC_metric)
                
                confidence = 1;
                ToC_sign = 1;

                ToC = yrs(nobs); 
            end

        else
            %Same as above but for upper limit less than 0
            yr_emerge = find(upperLimit < 0, 1); 
            
            if ~isempty(yr_emerge)
    
                % time-dependent fraction of future values greater than 0. 
                % each value, x, in the tmp vector indicates fraction of values
                % following x that is less than zero. 
                tmp = cumsum(flipud(upperLimit) < 0); 
                tmp = tmp'./[1:length(tmp)]; 
                tmp = flipud(tmp');
    
                ToC_metric = find(tmp(yr_emerge :end) > fraction_past_emergence, 1);
    
                if ~isempty(ToC_metric)
    
                    confidence = 1;
                    ToC_sign = -1; 
    
                    ToC = yrs(nobs); 
    
                end
             end

        end


    
        
% ToE definition if ToE is the first t_obs that is significant
    % Once ToC is confident, identify ToE
%         if confidence == 1
%             y_start = nobs; 
%             y_end = find(yrs == 2150);
% 
%             for nobs = y_start:y_end
% 
%                 [prior_mu, prior_cov, post_mu, post_cov] = GPR_fun(oos_gcm, nobs, GCM_data, sigma1, sigma2, sigma3, sigma_n);
%     
%                 post_noisy_cov = post_cov + sigma_n*eye(size(post_cov));
%                 
%                 upperLimit = post_mu(1:last_yr) + sig*diag(post_noisy_cov(1:last_yr, 1:last_yr).^(0.5)); 
%                 lowerLimit = post_mu(1:last_yr) - sig*diag(post_noisy_cov(1:last_yr, 1:last_yr).^(0.5)); 
%     
%                 % ToE requires that climate change is significanlty different
%                 % in time of observation
%                 if lowerLimit(nobs) > 0 && mean(lowerLimit(nobs:last_yr) > 0) > 0.5 
%                     ToE = yrs(nobs);
%                     ToE_sign = 1; 
%                     break;  
%                 elseif upperLimit(nobs) < 0 && mean(upperLimit(nobs:last_yr) < 0) > 0.5 
%                     ToE = yrs(nobs);
%                     ToE_sign = -1; 
%                     break;
%                 end
%             end
%             
% 
%         end

        if confidence == 0 && nobs == find(yrs == 2150)
            confidence = 1; 
        end
    end

    % ToE definition from paper, identified for posterior conditioned on
    % entire time series
    nobs = find(yrs == 2150);

    [prior_mu, prior_cov, post_mu, post_cov] = GPR_fun(oos_gcm, nobs, GCM_data, sigma1, sigma2, sigma3, sigma_n);
    
    post_noisy_cov = post_cov + sigma_n*eye(size(post_cov));

    upperLimit = post_mu(1:last_yr) + sig*diag(post_noisy_cov(1:last_yr, 1:last_yr).^(0.5)); 
    lowerLimit = post_mu(1:last_yr) - sig*diag(post_noisy_cov(1:last_yr, 1:last_yr).^(0.5)); 


    % if LoC has not yet been defined, then identify first year in
    % which lower limit is greater than zero (for wetter future)
    % Or upperLimit is less than zero
    yr_emerge = find(lowerLimit > 0, 1); 

    if ~isempty(yr_emerge)

        % time-dependent fraction of future values greater than 0. 
        % each value, x, in the tmp vector indicates fraction of values
        % following x that is greater than zero. 
        tmp = cumsum(flipud(lowerLimit) > 0); 
        tmp = tmp'./[1:length(tmp)]; 
        tmp = flipud(tmp');  

        % Finding initial window when at least 50% of values from then on are greater than 0 
        ToE_metric = find(tmp(yr_emerge : end) > fraction_past_emergence, 1);  

        if ~isempty(ToE_metric)
            
            ToE_sign = 1;
            ToE = yrs(yr_emerge + ToE_metric  - 1); 
        end

    else
        %Same as above but for upper limit less than 0
        yr_emerge = find(upperLimit < 0, 1); 
        
        if ~isempty(yr_emerge)

            % time-dependent fraction of future values greater than 0. 
            % each value, x, in the tmp vector indicates fraction of values
            % following x that is less than zero. 
            tmp = cumsum(flipud(upperLimit) < 0); 
            tmp = tmp'./[1:length(tmp)]; 
            tmp = flipud(tmp');

            ToE_metric = find(tmp(yr_emerge : end) > fraction_past_emergence, 1);

            if ~isempty(ToE_metric)

                ToE_sign = -1; 
                ToE = yrs(yr_emerge + ToE_metric  - 1); 

            end
        end


    end

end