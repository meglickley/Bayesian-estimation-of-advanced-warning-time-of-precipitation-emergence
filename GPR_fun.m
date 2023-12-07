function [prior_mu, prior_cov, post_mu, post_cov] = GPR_fun(oos_gcm, nobs, GCM_data, sigma1, sigma2, sigma3, sigma_n)

    tmp = GCM_data;
    tmp = [tmp(1:nobs,:); tmp];

    SmoothedCov = nancov(tmp');
    SmoothedMean = mean(tmp');

    time_diff = abs([1:nobs, 1:length(GCM_data)] - [1:nobs, 1:length(GCM_data)]'); % for the length scale in Kernel
    Noise_matrix = sigma1^2 * exp(-(1/(2*sigma2^2))* time_diff.^2) + sigma3^2*eye(length(SmoothedMean));

    Kernel = SmoothedCov + Noise_matrix;
    
    y = oos_gcm(1:nobs);
    mu = SmoothedMean(1:nobs); 
    mu_s = SmoothedMean(nobs+1:end); 
    Noise = sigma_n*eye(nobs, nobs);
    K_o = Kernel(1:nobs, 1:nobs); 
    K_s = Kernel(nobs+1:end, nobs+1:end); 
    K_os = Kernel(1:nobs, nobs+1:end);%Nobs+1:end); 

    diff = (y - mu'); 
    if size(diff,1) == 1
        diff = diff';
    end
    
    f = mu_s' + K_os'*inv(K_o+Noise)*diff;
    Sigma_s = K_s - K_os'*inv(K_o+Noise)*K_os;


    prior_mu = mu_s;
    prior_cov = K_s;
    post_mu = f; 
    post_cov = Sigma_s;
    
end