function preprocess_EFSubmission

% Processing GCM and CRU data to a 6x6 grid resolution and yearly timestep 
nModels = 16;
nLat = 90;
nLon = 180; 
nYrs = 450;

pr_rcp45 = NaN(nLon,nLat,nYrs*12,nModels);
HomeDir = pwd;

yrs = 1850.08:1/12:2300;
windowSize = 30;


% GCM name; 
GCM(1).name = 'bcc-csm1-1';  Yr1_Ind(1) = 1; 
GCM(2).name = 'CanESM2';     Yr1_Ind(2) = 1; 
GCM(3).name = 'CESM1-CAM5';  Yr1_Ind(3) = 1;
GCM(4).name = 'CSIRO-Mk3L-1-2'; Yr1_Ind(4) = 13;
GCM(5).name = 'CSIRO-Mk3-6-0';  Yr1_Ind(5) = 1; 
GCM(6).name = 'CNRM-CM5';       Yr1_Ind(6) = 1; 
GCM(7).name = 'FGOALS-g2';      Yr1_Ind(7) = 1; 
GCM(8).name = 'GFDL-CM3';       Yr1_Ind(8) = 121; 
GCM(9).name = 'GISS-E2-H';      Yr1_Ind(9) = 1; 
GCM(10).name = 'GISS-E2-R';     Yr1_Ind(10) = 1; 
GCM(11).name = 'HadGEM2-ES';    Yr1_Ind(11) = 120;
GCM(12).name = 'IPSL-CM5A-LR';  Yr1_Ind(12) = 1;
GCM(13).name = 'IPSL-CM5A-MR';  Yr1_Ind(13) = 1;
GCM(14).name = 'MIROC-ESM';     Yr1_Ind(14) = 1; 
GCM(15).name = 'MPI-ESM-LR';    Yr1_Ind(15) = 1; 
GCM(16).name = 'NorESM1-M';     Yr1_Ind(16) = 1; 

latgrid = flipud(linspace(-89,89,nLat)'); 
longrid = linspace(1,359,nLon);


cd GetData/cmip5/historical/Amon/pr


mmpersec_to_permonth = 60*60*24*365/12;


for gg = 1:16

    str = strcat('pr_Amon_', GCM(gg).name,'*');
    files = dir(str);
    filenum = 1;
    pr_tmp  = ncread(files(filenum).name, 'pr');
    lat_tmp = ncread(files(filenum).name, 'lat');
    lon_tmp = ncread(files(filenum).name, 'lon');
    time_tmp = ncread(files(filenum).name, 'time'); %Days since 1850


    t1 = Yr1_Ind(gg);
    t2 = t1 + size(pr_tmp,3) - 1;

    lat_tmp = ncread(files(filenum).name,'lat');
    lon_tmp = ncread(files(filenum).name,'lon');
    lon_tmp = [lon_tmp(end-5:end)-360;lon_tmp; lon_tmp(1:5)+360];
    pr_tmp = cat(1, pr_tmp(end-5:end,:,:), pr_tmp);
    pr_tmp = cat(1, pr_tmp, pr_tmp(1:5,:,:));

    [lon_in, lat_in, time_in] = meshgrid(lat_tmp,lon_tmp,1:1:size(pr_tmp,3));
    [lon_2,  lat_2, time_2] = meshgrid(linspace(-90, 90,90),linspace(1, 360,180),1:1:size(pr_tmp,3));

    pr_tmp2 = interp3(lon_in, lat_in, time_in, pr_tmp, lon_2, lat_2, time_2);
    pr_rcp45(:, :, t1:t2, gg) = mmpersec_to_permonth*pr_tmp2;
    
    if length(files) > 1
        tend = t2;
        
        for filenum = 2:length(files)
            pr_tmp  = ncread(files(filenum).name, 'pr');
            lat_tmp = ncread(files(filenum).name,'lat');
            lon_tmp = ncread(files(filenum).name,'lon');
        
        
            t1 = tend+1;
            t2 = size(pr_tmp,3);
            tend = tend+t2;

            % Ensuring that lon range is wider than interpolated space
            lon_tmp = [lon_tmp(end-5:end)-360;lon_tmp; lon_tmp(1:5)+360];
            pr_tmp = cat(1,pr_tmp(end-5:end,:,:),pr_tmp);
            pr_tmp = cat(1,pr_tmp,pr_tmp(1:5,:,:));

            [lon_in, lat_in,time_in] = meshgrid(lat_tmp,lon_tmp,1:1:size(pr_tmp,3));
            [lon_2, lat_2,time_2] = meshgrid(linspace(-90, 90,90),linspace(1, 360,180),1:1:size(pr_tmp,3));

            pr_tmp2 = interp3(lon_in, lat_in, time_in, pr_tmp, lon_2, lat_2, time_2);

            pr_rcp45(:, :, t1:tend, gg) = mmpersec_to_permonth*pr_tmp2;
        end
    end
end


cd(HomeDir)
cd GetData/cmip5/rcp45/Amon/pr
for gg = 1:16

    str = strcat('pr_Amon_', GCM(gg).name,'*');
    files = dir(str);
    filenum = 1;
    pr_tmp  = ncread(files(filenum).name, 'pr');
    lat_tmp = ncread(files(filenum).name, 'lat');
    lon_tmp = ncread(files(filenum).name, 'lon');
    
    if gg == 11
        t1 = 156*12;
    else
        t1 = 156*12+1;
    end
    
    t2 = t1 + size(pr_tmp,3) - 1;

    lat_tmp = ncread(files(filenum).name,'lat');
    lon_tmp = ncread(files(filenum).name,'lon');
    lon_tmp = [lon_tmp(end-5:end)-360;lon_tmp; lon_tmp(1:5)+360];
    pr_tmp = cat(1, pr_tmp(end-5:end,:,:), pr_tmp);
    pr_tmp = cat(1, pr_tmp, pr_tmp(1:5,:,:));

    [lon_in, lat_in, time_in] = meshgrid(lat_tmp,lon_tmp,1:1:size(pr_tmp,3));
    [lon_2,  lat_2, time_2] = meshgrid(linspace(-90, 90,90),linspace(1, 360,180),1:1:size(pr_tmp,3));

    pr_tmp2 = interp3(lon_in, lat_in, time_in, pr_tmp, lon_2, lat_2, time_2);
    pr_rcp45(:, :, t1:t2, gg) = mmpersec_to_permonth*pr_tmp2;
    
    if length(files) > 1
        tend = t2;
        
        for filenum = 2:length(files)
            pr_tmp  = ncread(files(filenum).name, 'pr');
            lat_tmp = ncread(files(filenum).name,'lat');
            lon_tmp = ncread(files(filenum).name,'lon');
        
        
            t1 = tend+1;
            t2 = size(pr_tmp,3);
            tend = tend+t2;

            % Ensuring that lon range is wider than interpolated space
            lon_tmp = [lon_tmp(end-5:end)-360;lon_tmp; lon_tmp(1:5)+360];
            pr_tmp = cat(1,pr_tmp(end-5:end,:,:),pr_tmp);
            pr_tmp = cat(1,pr_tmp,pr_tmp(1:5,:,:));

            [lon_in, lat_in,time_in] = meshgrid(lat_tmp,lon_tmp,1:1:size(pr_tmp,3));
            [lon_2, lat_2,time_2] = meshgrid(linspace(-90, 90,90),linspace(1, 360,180),1:1:size(pr_tmp,3));
            if gg == 11 && filenum == 13
                pr_tmp2 = interp2(lon_in, lat_in,  pr_tmp, lon_2, lat_2);
                pr_rcp45(:, :, t1, gg) = mmpersec_to_permonth*pr_tmp2;
            else
                pr_tmp2 = interp3(lon_in, lat_in, time_in, pr_tmp, lon_2, lat_2, time_2);
                pr_rcp45(:, :, t1:tend, gg) = mmpersec_to_permonth*pr_tmp2;
            end
        end
    end
end


cd(HomeDir)


pr_6x6 = nan(nLon, nLat-2, 451, nModels);
for lat_ii = 1 : 88

    tmp = nanmean(pr_rcp45(:, lat_ii:lat_ii+2, :,:), 2);
    tmp2 = cat(1, tmp, tmp(1,:,:,:)); 
    tmp2 = cat(1, tmp(end,:,:,:), tmp2); 

    for lon_ii = 1:180

        for yy = 1:451
            pr_6x6(lon_ii, lat_ii, yy,:) = nanmean(nanmean(tmp2(lon_ii:lon_ii+2, 1, 12*(yy-1)+1:12*yy, :),1), 3); 
        end

    end
    lat_ii
end

save('pr_6x6_grid.mat','pr_6x6', 'pr_rcp45', '-v7.3');


%%  loading CRU data
% Harris, I., Osborn, T.J., Jones, P. et al. Version 4 of the CRU TS monthly 
% high-resolution gridded multivariate climate dataset. Sci Data 7, 109 (2020). 
% https://doi.org/10.1038/s41597-020-0453-3

cd('CRU_v4.05')
files = dir('*.nc');

for ii = 1:12
    tmp = ncread(files(ii).name, 'pre'); 
    for yy = 1:10
        tmp1 = mean(tmp(:,:,12*(yy-1)+1:12*yy),3); 
        pr_highres(:,:,10*(ii-1)+yy) = tmp1; 
        for lon_ii = 1:180
            for lat_ii = 1:90
                yy2 = 10*(ii-1)+yy;
                pr_cru(lon_ii,lat_ii,yy2) = nanmean(nanmean(tmp1(4*(lon_ii-1)+1:4*lon_ii, 4*(lat_ii-1)+1:4*lat_ii),1),2);
            end
        end
    end
end
lat_cru = ncread(files(1).name, 'lat'); 
lon_cru = ncread(files(1).name, 'lon'); 
lat_cru = 0.5*(lat_cru(2:4:end) + lat_cru(3:4:end)); 
lon_cru = 0.5*(lon_cru(2:4:end) + lon_cru(3:4:end)); 
cd(HomeDir)

save('CRU_regridded.mat','pr_cru','lon_cru','lat_cru'); 
