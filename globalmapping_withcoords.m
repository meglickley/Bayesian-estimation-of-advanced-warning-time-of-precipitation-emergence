function globalmapping_withcoords(Field,Lat,Lon,landcovered,cmin,cmax,str,lonmin,lonmax,latmin,latmax, clrmp)

%load('colorscheme.mat','ColorMeg')
%clrmp = cbrewer('div','RdYlBu',60);
%clrmp = colormap(parula(50)); 
%clrmp = (1/max(max(clrmp)))*clrmp;
%clrmp = flipud(clrmp);
%clrmp = cbrewer('div','BrBG',45);
%clrmp = cbrewer('div','PRGn',60);
%clrmp = interp1([1, 45],[clrmp(1,:)', clrmp(45,:)']',[1:45]); 
%clrmp = clrmp(1:20,:);
% clrmp = cbrewer('seq','Purples',45);
% clrmp = clrmp(10:end,:);
if max(max(clrmp))>1
 clrmp = (1/max(max(clrmp)))*clrmp;
end
 clrmp(clrmp<0) = 0;
 %clrmp = flipud(clrmp);
%clrmp = clrmp(1:36,:);
%clrmp = colormap(jet);
% projection could be 'miller' or 'robinson
m_proj('robinson','long',[lonmin lonmax],'lat',[latmin latmax]); 
[lat, long] = meshgrid(linspace(lonmin, lonmax,360),linspace(latmin, latmax,360));
TempGrid = griddata(Lon,Lat, Field, lat, long); m_pcolor(lat, long, TempGrid); shading flat;  
if landcovered == 1
    m_coast('patch',[0.4 .4 0.4],'edgecolor',[0.2 .2 0.2]); 
    m_grid('linest','none','xticklabels',[],'yticklabels',[]);  
    colormap(clrmp);  
    c = colorbar;
    %c.FontSize = 12;
else
    m_coast('color',[0.2 .2 0.2]); 
    m_grid('linest','none','xticklabels',[],'yticklabels',[]);   colormap(clrmp);
    c = colorbar;
    %c.FontSize = 12;
end

if ~isnan(cmin)
    caxis([cmin, cmax]); 
end
title(str,'fontsize',12);

