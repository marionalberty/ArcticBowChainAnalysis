% Bchain_MMP1_spectra.m

%

clear all; close all; clc
warning off

set(0,'defaultaxesfontsize',18,'defaultaxeslinewidth',1,...
  'defaultlinelinewidth',2,'defaultpatchlinewidth',1,...
  'defaultFigureColor','white')


%% Set paths and directories
% Personnal paths:
driveName ='/Users/marionsofiaalberty/MATLAB/Arctic/ArcticMix/';
dataPath = 'bowchain/Data/Gridded/';
figBC = 'bowchain/Figures/';
figMMP = 'MMP/Figures/';


%% Load bowchain data
load([driveName dataPath 'MMPrepeat1.mat'])


%% Interp bowchain data onto regular grid
% Regular vertical grid
BC.z = (-3:-1:-20)';
% Initialize intermediate temperature grid
t_inter = nan(numel(BC.z),numel(data.dn));
% Interpolate in z
for i = 1:numel(data.dn)
  t_inter(:,i) = interp1(data.z(:,i),data.t(:,i),BC.z);
end


%% Get distance along the transect for bowchain
xxx = m_lldist(data.lon(1,:),data.lat(1,:));
distB = [0; cumsum(xxx)]';
dx_B = 0.001;
BC.x = 0:dx_B:max(distB);
% Interpolate onto regular distance grid
BC.t = interp1(distB,t_inter',BC.x)';
% Interpolate over nans
BC.t = naninterp(BC.t')';

clear t_inter
% Housekeeping
BC.data = data;


%% Load and subset MMP data
% Load file
load([driveName 'MMP/Data/MMPgrid.mat'])
% Get overlapping data
i_dep = MMP.dnum >= data.dn(1) & MMP.dnum <= data.dn(end);
N_dnm = numel(MMP.dnum);
% Go through MMP structure and truncate all variables
vlist = fieldnames(MMP);
for i = 1:numel(vlist)
  eval(sprintf('[~,n] = size(MMP.%s);',vlist{i}));
  if n == N_dnm
    % Truncate
    eval(sprintf('MMP.%s = MMP.%s(:,i_dep);',vlist{i},vlist{i}))
  end
end
data = MMP; clear MMP
MMP.data = data; clear data vlist N_dnm n i_dep

%% Get distance along transect for MMP data
dist = interp1(BC.data.dn,distB,MMP.data.dnum);
dx_M = mean(diff(dist));
MMP.x = 0:dx_M:max(dist);
MMP.z = (5:0.25:100)';
% Interpolate onto regular distance grid
MMP.t = interp2(dist,MMP.data.z',MMP.data.t,MMP.x,MMP.z);
MMP.s = interp2(dist,MMP.data.z',MMP.data.s,MMP.x,MMP.z);
MMP.sgth = interp2(dist,MMP.data.z',MMP.data.sgth,MMP.x,MMP.z);
MMP.eps = interp2(dist,MMP.data.z',MMP.data.eps,MMP.x,MMP.z,'nearest');
% Housekeeping
clear xxx dist
% Interpolate over nans
MMP.t = naninterp(MMP.t')';
MMP.s = naninterp(MMP.s')';
MMP.sgth = naninterp(MMP.sgth')';


%% Clean up eps (quick and dirty)
% Initialize mask
i_mask = ones(size(MMP.eps));
for i = 1:numel(MMP.x)
  [~,PP] = islocalmin(log10(MMP.eps(1:41,i)));
  [~,i_PP] = max(PP);
  i_mask(1:i_PP-1,i) = nan;
end

%% Calculte spectra for bowchain
% Number of points for fft
N_fft = 10/dx_B; % 10 km segments
% Initialize final spectra matrix
spec_bc = nan(numel(BC.z),N_fft/2);
for i = 1:numel(BC.z)
  % Make half-overlapping segments
  t_seg = segment(BC.t(i,:)',N_fft,0.5,0)';
  % Get number of segments
  N_seg = numel(t_seg(:,1));
  % Initialize segment spectra matrix
  spec_seg = nan(N_seg,N_fft/2);
  % Calc spectra for each segment
  for j = 1:N_seg
    [spec_seg(j,:),k_bc] = spectrum(t_seg(j,:),BC.x(1:N_fft));
  end
  % Average spectra for each depth
  spec_bc(i,:) = mean(spec_seg);
end


%% Calculte spectra for MMP
% Number of points for fft
N_fft = floor(15/dx_M); % 10 km segments
% Ensure round number
rm_extra = mod(N_fft,2);
N_fft = N_fft-rm_extra;
% Initialize final spectra matrix
spec_mmp = nan(numel(MMP.z),N_fft/2);
spec_S = spec_mmp;
spec_SGTH = spec_mmp;
for i = 1:numel(MMP.z)
  % Make half-overlapping segments
  t_seg = segment(MMP.t(i,:)',N_fft,0.5,0)';
  s_seg = segment(MMP.s(i,:)',N_fft,0.5,0)';
  sgth_seg = segment(MMP.sgth(i,:)',N_fft,0.5,0)';
  % Get number of segments
  N_seg = numel(t_seg(:,1));
  % Initialize segment spectra matrix
  spec_seg = nan(N_seg,N_fft/2);
  spec_S_seg = spec_seg;
  spec_SGTH_seg = spec_seg;
  % Calc spectra for each segment
  for j = 1:N_seg
    [spec_seg(j,:),k_mmp] = spectrum(t_seg(j,:),MMP.x(1:N_fft));
    [spec_S_seg(j,:),~] = spectrum(s_seg(j,:),MMP.x(1:N_fft));
    [spec_SGTH_seg(j,:),~] = spectrum(sgth_seg(j,:),MMP.x(1:N_fft));
  end
  % Average spectra for each depth
  spec_mmp(i,:) = mean(spec_seg);
  spec_S(i,:) = mean(spec_S_seg);
  spec_SGTH(i,:) = mean(spec_SGTH_seg);
end


%% Plot temperature from bc and mmp
figure('position',[0 0 1600 600])
pcolor(MMP.x,MMP.z,MMP.t)
hold on
pcolor(BC.x,-1*BC.z,BC.t)
shading flat
axis ij
xlim([0 max(MMP.x)])
ylim([min(-1*BC.z) 75])
colormap(jet)
freezeColors;
colorbar
[Cs,hs] = contour(MMP.x,MMP.z,MMP.s,24:32,'m','linewidth',2);
clabel(Cs,hs,'fontsize',14,'color','m','labelspacing',1000)
xlabel('Along Track Distance [km]')
ylabel('Depth [m]')
title('Bow Chain and MMP Temperature [^\circ C] Sep. 10, 2015')

print([driveName figMMP 'temperature_MMPrepeat1.png'],'-dpng')


%% Plot Dissipation from mmp
figure('position',[0 0 1600 600])
pcolor(MMP.x,MMP.z,log10(MMP.eps).*i_mask)
hold on
shading flat
axis ij
xlim([0 max(MMP.x)])
ylim([min(-1*BC.z) 75])
colormap(paruly(8))
caxis([-10 -2])
freezeColors;
colorbar
[Cs,hs] = contour(MMP.x,MMP.z,MMP.s,24:32,'m','linewidth',2);
clabel(Cs,hs,'fontsize',14,'color','m','labelspacing',1000)
xlabel('Along Track Distance [km]')
ylabel('Depth [m]')
title('MMP Dissipation [W kg^{-1}] Sep. 10, 2015')

print([driveName figMMP 'dissipation_MMPrepeat1.png'],'-dpng')


%% Get scaling for figure size so axis scale is consistent between plots
% bowchain units
y_bc = log10(max(spec_bc(:)))- log10(min(spec_bc(:)));
x_bc = log10(k_bc(end))-log10(k_bc(1));
% mmp temp units
y_mmp = log10(max(spec_SAT))- log10(min(spec_SAT));
x_mmp = log10(k_mmp(end))-log10(kr(1));
% mmp temp units
y_Smmp = log10(max(spec_S(:)))- log10(min(spec_S(:)));

% Set bc figure unit
left = 0.1;
titleHeight = 0.1;
funit = 0.6;
funit_y_mmp = funit*y_mmp/y_bc;
funit_x_mmp = funit*x_mmp/x_bc;
funit_y_Smmp = funit*y_Smmp/y_bc;


%% Load SST Spectra
% Mean
load('MUR_L4/CanadaBasin/Data/Spectra/MUR_mean_spectra.mat','spec','time')
spec_SATbar = spec(time == datenum(2017,9,10),:);
% 2015
load('MUR_L4/CanadaBasin/Data/Spectra/MUR_2015_spectra.mat','spec','kr',...
  'time')
spec_SAT = spec(time == datenum(2015,9,10),:);


%% Plot bowchain temperature line spectra
figure('position',[0 0 600 600])
zz = unique([-1*BC.z(1:4:end); MMP.z(61:32:end)]);
axes('position',[left 1-titleHeight-funit funit funit])
cb = paruly(numel(zz));
for i = 1:numel(zz)
  i_bc = find(BC.z == -1*zz(i));
  if ~isempty(i_bc)
    loglog(k_bc,spec_bc(i_bc,:),'color',cb(i,:))
    hold on
  end
end
colormap(cb)
xlabel('Wavenumber [1/km]')
title({'Bow Chain Temperature'; 'Wavenumber Spectra'})
c = colorbar;
c.Ticks = (1:numel(zz))./numel(zz);
c.TickLabels = num2str(zz);
set(c,'YDir','reverse' )
xlim([min(k_bc) max(k_bc)])
ylim([min(spec_bc(:)) max(spec_bc(:))])
loglog(k_bc,1e-3*k_bc.^-1,'--k','linewidth',3)
text(k_bc(100),5e-3*k_bc(100).^-1,'k^{-1}','fontsize',20)
loglog(k_bc,5e-4*k_bc.^-(5/3),'--k','linewidth',3)
text(k_bc(1000),1e-3*k_bc(1000).^-(5/3),'k^{-5/3}','fontsize',20)

print([driveName figBC 'Spectra/temperature_MMPrepeat1.png'],'-dpng')


%% Plot MMP temperature line spectra
figure('position',[0 0 600 600])
zz = unique([-1*BC.z(1:4:end); MMP.z(61:32:end)]);
kk = unique([kr k_mmp]);
axes('position',[left 1-titleHeight-funit_y_mmp funit_x_mmp funit_y_mmp])
cb = paruly(numel(zz));
loglog(kr,spec_SATbar,'m')
hold on
loglog(kr,spec_SAT,'r')
for i = 1:numel(zz)
  i_mmp = find(MMP.z == zz(i));
  if ~isempty(i_mmp)
    loglog(k_mmp,spec_mmp(i_mmp,:),'color',cb(i,:))
  end
end
colormap(cb)
xlabel('Wavenumber [1/km]')
title({'MMP Temperature'; 'Wavenumber Spectra'})
c = colorbar;
c.Ticks = (1:numel(zz))./numel(zz);
c.TickLabels = num2str(zz);
set(c,'YDir','reverse' )
xlim([min(kk) max(kk)])
ylim([min(spec_SAT) max(spec_SAT)])
loglog(kk,1e-4*kk.^-2,'--k','linewidth',3)
text(kk(2),1e-5*kk(2).^-2,'k^{-2}','fontsize',20)
legend('MUR SST Mean Sep. 10','MUR SST Sep. 10, 2015','location',...
  'southwest')

print([driveName figMMP 'Spectra/temperature_MMPrepeat1.png'],'-dpng')


%% Plot MMP salinity line spectra
figure('position',[0 0 600 600])
zz = unique([-1*BC.z(1:4:end); MMP.z(61:32:end)]);
axes('position',[left*2 1-titleHeight-funit_y_Smmp funit_x_mmp ...
  funit_y_Smmp])
cb = paruly(numel(zz));
for i = 1:numel(zz)
  i_mmp = find(MMP.z == zz(i));
  if ~isempty(i_mmp)
    loglog(k_mmp,spec_S(i_mmp,:),'color',cb(i,:))
    hold on
  end
end
colormap(cb)
xlabel('Wavenumber [1/km]')
title({'MMP Salinity'; 'Wavenumber Spectra'})
c = colorbar;
c.Ticks = (1:numel(zz))./numel(zz);
c.TickLabels = num2str(zz);
set(c,'YDir','reverse' )
xlim([min(k_mmp) max(k_mmp)])
ylim([min(spec_S(:)) max(spec_S(:))])
loglog(k_mmp,1e-2*k_mmp.^-1,'--k','linewidth',3)
text(k_mmp(2),1.5e-2*k_mmp(2).^-1,'k^{-1}','fontsize',20)
loglog(k_mmp,1e-4*k_mmp.^-2,'--k','linewidth',3)
text(k_mmp(2),2e-5*k_mmp(2).^-2,'k^{-2}','fontsize',20)

print([driveName figMMP 'Spectra/salinity_MMPrepeat1.png'],'-dpng')


%% Plot MMP density line spectra
figure('position',[0 0 600 600])
cb = paruly(numel(MMP.z));
for i = 1:20:numel(MMP.z)
  loglog(k_mmp,spec_SGTH(i,:),'color',cb(i,:))
  hold on
end
colormap(cb)
xlabel('Wavenumber [1/km]')
title('MMP Density Wavenumber Spectra')
c = colorbar;
c.Ticks = (1:20:numel(MMP.z))./numel(MMP.z);
c.TickLabels = num2str(MMP.z(1:20:numel(MMP.z)));
set(c,'YDir','reverse' )
xlim([min(k_mmp) max(k_mmp)])
ylim([min(spec_SGTH(:)) max(spec_SGTH(:))])
loglog(k_mmp,1e-1*k_mmp.^-1,'--k','linewidth',3)
text(k_mmp(10),5e-1*k_mmp(10).^-1,'k^{-1}','fontsize',20)

print([driveName figMMP 'Spectra/density_MMPrepeat1.png'],'-dpng')


%% Load data
load('MUR_L4/CanadaBasin/Data/MUR_20150909_11.mat')

%% Plot Sep 10, 2015
figure('position',[0 0 800 500])
m_proj('lambert','long',[-170 -135],'lat',[69 76])
m_pcolor(lon,lat,squeeze(data.SST(:,:,2))-273.15)
shading flat
colormap(paruly(25))
colorbar('location','eastoutside')
hold on
m_coast('patch',[.7 .7 .7]);
title('Sep. 10, 2015 MUR SST [^\circ C]')
m_grid('box','fancy','tickdir','in');

print('MUR_L4/CanadaBasin/Figures/SST_20150910.png','-dpng')



