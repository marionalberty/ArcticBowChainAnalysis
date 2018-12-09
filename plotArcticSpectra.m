function plotArcticSpectra(Cruise,Bcase)
% Loads bow chain temperature and chi for an individual deployment,
% calculates the horizontal wavenumber spectra for each depth and

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',1,...
  'defaultlinelinewidth',1,'defaultpatchlinewidth',1,...
  'defaultFigureColor','white')


%% Set path and params
% Personal path
driveName = '/Users/marionsofiaalberty/MATLAB/Arctic/';
% Inputs/Outputs directories:
pathin = [driveName Cruise '/bowchain/Data/AllDeployObs/'];
figOut = [driveName Cruise '/bowchain/Figures/MiniSummaries/' Bcase '/'];
dataout = [driveName Cruise '/bowchain/Data/Spectra/'];


%% Load deployment data
load([pathin Bcase '.mat'],'bchain')


%% Grid temperature to be uniform in space for spectra
% Get distances bewteen points
xx = m_lldist(bchain.ungridded.lon(1,:),bchain.ungridded.lat(1,:));
% Get distance along track
x = [0 cumsum(xx)']*1000;
[x,ix,~] = unique(x);
% Chose dx
dx = round(mean(diff(x)),4);
% Set up gridded structure
bchain.xgrid.x = 0:dx:max(x);
bchain.xgrid.z = bchain.z;
bchain.xgrid.dn = interp1(x,bchain.ungridded.dn(ix),bchain.xgrid.x);
bchain.xgrid.lat = interp1(bchain.ungridded.dn(ix),...
  bchain.ungridded.lat(ix),bchain.xgrid.dn);
bchain.xgrid.lon = interp1(bchain.ungridded.dn(ix),...
  bchain.ungridded.lon(ix),bchain.xgrid.dn);
if isfield(bchain,'MLD')
  bchain.xgrid.MLD = interp1(bchain.dn,bchain.MLD,bchain.xgrid.dn,...
    'nearest','extrap');
end
% Repmat x and add in positional offsets
x = bchain.ungridded.x(:,ix)+repmat(x,numel(bchain.ungridded.pos),1);
% z and t for interpolant input
z = -bchain.ungridded.z(:,ix);
t = bchain.ungridded.t(:,ix);
% Get meshgrids
[xq,zq] = meshgrid(bchain.xgrid.x,bchain.xgrid.z);
% Interpolate onto regular grid
F = scatteredInterpolant(x(:),z(:),t(:));
bchain.xgrid.t = F(xq,zq);
% Convert distance back to km
bchain.xgrid.x = bchain.xgrid.x/1000;
% Change dx units to km
dx = dx/1000;
% Calculate dTdz for xgrid
bchain.xgrid.dTdz = nan*bchain.xgrid.t;
for i = 1:numel(bchain.xgrid.dn)
  bchain.xgrid.dTdz(:,i) = gradient(bchain.xgrid.t(:,i),-bchain.xgrid.z);
end
% Calculate dx(dTdz) for xgrid
bchain.xgrid.dxdTdz = nan*bchain.xgrid.t;
for i = 1:numel(bchain.xgrid.z)
  bchain.xgrid.dxdTdz(i,:) = gradient(bchain.xgrid.dTdz(i,:),...
    bchain.xgrid.x);
end
% Housekeeping
clear F xq zq x ix xx z t


%% Calculate Spectra for bow chain temperature
% Initialize Spectra Structure
bchain.spectra.z = bchain.xgrid.z;
for i = 1:numel(bchain.xgrid.z)
  % Calc spectra for temperature
  [bchain.spectra.Pow_T(i,:),bchain.spectra.k] = ...
    spectrum(bchain.xgrid.t(i,:),bchain.xgrid.x);
  % Calc spectra for dTdz
  [bchain.spectra.Pow_Tz(i,:),~] = ...
    spectrum(bchain.xgrid.dxdTdz(i,:),bchain.xgrid.x);
end
if isfield(bchain,'MLD')
  % Calc spectra along MLD
  % Temperature
  [bchain.spectra.Pow_T_mld,~] = ...
    spectrum(bchain.xgrid.t(bchain.xgrid.z == bchain.xgrid.MLD),...
    bchain.xgrid.x);
  % dTdz
  [bchain.spectra.Pow_Tz_mld,~] = ...
    spectrum(bchain.xgrid.dxdTdz(bchain.xgrid.z == bchain.xgrid.MLD),...
    bchain.xgrid.x);
end


%% Calc Spectra for filtered and grided CTD observations
% Get indecies for CTDs observations
i_S = find(sum(~isnan(bchain.ungridded.s),2)/numel(bchain.ungridded.dn) ...
  > 0.9);
% Initialize depth
bchain.spectra.z_TS = round(mean(bchain.ungridded.p(i_S,:),2));
bchain.TSfiltered.z = round(mean(bchain.ungridded.p(i_S,:),2));
% Get 10 second window
win = hanning(round(10/mean(diff(bchain.ungridded.dn)*86400)));
% Get distances bewteen points [km]
xx = m_lldist(bchain.ungridded.lon(1,:),bchain.ungridded.lat(1,:));
% Get distance along track [m]
x = [0 cumsum(xx)']*1000;
[x,ix,~] = unique(x);
% Chose dx [km]
dx = (round(mean(diff(x)),4)*numel(win))/1000;
% Repmat x and add in positional offsets [km]
x = (bchain.ungridded.x(i_S,ix)+repmat(x,numel(i_S),1))/1000;
% Make x grid [km]
bchain.TSfiltered.x = 0:dx:max(x(end,:));
% Save dn, lat, lon
bchain.TSfiltered.dn = interp1(x(1,:),bchain.ungridded.dn(ix),...
  bchain.TSfiltered.x);
bchain.TSfiltered.lat = interp1(x(1,:),bchain.ungridded.lat(i_S(1),ix),...
  bchain.TSfiltered.x);
bchain.TSfiltered.lon = interp1(x(1,:),bchain.ungridded.lon(i_S(1),ix),...
  bchain.TSfiltered.x);
for i = 1:numel(bchain.spectra.z_TS)
  % Extract data to filter
  s_temp = bchain.ungridded.s(i_S(i),ix);
  t_temp = bchain.ungridded.t(i_S(i),ix);
  p_temp = bchain.ungridded.p(i_S(i),ix);
  % Check for nans
  iuse = ~isnan(s_temp);
  % Filter data
  s_temp = filt_ends(win,s_temp(iuse));
  t_temp = filt_ends(win,t_temp(iuse));
  p_temp = filt_ends(win,p_temp(iuse));
  % Interp onto regular xgrid
  bchain.TSfiltered.s(i,:) = interp1(x(i,iuse),s_temp,bchain.TSfiltered.x);
  bchain.TSfiltered.t(i,:) = interp1(x(i,iuse),t_temp,bchain.TSfiltered.x);
  bchain.TSfiltered.p(i,:) = interp1(x(i,iuse),p_temp,bchain.TSfiltered.x);
  % Calc density
  bchain.TSfiltered.sig0(i,:) = sw_pden(bchain.TSfiltered.s(i,:),...
    bchain.TSfiltered.t(i,:),bchain.TSfiltered.p(i,:),0);
  % Calc spectra
  [bchain.spectra.PowTS_T(i,:),bchain.spectra.k_TS] = ...
    spectrum(bchain.TSfiltered.t(i,:),bchain.TSfiltered.x);
  [bchain.spectra.PowTS_S(i,:),~] = spectrum(bchain.TSfiltered.s(i,:),...
    bchain.TSfiltered.x);
  [bchain.spectra.PowTS_Sigma(i,:),~] = spectrum(...
    bchain.TSfiltered.sig0(i,:),bchain.TSfiltered.x);
end
% Housekeeping
clear i_S win xx x ix dx *_temp


%% Plot all temperature spectra, colored by depth
cmap = paruly(numel(bchain.z));
figure
for i = numel(bchain.z):-1:1
  loglog(bchain.spectra.k,bchain.spectra.Pow_T(i,:),'color',cmap(i,:))
  hold on
end
colormap(cmap)
c = colorbar;
c.Ticks = (1:numel(bchain.z))./numel(bchain.z);
c.TickLabels = num2str(bchain.z);
set(c,'YDir','reverse' )
axis tight
xlabel('Wavenumber [1/km]')
title({'Bow Chain Temperature'; 'Wavenumber Spectra'})
if isfield(bchain,'MLD')
  % Wavenumber equal to mean mixed layer depth
  k_mld = 1000./mean(bchain.MLD);
  bchain.spectra.k_mld = k_mld;
  plot([k_mld k_mld],[min(bchain.spectra.Pow_T(:)) ...
    max(bchain.spectra.Pow_T(:))],':r')
  text(k_mld*1.2,min(bchain.spectra.Pow_T(:))*10,'k_{MLD}','fontsize',...
    20,'color','r')
  % Wavenumber of surface mixed layer Rossby wave
  k_mros = 1000*sw_f(mean(bchain.lat))/...
    (mean(bchain.N2(bchain.z < bchain.MLD))*mean(bchain.MLD)*2*pi);
  bchain.spectra.k_mros = k_mros;
  plot([k_mros k_mros],[min(bchain.spectra.Pow_T(:)) ...
    max(bchain.spectra.Pow_T(:))],':r')
  text(k_mros*1.2,min(bchain.spectra.Pow_T(:))*10,'k_{Rossby}',...
    'fontsize',20,'color','r')
  % Wavenumber scalings
  kL = [bchain.spectra.k(1) k_mros];
  kM = [k_mros bchain.spectra.k(end)];
  kH = [k_mld bchain.spectra.k(end)];
else
  % Wavenumber scalings
  kL = [bchain.spectra.k(1) 10];
  kM = [10 bchain.spectra.k(end)];
  kH = [100 bchain.spectra.k(end)];
end
% Wavenumber relationships
% -2
loglog(kL,1e-5*kL.^-2,':k','linewidth',3)
text(kL(1)*1.2,1e-4*kL(1).^-2,'k^{-2}','fontsize',20)
% -3
loglog(kL,1e-6*kL.^-3,':k','linewidth',3)
text(kL(end)*0.4,1e-7*kL(end).^-3,'k^{-3}','fontsize',20)
% -1
loglog(kM,1e-6*kM.^-1,'k','linewidth',3)
text(kM(1)*1.2,1e-9*kM(1).^-1,'k^{-1}','fontsize',20)
% -(5/3)
loglog(kH,1e-2*kH.^-(5/3),'k','linewidth',3)
text(kH(1)*1.2,kH(1).^-(5/3),'k^{-5/3}','fontsize',20)

% Print figure
print([figOut Bcase '_bchainTspectra'],'-dpng')
savefig([figOut Bcase '_bchainTspectra.fig'])

close


%% Plot all dTdz Spectra, colored by depth
cmap = paruly(numel(bchain.z));
figure
for i = numel(bchain.z):-1:1
  loglog(bchain.spectra.k,bchain.spectra.Pow_Tz(i,:),'color',cmap(i,:))
  hold on
end
colormap(cmap)
c = colorbar;
c.Ticks = (1:numel(bchain.z))./numel(bchain.z);
c.TickLabels = num2str(bchain.z);
set(c,'YDir','reverse' )
axis tight
xlabel('Wavenumber [1/km]')
title({'Bow Chain d/dx(dT/dz)'; 'Wavenumber Spectra'})
kk = [bchain.spectra.k(1) bchain.spectra.k(end)];
% Wavenumber relationships
% +1
loglog(kk,1e-2*kk,'k','linewidth',3)
text(kk(1)*1.2,5e-2*kk(1),'k','fontsize',20)
% +3/2
loglog(kk,1e-5*kk.^(3/2),'k','linewidth',3)
text(kk(1)*1.2,1e-6*kk(1).^(3/2),'k^{3/2}','fontsize',20)

% Print figure
print([figOut Bcase '_bchaindxdTdzspectra'],'-dpng')
savefig([figOut Bcase '_bchaindxdTdzspectra.fig'])

close


%% Plot all Sigma and Temp Spectra, colored by depth
cmap = winter(numel(bchain.spectra.z_TS));
figure
for i = numel(bchain.spectra.z_TS):-1:1
  loglog(bchain.spectra.k_TS,bchain.spectra.PowTS_T(i,:),':',...
    'color',cmap(i,:))
  hold on
  loglog(bchain.spectra.k_TS,bchain.spectra.PowTS_Sigma(i,:),...
    'color',cmap(i,:))
end
colormap(cmap)
c = colorbar;
c.Ticks = (1:numel(bchain.spectra.z_TS))./numel(bchain.spectra.z_TS);
c.TickLabels = num2str(bchain.spectra.z_TS);
set(c,'YDir','reverse' )
axis tight
xlabel('Wavenumber [1/km]')
title({'Bow Chain T(:) and \sigma (-)'; 'Wavenumber Spectra'})
% Add wavenumber limits
if isfield(bchain,'MLD')
  % Wavenumber of surface mixed layer Rossby wave
  plot([k_mros k_mros],[min(bchain.spectra.PowTS_T(:)) ...
    max(bchain.spectra.PowTS_T(:))],':r')
  text(k_mros*1.2,min(bchain.spectra.PowTS_T(:))*10,'k_{Rossby}',...
    'fontsize',20,'color','r')
  % Wavenumber scalings
  kL = [bchain.spectra.k_TS(1) k_mros];
  kM = [k_mros bchain.spectra.k_TS(end)];
else
  % Wavenumber scalings
  kL = [bchain.spectra.k_TS(1) 10];
  kM = [10 bchain.spectra.k_TS(end)];
end
% Wavenumber relationships
% -2
loglog(kL,1e-4*kL.^-2,':k','linewidth',3)
text(kL(1)*1.2,5e-6*kL(1).^-2,'k^{-2}','fontsize',20)
% -3
loglog(kL,1e-3*kL.^-3,':k','linewidth',3)
text(kL(1)*5,1e-4*kL(1).^-3,'k^{-3}','fontsize',20)
% -1
loglog(kM,1e-4*kM.^-1,'k','linewidth',3)
text(kM(1)*1.2,1e-2*kM(1).^-1,'k^{-1}','fontsize',20)

% Print figure
print([figOut Bcase '_bchainTSigmaspectra'],'-dpng')
savefig([figOut Bcase '_bchainTSigmaspectra.fig'])

close


%% Calculate Structure Function for bow chain temperature
% % Initialize Structure Function Structure
% bchain.strfunc.z = bchain.z;
% for i = 1:numel(bchain.z)
%   tic
%   % Calc spectra for segment
%   [bchain.strfunc.D(i,:),bchain.strfunc.r,bchain.strfunc.error(i,:)] = ...
%     structureFunction_SST(bchain.t(i,:),bchain.lat,bchain.lon,2*dx,100);
%   toc
% end
%
%
% %% Plot all structure functions, colored by depth
% cmap = paruly(numel(bchain.z));
% figure
% for i = numel(bchain.z):-1:1
%   loglog(1./bchain.strfunc.r,bchain.strfunc.D(i,:),'color',cmap(i,:))
%   hold on
% end
% colormap(cmap)
% c = colorbar;
% c.Ticks = (1:numel(bchain.z))./numel(bchain.z);
% c.TickLabels = num2str(bchain.z);
% set(c,'YDir','reverse' )
% axis tight
% xlabel('Wavenumber [1/km]')
% title({'Bow Chain Temperature'; 'Structure Function'})
% % loglog(bchain.spectra.k,1e-4*bchain.spectra.k.^-2,':k','linewidth',3)
% % text(bchain.spectra.k(5),1e-6*bchain.spectra.k(5).^-2,'k^{-2}',...
% %   'fontsize',20)
% % loglog(bchain.spectra.k,1e-3*bchain.spectra.k.^-3,'--k','linewidth',3)
% % text(300,1e-5*300.^-3,'k^{-3}','fontsize',20)
%
% % Print figure
% print([figOut Bcase '_bchainTstructfunc'],'-dpng')
% savefig([figOut Bcase '_bchainTstructfunc.fig'])
%
% close


%% Save data
% Make file name
fname_out = [dataout Bcase '.mat'];
% Save data
save(fname_out,'bchain')
end