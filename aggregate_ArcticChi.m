% aggregate_ArcticChi.m

% Aggregate all bowchain chi and mmp chi data into vectors containing the
% mmp chi, the mean bowchain chi within a minute of the mmp profile at the
% same depth, the depth of the observations, mmp time of observation, lat,
% lon, t, s, eps,

clear all; clc


set(0,'defaultaxesfontsize',16,'defaultaxeslinewidth',0.7,...
  'defaultlinelinewidth',1,'defaultpatchlinewidth',0.7,...
  'defaultFigureColor','white')


%% Set default paths
% Personnal paths:
driveName = '/Users/marionsofiaalberty/MATLAB/Arctic/';
% Inputs directories:
pathinA = [driveName 'ArcticMix/bowchain/Data/AllDeployObs/'];
flistA = dir(fullfile(pathinA,'*.mat'));
pathinS = [driveName 'SODA/bowchain/Data/AllDeployObs/'];
flistS = dir(fullfile(pathinS,'*.mat'));
flist = [flistA; flistS];
% Output directory:
pathout = [driveName 'Data/'];
% Figure out path:
figOut = [driveName 'Figures/PaperFigs/'];

% Housekeeping
clear pathin* flistA flistS


%% Initialize output structure
data.z = [];            data.dn = [];
data.lat = [];          data.lon = [];
data.t = [];            data.s = [];
data.eps_mmp = [];      data.chi_mmp = [];
data.chi_bar = [];      data.t_bar = [];


%% Set parameters
% Time range for averaging
dt = 1/(24*60);


for i = 1:numel(flist)
  %% Load the file with all the deployment observations
  % Clear previous observations
  clear bchain mmp
  % Load current deployment
  load([flist(i).folder '/' flist(i).name],'bchain','mmp')
  % Check if deployment had mmp
  if ~exist('mmp','var')
    continue
  end
  
  
  %% Reduce mmp record to only the depths that overlap with the bowchain
  % Find mmp indecies for common depths
  [~,i_mmpz,~] = intersect(mmp.z,bchain.chi.z);
  % Get all the fields that need to be reduced
  fname = fieldnames(mmp);
  % Reduce mmp size
  for ifield = 1: numel(fname)
    [rows,~] = size(mmp.(fname{ifield}));
    if rows > 1
      mmp.(fname{ifield}) = mmp.(fname{ifield})(i_mmpz,:);
    end
  end
  
  %% Nan out spurious chi
  bchain.chi.chi(bchain.chi.chi > 1e-3) = nan;
  bchain.chi.chi(bchain.chi.chi < 1e-11) = nan;
  mmp.chi(mmp.chi > 1e-3) = nan;
  if i == 8 || i > 9
    mmp.chi = mmp.chi*nan;
  end
  
  
  %% Get mean bow chain chi for each mmp profile, add to data structure
  % Get mean chi time
  time_chi = mean(bchain.chi.time);
  for i_time = 1:numel(mmp.dnum)
    % Indecies for averaging bow chain observations
    i_chi = find(time_chi >= mmp.dnum(i_time) - dt & ...
      time_chi <= mmp.dnum(i_time) + dt);
    i_temp = find(bchain.dn >= mmp.dnum(i_time) - dt & ...
      bchain.dn <= mmp.dnum(i_time) + dt);
    % Only use data above the MLD-3m
    iz = mmp.z < nanmean(bchain.MLD(i_temp) - 3);
    
    % Add vectors to data structures
    data.z = [data.z; mmp.z(iz)];
    data.dn = [data.dn; mmp.dnum(i_time)*ones(size(mmp.z(iz)))];
    data.lat = [data.lat; mmp.lat(i_time)*ones(size(mmp.z(iz)))];
    data.lon = [data.lon; mmp.lon(i_time)*ones(size(mmp.z(iz)))];
    data.t = [data.t; mmp.t(iz,i_time)];
    data.s = [data.s; mmp.s(iz,i_time)];
    data.eps_mmp = [data.eps_mmp; mmp.eps(iz,i_time)];
    data.chi_mmp = [data.chi_mmp; mmp.chi(iz,i_time)];
    % Calculate means
    data.chi_bar = [data.chi_bar; nanmean(bchain.chi.chi(iz,i_chi),2)];
    data.t_bar = [data.t_bar; nanmean(bchain.t(iz,i_temp),2)];
  end
end


%% Remove mising data
% Get indecies where there is neither mmp eps or chi
i_rm = (isnan(data.eps_mmp) & isnan(data.chi_mmp)) | isnan(data.chi_bar);
% Get all the fields that need to be reduced
fname = fieldnames(data);
% Remove data points
for ifield = 1: numel(fname)
  data.(fname{ifield})(i_rm) = [];
end

%% Save scattered data
save([pathout 'arcticChiEps.mat'],'data')


%% Plot data

figure
scatter(data.chi_bar,data.chi_mmp)
p = polyfit(log10(data.chi_bar(~isnan(data.chi_mmp))),...
  log10(data.chi_mmp(~isnan(data.chi_mmp))),1);
x = linspace(-11,-2,100);
y = polyval(p,x);
hold on
plot(10.^x,10.^y)
plot(10.^x,10.^x,'k')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel('Bow Chain')
ylabel('MMP')
axis equal

figure
scatter(data.chi_bar(data.z == 5),data.chi_mmp(data.z == 5))
p = polyfit(log10(data.chi_bar(~isnan(data.chi_mmp) & data.z == 5)),...
  log10(data.chi_mmp(~isnan(data.chi_mmp) & data.z == 5)),1);
x = linspace(-10,-3,100);
y = polyval(p,x);
hold on
plot(10.^x,10.^y)
plot(10.^x,10.^x,'k')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel('Bow Chain')
ylabel('MMP')
axis equal

