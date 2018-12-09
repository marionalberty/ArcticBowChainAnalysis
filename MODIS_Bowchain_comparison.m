% MODIS_Bowchain_comparison

% Get statistical comparison of how well the MODIS swath data agrees with
% the bowchain surface observations.
% Potential metrics:
% - mean and std of diff using closest available pass in time
% - diff as a function of time

clear all
close all
clc

set(0,'defaultaxesfontsize',18,'defaultaxeslinewidth',0.7,...
  'defaultlinelinewidth',0.7,'defaultpatchlinewidth',0.7)

%% Set cases, file lists and parameters

% Bowchain paths and files
path_bct = 'Arctic/ArcticMix/bowchain/gridded_data/';
flist_bct = dirs(fullfile(path_bct,'*.mat'));

% Initialize stat variables
dx = [];
dt = [];
dT = [];
cnum = [];

% Load Bathy
mlat = ncread('etopo1_CanadaBasin.nc','lat');
mlon = ncread('etopo1_CanadaBasin.nc','lon');
mbat = ncread('etopo1_CanadaBasin.nc','Band1');

% Satellite sensor
% A = AVHRR
% M = MODIS
% R = AMSR2
%
sensor_key = [];

% Set case
for case_bct = 1:length(flist_bct)
  
  fname_bct = flist_bct(case_bct).name;
  
  % Satellite paths and files
  path_sst = ['Arctic/ArcticMix/satellite/SST/' fname_bct(1:end-4) '/'];
  flist_sst = dirs(fullfile(path_sst,'*.nc'));
  
  
  %% Bowchain data
  
  % Load
  load([path_bct fname_bct],'grid')
  
  % Extract and simplify
  dnum_bct = grid.dnum(1,:)';
  lat_bct = grid.lat(1,:)';
  lon_bct = grid.lon(1,:)';
  tem_bct = [];
  for iiT = 1:length(grid.t(1,:))
    tem_bct(iiT,1) = interp1(grid.p(:,iiT),grid.t(:,iiT),0,'nearest',...
      'extrap');
  end
  clear grid
  
  % Get coordinates that encompass bowchain data
  minlat = floor(min(lat_bct));
  maxlat = ceil(max(lat_bct));
  minlon = floor(min(lon_bct));
  maxlon = ceil(max(lon_bct));
  
  %% Extract gridpoints of sst data that bowchain passed through
  %  Make vector of the sst value, time, lat extent, lon extent
  
  % Initialize satelitte SST vectors
  sec81_sat = [];
  lat_sat = [];
  lon_sat = [];
  tem_sat = [];
  
  % Index through each sst file
  for i = 1:length(flist_sst)
    % Load variables
    [~,SAT] = read_nc_file_struct([path_sst flist_sst(i).name]);
    
    % Initialize check to see if file has no data
    npts_pre = numel(tem_sat);
    
    % Extract points in bowchain box and append
    % seconds since 1981
    sec81_sat = [sec81_sat; SAT.sst_dtime(SAT.lat >= minlat & ...
      SAT.lat <= maxlat & SAT.lon >= minlon & SAT.lon <= maxlon & ...
      SAT.quality_level > 3) + double(SAT.time)];
    % latitude
    lat_sat = [lat_sat; SAT.lat(SAT.lat >= minlat & SAT.lat <= maxlat & ...
      SAT.lon >= minlon & SAT.lon <= maxlon & SAT.quality_level > 3)];
    % longitude
    lon_sat = [lon_sat; SAT.lon(SAT.lat >= minlat & SAT.lat <= maxlat & ...
      SAT.lon >= minlon & SAT.lon <= maxlon & SAT.quality_level > 3)];
    % temperature [degC]
    tem_sat = [tem_sat; SAT.sea_surface_temperature(SAT.lat >= minlat & ...
      SAT.lat <= maxlat & SAT.lon >= minlon & SAT.lon <= maxlon & ...
      SAT.quality_level > 3) - 273.15];
    
    % Number of points after file
    npts_post = numel(tem_sat);
    
    if npts_pre == npts_post
      delete([path_sst flist_sst(i).name])
    end
    
  end
  
  % Convert sec since 1981 to datenum
  dnum_sat = datenum(1981,1,1) + sec81_sat./86400;
  
  
  %% Run through bowchain and make list of best sst points to use
  
  % Initialize closest index
  i_sat4bct = NaN(size(lat_bct));
  dist_sat2bct = NaN(size(lat_bct));
  
  for i = 1:length(tem_bct)
    % Calculate distance between bowchain point and sat points
    dist = pos2dist(lat_bct(i),lon_bct(i),lat_sat,lon_sat);
    % Get the min
    [dist_sat2bct(i),i_sat4bct(i)] = min(dist);
  end
  
  % Get diffs
  diff_sstbct = tem_bct - tem_sat(i_sat4bct);
  diff_dnum = dnum_bct - dnum_sat(i_sat4bct);
  
  % Set aside case data
  dx = [dx; dist_sat2bct];
  dt = [dt; diff_dnum];
  dT = [dT; diff_sstbct];
  cnum = [cnum; case_bct*ones(size(diff_dnum))];
  
  %% Make some plots
  subplot(2,3,case_bct)
  contour(mlon,mlat,mbat',round(min(mbat(:)),-3):500:0,'color',...
    [0.7 0.7 0.7])
  hold on
  scatter(lon_sat,lat_sat,30,tem_sat,'fill')
  scatter(lon_bct,lat_bct,10,tem_bct,'fill')
  caxis([min([tem_sat; tem_bct]) max([tem_sat; tem_bct])])
  axis([minlon-1 maxlon+1 minlat-1 maxlat+1])
  colorbar
  title(fname_bct(1:end-4))
end

%% Calculate statistics of differences
% whole deployment, single sat grid points,...

% Mean diff of temp, dist, and time
