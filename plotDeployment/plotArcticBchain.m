function plotArcticBchain(Cruise,Bcase)
% Loads the gridded bow chian data for the given deployment case and then
% finds and loads all the available data from the concurrent SWIMS, MMP,
% SADCP, and met observations. All available data is then plotted and saved
% in 'AllDeployObs' folder.
%
% Inputs:
% Bcase - String that matches the filename for the desired bowchain
%         deployment. e.g. 'Test' or 'MMPrepeat1'


set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',1,...
  'defaultlinelinewidth',1,'defaultpatchlinewidth',1,...
  'defaultFigureColor','white')

%% Load Bowchain data
switch Cruise
  case 'ArcticMix'
    path.dir_base = 'Arctic/ArcticMix/';
  case 'SODA'
    path.dir_base = 'Arctic/SODA/';
end
% Bowchain data
path.dir_bchain = [path.dir_base 'bowchain/Data/Gridded/' Bcase '.mat'];
% Load and rename data variable
load(path.dir_bchain,'data')
bchain = data;
bchain.z = -bchain.z;
clear data
% Get deployment dnum range
dnum_range = [min(bchain.dn) max(bchain.dn)];


%% Load SADCP data
switch Cruise
  case 'ArcticMix'
    path.dir_sadcp = [path.dir_base 'wellADCP/Data/Processed/'];
    % get list of files
    flist = dir(fullfile(path.dir_sadcp,'*minute.mat'));
    for i = 1:numel(flist)
      % Load file
      load([path.dir_sadcp flist(i).name],'C')
      % Check if it overlaps with bow chain deployment
      if C.dnum(1) <= dnum_range(1) && C.dnum(end) >= dnum_range(2)
        break
      end
    end
    sadcp = C;
    clear C flist use_sec
  case 'SODA'
    path.dir_sadcp = [path.dir_base 'wellADCP/sadcp.mat'];
    load(path.dir_sadcp,'sadcp')
end

%% Load Surface fluxes
switch Cruise
  case 'ArcticMix'
    path.dir_metFLX = [path.dir_base ...
      'met_data_JM/AMX15_COARE_bulkFluxes.mat'];
    % Load data
    load(path.dir_metFLX,'COARE')
    % Extract data simultaneous with deployment
    i_flux = COARE.dtnum >= dnum_range(1) & COARE.dtnum <= dnum_range(2);
    % Variables to keep
    flds = {'tau','hnet','dtnum'};
    for i = 1:numel(flds)
      flux.(flds{i}) = COARE.(flds{i})(i_flux);
    end
    % Housekeeping
    clear COARE i_flux
  case 'SODA'
    path.dir_metFLX = [path.dir_base 'underway/fluxes.mat'];
    % Load data
    load(path.dir_metFLX,'fluxes')
    % Extract data simultaneous with deployment
    i_flux = fluxes.dnum >= dnum_range(1) & fluxes.dnum <= dnum_range(2);
    % Variables to keep
    flux.dtnum = fluxes.dnum(i_flux);
    flux.tau = fluxes.derived.stress_P(i_flux);
    flux.hnet = fluxes.derived.LH(i_flux) + fluxes.derived.SH(i_flux) + ...
      fluxes.derived.LW_net(i_flux) + fluxes.derived.SW_net(i_flux);
    % Housekeeping
    clear fluxes i_flux
end


%% Load MMP data
switch Cruise
  case 'ArcticMix'
    path.dir_MMP = [path.dir_base 'MMP/Data/MMPgridMLEpos.mat'];
  case 'SODA'
    path.dir_MMP = [path.dir_base 'MMP/MMPgrid.mat'];
end
load(path.dir_MMP,'MMP')
% Extract data simultaneous with deployment
i_mmp = MMP.dnum >= dnum_range(1) & MMP.dnum <= dnum_range(2);
% Variables to keep
flds = {'t','s','eps','chi','dnum','lat','lon'};
% Make mmp structure
mmp.z = MMP.p';
mmp.p = MMP.p';
for i = 1:numel(flds)
  mmp.(flds{i}) = real(MMP.(flds{i})(:,i_mmp));
end
% Find and remove bad data/profiles
mmp.eps(mmp.eps <= 1e-10) = nan;
mmp.chi(mmp.chi <= 1e-11) = nan;
mmp.chi(mmp.chi > 1) = nan;
% Remove bad profiles
i_rm = sum(isnan(mmp.eps))./numel(mmp.p) == 1;
for i = 1:numel(flds)
  mmp.(flds{i})(:,i_rm) = [];
end

% Housekeeping
clear MMP i_mmp flds i_rm


%% Fix time offset between mmp and bowchain
if ~isempty(mmp.eps)
  % Load lookup table for depth
  load([path.dir_base 'bowchain/Data/lagcov_lookup_mmp.mat'],'table')
  % set depths to use
  iz = table.z(strcmp(table.deployment,Bcase));
  % Get temperature
  T1 = bchain.t(bchain.z == iz,:);
  T2 = mmp.t(mmp.z == iz,:);
  % Find outliers
  i_out2 = isoutlier(T2);
  % Interp onto common time
  dt = 1/86400;
  t_q = max([bchain.dn(1) mmp.dnum(1)]):dt:min([bchain.dn(end) ...
    mmp.dnum(end)]);
  T1 = interp1(bchain.dn,T1,t_q);
  T2 = interp1(mmp.dnum(~i_out2),T2(~i_out2),t_q);
  % Get indicies for calc
  i_str = find(~isnan(T1.*T2) == 1,1,'first');
  i_end = find(~isnan(T1.*T2) == 1,1,'last');
  % Remove nans
  T1 = naninterp(T1);
  T2 = naninterp(T2);
  % Get time offset
  offset = lagcov_t_offset(t_q(i_str:i_end),T1(i_str:i_end),...
    T2(i_str:i_end));
  % Pause to let user check offset calculations
%   disp('Check MMP offset, then press any key')
%   pause;
  close all
  
  % Apply offset to mmp
  mmp.dnum = mmp.dnum + offset;
  % Interp bchain lat and lon onto mmp
  mmp.lon = interp1(bchain.dn,bchain.lon,mmp.dnum,'linear','extrap');
  mmp.lat = interp1(bchain.dn,bchain.lat,mmp.dnum,'linear','extrap');
  % Housekeeping
  clear offset T1 T2 i_end i_str t_q dt i_out2 z_off table
end


%% Clean up mmp epsilon and chi, then calculate heat fluxes
if ~isempty(mmp.eps)
  switch Cruise
    case 'ArcticMix'
      % Remove eps and chi above 10 m
      mmp.eps(mmp.z <= 10,:) = nan;
      mmp.chi(mmp.z <= 10,:) = nan;
    case 'SODA'
      %       % Remove eps and chi above 2 m
      %       mmp.eps(mmp.z <= 4,:) = nan;
      %       mmp.chi(mmp.z <= 4,:) = nan;
  end
  % Initialize mask
  i_mask = ones(size(mmp.eps));
  % Make mask for each profile
  for i = 1:numel(mmp.dnum)
    % Find all the local minima
    [TF,PP] = islocalmin(log10(mmp.eps(:,i)));
    if sum(TF) < 2
      continue
    end
    switch Cruise
      case 'ArcticMix'
        % Get the CDF of the local minima's prominence
        h = histogram(PP(TF),20,'normalization','probability');
        [cdf,i_cdf,~] = unique(cumsum(h.Values));
        p_cdf = h.BinEdges(1:end-1) + ((h.BinEdges(2)-h.BinEdges(2))/2);
        p_cdf = p_cdf(i_cdf);
        % Find the first prominent minimum in the top 20 m that is more
        % prominent than 66% of the located minima
        P_thres = interp1(cdf,p_cdf,0.66);
        i_PP = find(PP(mmp.z <= 20) > P_thres,1,'first');
      case 'SODA'
        i_PP = find(TF,1,'first');
    end
    % Adjust the nan mask
    if ~isempty(i_PP)
      i_mask(1:i_PP-1,i) = nan;
    end
  end
  % Apply mask
  mmp.eps = mmp.eps.*i_mask;
  mmp.chi = mmp.chi.*i_mask;
  % Housekeeping
  clear i_mask PP i_PP P_thres p_cdf cdf i_cdf h TF
  close all
  
  % Calculate density rho_0
  mmp.sgth = sw_pden(mmp.s,mmp.t,mmp.p,0);
  % Calculate heat fluxes
  % Prelocate dTdz
  mmp.dTdz = nan*mmp.t;
  % Calculate dTdz for each profile
  for i = 1:numel(mmp.dnum)
    % Get indecies of non-nans
    i_good = ~isnan(mmp.t(:,i));
    % skip if not enough data
    if sum(i_good) < 5
      continue
    end
    % dTdz from central differences
    mmp.dTdz(i_good,i) = cent_diff_n(mmp.t(i_good,i),mmp.z(2)-mmp.z(1),9)';
  end
  % Calc thermal diffusivity
  mmp.kappa_T = mmp.chi./(2*mmp.dTdz.^2);
  % Calc heat capacity
  Cp = sw_cp(mmp.s,mmp.t,mmp.p);
  mmp.J_q = -mmp.sgth.*Cp.*mmp.kappa_T.*mmp.dTdz;
  % Housekeeping
  clear Cp i_good
end


%% Load FastCTD
switch Cruise
  case 'SODA'
    path.dir_fctd = [path.dir_base 'fctd/'];
    % Load file
    load([path.dir_fctd 'FCTD_total.mat'],'F')
    % Variables to keep
    flds = {'time','temperature','salinity','lat','lon'};
    % Find data within deployment range
    i_fctd = F.time >= dnum_range(1) & F.time <= dnum_range(2);
    % Put data in structure
    for i = 1:numel(flds)
      if ~isempty(i_fctd)
        fctd.(flds{i}) = F.(flds{i})(:,i_fctd);
      else
        fctd.(flds{i}) = [];
      end
    end
    % Add z and p
    fctd.z = F.depth;
    % Quality control S
    fctd.salinity(fctd.salinity < 24) = nan;
    % Calculate density
    fctd.sgth = sw_pden(fctd.salinity,fctd.temperature,fctd.z,0);
    % Housekeeping
    clear F flds i_fctd i
    
    
    %% Fix time offset between FCTD and bowchain
    if ~isempty(fctd.salinity)
      % Load lookup table for depth
      load([path.dir_base 'bowchain/Data/lagcov_lookup_fctd.mat'],'table')
      % set depths to use
      iz = table.z(strcmp(table.deployment,Bcase));
      % Get temperature
      T1 = bchain.t(bchain.z == iz,:);
      T2 = fctd.temperature(fctd.z == iz,:);
      % Find outliers
      i_out2 = isoutlier(T2);
      % Interp onto common time
      dt = 1/86400;
      t_q = max([bchain.dn(1) fctd.time(1)]):dt:min([bchain.dn(end) ...
        fctd.time(end)]);
      T1 = interp1(bchain.dn,T1,t_q);
      T2 = interp1(fctd.time(~i_out2),T2(~i_out2),t_q);
      % Get indicies for calc
      % Need to
      i_str = find(~isnan(T2.*T1) == 1,1,'first');
      i_end = find(~isnan(T2.*T1) == 1,1,'last');
      % Remove nans
      T1 = naninterp(T1);
      T2 = naninterp(T2);
      % Get time offset
      offset = lagcov_t_offset(t_q(i_str:i_end),T1(i_str:i_end),...
        T2(i_str:i_end));
      % Pause to let user check offset calculations
%       disp('Check FastCTD offset, then press any key')
%       pause;
      close all
      % Apply offset to fctd
      fctd.time = fctd.time + offset;
      % Interp bchain lat and lon onto mmp
      fctd.lon = interp1(bchain.dn,bchain.lon,fctd.time,'linear','extrap');
      fctd.lat = interp1(bchain.dn,bchain.lat,fctd.time,'linear','extrap');
      % Housekeeping
      clear offset T1 T2 i_end i_str t_q dt i_out2 z_off table
    end
  case 'ArcticMix'
    % Variables to keep
    flds = {'time','temperature','salinity','lat','lon'};
    % Initialize FCTD
    for i = 1:numel(flds)
      fctd.(flds{i}) = [];
    end
    fctd.time = [];
end


%% Load Epsi CTD
switch Cruise
  case 'SODA'
    path.dir_fctd = [path.dir_base 'epsi/'];
    % Load file
    load([path.dir_fctd 'Gridded.mat'],'epsi')
    % Variables to keep
    flds = {'time','Temperature','Salinity'};
    % Find data within deployment range
    i_epsi = epsi.time >= dnum_range(1) & epsi.time <= dnum_range(2);
    % Put data in structure
    for i = 1:numel(flds)
      epsi.(flds{i})(:,~i_epsi) = [];
    end
    % Transpose pressure
    epsi.z = epsi.Pressure';
    % Quality control S
    epsi.Salinity(epsi.Salinity < 24) = nan;
    % Calculate density
    epsi.sgth = sw_pden(epsi.Salinity,epsi.Temperature,epsi.z,0);
    % Housekeeping
    clear flds i_epsi i
    
    
    %% Fix time offset between Epsi CTD and bowchain
    if ~isempty(epsi.Salinity)
      % Load lookup table for depth
      load([path.dir_base 'bowchain/Data/lagcov_lookup_epsi.mat'],'table')
      % set depths to use
      iz = table.z(strcmp(table.deployment,Bcase));
      % Get temperature
      T1 = bchain.t(bchain.z == iz,:);
      T2 = epsi.Temperature(epsi.z == iz,:);
      % Find outliers
      i_out2 = isoutlier(T2);
      % Interp onto common time
      dt = 1/86400;
      t_q = max([bchain.dn(1) epsi.time(1)]):dt:min([bchain.dn(end) ...
        epsi.time(end)]);
      T1 = interp1(bchain.dn,T1,t_q);
      T2 = interp1(epsi.time(~i_out2),T2(~i_out2),t_q);
      % Get indicies for calc
      % Need to
      i_str = find(~isnan(T2.*T1) == 1,1,'first');
      i_end = find(~isnan(T2.*T1) == 1,1,'last');
      % Remove nans
      T1 = naninterp(T1);
      T2 = naninterp(T2);
      % Get time offset
      offset = lagcov_t_offset(t_q(i_str:i_end),T1(i_str:i_end),...
        T2(i_str:i_end));
      % Pause to let user check offset calculations
%       disp('Check FastCTD offset, then press any key')
%       pause;
      close all
      % Apply offset to fctd
      epsi.time = epsi.time + offset;
      % Interp bchain lat and lon onto mmp
      epsi.lon = interp1(bchain.dn,bchain.lon,epsi.time,'linear','extrap');
      epsi.lat = interp1(bchain.dn,bchain.lat,epsi.time,'linear','extrap');
      % Housekeeping
      clear offset T1 T2 i_end i_str t_q dt i_out2 z_off table
    end
  case 'ArcticMix'
    % Variables to keep
    flds = {'time','Temperature','Salinity'};
    % Initialize Epsi CTD
    for i = 1:numel(flds)
      epsi.(flds{i}) = [];
    end
    epsi.time = [];
end


%% Load SWIMS
switch Cruise
  case 'ArcticMix'
    path.dir_SWIMS = [path.dir_base 'SWIMS/Data/'];
    % get list of files
    flist = dir(fullfile(path.dir_SWIMS,'SWIMSgrid-*.mat'));
    % Variables to keep
    flds = {'t1','t2','s1','s2','lat','lon'};
    % Initialize swims
    for i = 1:numel(flds)
      swims.(flds{i}) = [];
    end
    swims.dn = [];
    for i = 1:numel(flist)
      % Load file
      load([path.dir_SWIMS flist(i).name],'SWIMSgrid')
      % Check if data is within deployment range
      dn = datenum(2015,1,1) + SWIMSgrid.yday;
      i_swims = dn >= dnum_range(1) & dn <= dnum_range(2);
      if ~isempty(i_swims)
        swims.dn = [swims.dn dn(i_swims)];
        for k = 1:numel(flds)
          swims.(flds{k}) = [swims.(flds{k}) ...
            real(SWIMSgrid.(flds{k})(:,i_swims))];
        end
      end
    end
    % Add z and p
    swims.z = SWIMSgrid.z';
    swims.p = SWIMSgrid.p';
    % Quality control S
    swims.s1(swims.s1 < 24) = nan;
    swims.s2(swims.s2 < 24) = nan;
    % Calculate density
    swims.sgth1 = sw_pden(swims.s1,swims.t1,swims.p,0);
    swims.sgth2 = sw_pden(swims.s2,swims.t2,swims.p,0);
    % Housekeeping
    clear SWIMSgrid flist flds dn i_swims k i
    
    
    %% Fix time offset between swims and bowchain
    if ~isempty(swims.s1)
      % Load lookup table for depth
      load([path.dir_base 'bowchain/Data/lagcov_lookup_swims.mat'],'table')
      % set depths to use
      iz = table.z(strcmp(table.deployment,Bcase));
      % Get temperature
      T1 = bchain.t(bchain.z == iz,:);
      T2 = swims.t1(swims.z == iz,:);
      % Find outliers
      i_out2 = isoutlier(T2);
      % Interp onto common time
      dt = 1/86400;
      t_q = max([bchain.dn(1) swims.dn(1)]):dt:min([bchain.dn(end) ...
        swims.dn(end)]);
      T1 = interp1(bchain.dn,T1,t_q);
      T2 = interp1(swims.dn(~i_out2),T2(~i_out2),t_q);
      % Get indicies for calc
      % Need to
      i_str = find(~isnan(T2.*T1) == 1,1,'first');
      i_end = find(~isnan(T2.*T1) == 1,1,'last');
      % Remove nans
      T1 = naninterp(T1);
      T2 = naninterp(T2);
      % Get time offset
      offset = lagcov_t_offset(t_q(i_str:i_end),T1(i_str:i_end),...
        T2(i_str:i_end));
      % Pause to let user check offset calculations
%       disp('Check SWIMS offset, then press any key')
%       pause;
      close all
      % Apply offset to swims
      swims.dn = swims.dn + offset;
      % Interp bchain lat and lon onto mmp
      swims.lon = interp1(bchain.dn,bchain.lon,swims.dn,'linear','extrap');
      swims.lat = interp1(bchain.dn,bchain.lat,swims.dn,'linear','extrap');
      % Housekeeping
      clear offset T1 T2 i_end i_str t_q dt i_out2 z_off table
    end
    
    
    % If SWIMS is out of the water for 1 hour, add column of nans to improve
    % plotting
    dt = diff(swims.dn)*24;
    i_an = find(dt > 1);
    if ~isempty(i_an)
      flds = {'t1','t2','s1','s2','sgth1','sgth2'};
      for k = 1:numel(flds)
        swims.(flds{k}) = [swims.(flds{k})(:,1:i_an) nan(size(swims.z)) ...
          swims.(flds{k})(:,i_an+1:end)];
      end
      flds = {'lat','lon','dn'};
      for k = 1:numel(flds)
        swims.(flds{k}) = [swims.(flds{k})(1:i_an) ...
          mean(swims.(flds{k})(i_an:i_an+1)) swims.(flds{k})(i_an+1:end)];
      end
    end
    % Housekeeping
    clear dt i_an flds k
  case 'SODA'
    % Variables to keep
    flds = {'t1','t2','s1','s2','lat','lon'};
    % Initialize swims
    for i = 1:numel(flds)
      swims.(flds{i}) = [];
    end
    swims.dn = [];
end


%% Calculate chi from bowchain data
% Check that data is 1 Hz
if round(mean(diff(bchain.dn))*86400) == 1
  % Prelocate dT/dz
  bchain.dTdz = nan*bchain.t;
  % Average T in time to reduce noise in dTdz
  wi = hanning(300);
  T_4dz = bchain.t;
  for i = 1:numel(bchain.z)
    T_4dz(i,:) = filt_ends(wi,T_4dz(i,:));
  end
  % Calc dTdz for each profile
  for i = 1:numel(bchain.dn)
    % Get indecies of non-nans
    i_good = ~isnan(T_4dz(:,i));
    % skip if not enough data
    if sum(i_good) < 3
      continue
    end
    % Calculate dT/dz for profile
    bchain.dTdz(i_good,i) = gradient(T_4dz(i_good,i),...
      -1*bchain.z(i_good));
  end
  % Housekeeping
  clear i i_good wi T_4dz
  
  % Calculate the vertical velocity of the sensors
  % Get indecies where there is pressure
  i_w = find(sum(~isnan(bchain.ungridded.p),2)./...
    numel(bchain.ungridded.dn) > 0.5);
  % Prelocate w
  bchain.ungridded.w = nan*bchain.ungridded.t;
  w = nan*bchain.t;
  % Calculate w for each sensor
  for i = 1:numel(i_w)
    bchain.ungridded.w(i_w(i),:) = -1*...
      gradient(bchain.ungridded.p(i_w(i),:),bchain.ungridded.dn*86400);
  end
  % Interpolate onto gridded time
  w_p = interp1(bchain.ungridded.dn,bchain.ungridded.w(i_w,:)',bchain.dn)';
  z_wp = interp1(bchain.ungridded.dn,bchain.ungridded.z(i_w,:)',...
    bchain.dn)';
  for i = 1:numel(bchain.dn)
    % Interpolate onto bchain grid
    w(:,i) = interp1(-1*z_wp(:,i),w_p(:,i),bchain.z);
    if sum(~isnan(w(:,i))) > 1
      w(:,i) = interp1(bchain.z(~isnan(w(:,i))),w(~isnan(w(:,i)),i),...
        bchain.z,'nearest','extrap');
    end
  end
  % Housekeeping
  clear i_w i w_p z_wp
  
  % Calculate the ship speed
  % Get distances bewteen points
  dx = sign(diff(bchain.lon)).*m_lldist(bchain.lon,...
    mean(bchain.lat)*ones(size(bchain.lat)))';
  dy = sign(diff(bchain.lat)).*m_lldist(mean(bchain.lon)*...
    ones(size(bchain.lon)),bchain.lat)';
  % Get distance along track
  x = [0 cumsum(dx)];
  y = [0 cumsum(dy)];
  % Calc U = dx/dt [m/s]
  U = gradient(x*1000,bchain.dn*86400);
  V = gradient(y*1000,bchain.dn*86400);
  % Get mean ship speed
  bchain.Uship_bar = mean(abs(U + 1j*V));
  % Get SADCP U_20 and V_20
  U_s = interp1(sadcp.dnum,nanmean(sadcp.u(sadcp.z <= 25,:)),bchain.dn);
  V_s = interp1(sadcp.dnum,nanmean(sadcp.v(sadcp.z <= 25,:)),bchain.dn);
  % Replace nans with zeros
  U_s(isnan(U_s)) = 0;
  V_s(isnan(V_s)) = 0;
  % Combine
  U = complex(U+U_s,V+V_s);
  % Housekeeping
  clear x y dx dy U_s V_s V s_*
  
  % Estimate N2
  % Make z grid
  z_sgth = [0:50]';
  % Prelocate density matricies
  % mmp
  sgth_mmp = [];
  t_mmp = [];
  s_mmp = [];
  dnum_mmp = [];
  % fctd
  sgth_fctd = [];
  t_fctd = [];
  s_fctd = [];
  dnum_fctd = [];
  % epsi
  sgth_epsi = [];
  t_epsi = [];
  s_epsi = [];
  dnum_epsi = [];
  % swims
  sgth_swims = [];
  t_swims = [];
  s_swims = [];
  dnum_swims = [];
  % Interp MMP sgth onto grid, if available
  if ~isempty(mmp.t)
    sgth_mmp = interp1(mmp.z,mmp.sgth,z_sgth);
    t_mmp = interp1(mmp.z,mmp.t,z_sgth);
    s_mmp = interp1(mmp.z,mmp.s,z_sgth);
    dnum_mmp = mmp.dnum;
  end
  % Interp FastCTD sgth onto grid, if available
  if ~isempty(fctd.temperature)
    sgth_fctd = interp1(fctd.z,fctd.sgth,z_sgth);
    t_fctd = interp1(fctd.z,fctd.temperature,z_sgth);
    s_fctd = interp1(fctd.z,fctd.salinity,z_sgth);
    dnum_fctd = fctd.time;
  end
  % Interp Epsi sgth onto grid, if available
  if ~isempty(epsi.Temperature)
    sgth_epsi = interp1(epsi.z,epsi.sgth,z_sgth);
    t_epsi = interp1(epsi.z,epsi.Temperature,z_sgth);
    s_epsi = interp1(epsi.z,epsi.Salinity,z_sgth);
    dnum_epsi = epsi.time;
  end
  % Interp SWIMS sgth onto grid, if available
  if ~isempty(swims.t1)
    sgth_swims = interp1(swims.z,swims.sgth1,z_sgth);
    t_swims = interp1(swims.z,swims.t1,z_sgth);
    s_swims = interp1(swims.z,swims.s1,z_sgth);
    dnum_swims = swims.dn;
  end
  % Concatinate sgth observations
  sgth = [sgth_mmp sgth_fctd sgth_epsi sgth_swims];
  t_sgth = [t_mmp t_fctd t_epsi t_swims];
  s_sgth = [s_mmp s_fctd s_epsi s_swims];
  dnum_sgth = [dnum_mmp dnum_fctd dnum_epsi dnum_swims];
  % Housekeeping
  clear sgth_* t_mmp t_fctd t_epsi t_swims dnum_mmp dnum_fctd dnum_epsi ...
    dnum_swims s_mmp s_fctd s_epsi s_swims
  
  % Calc chi if sgth is not empty
  if ~isempty(sgth)
    % Sort by time
    [dnum_sgth,i_dn] = sort(dnum_sgth);
    sgth = sgth(:,i_dn);
    t_sgth = t_sgth(:,i_dn);
    s_sgth = s_sgth(:,i_dn);
    % Find profiles that don't have observations in the top 8 m
    i_bad = isnan(sgth(z_sgth == 8,:));
    % Remove profiles without top 8 m
    sgth(:,i_bad) = [];
    t_sgth(:,i_bad) = [];
    s_sgth(:,i_bad) = [];
    dnum_sgth(i_bad) = [];
    % Housekeeping
    clear i_dn i_bad
    
    % Sort to be statically stable and find MLD relationship
    % Mixed layer delta sgth from Spicy Seas paper (MacKinnon et al. 2016)
    MLS = 0.5;  % [kg/m^{-3}]
    N2_fill = 1e-5;
    dz = z_sgth(2)-z_sgth(1);
    % Prelocate mld variables
    MLD = dnum_sgth*nan;
    dT_ML = dnum_sgth*nan;
    sgth_ML = dnum_sgth*nan;
    t_ML = dnum_sgth*nan;
    s_ML = dnum_sgth*nan;
    N2_sgth = sgth(2:end,:)*nan;
    z_N2 = z_sgth(2:end) - (dz/2);
    for i = 1:numel(dnum_sgth)
      % Find non nans in profile
      i_good = ~isnan(sgth(:,i));
      % Get temporary temp
      t_temp = t_sgth(i_good,i);
      s_temp = s_sgth(i_good,i);
      % Sort to be stable
      [sgth(i_good,i),i_srtd] = sort(sgth(i_good,i));
      % Sort temp to match sorted sgth
      t_sgth(i_good,i) = t_temp(i_srtd);
      s_sgth(i_good,i) = s_temp(i_srtd);
      % Find min sgth
      sgth_min = min(sgth(:,i));
      % Get index of MLD
      iBase = find(sgth(:,i) > sgth_min + MLS,1,'first');
      if isempty(iBase)
        % Use bottom if no mixed layer found
        iBase = numel(z_sgth);
      end
      % Find top of the water column
      i_top = find(i_good,1,'first');
      % Save MLD
      MLD(i) = z_sgth(iBase);
      % Find the equivalent temperature gradient
      dT_ML(i) = t_sgth(i_top,i) - t_sgth(iBase,i);
      % Save the mean mixed layer properties
      sgth_ML(i) = nanmean(sgth(1:iBase,i));
      t_ML(i) = nanmean(t_sgth(1:iBase,i));
      s_ML(i) = nanmean(s_sgth(1:iBase,i));
      % Calculate N2 for each profile
      N2_sgth(:,i) = 9.8.*diff(sgth(:,i))./(dz*nanmean(sgth(:,i)));
      if z_sgth(i_top) < MLD(i)
        if N2_sgth(i_top,i) > 1e-3
          % If the topmost N2 is in the mixed layer
          N2_sgth(1:i_top-1,i) = N2_fill;
        else
          % Slab interpolate the profile to the surface
          N2_sgth(1:i_top,i) = N2_sgth(i_top,i);
        end
      end
    end
    % Get the depth of each N2 profile relative to the mixed layer
    % Negative is above the mixed layer, positive below
    z_rel2ML = z_N2*ones(1,numel(dnum_sgth)) - ones(numel(z_N2),1)*MLD;
    % Average N2 as a function of distance above/below mixed layer
    % Get z to query over
    z_q = min(z_rel2ML(:)):max(z_rel2ML(:));
    % Prelocate N2 bar
    N2_bar = nan*z_q;
    for i = 1:numel(z_q)
      N2_bar(i) = nanmean(N2_sgth(z_rel2ML == z_q(i)));
    end
    % Calc dT_bar
    dT_bar = nanmean(dT_ML)*2;
    % Shift z_q
    z_q = z_q + 0.5;
    % Interp mean ML sgth onto bchain time
    bchain.rho_ML = interp1(dnum_sgth,sgth_ML,bchain.dn,'nearest',...
      'extrap');
    bchain.t_ML = interp1(dnum_sgth,t_ML,bchain.dn,'nearest',...
      'extrap');
    bchain.s_ML = interp1(dnum_sgth,s_ML,bchain.dn,'nearest',...
      'extrap');
    % Housekeeping
    clear MLS N2_fill dz MLD dT_ML N2_sgth z_N2 i i_good temp i_srtd ...
      sgth t_sgth sgth_min iBase i_top z_sgth z_rel2ML dnum_sgth ...
      sgth_ML s_sgth t_sgth s_sgth
    
    % Find mixed layer depth from bchain record using observed dT
    bchain.N2 = nan*bchain.t;
    i_z = find(bchain.z > 7,1,'first')-1;
    bchain.MLD = nan*bchain.dn;
    for i = 1:numel(bchain.dn)
      % Get topmost T
      i_top = find(~isnan(bchain.t(:,i)),1,'first');
      % Get index of MLD
      [~,iBase] = min(abs((bchain.t(i_top,i) - bchain.t(bchain.z > 7,i))...
        - dT_bar));
      bchain.MLD(i) = bchain.z(iBase+i_z);
      % Interp N2 model onto data
      bchain.N2(:,i) = interp1(z_q,N2_bar,bchain.z - bchain.z(iBase+i_z));
      % Extrap to top if there are near-surface nans
      i_top = find(~isnan(bchain.N2(:,i)),1,'first');
      bchain.N2(1:i_top,i) = bchain.N2(i_top,i);
    end
    % Housekeeping
    clear i_top iBase N2_bar z_q i_z
    
    % Prelocate Eps, Chi, and time_chi
    bchain.chi.eps = [];
    bchain.chi.chi = [];
    bchain.chi.time = [];
    bchain.chi.z = [];
    % Calculate chi for each depth
    for i = 1:numel(bchain.z)
      [bchain.chi.eps(i,:),bchain.chi.chi(i,:),U_chi,bchain.chi.time] = ...
        sbe56EpsiChi(bchain.t(i,:),filt_ends(hanning(10),abs(U))',...
        bchain.dTdz(i,:),bchain.N2(i,:),w(i,:),bchain.dn);
      bchain.chi.z(i,1) = bchain.z(i);
    end
    
    %% Calculate MLD, dTdz_ML, and heat flux for chi profiles
    % Initialize matrices
    bchain.chi.MLD = nan*bchain.chi.chi(1,:);
    bchain.chi.T = nan*bchain.chi.chi;
    bchain.chi.N2 = nan*bchain.chi.chi;
    bchain.chi.dTdz = nan*bchain.chi.chi;
    bchain.chi.Jq = nan*bchain.chi.chi;
    % Calculate Cp for each time [J kg^-1 C^-1]
    Cp = sw_cp(bchain.s_ML,bchain.t_ML,bchain.MLD);
    % Work by chi profile
    for i = 1:numel(bchain.chi.MLD)
      % Get indecies of bow chain profiles within chi time scale
      i_chi = bchain.dn >= bchain.chi.time(1,i) & ...
        bchain.dn <= bchain.chi.time(2,i);
      % Average MLD
      bchain.chi.MLD(i) = nanmean(bchain.MLD(i_chi));
      % Average T
      bchain.chi.T(:,i) = nanmean(bchain.t(:,i_chi),2);
      % Average N2
      bchain.chi.N2(:,i) = nanmean(bchain.N2(:,i_chi),2);
      % Calculate dTdz from mean T profiles [C/m]
      bchain.chi.dTdz(:,i) = gradient(bchain.chi.T(:,i),-1*bchain.chi.z);
      % Calculate heat flux [W/m^2]
      bchain.chi.Jq(:,i) = -mean(bchain.rho_ML(i_chi))*mean(Cp(i_chi))*...
        bchain.chi.chi(:,i)./(2*bchain.chi.dTdz(:,i));
    end
    
    
    % Plot chi and inputs for chi
    % Set up pannel sizing
    titleheight = 0.04;
    labelheight = 0.02;
    n_fig_units = 9;
    n_figs = 5;
    totalH = 1 - (titleheight*n_figs) - (labelheight*(n_figs+1));
    HperUnit = totalH/n_fig_units;
    left = 0.05;
    width = 0.88;
    % Make figure
    figure('position',[0 0 1100 1000])
    % Bow chain T
    ah(1) = axes('position',[left 1-titleheight-HperUnit*2 width ...
      HperUnit*2]);
    pcolor(bchain.dn,bchain.z,bchain.t)
    hold on
    plot(mean(bchain.chi.time),bchain.chi.MLD,'k')
    shading interp
    xlim(dnum_range)
    ylim([0 bchain.z(end)])
    axis ij
    datetick('x','keeplimits')
    ylabel('Depth [m]')
    colorbar(ah(1),'position',[left+width+0.01 1-titleheight-HperUnit*2 ...
      0.02 HperUnit*2])
    colormap(ah(1),parula)
    set(ah(1),'Layer','top','ytick',0:5:20)
    title(['Bowchain Temperature [^\circ C], ' strrep(Bcase,'_','.') ...
      ' on ' datestr(dnum_range(1),1) ' to ' datestr(dnum_range(2),1)])
    % U_c
    ah(2) = axes('position',[left 1-titleheight*2-HperUnit*3-...
      labelheight width HperUnit]);
    plot(mean(bchain.chi.time),U_chi)
    xlim(dnum_range)
    datetick('x','keeplimits')
    ylabel('Speed [m/s]')
    title('Speed of flow past sensor [m/s]')
    % N2
    ah(3) = axes('position',[left 1-titleheight*3-HperUnit*5-...
      labelheight*2 width HperUnit*2]);
    pcolor(bchain.dn,bchain.z,bchain.N2)
    shading interp
    xlim(dnum_range)
    ylim([0 bchain.z(end)])
    axis ij
    datetick('x','keeplimits')
    ylabel('Depth [m]')
    colorbar(ah(3),'position',[left+width+0.01 1-titleheight*3-...
      HperUnit*5-labelheight*2 0.02 HperUnit*2])
    colormap(ah(3),parula)
    set(ah(3),'Layer','top','ytick',0:5:20)
    title('N^2 [s^{-2}]')
    % w
    ah(4) = axes('position',[left 1-titleheight*4-HperUnit*7-...
      labelheight*3 width HperUnit*2]);
    pcolor(bchain.dn,bchain.z,w)
    shading interp
    xlim(dnum_range)
    ylim([0 bchain.z(end)])
    axis ij
    datetick('x','keeplimits')
    ylabel('Depth [m]')
    colorbar(ah(4),'position',[left+width+0.01 1-titleheight*4-...
      HperUnit*7-labelheight*3 0.02 HperUnit*2])
    colormap(ah(4),redblue(26))
    caxis([-3 3])
    set(ah(4),'Layer','top','ytick',0:5:20)
    title('Vertical Velocity of Sensors [m/s]')
    % chi
    ah(5) = axes('position',[left 1-titleheight*5-HperUnit*9-...
      labelheight*4 width HperUnit*2]);
    pcolor(mean(bchain.chi.time),bchain.chi.z,log10(bchain.chi.chi))
    shading interp
    xlim(dnum_range)
    ylim([0 bchain.z(end)])
    axis ij
    datetick('x','keeplimits')
    caxis([-10 max(log10(bchain.chi.chi(:)))])
    ylabel('Depth [m]')
    colorbar(ah(5),'position',[left+width+0.01 1-titleheight*5-...
      HperUnit*9-labelheight*4 0.02 HperUnit*2])
    colormap(ah(5),parula)
    set(ah(5),'Layer','top','ytick',0:5:20)
    title('\chi_T log_{10}([K^2/s])')
    
    % link axes
    linkaxes(ah,'x')
    
    % Print figure
    print([path.dir_base 'bowchain/Figures/MiniSummaries/' Bcase '/' ...
      Bcase '_chi'],'-dpng')
    savefig([path.dir_base 'bowchain/Figures/MiniSummaries/' Bcase '/' ...
      Bcase '_chi.fig'])
    
    close
    
    
    % Plot bow chain t, chi and heat flux
    % Set up pannel sizing
    titleheight = 0.05;
    labelheight = 0.03;
    n_figs = 3;
    totalH = 1 - (titleheight*n_figs) - (labelheight*(n_figs+1));
    HperUnit = totalH/n_figs;
    left = 0.05;
    width = 0.86;
    % Prep for Flux Figures
    % Pick countours
    J_levels = [-1000 -500 -100 -50 -10 -5 0 5 10 50 100 500 1000];
    % Prep a cleaned up flux
    J_4bin = bchain.chi.Jq;
    J_4bin(J_4bin < J_levels(1)) = J_levels(1);
    J_4bin(J_4bin > J_levels(end)) = J_levels(end);
    % Get binning
    [~,~,i_bin] = histcounts(J_4bin,J_levels);
    i_bin(i_bin == 0) = nan;
        
    % Make figure
    figure('position',[0 0 1100 600])
    % Bow chain T
    ah(1) = axes('position',[left 1-titleheight-HperUnit width ...
      HperUnit]);
    pcolor(bchain.dn,bchain.z,bchain.t)
    hold on
    plot(mean(bchain.chi.time),bchain.chi.MLD,'k')
    shading interp
    xlim(dnum_range)
    ylim([0 bchain.z(end)])
    axis ij
    datetick('x','keeplimits')
    ylabel('Depth [m]')
    colorbar(ah(1),'position',[left+width+0.01 1-titleheight-HperUnit ...
      0.02 HperUnit])
    colormap(ah(1),parula)
    set(ah(1),'Layer','top','ytick',0:5:20)
    set(ah(1),'XTickLabel','')
    title(['Bow Chain Temperature [^\circ C], ' strrep(Bcase,'_','.') ...
      ' on ' datestr(dnum_range(1),1) ' to ' ...
      datestr(dnum_range(2),1)])
    % Bow chain chi
    ah(2) = axes('position',[left 1-titleheight*2-HperUnit*2-...
      labelheight width HperUnit]);
    pcolor(mean(bchain.chi.time),bchain.chi.z,log10(bchain.chi.chi))
    hold on
    plot(mean(bchain.chi.time),bchain.chi.MLD,'k')
    shading interp
    xlim(dnum_range)
    ylim([0 bchain.z(end)])
    axis ij
    datetick('x','keeplimits')
    caxis([-10 max(log10(bchain.chi.chi(:)))])
    ylabel('Depth [m]')
    colorbar(ah(2),'position',[left+width+0.01 1-titleheight*2-...
      HperUnit*2-labelheight 0.02 HperUnit])
    colormap(ah(2),parula)
    set(ah(2),'Layer','top','ytick',0:5:20)
    set(ah(2),'XTickLabel','')
    title('Bow Chain \chi_T log_{10}([K^2/s])')
    % Bow chain heat flux
    ah(3) = axes('position',[left 1-titleheight*3-HperUnit*3-...
      labelheight*2 width HperUnit]);
    contourf(mean(bchain.chi.time),bchain.chi.z,i_bin,...
      0.5:(max(i_bin(:))+0.5),'LineStyle','none')
    hold on
    plot(mean(bchain.chi.time),bchain.chi.MLD,'k')
    xlim(dnum_range)
    ylim([0 bchain.z(end)])
    axis ij
    datetick('x','keeplimits')
    caxis([0.5 max(i_bin(:))+0.5])
    ylabel('Depth [m]')
    cbh = colorbar(ah(3),'position',[left+width+0.01 1-titleheight*3-...
      HperUnit*3-labelheight*2 0.02 HperUnit]);
    colormap(ah(3),redblue(numel(J_levels)-1))
    set(cbh,'XTick',0.5:(max(i_bin(:))+0.5))
    set(cbh,'XTickLabel',sprintfc('%d',J_levels'))
    set(ah(3),'Layer','top','ytick',0:5:20)
    title('Bow Chain J_q [W/m^2]')
    
    % link axes
    linkaxes(ah,'x')
    
    % Print figure
    print([path.dir_base 'bowchain/Figures/MiniSummaries/' Bcase '/' ...
      Bcase '_heatflux'],'-dpng')
    savefig([path.dir_base 'bowchain/Figures/MiniSummaries/' Bcase '/' ...
      Bcase '_heatflux.fig'])
    
    close
    
    
    if ~isempty(mmp.eps)
      % Plot bowchain t and chi comapred with MMP t, chi, and eps
      % Set up pannel sizing
      titleheight = 0.04;
      labelheight = 0.02;
      n_figs = 5;
      totalH = 1 - (titleheight*n_figs) - (labelheight*(n_figs+1));
      HperUnit = totalH/n_figs;
      left = 0.05;
      width = 0.88;
      % Get salinity range for contouring
      slim = [min([mmp.s(:); swims.s1(:)]) max([mmp.s(:); swims.s1(:)])];
      % Get temperature range for contouring
      tlim = [min(bchain.t(:)) max(bchain.t(:))];
      % Make figure
      figure('position',[0 0 1100 1000])
      % Bow chain T
      ah(1) = axes('position',[left 1-titleheight-HperUnit width ...
        HperUnit]);
      pcolor(bchain.dn,bchain.z,bchain.t)
      shading interp
      xlim(dnum_range)
      ylim([0 bchain.z(end)])
      axis ij
      datetick('x','keeplimits')
      ylabel('Depth [m]')
      colorbar(ah(1),'position',[left+width+0.01 1-titleheight-HperUnit ...
        0.02 HperUnit])
      colormap(ah(1),parula)
      set(ah(1),'Layer','top','ytick',0:5:20)
      freezeColors;
      hold on
      contour(mmp.dnum,mmp.z,real(mmp.s),linspace(slim(1),slim(2),20),'k')
      title(['Bow Chain Temperature [^\circ C], ' strrep(Bcase,'_','.') ...
        ' on ' datestr(dnum_range(1),1) ' to ' ...
        datestr(dnum_range(2),1)])
      % Bow chain chi
      ah(2) = axes('position',[left 1-titleheight*2-HperUnit*2-...
        labelheight width HperUnit]);
      pcolor(mean(bchain.chi.time),bchain.chi.z,log10(bchain.chi.chi))
      shading interp
      xlim(dnum_range)
      ylim([0 bchain.z(end)])
      axis ij
      datetick('x','keeplimits')
      caxis([-10 max(log10(bchain.chi.chi(:)))])
      ylabel('Depth [m]')
      colorbar(ah(2),'position',[left+width+0.01 1-titleheight*2-...
        HperUnit*2-labelheight 0.02 HperUnit])
      colormap(ah(2),parula)
      set(ah(2),'Layer','top','ytick',0:5:20)
      freezeColors;
      hold on
      contour(mmp.dnum,mmp.z,real(mmp.s),linspace(slim(1),slim(2),20),'k')
      title('Bow Chain \chi_T log_{10}([K^2/s])')
      % MMP T with S contours
      ah(3) = axes('position',[left 1-titleheight*3-HperUnit*3-...
        labelheight*2 width HperUnit]);
      pcolor(mmp.dnum,mmp.z,mmp.t)
      shading interp
      caxis(tlim)
      xlim(dnum_range)
      ylim([0 bchain.z(end)])
      axis ij
      datetick('x','keeplimits')
      ylabel('Depth [m]')
      colorbar(ah(3),'position',[left+width+0.01 1-titleheight*3- ...
        HperUnit*3-labelheight*2 0.02 HperUnit])
      colormap(ah(3),parula)
      set(ah(3),'Layer','top','ytick',0:5:20)
      freezeColors;
      hold on
      contour(mmp.dnum,mmp.z,mmp.s,linspace(slim(1),slim(2),20),'k')
      title('MMP Temperature [^\circ C]')
      % MMP chi with S contours
      ah(4) = axes('position',[left 1-titleheight*4-HperUnit*4- ...
        labelheight*3 width HperUnit]);
      pcolor(mmp.dnum,mmp.z,log10(mmp.chi))
      shading interp
      caxis([-10 max(log10(bchain.chi.chi(:)))])
      xlim(dnum_range)
      ylim([0 bchain.z(end)])
      axis ij
      datetick('x','keeplimits')
      ylabel('Depth [m]')
      colorbar(ah(4),'position',[left+width+0.01 1-titleheight*4-...
        HperUnit*4-labelheight*3 0.02 HperUnit])
      colormap(ah(4),parula)
      set(ah(4),'Layer','top','ytick',0:5:20)
      freezeColors;
      hold on
      contour(mmp.dnum,mmp.z,real(mmp.s),linspace(slim(1),slim(2),20),'k')
      title('MMP \chi_T log_{10}[K^2/s]')
      % MMP eps with S contours
      ah(5) = axes('position',[left 1-titleheight*5-HperUnit*5-...
        labelheight*4 width HperUnit]);
      pcolor(mmp.dnum,mmp.z,log10(mmp.eps))
      shading interp
      clim = [min(log10(mmp.eps(:))) max(log10(mmp.eps(:)))];
      caxis(clim)
      xlim(dnum_range)
      ylim([0 bchain.z(end)])
      axis ij
      datetick('x','keeplimits')
      ylabel('Depth [m]')
      colorbar(ah(5),'position',[left+width+0.01 1-titleheight*5-...
        HperUnit*5-labelheight*4 0.02 HperUnit])
      colormap(ah(5),parula)
      set(ah(5),'Layer','top','ytick',0:5:20)
      freezeColors;
      hold on
      contour(mmp.dnum,mmp.z,real(mmp.s),linspace(slim(1),slim(2),20),'k')
      title('MMP \epsilon log_{10}[W kg^{-1}]')
      
      % link axes
      linkaxes(ah,'x')
      
      % Print figure
      print([path.dir_base 'bowchain/Figures/MiniSummaries/' Bcase '/' ...
        Bcase '_bchainMMPchiComp'],'-dpng')
      savefig([path.dir_base 'bowchain/Figures/MiniSummaries/' Bcase ...
        '/' Bcase '_bchainMMPchiComp.fig'])
      
      close
    end
  end
end
% Housekeeping
clear U z_sgth t_sgth dnum_sgth sgth w


%% Plot data for deployment
% Set up pannel sizing
titleheight = 0.04;
labelheight = 0.02;
n_fig_units = 10;
n_figs = 6;
totalH = 1 - (titleheight*n_figs) - (labelheight*(n_figs+1));
HperUnit = totalH/n_fig_units;
left = 0.05;
width = 0.88;
% Get salinity range for contouring
slim = [min([mmp.s(:); swims.s1(:); fctd.salinity(:); epsi.Salinity(:)])...
  max([mmp.s(:); swims.s1(:); fctd.salinity(:); epsi.Salinity(:)])];
% Get temperature range for contouring
tlim = [min(bchain.t(:)) max(bchain.t(:))];
% Make shear^2 matrix
dz = sadcp.z(2)-sadcp.z(1);
shear = (diff(sadcp.u)./dz).^2 + (diff(sadcp.v)./dz).^2;

% Make figure
figure('position',[0 0 1100 1000])
% Bowchain
ah(1) = axes('position',[left 1-titleheight-HperUnit*2 width HperUnit*2]);
pcolor(bchain.dn,bchain.z,bchain.t)
shading interp
caxis(tlim)
xlim(dnum_range)
ylim([0 20])
axis ij
datetick('x','keeplimits')
ylabel('Depth [m]')
colorbar(ah(1),'position',[left+width+0.01 1-titleheight-HperUnit*2 ...
  0.02 HperUnit*2])
colormap(jet)
set(ah(1),'Layer','top','ytick',0:5:20)
title(['Bowchain Temperature [^\circ C], ' strrep(Bcase,'_','.') ' on '...
  datestr(dnum_range(1),1) ' to ' datestr(dnum_range(2),1)])
if ~isempty(mmp.eps)
  freezeColors;
  hold on
  contour(mmp.dnum,mmp.z,real(mmp.s),linspace(slim(1),slim(2),20),'k')
end
if ~isempty(swims.t1)
  freezeColors;
  hold on
  contour(swims.dn,swims.z,swims.s1,linspace(slim(1),slim(2),20),'k')
end
if ~isempty(fctd.temperature)
  freezeColors;
  hold on
  contour(fctd.time,fctd.z,fctd.salinity,linspace(slim(1),slim(2),20),'k')
end
if ~isempty(epsi.Temperature)
  freezeColors;
  hold on
  contour(epsi.time,epsi.z,epsi.Salinity,linspace(slim(1),slim(2),20),'k')
end

% Temp from SWIMS,FCTD, and MMP
if ~isempty(swims.t1) && ~isempty(mmp.t)
  % Both available
  ah(2) = axes('position',[left 1-titleheight*2-labelheight-HperUnit*4 ...
    width HperUnit*2]);
  pcolor(swims.dn,swims.z,swims.t1)
  hold on
  pcolor(mmp.dnum,mmp.z,mmp.t)
  shading interp
  caxis(tlim)
  xlim(dnum_range)
  ylim([0 50])
  axis ij
  datetick('x','keeplimits')
  ylabel('Depth [m]')
  colorbar(ah(2),'position',[left+width+0.01 1-titleheight*2-HperUnit*4-...
    labelheight 0.02 HperUnit*2])
  colormap(jet)
  set(ah(2),'Layer','top','ytick',0:10:50)
  freezeColors;
  contour(swims.dn,swims.z,swims.s1,linspace(slim(1),slim(2),20),'k')
  contour(mmp.dnum,mmp.z,mmp.s,linspace(slim(1),slim(2),20),'k')
  title('SWIMS and MMP Temperature [^\circ C]')
  
elseif ~isempty(swims.t1) && isempty(mmp.t)
  % Just SWIMS
  ah(2) = axes('position',[left 1-titleheight*2-labelheight-HperUnit*4 ...
    width HperUnit*2]);
  pcolor(swims.dn,swims.z,swims.t1)
  shading interp
  caxis(tlim)
  xlim(dnum_range)
  ylim([0 50])
  axis ij
  datetick('x','keeplimits')
  ylabel('Depth [m]')
  colorbar(ah(2),'position',[left+width+0.01 1-titleheight*2-HperUnit*4-...
    labelheight 0.02 HperUnit*2])
  colormap(jet)
  set(ah(2),'Layer','top','ytick',0:10:50)
  freezeColors;
  hold on
  contour(swims.dn,swims.z,swims.s1,linspace(slim(1),slim(2),20),'k')
  title('SWIMS Temperature [^\circ C]')
  
elseif ~isempty(epsi.Temperature)
  % Just Epsi
  ah(2) = axes('position',[left 1-titleheight*2-labelheight-HperUnit*4 ...
    width HperUnit*2]);
  pcolor(epsi.time,epsi.z,epsi.Temperature)
  shading interp
  caxis(tlim)
  xlim(dnum_range)
  ylim([0 50])
  axis ij
  datetick('x','keeplimits')
  ylabel('Depth [m]')
  colorbar(ah(2),'position',[left+width+0.01 1-titleheight*2-HperUnit*4-...
    labelheight 0.02 HperUnit*2])
  colormap(jet)
  set(ah(2),'Layer','top','ytick',0:10:50)
  freezeColors;
  hold on
  contour(epsi.time,epsi.z,epsi.Salinity,linspace(slim(1),slim(2),20),'k')
  title('Epsi Temperature [^\circ C]')
  
elseif ~isempty(fctd.temperature) && ~isempty(mmp.t)
  % FCTD and MMP available
  ah(2) = axes('position',[left 1-titleheight*2-labelheight-HperUnit*4 ...
    width HperUnit*2]);
  pcolor(fctd.time,fctd.z,fctd.temperature)
  hold on
  pcolor(mmp.dnum,mmp.z,mmp.t)
  shading interp
  caxis(tlim)
  xlim(dnum_range)
  ylim([0 50])
  axis ij
  datetick('x','keeplimits')
  ylabel('Depth [m]')
  colorbar(ah(2),'position',[left+width+0.01 1-titleheight*2-HperUnit*4-...
    labelheight 0.02 HperUnit*2])
  colormap(jet)
  set(ah(2),'Layer','top','ytick',0:10:50)
  freezeColors;
  contour(fctd.time,fctd.z,fctd.salinity,linspace(slim(1),slim(2),20),'k')
  title('FastCTD and MMP Temperature [^\circ C]')
  
elseif isempty(swims.t1) && isempty(fctd.temperature) && ~isempty(mmp.eps)
  % Just MMP
  ah(2) = axes('position',[left 1-titleheight*2-labelheight-HperUnit*4 ...
    width HperUnit*2]);
  pcolor(mmp.dnum,mmp.z,mmp.t)
  shading interp
  caxis(tlim)
  xlim(dnum_range)
  ylim([0 50])
  axis ij
  datetick('x','keeplimits')
  ylabel('Depth [m]')
  colorbar(ah(2),'position',[left+width+0.01 1-titleheight*2-HperUnit*4-...
    labelheight 0.02 HperUnit*2])
  colormap(jet)
  set(ah(2),'Layer','top','ytick',0:10:50)
  freezeColors;
  hold on
  contour(mmp.dnum,mmp.z,mmp.s,linspace(slim(1),slim(2),20),'k')
  title('MMP Temperature [^\circ C]')
  
elseif isempty(swims.t1) && ~isempty(fctd.temperature) && isempty(mmp.eps)
  % Just FCTD
  ah(2) = axes('position',[left 1-titleheight*2-labelheight-HperUnit*4 ...
    width HperUnit*2]);
  pcolor(fctd.time,fctd.z,fctd.temperature)
  shading interp
  caxis(tlim)
  xlim(dnum_range)
  ylim([0 50])
  axis ij
  datetick('x','keeplimits')
  ylabel('Depth [m]')
  colorbar(ah(2),'position',[left+width+0.01 1-titleheight*2-HperUnit*4-...
    labelheight 0.02 HperUnit*2])
  colormap(jet)
  set(ah(2),'Layer','top','ytick',0:10:50)
  freezeColors;
  hold on
  contour(fctd.time,fctd.z,fctd.salinity,linspace(slim(1),slim(2),20),'k')
  title('FastCTD Temperature [^\circ C]')
end

% MMP
if ~isempty(mmp.eps)
  ah(3) = axes('position',[left 1-titleheight*3-labelheight*2-HperUnit*6 ...
    width HperUnit*2]);
  pcolor(mmp.dnum,mmp.z,log10(mmp.eps))
  shading interp
  clim = [min(log10(mmp.eps(:))) max(log10(mmp.eps(:)))];
  caxis(clim)
  xlim(dnum_range)
  ylim([0 50])
  axis ij
  datetick('x','keeplimits')
  ylabel('Depth [m]')
  colorbar(ah(3),'position',[left+width+0.01 1-titleheight*3-HperUnit*6-...
    labelheight*2 0.02 HperUnit*2])
  colormap(jet)
  set(ah(3),'Layer','top','ytick',0:10:50)
  freezeColors;
  hold on
  contour(mmp.dnum,mmp.z,real(mmp.s),linspace(slim(1),slim(2),20),'k')
  title('MMP \epsilon log_{10}[W kg^{-1}]')
end

% SADCP
ah(4) = axes('position',[left 1-titleheight*4-labelheight*3-HperUnit*8 ...
  width HperUnit*2]);
pcolor(sadcp.dnum,sadcp.z(2:end)-(dz/2),shear)
shading interp
xlim(dnum_range)
caxis([0 0.001])
ylim([0 20])
axis ij
datetick('x','keeplimits')
ylabel('Depth [m]')
colorbar(ah(4),'position',[left+width+0.01 1-titleheight*4-HperUnit*8-...
  labelheight*3 0.02 HperUnit*2])
colormap(jet)
set(ah(4),'Layer','top','ytick',0:5:20)
title('Well ADCP Vertical Shear^2 [s^{-2}]')
if ~isempty(mmp.eps)
  freezeColors;
  hold on
  contour(mmp.dnum,mmp.z,real(mmp.s),linspace(slim(1),slim(2),20),'k')
end
if ~isempty(swims.t1)
  freezeColors;
  hold on
  contour(swims.dn,swims.z,swims.s1,linspace(slim(1),slim(2),20),'k')
end
if ~isempty(fctd.temperature)
  freezeColors;
  hold on
  contour(fctd.time,fctd.z,fctd.salinity,linspace(slim(1),slim(2),20),'k')
end
if ~isempty(epsi.Temperature)
  freezeColors;
  hold on
  contour(epsi.time,epsi.z,epsi.Salinity,linspace(slim(1),slim(2),20),'k')
end

% Wind stress
ah(5) = axes('position',[left 1-titleheight*5-labelheight*4-HperUnit*9 ...
  width HperUnit]);
plot(flux.dtnum,flux.tau)
grid on
grid minor
axis tight
ylim([0 round(max(flux.tau),3)])
xlim(dnum_range)
datetick('x','keeplimits')
title('Wind Stress [N/m]')


% Heat flux
ah(6) = axes('position',[left 1-titleheight*6-labelheight*5-HperUnit*10 ...
  width HperUnit]);
plot(flux.dtnum,flux.hnet)
grid on
grid minor
axis tight
xlim(dnum_range)
datetick('x','keeplimits')
ylim([min(flux.hnet) max(flux.hnet)])
title('Net Heat Flux [W m^{-2}]')

% link axes
linkaxes(ah,'x')

%% Print figure
print([path.dir_base 'bowchain/Figures/MiniSummaries/' Bcase '/' Bcase ...
  '_allData'],'-dpng')
savefig([path.dir_base 'bowchain/Figures/MiniSummaries/' Bcase '/' ...
  Bcase '_allData.fig'])

close

%% Save all the deployment's data
% Make file name
fname_out = [path.dir_base  'bowchain/Data/AllDeployObs/' Bcase '.mat'];
% Save data
if ~isempty(swims.t1) && ~isempty(mmp.eps)
  save(fname_out,'bchain','sadcp','flux','mmp','swims')
elseif isempty(swims.t1) && isempty(fctd.temperature) && ~isempty(mmp.eps)
  save(fname_out,'bchain','sadcp','flux','mmp')
elseif ~isempty(fctd.temperature) && ~isempty(mmp.eps)
  save(fname_out,'bchain','sadcp','flux','mmp','fctd')
elseif ~isempty(swims.t1) && isempty(mmp.eps)
  save(fname_out,'bchain','sadcp','flux','swims')
elseif ~isempty(fctd.temperature) && isempty(mmp.eps)
  save(fname_out,'bchain','sadcp','flux','fctd')
elseif ~isempty(epsi.Temperature) && isempty(mmp.eps)
  save(fname_out,'bchain','sadcp','flux','epsi')
else
  save(fname_out,'bchain','sadcp','flux')
end
end

function offset = lagcov_t_offset(time,T1,T2)
% Smooth signals
win = hanning(60*5);
T1 = filt_ends(win,T1);
T2 = filt_ends(win,T2);
% Get skill of lagged covariance
[lag,skill]=cov_lag(T1,T2,round(numel(T1)/2));
% Find max skill
[~,i_lag] = max(skill);
% Get offset
offset = lag(i_lag)*mean(diff(time));
% Plot resutls
figure
plot(time,T1,time,T2,time+offset,T2)
datetick('x','keeplimits')
legend('Singal 1','Singal 2','Signal 2 + offset')
title(['Offset = ' num2str(offset*86400) ' sec'])
end