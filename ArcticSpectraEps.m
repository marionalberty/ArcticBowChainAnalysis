function ArcticSpectraEps
% Load all bow chain combined data and calculate the horizontal wavenumber
% spectra of the segment nearest the surface. Calculate the mixed layer
% depth (MLD) for all the mmp profiles within the segment and integrate
% the rate of dissipation of kinetic energy (epsilon) over the mixed layer
% from the available eps observations. Use observed chi to calculate the
% heat fluxed through the base of the mixed layer and save the mean
% temperature within the spectral segment, dT/dz, MLD, lat, lon, and dnum.


set(0,'defaultaxesfontsize',16,'defaultaxeslinewidth',0.7,...
  'defaultlinelinewidth',1,'defaultpatchlinewidth',0.7,...
  'defaultFigureColor','white')


%% Set path and params
% Personnal paths:
driveName = '/Users/marionsofiaalberty/MATLAB/Arctic/ArcticMix/bowchain/';
% Inputs/Outputs directories:
pathin = [driveName 'Data/AllDeployObs/'];
dataOut = [driveName 'Data/Spectra/'];

% Spectral segement length
L = 30;   % [km]

% Mixed layer delta sgth from Spicy Seas paper (MacKinnon et al. 2016)
MLS = 0.5;  % [kg/m^{-3}]
% Second option for well mixed water column that fails first test
MLS2 = 0.03;  % [kg/m^{-3}]

% Bootstrap constants for confidence interval
alpha = 0.1;
nboot = 1000;


%% Initialize data structures
data.deployment = {};   data.z = [];            data.Tbar = [];
data.lat = [];          data.lon = [];          data.dn_range = [];
data.slope_lowK = [];   data.slope_midK = [];   data.slope_highK = [];
data.var_lowK = [];     data.var_midK = [];     data.var_highK = [];
data.MLD = [];          data.dTdz_ML = [];         data.heat_flux = [];
data.eps_bar = [];      data.eps_ci = [];
data.chi_bar = [];      data.chi_ci = [];


%% Get list of files
flist = dir(fullfile(pathin,'*.mat'));


for i_dep = 1:numel(flist)
  %% Load data, check for mmp observations
  % Clear previous observations
  clear bchain mmp
  % Load current deployment
  load([pathin flist(i_dep).name],'bchain','mmp')
  % Check if deployment had mmp
  if ~exist('mmp','var')
    continue
  end
  
  
  %% Interpolate bowchain onto regular grid
  % Get distances bewteen points
  xx = m_lldist(bchain.lon,bchain.lat);
  % Get distance along track
  x = [0 cumsum(xx)'];
  % Chose dx
  dx = round(mean(diff(x)),4);
  if dx > 0
    % Get list of flds
    flds = fieldnames(bchain);
    % Set along track distance
    bchain.x = 0:dx:max(x);
    % Grid in space
    for i_f = 1:numel(flds)
      [~,nc] = size(bchain.(flds{i_f}));
      if nc > 1
        bchain.(flds{i_f}) = interp1(x,bchain.(flds{i_f})',bchain.x)';
      end
    end
  end
  % Housekeeping
  clear xx x dx nc flds
  
  
  %% Calculate MLD, dTdz_ML, chi, and heat flux for MMP profiles
  % Initialize matrices for mmp quantities
  mmp.sgth_srtd = mmp.sgth*nan;
  mmp.MLD = mmp.dnum*nan;
  mmp.dTdz_ML = mmp.dnum*nan;
  mmp.chi_bar = mmp.dnum*nan;
  mmp.rhoML_bar = mmp.dnum*nan;
  mmp.heat_flux = mmp.dnum*nan;
  % Sort each profile and get MLD
  for i_dn = 1:numel(mmp.dnum)
    % Find non-nans
    i_good = ~isnan(mmp.sgth(:,i_dn));
    % Sort sgth, preserving depths
    mmp.sgth_srtd(i_good,i_dn) = sort(mmp.sgth(i_good,i_dn));
    % Find min sgth
    sgth_min = min(mmp.sgth_srtd(:,i_dn));
    % Get index of MLD
    iBase = find(mmp.sgth_srtd(:,i_dn) > sgth_min+MLS,1,'first');
    % If well mixed, try deBoyer et al. 2014 criterion
    if isempty(iBase)
      iBase = find(mmp.sgth_srtd(:,i_dn) > sgth_min+MLS2,1,'first');
    end
    % If still no mixed layer detected
    if isempty(iBase)
      if sum(~isnan(mmp.sgth_srtd(:,i_dn))) > 10
        % There is a profile but very well mixed to the bottom
        [~,iBase] = max(mmp.sgth_srtd(:,i_dn));
      else
        % Not enough data
        continue
      end
    end
    % Save MLD
    mmp.MLD(i_dn) = mmp.p(iBase);
    % Mean density in mixed layer
    mmp.rhoML_bar(i_dn) = nanmean(mmp.sgth(1:iBase,i_dn));
    % Calculate dTdz_ML over 10 m window
    dT = mean(mmp.t(mmp.p >= mmp.MLD(i_dn)+4 & ...
      mmp.p <= mmp.MLD(i_dn)+6,i_dn)) - ...
      mean(mmp.t(mmp.p >= mmp.MLD(i_dn)-6 & ...
      mmp.p <= mmp.MLD(i_dn)-4,i_dn));
    mmp.dTdz_ML(i_dn) = -dT/10;
    % Get mean chi over the base of the mixed layer
    chi = mmp.chi(mmp.p >= mmp.MLD(i_dn)-5 & ...
      mmp.p <= mmp.MLD(i_dn)+5,i_dn);
    chi(isnan(chi)) = [];
    % Get statistics for chi
    if numel(chi) == 0
      % if empty set quantities to nan
      mmp.chi_bar(i_dn) = nan;
    elseif numel(chi) == 1
      % if only one chi, save as mean and nan as ci
      mmp.chi_bar(i_dn) = chi;
    else
      % Two or more chi
      mmp.chi_bar(i_dn) =  mean(bootstrp(nboot,@mean,chi));
    end
    % Calculate heat flux [W/m^2]
    Cp = nanmean(sw_cp(mmp.s(1:iBase,i_dn),mmp.t(1:iBase,i_dn),...
      mmp.p(1:iBase)));
    mmp.heat_flux(i_dn) = -mmp.rhoML_bar(i_dn)*Cp*mmp.chi_bar(i_dn)/...
      (2*mmp.dTdz_ML(i_dn));
  end
  
  
  %% Calculate spectra for segments
  % Number of points for fft
  dx = bchain.x(2) - bchain.x(1);
  N_fft = floor(L/dx/2)*2; % 30 km segments
  % Get segment time stamps
  dn_seg = mk_segment(bchain.dn',N_fft,0.5,0)';
  
  % Skip if deployment if it is too short
  if isempty(dn_seg)
    continue
  end
  % Get segment lats and lons
  lat_seg = mk_segment(bchain.lat',N_fft,0.5,0)';
  lon_seg = mk_segment(bchain.lon',N_fft,0.5,0)';
  % Get number of segments
  N_seg = numel(dn_seg(:,1));
  % Calculate spectra and spectral quantities for each segment
  for i_seg = 1:N_seg
    % Save deployment name
    data.deployment{end+1,1} = flist(i_dep).name(1:end-4);
    % Save spectra dnum range
    data.dn_range(end+1,:) = [dn_seg(i_seg,1) dn_seg(i_seg,end)];
    % Save spectra mean lat and lon
    data.lat(end+1) = mean(lat_seg(i_seg,:));
    data.lon(end+1) = mean(lon_seg(i_seg,:));
    
    % Set indicator for while loop
    i_z = 1;
    % Use while loop to get shallowest full temperature segment
    while i_z > 0
      % Make half-overlapping segments
      t_seg = mk_segment(bchain.t(i_z,:)',N_fft,0.5,0)';
      t_seg = t_seg(i_seg,:);
      % Calculate fraction of timeseries that is nan
      nanFrac = sum(isnan(t_seg))./N_fft;
      % If less than 5% of the vector is nan, interp and move on to spectra
      if nanFrac < 0.05
        % Interp over nans
        t_seg = naninterp(t_seg);
        % Save spectra depth
        data.z(end+1,1) = bchain.z(i_z);
        % Save mean segment temperature
        data.Tbar(end+1,1) = mean(t_seg);
        % Update index to break loop
        i_z = 0;
      else
        % Else, try the next depth
        i_z = i_z +1;
      end
    end
    
    %% Get mean MLD, dTdz, heat flux, eps, and chi for segment
    % Get indecies of profiles during segment
    i_profiles = mmp.dnum >= dn_seg(i_seg,1) & ...
      mmp.dnum <= dn_seg(i_seg,end);
    % Mean MLD
    data.MLD(end+1,1) = nanmean(mmp.MLD(i_profiles));
    % Mean dTdz
    data.dTdz_ML(end+1,1) = nanmean(mmp.dTdz_ML(i_profiles));
    % Mean heat flux
    data.heat_flux(end+1,1) = nanmean(mmp.heat_flux(i_profiles));
    
    % Find all posible eps
    eps = mmp.eps(mmp.z <= data.MLD(end),i_profiles);
    % Find all posible chi
    chi = mmp.chi(mmp.z <= data.MLD(end),i_profiles);
    % Remove nans
    eps(isnan(eps)) = [];
    chi(isnan(chi)) = [];
    
    % Get statistics for eps
    if numel(eps) == 0
      % if empty set quantities to nan
      data.eps_bar(end+1,1) = nan;
      data.eps_ci(end+1,:) = [nan nan];
    elseif numel(eps) == 1
      % if only one eps, save as mean and nan as ci
      data.eps_bar(end+1,1) = eps;
      data.eps_ci(end+1,:) = [nan nan];
    else
      % Two or more eps
      data.eps_bar(end+1,1) =  mean(bootstrp(nboot,@mean,eps));
      CI = bootci(nboot,{@mean,eps},'type','per','alpha',alpha)';
      data.eps_ci(end+1,:) = [data.eps_bar(end) - CI(1) ...
        CI(2) - data.eps_bar(end)];
    end
    
    % Get statistics for chi
    if numel(chi) == 0
      % if empty set quantities to nan
      data.chi_bar(end+1,1) = nan;
      data.chi_ci(end+1,:) = [nan nan];
    elseif numel(chi) == 1
      % if only one chi, save as mean and nan as ci
      data.chi_bar(end+1,1) = chi;
      data.chi_ci(end+1,:) = [nan nan];
    else
      % Two or more chi
      data.chi_bar(end+1,1) =  mean(bootstrp(nboot,@mean,chi));
      CI = bootci(nboot,{@mean,chi},'type','per','alpha',alpha)';
      data.chi_ci(end+1,:) = [data.chi_bar(end) - CI(1) ...
        CI(2) - data.chi_bar(end)];
    end
    
    
    %% Calculate spectral quantities
    % Calc spectra for segment
    [spec,k] = spectrum(t_seg,bchain.x(1:N_fft));
    % Low wavenumber range (wavelengths >= 1 km)
    ik = k <= 1;
    % Estimate slope
    fobj = fit(log(k(ik))',log(spec(ik))','poly1');
    p = coeffvalues(fobj);
    q = confint(fobj);
    data.slope_lowK(end+1,:) = [p(1) p(1)-q(1,1)];
    % Quantify spectral energy
    data.var_lowK(end+1,1) = trapz(k(ik),spec(ik));
    
    % Mid wavenumber range (wavelengths < 500 m & >= 50 m)
    ik = k > 1/0.5 & k <= 1/0.05;
    fobj = fit(log(k(ik))',log(spec(ik))','poly1');
    p = coeffvalues(fobj);
    q = confint(fobj);
    data.slope_midK(end+1,:) = [p(1) p(1)-q(1,1)];
    % Quantify spectral energy
    data.var_midK(end+1,1) = trapz(k(ik),spec(ik));
    
    % High wavenumber range (wavelengths < max MLD)
    mmld = min([0.05 max(mmp.MLD(i_profiles))/1000]);
    mmld = max([mmld 2.1*dx]);
    ik = k >= 1/mmld;
    fobj = fit(log(k(ik))',log(spec(ik))','poly1');
    p = coeffvalues(fobj);
    q = confint(fobj);
    data.slope_highK(end+1,:) = [p(1) p(1)-q(1,1)];
    % Quantify spectral energy
    data.var_highK(end+1,1) = trapz(k(ik),spec(ik));
  end
  % End deployment loop
end


%% Save data
save([dataOut 'allMMPdeployments.mat'],'data')
end