function ArcticSigmaSpectraChi
% Load all bow chain data with chi estimates and surface density horizontal
% wavenumber spectra

warning off

set(0,'defaultaxesfontsize',16,'defaultaxeslinewidth',0.7,...
  'defaultlinelinewidth',1,'defaultpatchlinewidth',0.7,...
  'defaultFigureColor','white')


%% Set path and params
% Personnal paths:
driveName = '/Users/marionsofiaalberty/MATLAB/Arctic/';
% Inputs/Outputs directories:
pathinA = [driveName 'ArcticMix/bowchain/Data/Spectra/'];
pathinS = [driveName 'SODA/bowchain/Data/Spectra/'];
dataOut = [driveName 'Data/Spectra/'];
figOut = [driveName 'Figures/Spectra/'];

% Spectral segement length
L = 30;     % [km]
z_max = 8;  % [m]

% Bootstrap constants for confidence interval
alpha = 0.05;
nboot = 1000;


%% Get list of files
flistA = dir(fullfile(pathinA,'*.mat'));
flistS = dir(fullfile(pathinS,'*.mat'));
flist = [flistA(3:end); flistS];
% Housekeeping
clear flistA flistS


%% Initialize data structures
data.deployment = {};   data.cruise = {};       data.z = [];
data.dn_range = [];     data.lat = [];          data.lon = [];
data.k = {};            data.Pow_T = {};        data.Pow_Sig = {};
data.k_mros = [];       data.k_mld = [];        data.k_min = [];
data.k_ship = [];
data.slopeT_lowK = [];  data.varT_lowK = [];
data.slopeSig_lowK = [];data.varSig_lowK = [];
data.slopeT_smK = [];   data.varT_smK = [];
data.slopeSig_smK = []; data.varSig_smK = [];
data.chi_bar = [];      data.chi_ci = [];


for i_dep = 1:numel(flist)
  %% Load data
  % Clear previous observations
  clear bchain
  % Load current deployment
  load([flist(i_dep).folder '/' flist(i_dep).name],'bchain')
  % Check if deployment has bchain chi
  if ~isfield(bchain,'chi')
    continue
  end
  % Get Cruise Name
  strPath = strsplit(flist(i_dep).folder,'/');
  Cruise = strPath{6};
  
  
  %% Extract deployment data and calculate number of segments
  % Check if shallowest observations are in the mixed layer
  if bchain.TSfiltered.z(1) > z_max
    continue
  end
  % Extract temperature and sigma data to calc spectra
  T = bchain.TSfiltered.t(1,:);
  Sigma = bchain.TSfiltered.sig0(1,:);
  % Extract x, z and dn
  x = bchain.TSfiltered.x;
  dn_T = bchain.TSfiltered.dn;
  z = bchain.TSfiltered.z(1);
  % Extract Chi for averaging
  chi = bchain.chi.chi(bchain.chi.z == z,:);
  dn_chi = mean(bchain.chi.time);
  % Remove nans
  dn_chi(isnan(chi)) = [];
  chi(isnan(chi)) = [];
  % Number of points for fft
  dx = mean(diff(x));
  N_fft = floor(L/dx/2)*2;    % always even
  % Number of points in temp reccord
  N_T = numel(T);
  % Number of spectral segments
  n_segs = floor(2*(N_T/N_fft)-1);
  
  
  %% Calculate spectra
  if n_segs == 0
    % Segement is shorter than desired segement length, calculate spectra
    % anyways
    % Save deployment information
    % Save deployment name
    data.deployment{end+1,1} = flist(i_dep).name(1:end-4);
    % Save cruise name
    data.cruise{end+1,1} = Cruise;
    % Save spectra mean lat and lon
    data.lat(end+1,1) = mean(bchain.xgrid.lat);
    data.lon(end+1,1) = mean(bchain.xgrid.lon);
    % Save spectra dnum range
    data.dn_range(end+1,:) = [dn_T(1) dn_T(end)];
    % Save spectra depth
    data.z(end+1,1) = z;
    
    % Calc spectra for segment
    [specT,k] = spectrum(T,x);
    [specS,~] = spectrum(Sigma,x);
    data.Pow_T{end+1,1} = specT;
    data.Pow_Sig{end+1,1} = specS;
    data.k{end+1,1} = k;
    data.k_mros(end+1,1) = bchain.spectra.k_mros;
    data.k_mld(end+1,1) = bchain.spectra.k_mld;
    data.k_min(end+1,1) = k(1);
    % Get integration limits to separate submesoscale wavenumber band from
    % chi integration limits
    % k_ship = 1/[(max period from chi)(ship speed)(factor of 3)] [1/km]
    k_ship = 1/(20*bchain.Uship_bar*3/1000);
    data.k_ship(end+1,1) = k_ship;
    
    
    % Estimate low wavenumber slope
    % Low wavenumber range (wavelengths greater than Rossby)
    ik = k <= bchain.spectra.k_mros;
    % Estimate slope for temperature spectra
    fobj = fit(log(k(ik))',log(specT(ik))','poly1');
    p = coeffvalues(fobj);
    q = confint(fobj);
    data.slopeT_lowK(end+1,:) = [p(1) p(1)-q(1,1)];
    % Estimate slope for density spectra
    fobj = fit(log(k(ik))',log(specS(ik))','poly1');
    p = coeffvalues(fobj);
    q = confint(fobj);
    data.slopeSig_lowK(end+1,:) = [p(1) p(1)-q(1,1)];
    % Quantify spectral energy
    data.varT_lowK(end+1,1) = trapz(k(ik),specT(ik));
    data.varSig_lowK(end+1,1) = trapz(k(ik),specS(ik));
    
    % Estimate submesoscale wavenumber slope
    if bchain.spectra.k_mros < k_ship
      % wavelengths less than Rossby and greater than k_ship
      ik = k > bchain.spectra.k_mros & k <= k_ship;
      % Estimate slope for temperature spectra
      fobj = fit(log(k(ik))',log(specT(ik))','poly1');
      p = coeffvalues(fobj);
      q = confint(fobj);
      data.slopeT_smK(end+1,:) = [p(1) p(1)-q(1,1)];
      % Estimate slope for temperature spectra
      fobj = fit(log(k(ik))',log(specS(ik))','poly1');
      p = coeffvalues(fobj);
      q = confint(fobj);
      data.slopeSig_smK(end+1,:) = [p(1) p(1)-q(1,1)];
      % Quantify spectral energy
      data.varT_smK(end+1,1) = trapz(k(ik),specT(ik));
      data.varSig_smK(end+1,1) = trapz(k(ik),specS(ik));
    else
      data.slopeT_smK(end+1,:) = [nan nan];
      data.varT_smK(end+1,1) = nan;
      data.slopeSig_smK(end+1,:) = [nan nan];
      data.varSig_smK(end+1,1) = nan;
    end
    
    % Get mean chi
    data.chi_bar(end+1,1) =  median(bootstrp(nboot,@median,chi));
    CI = bootci(nboot,{@median,chi},'type','per','alpha',alpha)';
    data.chi_ci(end+1,:) = [data.chi_bar(end) - CI(1) ...
      CI(2) - data.chi_bar(end)];
    
  else
    % While loop counter
    i_seg = 1;
    % Set window limits
    i_str = 1;
    i_end = N_fft;
    % Transect is at least one segment length
    while i_seg <= n_segs
      % Save deployment information
      % Save deployment name
      data.deployment{end+1,1} = flist(i_dep).name(1:end-4);
      % Save cruise name
      data.cruise{end+1,1} = Cruise;
      % Save spectra mean lat and lon
      data.lat(end+1,1) = mean(bchain.xgrid.lat(i_str:i_end));
      data.lon(end+1,1) = mean(bchain.xgrid.lon(i_str:i_end));
      % Save spectra dnum range
      data.dn_range(end+1,:) = [dn_T(i_str) dn_T(i_end)];
      % Save spectra depth
      data.z(end+1,1) = z;
      
      % Calc spectra for segment
      [specT,k] = spectrum(T(i_str:i_end),x(i_str:i_end));
      [specS,~] = spectrum(Sigma(i_str:i_end),x(i_str:i_end));
      data.Pow_T{end+1,1} = specT;
      data.Pow_Sig{end+1,1} = specS;
      data.k{end+1,1} = k;
      data.k_mros(end+1,1) = bchain.spectra.k_mros;
      data.k_mld(end+1,1) = bchain.spectra.k_mld;
      data.k_min(end+1,1) = k(1);
      % Get integration limits to separate submesoscale wavenumber band
      % chi integration limits
      % k_ship = 1/[(max period from chi)(ship speed)(factor of 3)] [1/km]
      k_ship = 1/(20*bchain.Uship_bar*3/1000);
      data.k_ship(end+1,1) = k_ship;
      
      % Estimate low wavenumber slope
      % Low wavenumber range (wavelengths greater than Rossby)
      ik = k <= bchain.spectra.k_mros;
      % Estimate slope for temperature spectra
      fobj = fit(log(k(ik))',log(specT(ik))','poly1');
      p = coeffvalues(fobj);
      q = confint(fobj);
      data.slopeT_lowK(end+1,:) = [p(1) p(1)-q(1,1)];
      % Estimate slope for density spectra
      fobj = fit(log(k(ik))',log(specS(ik))','poly1');
      p = coeffvalues(fobj);
      q = confint(fobj);
      data.slopeSig_lowK(end+1,:) = [p(1) p(1)-q(1,1)];
      % Quantify spectral energy
      data.varT_lowK(end+1,1) = trapz(k(ik),specT(ik));
      data.varSig_lowK(end+1,1) = trapz(k(ik),specS(ik));
      
      % Estimate submesoscale wavenumber slope
      if bchain.spectra.k_mros < k_ship
        % wavelengths less than Rossby and greater than k_ship
        ik = k > bchain.spectra.k_mros & k <= k_ship;
        % Estimate slope for temperature spectra
        fobj = fit(log(k(ik))',log(specT(ik))','poly1');
        p = coeffvalues(fobj);
        q = confint(fobj);
        data.slopeT_smK(end+1,:) = [p(1) p(1)-q(1,1)];
        % Estimate slope for temperature spectra
        fobj = fit(log(k(ik))',log(specS(ik))','poly1');
        p = coeffvalues(fobj);
        q = confint(fobj);
        data.slopeSig_smK(end+1,:) = [p(1) p(1)-q(1,1)];
        % Quantify spectral energy
        data.varT_smK(end+1,1) = trapz(k(ik),specT(ik));
        data.varSig_smK(end+1,1) = trapz(k(ik),specS(ik));
      else
        data.slopeT_smK(end+1,:) = [nan nan];
        data.varT_smK(end+1,1) = nan;
        data.slopeSig_smK(end+1,:) = [nan nan];
        data.varSig_smK(end+1,1) = nan;
      end
      
      % Get mean chi
      % Get indecies for chi in segment
      i_chi = dn_chi >= dn_T(i_str) & dn_chi <= dn_T(i_end);
      data.chi_bar(end+1,1) =  median(bootstrp(nboot,@median,chi(i_chi)));
      CI = bootci(nboot,{@median,chi(i_chi)},'type','per','alpha',alpha)';
      data.chi_ci(end+1,:) = [data.chi_bar(end) - CI(1) ...
        CI(2) - data.chi_bar(end)];
      
      % Update segment count
      i_seg = i_seg+1;
      % Update segment start and end
      i_str = i_str + (N_fft/2);
      i_end = i_end + (N_fft/2);
    end
    
    % Check if remaining Segement is long enough to also have spectra
    % calculated
    i_end = min([i_end numel(T)]);
    if x(i_end)-x(i_str) >= 1/bchain.spectra.k_mros
      % Save deployment information
      % Save deployment name
      data.deployment{end+1,1} = flist(i_dep).name(1:end-4);
      % Save cruise name
      data.cruise{end+1,1} = Cruise;
      % Save spectra mean lat and lon
      data.lat(end+1,1) = mean(bchain.xgrid.lat(i_str:i_end));
      data.lon(end+1,1) = mean(bchain.xgrid.lon(i_str:i_end));
      % Save spectra dnum range
      data.dn_range(end+1,:) = [dn_T(i_str) dn_T(i_end)];
      % Save spectra depth
      data.z(end+1,1) = z;
      
      % Calc spectra for segment
      [specT,k] = spectrum(T(i_str:i_end),x(i_str:i_end));
      [specS,~] = spectrum(Sigma(i_str:i_end),x(i_str:i_end));
      data.Pow_T{end+1,1} = specT;
      data.Pow_Sig{end+1,1} = specS;
      data.k{end+1,1} = k;
      data.k_mros(end+1,1) = bchain.spectra.k_mros;
      data.k_mld(end+1,1) = bchain.spectra.k_mld;
      data.k_min(end+1,1) = k(1);
      % Get integration limits to separate submesoscale wavenumber band
      % chi integration limits
      % k_ship = 1/[(max period from chi)(ship speed)(factor of 3)] [1/km]
      k_ship = 1/(20*bchain.Uship_bar*3/1000);
      data.k_ship(end+1,1) = k_ship;
      
      % Estimate low wavenumber slope
      % Low wavenumber range (wavelengths greater than Rossby)
      ik = k <= bchain.spectra.k_mros;
      % Estimate slope for temperature spectra
      fobj = fit(log(k(ik))',log(specT(ik))','poly1');
      p = coeffvalues(fobj);
      q = confint(fobj);
      data.slopeT_lowK(end+1,:) = [p(1) p(1)-q(1,1)];
      % Estimate slope for density spectra
      fobj = fit(log(k(ik))',log(specS(ik))','poly1');
      p = coeffvalues(fobj);
      q = confint(fobj);
      data.slopeSig_lowK(end+1,:) = [p(1) p(1)-q(1,1)];
      % Quantify spectral energy
      data.varT_lowK(end+1,1) = trapz(k(ik),specT(ik));
      data.varSig_lowK(end+1,1) = trapz(k(ik),specS(ik));
      
      % Estimate submesoscale wavenumber slope
      if bchain.spectra.k_mros < k_ship
        % wavelengths less than Rossby and greater than k_ship
        ik = k > bchain.spectra.k_mros & k <= k_ship;
        % Estimate slope for temperature spectra
        fobj = fit(log(k(ik))',log(specT(ik))','poly1');
        p = coeffvalues(fobj);
        q = confint(fobj);
        data.slopeT_smK(end+1,:) = [p(1) p(1)-q(1,1)];
        % Estimate slope for temperature spectra
        fobj = fit(log(k(ik))',log(specS(ik))','poly1');
        p = coeffvalues(fobj);
        q = confint(fobj);
        data.slopeSig_smK(end+1,:) = [p(1) p(1)-q(1,1)];
        % Quantify spectral energy
        data.varT_smK(end+1,1) = trapz(k(ik),specT(ik));
        data.varSig_smK(end+1,1) = trapz(k(ik),specS(ik));
      else
        data.slopeT_smK(end+1,:) = [nan nan];
        data.varT_smK(end+1,1) = nan;
        data.slopeSig_smK(end+1,:) = [nan nan];
        data.varSig_smK(end+1,1) = nan;
      end
      
      % Get mean chi
      % Get indecies for chi in segment
      i_chi = dn_chi >= dn_T(i_str) & dn_chi <= dn_T(i_end);
      data.chi_bar(end+1,1) =  median(bootstrp(nboot,@median,chi(i_chi)));
      CI = bootci(nboot,{@median,chi(i_chi)},'type','per','alpha',alpha)';
      data.chi_ci(end+1,:) = [data.chi_bar(end) - CI(1) ...
        CI(2) - data.chi_bar(end)];
    end
  end
  % End deployment loop
end


%% Save data
save([dataOut 'allDeployments_SigmaSpectraChi.mat'],'data')


%% Scatter plot of chi and temperature low wavenumber slope
% Get colorbar limits for lowK normalized variance
clim = [floor(log10(min(data.varT_lowK./(data.k_mros-data.k_min)))) ...
  ceil(log10(max(data.varT_lowK./(data.k_mros-data.k_min))))];
% Use integer tick marks
cticks = clim(1):clim(2);
% Get colormap for plotting
cmap = jet(numel(cticks)-1);
% Get colorbin for each Spectra
[~,~,cbin] = histcounts(log10(data.varT_lowK./(data.k_mros-data.k_min)),...
  cticks);

figure
for i = 1:numel(data.varT_lowK)
  if ~isnan(data.varT_lowK(i)) && strcmp(data.cruise{i},'SODA')
    errorbar(data.slopeT_lowK(i,1),data.chi_bar(i),data.chi_ci(i,1),...
      data.chi_ci(i,2),data.slopeT_lowK(i,2),data.slopeT_lowK(i,2),'o',...
      'MarkerSize',15,'MarkerEdgeColor',cmap(cbin(i),:),...
      'MarkerFaceColor',cmap(cbin(i),:),'color',cmap(cbin(i),:))
    hold on
  elseif ~isnan(data.varT_lowK(i)) && strcmp(data.cruise{i},'ArcticMix')
    errorbar(data.slopeT_lowK(i,1),data.chi_bar(i),data.chi_ci(i,1),...
      data.chi_ci(i,2),data.slopeT_lowK(i,2),data.slopeT_lowK(i,2),'d',...
      'MarkerSize',15,'MarkerEdgeColor',cmap(cbin(i),:),...
      'MarkerFaceColor',cmap(cbin(i),:),'color',cmap(cbin(i),:))
  end
end
set(gca,'YScale','log')
colormap(cmap)
c = colorbar;
c.Ticks = (0:numel(cticks)-1)./(numel(cticks)-1);
c.TickLabels = num2str(cticks');
xlabel('Mesoscale Spectral Slope (k < k_{Rossby})')
ylabel('\chi_T [K^2/s]')
title('Normalized Mesoscale Temperature Spectral Variance')
ylabel(c,'log_{10}[K^2 km]')

print([figOut 'allDeploys_LKslopeT_Chi'],'-dpng')
savefig([figOut 'allDeploys_LKslopeT_Chi.fig'])


%% Scatter plot of chi and density low wavenumber slope
% Get colorbar limits for lowK normalized variance
clim = [floor(log10(min(data.varSig_lowK./(data.k_mros-data.k_min)))) ...
  ceil(log10(max(data.varSig_lowK./(data.k_mros-data.k_min))))];
% Use integer tick marks
cticks = clim(1):clim(2);
% Get colormap for plotting
cmap = jet(numel(cticks)-1);
% Get colorbin for each Spectra
[~,~,cbin] = histcounts(log10(data.varSig_lowK./...
  (data.k_mros-data.k_min)),cticks);

figure
for i = 1:numel(data.varSig_lowK)
  if ~isnan(data.varSig_lowK(i)) && strcmp(data.cruise{i},'SODA')
    errorbar(data.slopeSig_lowK(i,1),data.chi_bar(i),data.chi_ci(i,1),...
      data.chi_ci(i,2),data.slopeSig_lowK(i,2),...
      data.slopeSig_lowK(i,2),'o','MarkerSize',15,'MarkerEdgeColor',...
      cmap(cbin(i),:),'MarkerFaceColor',cmap(cbin(i),:),'color',...
      cmap(cbin(i),:))
    hold on
  elseif ~isnan(data.varSig_lowK(i)) && strcmp(data.cruise{i},'ArcticMix')
    errorbar(data.slopeSig_lowK(i,1),data.chi_bar(i),data.chi_ci(i,1),...
      data.chi_ci(i,2),data.slopeSig_lowK(i,2),...
      data.slopeSig_lowK(i,2),'d','MarkerSize',15,'MarkerEdgeColor',...
      cmap(cbin(i),:),'MarkerFaceColor',cmap(cbin(i),:),'color',...
      cmap(cbin(i),:))
  end
end
set(gca,'YScale','log')
colormap(cmap)
c = colorbar;
c.Ticks = (0:numel(cticks)-1)./(numel(cticks)-1);
c.TickLabels = num2str(cticks');
xlabel('Mesoscale Spectral Slope (k < k_{Rossby})')
ylabel('\chi_T [K^2/s]')
title('Normalized Mesoscale Density Spectral Variance')
ylabel(c,'log_{10}[(kg/m^3)^2 km]')

print([figOut 'allDeploys_LKslopeSigma_Chi'],'-dpng')
savefig([figOut 'allDeploys_LKslopeSigma_Chi.fig'])


%% Scatter plot of chi and submesoscale Temperature wavenumber slope
% Get colorbar limits for smK variance
clim = [floor(log10(min(data.varT_smK./(data.k_ship-data.k_mros)))) ...
  ceil(log10(max(data.varT_smK./(data.k_ship-data.k_mros))))];
% Use integer tick marks
cticks = clim(1):clim(2);
% Get colormap for plotting
cmap = jet(numel(cticks)-1);
% Get colorbin for each Spectra
[~,~,cbin] = histcounts(log10(data.varT_smK./(data.k_ship-data.k_mros)),...
  cticks);

figure
for i = 1:numel(data.varT_smK)
  if ~isnan(data.slopeT_smK(i)) && strcmp(data.cruise{i},'SODA')
    errorbar(data.slopeT_smK(i,1),data.chi_bar(i),data.chi_ci(i,1),...
      data.chi_ci(i,2),data.slopeT_smK(i,2),data.slopeT_smK(i,2),'o',...
      'MarkerSize',15,'MarkerEdgeColor',cmap(cbin(i),:),...
      'MarkerFaceColor',cmap(cbin(i),:),'color',cmap(cbin(i),:))
    hold on
  elseif ~isnan(data.slopeT_smK(i)) && strcmp(data.cruise{i},'ArcticMix')
    errorbar(data.slopeT_smK(i,1),data.chi_bar(i),data.chi_ci(i,1),...
      data.chi_ci(i,2),data.slopeT_smK(i,2),data.slopeT_smK(i,2),'d',...
      'MarkerSize',15,'MarkerEdgeColor',cmap(cbin(i),:),...
      'MarkerFaceColor',cmap(cbin(i),:),'color',cmap(cbin(i),:))
    hold on
  end
end
set(gca,'YScale','log')
colormap(cmap)
c = colorbar;
c.Ticks = (0:numel(cticks)-1)./(numel(cticks)-1);
c.TickLabels = num2str(cticks');
xlabel('Submesoscale Spectral Slope (k_{Rossby} < k < k_{MLD})')
ylabel('\chi_T [K^2/s]')
title('Normalized Submesoscale Temperature Spectral Variance')
ylabel(c,'log_{10}[K^2 km]')

print([figOut 'allDeploys_SmKslopeT_Chi'],'-dpng')
savefig([figOut 'allDeploys_SmKslopeT_Chi.fig'])


%% Scatter plot of chi and submesoscale density wavenumber slope
% Get colorbar limits for smK variance
clim = [floor(log10(min(data.varSig_smK./(data.k_ship-data.k_mros)))) ...
  ceil(log10(max(data.varSig_smK./(data.k_ship-data.k_mros))))];
% Use integer tick marks
cticks = clim(1):clim(2);
% Get colormap for plotting
cmap = jet(numel(cticks)-1);
% Get colorbin for each Spectra
[~,~,cbin] = histcounts(log10(data.varSig_smK./...
  (data.k_ship-data.k_mros)),cticks);

figure
for i = 1:numel(data.varSig_smK)
  if ~isnan(data.slopeSig_smK(i)) && strcmp(data.cruise{i},'SODA')
    errorbar(data.slopeSig_smK(i,1),data.chi_bar(i),data.chi_ci(i,1),...
      data.chi_ci(i,2),data.slopeSig_smK(i,2),data.slopeSig_smK(i,2),...
      'o','MarkerSize',15,'MarkerEdgeColor',cmap(cbin(i),:),...
      'MarkerFaceColor',cmap(cbin(i),:),'color',cmap(cbin(i),:))
    hold on
  elseif ~isnan(data.slopeSig_smK(i)) && strcmp(data.cruise{i},'ArcticMix')
    errorbar(data.slopeSig_smK(i,1),data.chi_bar(i),data.chi_ci(i,1),...
      data.chi_ci(i,2),data.slopeSig_smK(i,2),data.slopeSig_smK(i,2),...
      'd','MarkerSize',15,'MarkerEdgeColor',cmap(cbin(i),:),...
      'MarkerFaceColor',cmap(cbin(i),:),'color',cmap(cbin(i),:))
    hold on
  end
end
set(gca,'YScale','log')
colormap(cmap)
c = colorbar;
c.Ticks = (0:numel(cticks)-1)./(numel(cticks)-1);
c.TickLabels = num2str(cticks');
xlabel('Submesoscale Spectral Slope (k_{Rossby} < k < k_{MLD})')
ylabel('\chi_T [K^2/s]')
title('Normalized Submesoscale Density Spectral Variance')
ylabel(c,'log_{10}[(k/m^3)^2 km]')

print([figOut 'allDeploys_SmKslopeSigma_Chi'],'-dpng')
savefig([figOut 'allDeploys_SmKslopeSigma_Chi.fig'])


%% Scatter plot of chi and submesoscale temperature wavenumber variance
%  normalized by submesoscale wavenumber range
% Get colorbar limits for lowK variance
clim = [floor(log10(min(data.varT_lowK./(data.k_mros-data.k_min)))) ...
  ceil(log10(max(data.varT_lowK./(data.k_mros-data.k_min))))];
% Use integer tick marks
cticks = clim(1):clim(2);
% Get colormap for plotting
cmap = jet(numel(cticks)-1);
% Get colorbin for each Spectra
[~,~,cbin] = histcounts(log10(data.varT_lowK./(data.k_mros-data.k_min)),...
  cticks);

figure
for i = 1:numel(data.varT_smK)
  if ~isnan(data.varT_smK(i)) && strcmp(data.cruise{i},'SODA')
    errorbar(data.varT_smK(i)/(data.k_ship(i)-data.k_mros(i)),...
      data.chi_bar(i),data.chi_ci(i,1),data.chi_ci(i,2),'o',...
      'MarkerSize',15,'MarkerEdgeColor',cmap(cbin(i),:),...
      'MarkerFaceColor',cmap(cbin(i),:),'color',cmap(cbin(i),:))
    hold on
  elseif ~isnan(data.varT_smK(i)) && strcmp(data.cruise{i},'ArcticMix')
    errorbar(data.varT_smK(i)/(data.k_ship(i)-data.k_mros(i)),...
      data.chi_bar(i),data.chi_ci(i,1),data.chi_ci(i,2),'d',...
      'MarkerSize',15,'MarkerEdgeColor',cmap(cbin(i),:),...
      'MarkerFaceColor',cmap(cbin(i),:),'color',cmap(cbin(i),:))
    hold on
  end
end
set(gca,'YScale','log')
set(gca,'XScale','log')
colormap(cmap)
c = colorbar;
c.Ticks = (0:numel(cticks)-1)./(numel(cticks)-1);
c.TickLabels = num2str(cticks');
title({'Normalized Temperature Variance';...
  'Submesoscale (k_{Rossby} < k < k_{ship})'})
xlabel('[K^2 km]')
ylabel('\chi_T [K^2/s]')
ylabel(c,['Normalized Mesoscale (k < k_{Rossby}) Variance ' ...
  'log_{10}[K^2 km]'])

print([figOut 'allDeploys_SmKvarT_Chi'],'-dpng')
savefig([figOut 'allDeploys_SmKvarT_Chi.fig'])


%% Scatter plot of chi and submesoscale density wavenumber variance
%  normalized by submesoscale wavenumber range
% Get colorbar limits for lowK variance
clim = [floor(log10(min(data.varSig_lowK./(data.k_mros-data.k_min)))) ...
  ceil(log10(max(data.varSig_lowK./(data.k_mros-data.k_min))))];
% Use integer tick marks
cticks = clim(1):clim(2);
% Get colormap for plotting
cmap = jet(numel(cticks)-1);
% Get colorbin for each Spectra
[~,~,cbin] = histcounts(log10(data.varSig_lowK./...
  (data.k_mros-data.k_min)),cticks);

figure
for i = 1:numel(data.varSig_smK)
  if ~isnan(data.varSig_smK(i)) && strcmp(data.cruise{i},'SODA')
    errorbar(data.varSig_smK(i)/(data.k_ship(i)-data.k_mros(i)),...
      data.chi_bar(i),data.chi_ci(i,1),data.chi_ci(i,2),'o',...
      'MarkerSize',15,'MarkerEdgeColor',cmap(cbin(i),:),...
      'MarkerFaceColor',cmap(cbin(i),:),'color',cmap(cbin(i),:))
    hold on
  elseif ~isnan(data.varSig_smK(i)) && strcmp(data.cruise{i},'ArcticMix')
    errorbar(data.varSig_smK(i)/(data.k_ship(i)-data.k_mros(i)),...
      data.chi_bar(i),data.chi_ci(i,1),data.chi_ci(i,2),'d',...
      'MarkerSize',15,'MarkerEdgeColor',cmap(cbin(i),:),...
      'MarkerFaceColor',cmap(cbin(i),:),'color',cmap(cbin(i),:))
    hold on
  end
end
set(gca,'YScale','log')
set(gca,'XScale','log')
colormap(cmap)
c = colorbar;
c.Ticks = (0:numel(cticks)-1)./(numel(cticks)-1);
c.TickLabels = num2str(cticks');
title({'Normalized Density Variance';...
  'Submesoscale (k_{Rossby} < k < k_{ship})'})
xlabel('[(kg/m^3)^2 km]')
ylabel('\chi_T [K^2/s]')
ylabel(c,['Normalized Mesoscale (k < k_{Rossby}) Variance ' ...
  'log_{10}[(kg/m^3)^2 km]'])

print([figOut 'allDeploys_SmKvarSigma_Chi'],'-dpng')
savefig([figOut 'allDeploys_SmKvarSigma_Chi.fig'])


%% Submesoscale vs Mesoscale Normalized Temperature Spectral Variance
% Use integer tick marks
cticks = -3.5:-0.5;
% Get colormap for plotting
cmap = jet(numel(cticks)-1);
% Get colorbin for each Spectra
[~,~,cbin] = histcounts(data.slopeT_lowK(:,1),cticks);

figure
for i = 1:numel(data.varT_smK)
  if ~isnan(data.varT_smK(i)) && strcmp(data.cruise{i},'SODA')
    scatter(data.varT_lowK(i)/(data.k_mros(i)-data.k_min(i)),...
      data.varT_smK(i)/(data.k_ship(i)-data.k_mros(i)),200,'o',...
      'MarkerEdgeColor',cmap(cbin(i),:),'MarkerFaceColor',...
      cmap(cbin(i),:))
    hold on
  elseif ~isnan(data.varT_smK(i)) && strcmp(data.cruise{i},'ArcticMix')
    scatter(data.varT_lowK(i)/(data.k_mros(i)-data.k_min(i)),...
      data.varT_smK(i)/(data.k_ship(i)-data.k_mros(i)),200,'d',...
      'MarkerEdgeColor',cmap(cbin(i),:),'MarkerFaceColor',...
      cmap(cbin(i),:))
    hold on
  end
end
set(gca,'YScale','log')
set(gca,'XScale','log')
c = colorbar;
c.Ticks = (0:numel(cticks)-1)./(numel(cticks)-1);
c.TickLabels = num2str(cticks');
xlabel('Mesoscale (k < k_{Rossby})')
ylabel('Submesoscale (k_{Rossby} < k < k_{MLD})')
title('Normalized Temperature Spectral Variance [K^2 km]')
ylabel(c,'Mesoscale Spectral Slope')
colormap(cmap)

print([figOut 'allDeploys_SmKvarT_LKvarT'],'-dpng')
savefig([figOut 'allDeploys_SmKvarT_LKvarT.fig'])


%% Submesoscale vs Mesoscale Normalized Density Spectral Variance
% Use integer tick marks
cticks = [-3.5:-0.5 0];
% Get colormap for plotting
cmap = jet(numel(cticks)-1);
% Get colorbin for each Spectra
[~,~,cbin] = histcounts(data.slopeSig_lowK(:,1),cticks);

figure
for i = 1:numel(data.varSig_smK)
  if ~isnan(data.varSig_smK(i)) && strcmp(data.cruise{i},'SODA')
    scatter(data.varSig_lowK(i)/(data.k_mros(i)-data.k_min(i)),...
      data.varSig_smK(i)/(data.k_ship(i)-data.k_mros(i)),200,'o',...
      'MarkerEdgeColor',cmap(cbin(i),:),'MarkerFaceColor',...
      cmap(cbin(i),:))
    hold on
  elseif ~isnan(data.varSig_smK(i)) && strcmp(data.cruise{i},'ArcticMix')
    scatter(data.varSig_lowK(i)/(data.k_mros(i)-data.k_min(i)),...
      data.varSig_smK(i)/(data.k_ship(i)-data.k_mros(i)),200,'d',...
      'MarkerEdgeColor',cmap(cbin(i),:),'MarkerFaceColor',...
      cmap(cbin(i),:))
    hold on
  end
end
set(gca,'YScale','log')
set(gca,'XScale','log')
c = colorbar;
c.Ticks = (0:numel(cticks)-1)./(numel(cticks)-1);
c.TickLabels = num2str(cticks');
xlabel('Mesoscale (k < k_{Rossby})')
ylabel('Submesoscale (k_{Rossby} < k < k_{ship})')
title('Normalized Density Spectral Variance [(kg/m^3)^2 km]')
ylabel(c,'Mesoscale Spectral Slope')
colormap(cmap)

print([figOut 'allDeploys_SmKvarSigma_LKvarSigma'],'-dpng')
savefig([figOut 'allDeploys_SmKvarSigma_LKvarSigma.fig'])


%% Submesoscale vs. Mesoscale Temperature Spectral Slope
% Get colorbar limits for lowK normalized variance
clim = [floor(log10(min(data.varT_lowK./(data.k_mros-data.k_min)))) ...
  ceil(log10(max(data.varT_lowK./(data.k_mros-data.k_min))))];
% Use integer tick marks
cticks = clim(1):clim(2);
% Get colormap for plotting
cmap = jet(numel(cticks)-1);
% Get colorbin for each Spectra
[~,~,cbin] = histcounts(log10(data.varT_lowK./(data.k_mros-data.k_min)),...
  cticks);

figure
for i = 1:numel(data.slopeT_smK(:,1))
  if ~isnan(data.slopeT_smK(i)) && strcmp(data.cruise{i},'SODA')
    errorbar(data.slopeT_lowK(i,1),data.slopeT_smK(i,1),...
      data.slopeT_smK(i,2),data.slopeT_smK(i,2),...
      data.slopeT_lowK(i,2),data.slopeT_lowK(i,2),'o',...
      'MarkerSize',15,'MarkerEdgeColor',cmap(cbin(i),:),...
      'MarkerFaceColor',cmap(cbin(i),:),'color',cmap(cbin(i),:))
    hold on
  elseif ~isnan(data.slopeT_smK(i)) && strcmp(data.cruise{i},'ArcticMix')
    errorbar(data.slopeT_lowK(i,1),data.slopeT_smK(i,1),...
      data.slopeT_smK(i,2),data.slopeT_smK(i,2),...
      data.slopeT_lowK(i,2),data.slopeT_lowK(i,2),'d',...
      'MarkerSize',15,'MarkerEdgeColor',cmap(cbin(i),:),...
      'MarkerFaceColor',cmap(cbin(i),:),'color',cmap(cbin(i),:))
    hold on
  end
end
c = colorbar;
c.Ticks = (0:numel(cticks)-1)./(numel(cticks)-1);
c.TickLabels = num2str(cticks');
xlabel('Mesoscale (k < k_{Rossby})')
ylabel('Submesoscale (k_{Rossby} < k < k_{MLD})')
title('Temperature Spectral Slope')
ylabel(c,'Mesoscale Variance log_{10}[K^2 km]')
colormap(cmap)

print([figOut 'allDeploys_SmKslopeT_LKslopeT'],'-dpng')
savefig([figOut 'allDeploys_SmKslopeT_LKslopeT.fig'])


%% Submesoscale vs. Mesoscale Density Spectral Slope
% Get colorbar limits for lowK normalized variance
clim = [floor(log10(min(data.varSig_lowK./(data.k_mros-data.k_min)))) ...
  ceil(log10(max(data.varSig_lowK./(data.k_mros-data.k_min))))];
% Use integer tick marks
cticks = clim(1):clim(2);
% Get colormap for plotting
cmap = jet(numel(cticks)-1);
% Get colorbin for each Spectra
[~,~,cbin] = histcounts(log10(data.varSig_lowK./...
  (data.k_mros-data.k_min)),cticks);

figure
for i = 1:numel(data.slopeSig_smK(:,1))
  if ~isnan(data.slopeSig_smK(i)) && strcmp(data.cruise{i},'SODA')
    errorbar(data.slopeSig_lowK(i,1),data.slopeSig_smK(i,1),...
      data.slopeSig_smK(i,2),data.slopeSig_smK(i,2),...
      data.slopeSig_lowK(i,2),data.slopeSig_lowK(i,2),'o',...
      'MarkerSize',15,'MarkerEdgeColor',cmap(cbin(i),:),...
      'MarkerFaceColor',cmap(cbin(i),:),'color',cmap(cbin(i),:))
    hold on
  elseif ~isnan(data.slopeSig_smK(i)) && strcmp(data.cruise{i},'ArcticMix')
    errorbar(data.slopeSig_lowK(i,1),data.slopeSig_smK(i,1),...
      data.slopeSig_smK(i,2),data.slopeSig_smK(i,2),...
      data.slopeSig_lowK(i,2),data.slopeSig_lowK(i,2),'d',...
      'MarkerSize',15,'MarkerEdgeColor',cmap(cbin(i),:),...
      'MarkerFaceColor',cmap(cbin(i),:),'color',cmap(cbin(i),:))
    hold on
  end
end
c = colorbar;
c.Ticks = (0:numel(cticks)-1)./(numel(cticks)-1);
c.TickLabels = num2str(cticks');
xlabel('Mesoscale (k < k_{Rossby})')
ylabel('Submesoscale (k_{Rossby} < k < k_{MLD})')
title('Density Spectral Slope')
ylabel(c,'Mesoscale Variance log_{10}[(kg/m^3)^2 km]')
colormap(cmap)

print([figOut 'allDeploys_SmKslopeSig_LKslopeSig'],'-dpng')
savefig([figOut 'allDeploys_SmKslopeSig_LKslopeSig.fig'])


%% Density vs. Temperature Mesoscale Spectral Slopes
figure
for i = 1:numel(data.slopeSig_lowK(:,1))
  if ~isnan(data.slopeSig_lowK(i)) && strcmp(data.cruise{i},'SODA')
    errorbar(data.slopeSig_lowK(i,1),data.slopeT_lowK(i,1),...
      data.slopeT_lowK(i,2),data.slopeT_lowK(i,2),...
      data.slopeSig_lowK(i,2),data.slopeSig_lowK(i,2),'o',...
      'MarkerSize',15,'MarkerEdgeColor','r','MarkerFaceColor','r',...
      'color','r')
    hold on
  elseif ~isnan(data.slopeSig_lowK(i)) && strcmp(data.cruise{i},'ArcticMix')
    errorbar(data.slopeSig_lowK(i,1),data.slopeT_lowK(i,1),...
      data.slopeT_lowK(i,2),data.slopeT_lowK(i,2),...
      data.slopeSig_lowK(i,2),data.slopeSig_lowK(i,2),'d',...
      'MarkerSize',15,'MarkerEdgeColor','b',...
      'MarkerFaceColor','b','color','b')
    hold on
  end
end
xlabel('Density')
ylabel('Temperature')
plot([-5 1],[-5 1],':k')
title('Mesoscale (k < k_{Rossby}) Spectral Slope')

print([figOut 'allDeploys_LKslopeT_LKslopeSig'],'-dpng')
savefig([figOut 'allDeploys_LKslopeT_LKslopeSig.fig'])


%% Density vs. Temperature Submesoscale Spectral Slopes
figure
for i = 1:numel(data.slopeSig_smK(:,1))
  if ~isnan(data.slopeSig_smK(i)) && strcmp(data.cruise{i},'SODA')
    errorbar(data.slopeSig_smK(i,1),data.slopeT_smK(i,1),...
      data.slopeT_smK(i,2),data.slopeT_smK(i,2),...
      data.slopeSig_smK(i,2),data.slopeSig_smK(i,2),'o',...
      'MarkerSize',15,'MarkerEdgeColor','r','MarkerFaceColor','r',...
      'color','r')
    hold on
  elseif ~isnan(data.slopeSig_smK(i)) && strcmp(data.cruise{i},'ArcticMix')
    errorbar(data.slopeSig_smK(i,1),data.slopeT_smK(i,1),...
      data.slopeT_smK(i,2),data.slopeT_smK(i,2),...
      data.slopeSig_smK(i,2),data.slopeSig_smK(i,2),'d',...
      'MarkerSize',15,'MarkerEdgeColor','b',...
      'MarkerFaceColor','b','color','b')
    hold on
  end
end
xlabel('Density')
ylabel('Temperature')
plot([-6 2],[-6 2],':k')
title('Submesoscale (k_{Rossby} < k < k_{ship}) Spectral Slope')

print([figOut 'allDeploys_SmKslopeT_SmKslopeSig'],'-dpng')
savefig([figOut 'allDeploys_SmKslopeT_SmKslopeSig.fig'])


%% Density vs. Temperature Mesoscale Spectral Variance
figure
for i = 1:numel(data.varSig_lowK)
  if ~isnan(data.varSig_lowK(i)) && strcmp(data.cruise{i},'SODA')
    scatter(data.varSig_lowK(i),data.varT_lowK(i),200,'o',...
      'MarkerEdgeColor','r','MarkerFaceColor','r')
    hold on
  elseif ~isnan(data.varSig_lowK(i)) && strcmp(data.cruise{i},'ArcticMix')
    scatter(data.varSig_lowK(i),data.varT_lowK(i),200,'d',...
      'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
  end
end
set(gca,'YScale','log')
set(gca,'XScale','log')
plot([1e-6 1e0],[1e-6 1e0],':k')
xlabel('Density [(kg/m^3)^2 km]')
ylabel('Temperature [K^2 km]')
title('Mesoscale (k < k_{Rossby}) Spectral Variance')

print([figOut 'allDeploys_LKvarSigma_LKvarT'],'-dpng')
savefig([figOut 'allDeploys_LKvarSigma_LKvarT.fig'])


%% Density vs. Temperature Submesoscale Spectral Variance
figure
for i = 1:numel(data.varSig_smK)
  if ~isnan(data.varSig_smK(i)) && strcmp(data.cruise{i},'SODA')
    scatter(data.varSig_smK(i),data.varT_smK(i),200,'o',...
      'MarkerEdgeColor','r','MarkerFaceColor','r')
    hold on
  elseif ~isnan(data.varSig_smK(i)) && strcmp(data.cruise{i},'ArcticMix')
    scatter(data.varSig_smK(i),data.varT_smK(i),200,'d',...
      'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
  end
end
set(gca,'YScale','log')
set(gca,'XScale','log')
plot([1e-7 1e-1],[1e-7 1e-1],':k')
xlabel('Density [(kg/m^3)^2 km]')
ylabel('Temperature [K^2 km]')
title('Submesoscale (k_{Rossby} < k < k_{ship}) Spectral Variance')

print([figOut 'allDeploys_SmKvarSigma_SmKvarT'],'-dpng')
savefig([figOut 'allDeploys_SmKvarSigma_SmKvarT.fig'])

end