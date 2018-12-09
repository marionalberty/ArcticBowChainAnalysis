function ArcticSpectraChi
% Load all bow chain data with chi estimates and surface horizontal
% wavenumber spectra.

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
z_comp = 5; % [m]

% Bootstrap constants for confidence interval
alpha = 0.1;
nboot = 1000;


%% Get list of files
flistA = dir(fullfile(pathinA,'*.mat'));
flistS = dir(fullfile(pathinS,'*.mat'));
flist = [flistA; flistS];
% Housekeeping
clear flistA flistS


%% Initialize data structures
data.deployment = {};   data.cruise = {};
data.dn_range = [];     data.lat = [];          data.lon = [];
data.Pow_T = {};        data.k = {};
data.k_mros = [];       data.k_mld = [];        data.k_min = [];
data.slope_lowK = [];   data.var_lowK = [];
data.slope_smK = [];    data.var_smK = [];
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
  % Extract temperature data to calc spectra for
  T = bchain.xgrid.t(bchain.xgrid.z == z_comp,:);
  x = bchain.xgrid.x;
  dn_T = bchain.xgrid.dn;
  % Extract Chi for averaging
  chi = bchain.chi.chi(bchain.chi.z == z_comp,:);
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
    
    % Calc spectra for segment
    [spec,k] = spectrum(T,x);
    data.Pow_T{end+1,1} = spec;
    data.k{end+1,1} = k;
    data.k_mros(end+1,1) = bchain.spectra.k_mros;
    data.k_mld(end+1,1) = bchain.spectra.k_mld;
    data.k_min(end+1,1) = k(1);
    
    % Estimate low wavenumber slope
    % Low wavenumber range (wavelengths greater than Rossby)
    ik = k <= bchain.spectra.k_mros;
    % Estimate slope
    fobj = fit(log(k(ik))',log(spec(ik))','poly1');
    p = coeffvalues(fobj);
    q = confint(fobj);
    data.slope_lowK(end+1,:) = [p(1) p(1)-q(1,1)];
    % Quantify spectral energy
    data.var_lowK(end+1,1) = trapz(k(ik),spec(ik));
    
    % Estimate submesoscale wavenumber slope
    if bchain.spectra.k_mros < bchain.spectra.k_mld
      % wavelengths less than Rossby and greater than MLD
      ik = k > bchain.spectra.k_mros & k <= bchain.spectra.k_mld;
      % Estimate slope
      fobj = fit(log(k(ik))',log(spec(ik))','poly1');
      p = coeffvalues(fobj);
      q = confint(fobj);
      data.slope_smK(end+1,:) = [p(1) p(1)-q(1,1)];
      % Quantify spectral energy
      data.var_smK(end+1,1) = trapz(k(ik),spec(ik));
    else
      data.slope_smK(end+1,:) = [nan nan];
      data.var_smK(end+1,1) = nan;
    end
    
    % Get mean chi
    data.chi_bar(end+1,1) =  mean(bootstrp(nboot,@mean,chi));
    CI = bootci(nboot,{@mean,chi},'type','per','alpha',alpha)';
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
      
      % Calc spectra for segment
      [spec,k] = spectrum(T(i_str:i_end),x(i_str:i_end));
      data.Pow_T{end+1,1} = spec;
      data.k{end+1,1} = k;
      data.k_mros(end+1,1) = bchain.spectra.k_mros;
      data.k_mld(end+1,1) = bchain.spectra.k_mld;
      data.k_min(end+1,1) = k(1);
      
      % Estimate low wavenumber slope
      % Low wavenumber range (wavelengths greater than Rossby)
      ik = k <= bchain.spectra.k_mros;
      % Estimate slope
      fobj = fit(log(k(ik))',log(spec(ik))','poly1');
      p = coeffvalues(fobj);
      q = confint(fobj);
      data.slope_lowK(end+1,:) = [p(1) p(1)-q(1,1)];
      % Quantify spectral energy
      data.var_lowK(end+1,1) = trapz(k(ik),spec(ik));
      
      % Estimate submesoscale wavenumber slope
      if bchain.spectra.k_mros < bchain.spectra.k_mld
        % wavelengths less than Rossby and greater than MLD
        ik = k > bchain.spectra.k_mros & k <= bchain.spectra.k_mld;
        % Estimate slope
        fobj = fit(log(k(ik))',log(spec(ik))','poly1');
        p = coeffvalues(fobj);
        q = confint(fobj);
        data.slope_smK(end+1,:) = [p(1) p(1)-q(1,1)];
        % Quantify spectral energy
        data.var_smK(end+1,1) = trapz(k(ik),spec(ik));
      else
        data.slope_smK(end+1,:) = [nan nan];
        data.var_smK(end+1,1) = nan;
      end
      
      % Get mean chi
      % Get indecies for chi in segment
      i_chi = dn_chi >= dn_T(i_str) & dn_chi <= dn_T(i_end);
      data.chi_bar(end+1,1) =  mean(bootstrp(nboot,@mean,chi(i_chi)));
      CI = bootci(nboot,{@mean,chi(i_chi)},'type','per','alpha',alpha)';
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
      
      % Calc spectra for segment
      [spec,k] = spectrum(T(i_str:i_end),x(i_str:i_end));
      data.Pow_T{end+1,1} = spec;
      data.k{end+1,1} = k;
      data.k_mros(end+1,1) = bchain.spectra.k_mros;
      data.k_mld(end+1,1) = bchain.spectra.k_mld;
      data.k_min(end+1,1) = k(1);
      
      % Estimate low wavenumber slope
      % Low wavenumber range (wavelengths greater than Rossby)
      ik = k <= bchain.spectra.k_mros;
      % Estimate slope
      fobj = fit(log(k(ik))',log(spec(ik))','poly1');
      p = coeffvalues(fobj);
      q = confint(fobj);
      data.slope_lowK(end+1,:) = [p(1) p(1)-q(1,1)];
      % Quantify spectral energy
      data.var_lowK(end+1,1) = trapz(k(ik),spec(ik));
      
      % Estimate submesoscale wavenumber slope
      if bchain.spectra.k_mros < bchain.spectra.k_mld
        % wavelengths less than Rossby and greater than MLD
        ik = k > bchain.spectra.k_mros & k <= bchain.spectra.k_mld;
        % Estimate slope
        fobj = fit(log(k(ik))',log(spec(ik))','poly1');
        p = coeffvalues(fobj);
        q = confint(fobj);
        data.slope_smK(end+1,:) = [p(1) p(1)-q(1,1)];
        % Quantify spectral energy
        data.var_smK(end+1,1) = trapz(k(ik),spec(ik));
      else
        data.slope_smK(end+1,:) = [nan nan];
        data.var_smK(end+1,1) = nan;
      end
      
      % Get mean chi
      % Get indecies for chi in segment
      i_chi = dn_chi >= dn_T(i_str) & dn_chi <= dn_T(i_end);
      data.chi_bar(end+1,1) =  mean(bootstrp(nboot,@mean,chi(i_chi)));
      CI = bootci(nboot,{@mean,chi(i_chi)},'type','per','alpha',alpha)';
      data.chi_ci(end+1,:) = [data.chi_bar(end) - CI(1) ...
        CI(2) - data.chi_bar(end)];
    end
  end
  % End deployment loop
end


%% Save data
save([dataOut 'allDeployments_SpectraChi.mat'],'data')


%% Loglog plot of 5m temp specrta colored by mean chi
% Get colorbar limits from Chi_bar
clim = [floor(log10(min(data.chi_bar))) ceil(log10(max(data.chi_bar)))];
% Use integer tick marks
cticks = clim(1):clim(2);
% Get colormap for plotting
cmap = jet(numel(cticks)-1);
% Get colorbin for each Spectra
[~,~,cbin] = histcounts(log10(data.chi_bar),cticks);

% Make figure
figure
for i = 1:numel(data.chi_bar)
  loglog(data.k{i},data.Pow_T{i},'color',cmap(cbin(i),:))
  hold on
end
colormap(cmap)
c = colorbar;
c.Ticks = (0:numel(cticks)-1)./(numel(cticks)-1);
c.TickLabels = num2str(cticks');
axis tight
xlim([1/L mean(data.k_mros)])
xlabel('Horizontal Wavenumber [1/km]')
title('Temperature Spectra [K^2 km]')
ylabel(c,'Segement Mean \chi_T log_{10}[K^2/s]')

print([figOut 'allDeploys_bchainTspectra'],'-dpng')
savefig([figOut 'allDeploys_bchainTspectra.fig'])


%% Scatter plot of chi and low wavenumber slope
% Get colorbar limits for lowK normalized variance
clim = [floor(log10(min(data.var_lowK./(data.k_mros-data.k_min)))) ...
  ceil(log10(max(data.var_lowK./(data.k_mros-data.k_min))))];
% Use integer tick marks
cticks = clim(1):clim(2);
% Get colormap for plotting
cmap = jet(numel(cticks)-1);
% Get colorbin for each Spectra
[~,~,cbin] = histcounts(log10(data.var_lowK./(data.k_mros-data.k_min)),...
  cticks);

figure
for i = 1:numel(data.var_lowK)
  if ~isnan(data.var_lowK(i)) && strcmp(data.cruise{i},'SODA')
    errorbar(data.slope_lowK(i,1),data.chi_bar(i),data.chi_ci(i,1),...
      data.chi_ci(i,2),data.slope_lowK(i,2),data.slope_lowK(i,2),'o',...
      'MarkerSize',15,'MarkerEdgeColor',cmap(cbin(i),:),...
      'MarkerFaceColor',cmap(cbin(i),:),'color',cmap(cbin(i),:))
    hold on
  elseif ~isnan(data.var_lowK(i)) && strcmp(data.cruise{i},'ArcticMix')
    errorbar(data.slope_lowK(i,1),data.chi_bar(i),data.chi_ci(i,1),...
      data.chi_ci(i,2),data.slope_lowK(i,2),data.slope_lowK(i,2),'d',...
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
title('Normalized Mesoscale Spectral Variance')
ylabel(c,'log_{10}[K^2 km]')

print([figOut 'allDeploys_LKslope_Chi'],'-dpng')
savefig([figOut 'allDeploys_LKslope_Chi.fig'])


%% Scatter plot of chi and submesoscale wavenumber slope
% Get colorbar limits for smK variance
clim = [floor(log10(min(data.var_smK./(data.k_mld-data.k_mros)))) ...
  ceil(log10(max(data.var_smK./(data.k_mld-data.k_mros))))];
% Use integer tick marks
cticks = clim(1):clim(2);
% Get colormap for plotting
cmap = jet(numel(cticks)-1);
% Get colorbin for each Spectra
[~,~,cbin] = histcounts(log10(data.var_smK./(data.k_mld-data.k_mros)),...
  cticks);

figure
for i = 1:numel(data.var_smK)
  if ~isnan(data.slope_smK(i)) && strcmp(data.cruise{i},'SODA')
    errorbar(data.slope_smK(i,1),data.chi_bar(i),data.chi_ci(i,1),...
      data.chi_ci(i,2),data.slope_smK(i,2),data.slope_smK(i,2),'o',...
      'MarkerSize',15,'MarkerEdgeColor',cmap(cbin(i),:),...
      'MarkerFaceColor',cmap(cbin(i),:),'color',cmap(cbin(i),:))
    hold on
  elseif ~isnan(data.slope_smK(i)) && strcmp(data.cruise{i},'ArcticMix')
    errorbar(data.slope_smK(i,1),data.chi_bar(i),data.chi_ci(i,1),...
      data.chi_ci(i,2),data.slope_smK(i,2),data.slope_smK(i,2),'d',...
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
title('Normalized Submesoscale Spectral Variance')
ylabel(c,'log_{10}[K^2 km]')

print([figOut 'allDeploys_SmKslope_Chi'],'-dpng')
savefig([figOut 'allDeploys_SmKslope_Chi.fig'])


%% Scatter plot of chi and submesoscale wavenumber variance normalized by
%  submesoscale wavenumber range
% Get colorbar limits for lowK variance
clim = [floor(log10(min(data.var_lowK./(data.k_mros-data.k_min)))) ...
  ceil(log10(max(data.var_lowK./(data.k_mros-data.k_min))))];
% Use integer tick marks
cticks = clim(1):clim(2);
% Get colormap for plotting
cmap = jet(numel(cticks)-1);
% Get colorbin for each Spectra
[~,~,cbin] = histcounts(log10(data.var_lowK./(data.k_mros-data.k_min)),...
  cticks);

figure
for i = 1:numel(data.var_smK)
  if ~isnan(data.var_smK(i)) && strcmp(data.cruise{i},'SODA')
    errorbar(data.var_smK(i)/(data.k_mld(i)-data.k_mros(i)),...
      data.chi_bar(i),data.chi_ci(i,1),...
      data.chi_ci(i,2),'o','MarkerSize',15,'MarkerEdgeColor',...
      cmap(cbin(i),:),'MarkerFaceColor',cmap(cbin(i),:),'color',...
      cmap(cbin(i),:))
    hold on
  elseif ~isnan(data.var_smK(i)) && strcmp(data.cruise{i},'ArcticMix')
    errorbar(data.var_smK(i)/(data.k_mld(i)-data.k_mros(i)),...
      data.chi_bar(i),data.chi_ci(i,1),...
      data.chi_ci(i,2),'d','MarkerSize',15,'MarkerEdgeColor',...
      cmap(cbin(i),:),'MarkerFaceColor',cmap(cbin(i),:),'color',...
      cmap(cbin(i),:))
    hold on
  end
end
set(gca,'YScale','log')
set(gca,'XScale','log')
colormap(cmap)
c = colorbar;
c.Ticks = (0:numel(cticks)-1)./(numel(cticks)-1);
c.TickLabels = num2str(cticks');
title('Normalized Submesoscale (k_{Rossby} < k < k_{MLD}) Variance')
xlabel('log_{10}[K^2 km]')
ylabel('\chi_T [K^2/s]')
ylabel(c,{'Normalized Mesoscale (k < k_{Rossby}) Variance';...
  'log_{10}[K^2 km]'})

print([figOut 'allDeploys_SmKvar_Chi'],'-dpng')
savefig([figOut 'allDeploys_SmKvar_Chi.fig'])


%% Submesoscale vs Mesoscale Normalized Spectral Variance
% Use integer tick marks
cticks = -3.5:-0.5;
% Get colormap for plotting
cmap = jet(numel(cticks)-1);
% Get colorbin for each Spectra
[~,~,cbin] = histcounts(data.slope_lowK(:,1),cticks);

figure
for i = 1:numel(data.var_smK)
  if ~isnan(data.var_smK(i)) && strcmp(data.cruise{i},'SODA')
    scatter(data.var_lowK(i),data.var_smK(i),200,'o',...
      'MarkerEdgeColor',cmap(cbin(i),:),'MarkerFaceColor',...
      cmap(cbin(i),:))
    hold on
  elseif ~isnan(data.var_smK(i)) && strcmp(data.cruise{i},'ArcticMix')
    scatter(data.var_lowK(i),data.var_smK(i),200,'d',...
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
title('Normalized Spectral Variance log_{10}[K^2 km]')
ylabel(c,'Mesoscale Spectral Slope')
colormap(cmap)

print([figOut 'allDeploys_SmKvar_LKvar'],'-dpng')
savefig([figOut 'allDeploys_SmKvar_LKvar.fig'])


%% Submesoscale vs. Mesoscale Spectral Slope
% Get colorbar limits for lowK normalized variance
clim = [floor(log10(min(data.var_lowK./(data.k_mros-data.k_min)))) ...
  ceil(log10(max(data.var_lowK./(data.k_mros-data.k_min))))];
% Use integer tick marks
cticks = clim(1):clim(2);
% Get colormap for plotting
cmap = jet(numel(cticks)-1);
% Get colorbin for each Spectra
[~,~,cbin] = histcounts(log10(data.var_lowK./(data.k_mros-data.k_min)),...
  cticks);

figure
for i = 1:numel(data.slope_smK(:,1))
  if ~isnan(data.slope_smK(i)) && strcmp(data.cruise{i},'SODA')
    errorbar(data.slope_lowK(i,1),data.slope_smK(i,1),...
      data.slope_smK(i,2),data.slope_smK(i,2),...
      data.slope_lowK(i,2),data.slope_lowK(i,2),'o',...
      'MarkerSize',15,'MarkerEdgeColor',cmap(cbin(i),:),...
      'MarkerFaceColor',cmap(cbin(i),:),'color',cmap(cbin(i),:))
    hold on
  elseif ~isnan(data.slope_smK(i)) && strcmp(data.cruise{i},'ArcticMix')
    errorbar(data.slope_lowK(i,1),data.slope_smK(i,1),...
      data.slope_smK(i,2),data.slope_smK(i,2),...
      data.slope_lowK(i,2),data.slope_lowK(i,2),'d',...
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
title('Spectral Slope')
ylabel(c,'Mesoscale Variance log_{10}[K^2 km]')
colormap(cmap)

print([figOut 'allDeploys_SmKslope_LKslope'],'-dpng')
savefig([figOut 'allDeploys_SmKslope_LKslope.fig'])

end