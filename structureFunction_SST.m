function [D,r,error] = structureFunction_SST(SST,lat,lon,r_res,r_max)
% Calculate the structure function, D, of the field, SST, as a function of
% the horizontal separation of observations, r, using the lat,lon position
% of each observation.
%
% SST should be in Kelvin or Celcius
%
% D(r) = (2/N(r)) sum_{i=1} ^{N(r)} ([SST(x+r) - SST(x)]^2)
%
% Ndifs = (Npts-1)*Npts/2;


%% Make sure vectors are columns
SST = SST(:);
lat = lat(:);
lon = lon(:);


%% Prelocate Variables
% Make bin edges
dr_edges = 0:r_res:r_max;
% Calculate r
r = dr_edges(1:end-1) + r_res/2;
% Prelocate sum_dSST2
sum_dSST2 = 0*r;
% Prelocate total bin count N(r)
N_tot = 0*r;
% Get number of observations
Npts = numel(SST);
% Max segment length
lseg = 1e7;
% Turn while loop on until all differences have been calculated
keepRunning = 1;
% Start i at 1
istrt = 1;
% Start k at 2
kstrt = istrt + 1;


while keepRunning
  %% Initialize Variables
  % Initialize count of differences calculated
  i_dif = 1;
  % Prelocate dSST^2 and dr
  dSST2 = nan(1,lseg);
  dr = nan(1,lseg);
  
  
  %% Calculate dSST and dr
  % Must use loops because some vectors are too large to do this with
  % symmetric matricies
  % Loop through calculation
  for i = istrt:Npts-1
    if i_dif <= lseg
      kstrt = max([i+1 kstrt]);
      for k = kstrt:Npts
        if i_dif <= lseg
          dr_temp = m_lldist([lon(i) lon(k)],[lat(i) lat(k)]);
          if dr_temp <= r_max
            % Find the squared difference of SST
            dSST2(i_dif) = (SST(i) - SST(k))^2;
            % Find the distance between the observations
            dr(i_dif) = dr_temp;
            % Add to difference count
            i_dif = i_dif +1;
          end
        else
          % Break out of k
          break
        end
      end
    else
      % Break out of i
      break
    end
    % Reset kstart
    kstrt = i+1;
  end
  % Shorten dSST2 and dr to remove nans
  dr(isnan(dSST2)) = [];
  dSST2(isnan(dSST2)) = [];
  
  
  %% Make D(r)
  % Get distribution of distances between points
  [N,~,bin] = histcounts(dr,dr_edges);
  % Sum up N's
  N_tot = N_tot + N;
  % Sum up dSST^2
  for ir = 1:numel(r)
    sum_dSST2(ir) = sum_dSST2(ir) + sum(dSST2(bin == ir));
  end
  
  
  %% Reset Variables and update counts
  % Update starting count for i and k
  istrt = i;
  kstrt = k;
  
  % Break keeprunning while loop if all differences have been calculated
  if i == Npts-1 && k == Npts
    keepRunning = 0;
  end
  % End of keepRunning while loop
end

% Calculate D(r)
D = 2*sum_dSST2./N_tot;
% Set D(r) equal to nan when no points were available with r separation
D(N_tot == 0) = nan;
% Get errorbars
error = [N_tot./chi2inv(.05/2,N_tot); N_tot./chi2inv(1-.05/2,N_tot)];
end