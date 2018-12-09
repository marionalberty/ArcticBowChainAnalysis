function [epsi,chi,U_S,time_S] = sbe56EpsiChi(T,U,Tz,N2,w,time)
% Following Zhang and Moum 2010, Section 3b.
% 1. Computes nsegs spectra of sections of timeseries Ns elements long
% 2. Corrects for roll off using transfer function calculated by Byungho
% 3. Least squares fits to eq. 8 in Z&M over range defined by l0,l1
% 4. Solves for epsilon
% T, U, Tz, N2 and time are vectors of the same length, at 1Hz
% (might need to up-interpolate for U).
% T: temperature in C
% U: flow speed past sensor in m/s
% Tz: dTdz in C/m
% N2: in (rad/s)^2
% w: vertical velocity of the sensor in m/s
% time: in days (datenum is fine)

% RCM Jul 2014
% Edits by MSA August 2018


%% Set constants
% no of seconds over which to fit
Nfft = 300;
% Obukhov-Corrsin constant
Ct = 0.4;
% Mixing efficiency
Gam = 0.2;
% Transfer function constant
fc = 0.255;

%% Reduce length of vectors to maximize N*Nfft
npts = floor(numel(T)/Nfft)*Nfft;
T = T(1:npts);
U = U(1:npts);
Tz = Tz(1:npts);
N2 = N2(1:npts);
w = w(1:npts);
time = time(1:npts);


%% Calculate measured dTdx
% Convert time to seconds
times = (time - time(1))*86400;     % in [s] since start of record
dTdx_M = gradient(T,times)./U;
% Make half overlapping segments of dTdx
dTdx_MS = reshape(dTdx_M,Nfft,[]);


%% Calculate wave induced dTdx
% Calculate dTdx_wave
dTdx_W = w.*Tz./U;
% Make half overlapping segments of dTdx
dTdx_WS = reshape(dTdx_W,Nfft,[]);


%% Compute spectra
% Get number of segments
n_seg = numel(dTdx_MS(1,:));
% Prelocate spectra for speed
spec_M = nan(Nfft/2,n_seg);
spec_W = spec_M;
Xspec = spec_M;
gamma = spec_M;
% Calculate spectra
for i = 1:n_seg
  % Compute spectra for measured dTdx
  [spec_M(:,i),f] = spectrum(dTdx_MS(:,i),times(1:Nfft));
  % Compute spectra for wave induced dTdx
  [spec_W(:,i),~] = spectrum(dTdx_WS(:,i),times(1:Nfft));
  % Computre cross-spectra of the two signals
  [Xspec(:,i),~] = cross_spectrum(dTdx_MS(:,i),dTdx_WS(:,i),times(1:Nfft));
end
% Evaluate transfer function
tf = repmat(1./(1+(f'/fc).^2),1,n_seg);
% Correct spectra with transfer function for high wave number attenuation
spec_M = spec_M./tf;
spec_W = spec_W./tf;
Xspec = Xspec./tf;

% Calculate gamma with 2-hr spectral averages
% Get time ranges to average over
t_gam = unique([0:7200:times(end) times(end)]);
% Mean elapsed time in secs
times_S = nanmean(reshape(times,Nfft,[]));
for i = 1:numel(t_gam)-1
  % Get indecies of spectra within 2hr window
  i_bar = find(times_S > t_gam(i) & times_S < t_gam(i+1));
  % Calculate gamma
  gamma(:,i_bar) = ((abs(nanmean(Xspec(:,i_bar),2)).^2)./...
    (nanmean(spec_M(:,i_bar),2).*nanmean(spec_W(:,i_bar),2)))*...
    ones(1,numel(i_bar));
end
%   gamma = (abs(Xspec).^2)./(spec_M.*spec_W);
% Calculate environmental spectra
spec_E = spec_M.*(1-gamma);


%% Least squares fit to f^(1/3) spectrum
% Get the segment means for dTdz, N2, U and time
Tz_S = nanmean(reshape(Tz,Nfft,[]));
N2_S = nanmean(reshape(N2,Nfft,[]));
U_S = nanmean(reshape(U,Nfft,[]));
time_S = reshape(time,Nfft,[]);
time_S = [time_S(1,:); time_S(end,:)];
% Set statically unstable segments to nan
N2_S(N2_S < 0) = NaN;
N = sqrt(N2_S);
% fit range: somewhat arbitrarily chosen
% [~,l0] = min(abs(f-3e-2));
% [~,l1] = min(abs(f-2e-1));
[~,l0] = min(abs(f-5e-2));
[~,l1] = min(abs(f-5e-1));


%% Calculate Epsilon and Chi
% Prelocate epsilon and chi
epsi = NaN(1,n_seg);
chi = NaN(1,n_seg);
for jj = 1:n_seg
  %% 1 param fit for m in phi = m*f^{1/3}
  % LS fit minimizing normally distributed error
  % G = f(l0:l1).^(1/3);
  % D = sp_Tm(l0:l1,jj)';
  % lhs = 2*G*G';
  % rhs = 2*G*D';
  % m = lhs\rhs;
  % fi = m*G;
  
  % LS fit minimizing log-normally distributed error
  G = ones(1,l1-l0+1);
  D = log10(spec_E(l0:l1,jj)'./(f(l0:l1).^(1/3)));
  lhs = 2*G*G';
  rhs = 2*G*D';
  m = 10^(lhs\rhs);
  %m = mean(sp_Tm(l0:l1,jj)'./(f(l0:l1).^(1/3))); % this is the same
  
  
  %% 2 param fit for m2 and p2 in phi = m2*f^p2
  % G = [log10(f(l0:l1)); ones(1,l1-l0+1)];
  % D = log10(sp_Tm(l0:l1,jj))';
  % lhs = 2*G*G';
  % rhs = 2*G*D';
  % a(jj,:) = lhs\rhs;
  
  epsi(jj) = ((U_S(jj)^2)*(N(jj)^3)*(m^(3/2)))/...
    (((2*pi)^2)*((2*Ct*Gam*(Tz_S(jj)^2))^(3/2)));
  
  chi(jj) = ((U_S(jj)^2)*N(jj)*(m^(3/2)))/...
    (((2*pi)^2)*(Ct^(3/2))*((2*Gam*(Tz_S(jj)^2))^(1/2)));
end
% Ensure epsilon and chi are all real
epsi = real(epsi);
chi = real(chi);
%  and positive
epsi(epsi < 0) = nan;
chi(chi < 0) = nan;

% NAN eps/chi computed over segments where u bar was less than 0.25 m/s
epsi(U_S <= 0.25) = nan;
chi(U_S <= 0.25) = nan;
end

function [spec,f] = spectrum(qq,x1)
% function [Spec,f] = spectrum(qq,x1)
% computes the 1-dimensional spectrum of the given field, qq, with a
% resolution of dx1. Spectra are detrended, pre-whitened with a hanning
% window and normalized to preserve variance.

%% Detrend qq
qq = qq(:);   x1 = x1(:);
qq = detrend(qq);

%% Window qq
N = numel(qq);
win = hann(N)*sqrt(N./sum(hann(N).^2));
qq = qq.*win;

%% Calculate frequencies/wavenumbers
fs = 1/mean(diff(x1));
f = linspace(fs/N,fs/2,N/2);

%% Compute spectra
QQ = fft(qq);
phi = 2*QQ(2:N/2+1).*conj(QQ(2:N/2+1))/N;
spec = phi'./fs;
end

function [spec,f] = cross_spectrum(q1,q2,x)
% function [Spec,f] = spectrum(q1,q2,x)
% computes the 1-dimensional cross spectrum of the given fields, q1 and q2,
% which have a resolution of dx. The vectors are detrened, pre-whitened
% with a hanning window, and normalized to preserve variance.

%% Detrend q1
q1 = q1(:);   x = x(:);
q1 = detrend(q1);

%% Detrend q2
q2 = q2(:);
q2 = detrend(q2);

%% Window qq
N = numel(q1);
win = hann(N)*sqrt(N./sum(hann(N).^2));
q1 = q1.*win;
q2 = q2.*win;

%% Calculate frequencies/wavenumbers
fs = 1/mean(diff(x));
f = linspace(fs/N,fs/2,N/2);

%% Compute spectra
Q1 = fft(q1);
Q2 = fft(q2);
phi = 2*Q2(2:N/2+1).*conj(Q1(2:N/2+1))/N;
spec = phi'./fs;
end