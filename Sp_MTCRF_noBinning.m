function [time, radRF, transRF, varR, varT] = ...
    Sp_MTCRF_noBinning(RR, TT, ZZ, NN, Dt, Fcut)

% Get dimensions
[N, nwaves] = size(RR);

% Initialize output matrices for frequency-domain RFs and variance
radRF_freq = zeros(N, nwaves);
transRF_freq = zeros(N, nwaves);
varR = zeros(N, nwaves); % Returning freq-domain variance just in case
varT = zeros(N, nwaves);

% Pre-define slepian taper (Whole-trace tapering)
P = 2.5;
K = 2;
[PSI,~] = sleptap(N, P, K);

% --- Loop 1: Calculate RFs in the Frequency Domain ---
for iwave = 1:nwaves
    
    R = RR(:, iwave);
    T = TT(:, iwave);
    Z = ZZ(:, iwave);
    dnoi = NN(:, iwave);
    
    % Call the standard MT core function (NOT the ET version)
    [~, radRF_freq(:, iwave), transRF_freq(:,iwave), varR(:, iwave), varT(:,iwave)] = ...
        Sp_MTCRF_se(R, T, Z, dnoi, PSI, K, Dt, Fcut);
end

% --- Loop 2: Transform back to Time Domain ---
radRF = zeros(N, nwaves);
transRF = zeros(N, nwaves);
fmax = 1/(2.0*Dt);

for iwave = 1:nwaves
    
    HRin = radRF_freq(:, iwave);
    HTin = transRF_freq(:, iwave);
    
    % Apply cosine taper (low-pass filter)
    [HR, ~] = costaper(HRin, Fcut, Dt);
    [HT, ~] = costaper(HTin, Fcut, Dt);
    
    % Inverse FFT to time domain (shifts -tlag to +tlag)
    radRF(:, iwave) = (2*fmax/Fcut) .* fftshift(real(ifft(HR)));  
    transRF(:, iwave) = (2*fmax/Fcut) .* fftshift(real(ifft(HT)));
    
end

% --- Step 3: Define Time Vector ---
time = ((0:N-1)-floor(N/2)) * Dt;
time = time(:); % Ensure column vector

end