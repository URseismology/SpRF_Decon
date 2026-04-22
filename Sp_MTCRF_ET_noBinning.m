function [time, radRF, transRF] = ...
    Sp_MTCRF_ET_noBinning(RR, TT, ZZ, NN, Dt, Fcut, taperWin, olap, p, k)

% --- STEP 1: AUTO-PADDING ---
%pad to next power of 2 that gives at least 800s total duration.
[N_orig, nwaves] = size(RR);
target_time = 800; % Seconds (Adjust if needed)
N_pad = 2^nextpow2(target_time / Dt); 

if N_pad > N_orig
    RR(N_pad, 1) = 0; 
    TT(N_pad, 1) = 0; 
    ZZ(N_pad, 1) = 0; 
    NN(N_pad, 1) = 0; 
    N = N_pad; 
else
    N = N_orig;
end

% Initialize frequency domain matrices
radRF_freq = zeros(N, nwaves);
transRF_freq = zeros(N, nwaves);

% Pre-define slepian taper
Ntaper = round(taperWin / Dt);
[PSI, ~] = sleptap(Ntaper, p, k); %generates the slepian tapers

% --- Loop 1: Calculate RFs (Frequency Domain) ---
for iwave = 1: nwaves
    R = RR(:, iwave);
    T = TT(:, iwave);
    Z = ZZ(:, iwave);
    dnoi = NN(:, iwave);
    
    % Call the ET core function 
    [~, radRF_freq(:, iwave), transRF_freq(:,iwave), ~, ~] = ...
        Sp_MTCRF_ET_se_noBinning(R, T, Z, dnoi, PSI, k, Dt, taperWin, olap);
end

% --- Loop 2: Transform to Time Domain ---
radRF = zeros(N, nwaves);
transRF = zeros(N, nwaves);
fmax = 1/(2.0*Dt);

for iwave = 1: nwaves
    HRin = radRF_freq(:, iwave);
    HTin = transRF_freq(:, iwave);
    
    % Apply cosine taper (low-pass filter)
    [HR, ~] = costaper(HRin, Fcut, Dt);
    [HT, ~] = costaper(HTin, Fcut, Dt);
    
    % Inverse FFT to time domain
    radRF(:, iwave) = (2*fmax/Fcut) .* fftshift(real(ifft(HR))); 
    transRF(:, iwave) = (2*fmax/Fcut) .* fftshift(real(ifft(HT)));
end

% --- STEP 3: Time Vector ---
time = ((0:N-1)-floor(N/2)) * Dt;

end
