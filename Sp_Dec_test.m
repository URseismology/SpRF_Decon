%Load data
data_all = load('S_waveforms.mat'); 
data = data_all.plot_data_matrix_R;            
time = data_all.reduced_time_matrix;
epidist = data_all.epicentral_distances;
phases = data_all.theoretical_phases;

rayp = zeros(length(epidist), 1);
for i = 1:length(epidist)
    % direction = epiDist -> rayP
    [~, rp_val] = raypToEpiDist(epidist(i), true, false, '.'); 
    rayp(i) = rp_val; 
end
rayp = rayp / 111.19;

%Plot raw data
figure(1); clf;
hold on
plot(time', data', 'k', 'LineWidth', 0.6);

hS   = plot(phases.S_times, epidist, 'r--', 'LineWidth', 1.5, 'DisplayName', 'S');
hScS = plot(phases.ScS_times, epidist, 'b--', 'LineWidth', 1.5, 'DisplayName', 'ScS');
hSKS = plot(phases.SKS_times, epidist, 'g--', 'LineWidth', 1.5, 'DisplayName', 'SKS');
hSmp  = plot(phases.Sp_times, epidist, 'm--', 'LineWidth', 1.5, 'DisplayName', 'Sp');
hSLP = plot(phases.SLP_times, epidist, 'y--', 'LineWidth', 1.5, 'DisplayName', 'SLP');
h410 = plot(phases.S410P_times, epidist, '--', 'Color', [1 0.5 0], 'LineWidth', 1.5, 'DisplayName', 'S410P'); 
h660 = plot(phases.S660P_times, epidist, 'c--', 'LineWidth', 1.5, 'DisplayName', 'S660P');

xlim([400 700]);
ylim([min(epidist) max(epidist)]); 
xlabel('Time - 10 * \Delta (s)');
ylabel('Epicentral Distance (\circ)');
title('Raw AxiSEM S-Wave Data (Radial Component)');
legend([hS, hScS, hSKS, hSmp, hSLP, h410, h660], 'Location', 'northeast');
hold off
%% -----------------------------------
% 3. PREPARE data for deconvolution
% ------------------------------------
RR = data_all.raw_data_matrix_R.';  
ZZ = data_all.raw_data_matrix_Z.';  
TT = zeros(size(RR));

dt  = data_all.original_time(2) - data_all.original_time(1);
F_cut   = 0.4;

% ETMT Setup paramters
% ------------------------------
taperWin = 60;   % Window size 
olap     = 0.75; % Overlap                    
p        = 3;    % Time-bandwidth product                                                                                     
k        = 4;    % Number of tapers  

%%
%Rotate data
pvel = 7.5;  
svel = 4.687;

% fprintf('Rotating Z and R into upgoing P and S...\n');
[PP, SS, ~] = psvrotate(ZZ, RR, TT, rayp, pvel, svel);
%[PP, SS, ~] = lqtrotate(ZZ, RR, TT, rayp, svel);

%% ========================
% --- ASYMMETRIC WINDOWING 
% =========================
fprintf('Applying asymmetric time windows...\n');

% 1. DEFINE WINDOW PARAMETERS (Relative to S-arrival, in seconds)
% P-Wave (Numerator) Window: Captures the deep precursors
p_win_pre  = 100;  % Seconds BEFORE S-wave
p_win_post = 20;   % Seconds AFTER S-wave

% S-Wave (Denominator) Window: Isolates the pure source pulse
s_win_pre  = 7;   % Seconds BEFORE S-wave
s_win_post = 7;   % Seconds AFTER S-wave

PP_win = zeros(size(PP));
SS_win = zeros(size(SS));
T_matrix = data_all.reduced_time_matrix.'; % Time x Traces [Size: nn x nwaves]

for i = 1:size(PP, 2)
    t_S = phases.S_times(i);
    
    p_mask = (T_matrix(:, i) >= t_S - p_win_pre) & (T_matrix(:, i) <= t_S + p_win_post);
    s_mask = (T_matrix(:, i) >= t_S - s_win_pre) & (T_matrix(:, i) <= t_S + s_win_post);
    
    % Apply masks (data outside these windows remains exactly zero)
    PP_win(p_mask, i) = PP(p_mask, i);
    SS_win(s_mask, i) = SS(s_mask, i);
end

%% ========================
% --- NOISE EXTRACTION 
% =========================
% Extract background noise from the UN-WINDOWED rotated SS matrix
% T_red = data_all.reduced_time_matrix.'; 
% noise_window_2D = (T_red > 350) & (T_red < 400);
% 
% NN = zeros(size(SS));
% NN(noise_window_2D) = SS(noise_window_2D); 
% 
% % Add a tiny computational noise floor to prevent divide-by-zero in synthetics
% NN = NN + (randn(size(SS)) * 1e-12);

noise_level = 1e-12; 
NN = randn(size(SS)) * noise_level;

%%
% % Run ETMT deconvolution
% fprintf('Running ETMT Kernel for Fcut = %.2f Hz...\n', F_cut); 
% [time_rf, rf_result, ~] = Sp_MTCRF_ET_noBinning(PP_win, TT, SS_win, NN,dt, F_cut, taperWin, olap, p, k);

%Run standard deconvolution
fprintf('Running Standard MT Kernel on Rotated Data (Fcut = %.2f Hz)...\n', F_cut);
[time_rf, rf_result, ~] = Sp_MTCRF_noBinning(PP_win, TT, SS_win, NN, dt, F_cut);

radRF = rf_result * -1; % Flip polarity so Moho is a positive peak
time_rf = -time_rf(:);  % Convert to positive precursor time

%%
% --- 3. PLOT DECONVOLVED DATA ---
figure(3); clf; hold on;
rf_scale = 1.2; 
epidist = data_all.epicentral_distances; % Use the original distance vector

for i = 1:length(epidist)
    x = time_rf';             
    y_raw = radRF(:, i)';     
    
    norm_val = max(abs(y_raw));
    
    % Safety check: avoid division by zero if a trace is empty
    if norm_val == 0, norm_val = 1; end
    
    % Apply scaling and shift to the baseline epicentral distance
    y = (y_raw / norm_val * rf_scale) + epidist(i);
    baseline = epidist(i);
    
    % --- Shade POSITIVE peaks (Black) ---
    y_pos = max(y, baseline);
    patch([x, fliplr(x)], [y_pos, baseline*ones(1,length(x))], 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    
    % --- Shade NEGATIVE troughs (Grey) ---
    % y_neg = min(y, baseline);
    % patch([x, fliplr(x)], [y_neg, baseline*ones(1,length(x))], [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    
    % Plot the wiggle line
    plot(x, y, 'k', 'LineWidth', 0.05);
end

% --- 4. PLOT THEORETICAL PHASES ---
% Precursor Time = (S_time - Phase_time)
t_S_ref = phases.S_times;
hS   = plot(t_S_ref - phases.S_times, epidist, 'r--', 'LineWidth', 1.5, 'DisplayName', 'S (t=0)');
hSmp  = plot(t_S_ref - phases.Sp_times, epidist, 'm--', 'LineWidth', 1.5, 'DisplayName', 'Moho (Sp)');
hSLP = plot(t_S_ref - phases.SLP_times, epidist, 'y--', 'LineWidth', 1.5, 'DisplayName', 'LAB');
h410 = plot(t_S_ref - phases.S410P_times, epidist, '--', 'Color', [1 0.5 0], 'LineWidth', 1, 'DisplayName', '410'); 
h660 = plot(t_S_ref - phases.S660P_times, epidist, 'c--', 'LineWidth', 1, 'DisplayName', '660');

xlim([-20 100]); 
ylim([min(epidist)-1 105]); 
%set(gca, 'XDir', 'reverse'); % S-wave on right, deep conversions on left

xlabel('Precursor Time (s)');
ylabel('Epicentral Distance (\circ)');
title('Final S-to-P Receiver Functions (Full Dataset)');
legend([hS, hSmp, hSLP, h410, h660], 'Location', 'northeast');
grid on; box on; hold off;
