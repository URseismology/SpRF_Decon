function [time, radRF, transRF, varR, varT] = ...
    Sp_MTCRF_ET_se_noBinning(R, T, Z, dnoi, PSI, K, Dt, taper_win, olap)

N = length(Z); % length of input/output
Ntaper = round(taper_win / Dt); % converts to samples
ofac   = 1 / (1 - olap); % overlap factor
Nwin   = round(ofac * N / Ntaper); % number of moving windows

YZ = repmat(complex(0), N, K);
YR = repmat(complex(0), N, K);
YT = repmat(complex(0), N, K);
So = repmat(complex(0), N, K);

for j = 1:Nwin  %Window of size taperWin across the trace

    %js = round(j * Ntaper / ofac);
    js = round((j-1) * Ntaper / ofac) + 1;
    je = round(js + Ntaper); % start and end indices for current window

    if je > N
        break;
    end

    for k = 1:K  %Calculate the spectrum for this windows and sum it up

        w = zeros(N, 1);
        w(js:je - 1) = PSI(:, k); % multiply the data by the k-th slepian taper

        % taper the data and do FFT
        Dat = w .* [Z, R, T, dnoi];
        fDat = fft(Dat);      %spectrum of the current chunk

        % k-th multi-taper spectral estimates and add to global total
        YZ(:,k) = YZ(:,k) + fDat(:,1);
        YR(:,k) = YR(:,k) + fDat(:,2);
        YT(:,k) = YT(:,k) + fDat(:,3);
        So(:,k) = So(:,k) + fDat(:,4);

    end
end

% cross-spectrum between components
SZT = repmat(complex(0), N, 1);
SZR = repmat(complex(0), N, 1);
SZZ = repmat(complex(0), N, 1);
SRR = repmat(complex(0), N, 1);
STT = repmat(complex(0), N, 1);
S0  = repmat(complex(0), N, 1);   % rename to S0 to avoid clash

for kk = 1:K
    SZR = SZR + conj(YZ(:,kk)) .* YR(:,kk);
    SZT = SZT + conj(YZ(:,kk)) .* YT(:,kk);
    SZZ = SZZ + conj(YZ(:,kk)) .* YZ(:,kk);
    SRR = SRR + conj(YR(:,kk)) .* YR(:,kk);
    STT = STT + conj(YT(:,kk)) .* YT(:,kk);
    S0  = S0  + conj(So(:,kk)) .* So(:,kk);   % use accumulated noise
end

% average noise power and inflate a bit
S0    = S0 / K;
alpha = 3;
S0    = alpha * S0;

% Perform deconvolution frequency domain RFs
HR = SZR ./ (SZZ + S0);
HT = SZT ./ (SZZ + S0);

% Coherence estimates
CR = abs(SZR) ./ (sqrt(SZZ) .* sqrt(SRR) );
CT = abs(SZT) ./ (sqrt(SZZ) .* sqrt(STT) );

% variance

varR = ((1 - CR.^2) ./ ( (K-1) .*  (CR.^2))) .* (abs(HR).^2) ;
varT = ((1 - CT.^2) ./ ( (K-1) .*  (CT.^2))) .* (abs(HT).^2) ;

% return frequency domain RFs

radRF = HR;
transRF = HT;
time = ((0:N-1)-floor(N/2)) * Dt;