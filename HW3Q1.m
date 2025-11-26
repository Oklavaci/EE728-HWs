% =========================================================================
%  Channel Capacity C* for multiple constellations based on Ungerboeck 1982
%  Discrete-Input / Continuous-Output AWGN Channel (Soft ML)
%  Constellations Included:
%    - 4-PSK   (QPSK)
%    - 8-PSK
%    - 16-QASK
%    - 16-QAM
%    - 64-QAM
% =========================================================================

clear; close all; clc;

%% === CONSTELLATIONS (unit average power) ===
% 4-PSK (QPSK)
PSK4 = exp(1j*(0:3)*pi/2);
PSK4 = PSK4/sqrt(mean(abs(PSK4).^2));

% 8-PSK
PSK8 = exp(1j*(0:7)*2*pi/8);
PSK8 = PSK8/sqrt(mean(abs(PSK8).^2));

% 16-QASK (same as 16-QAM with unit average power)
QASK16 = qammod(0:15,16,'UnitAveragePower',true);

% 16-QAM
QAM16 = qammod(0:15,16,'UnitAveragePower',true);

% 64-QAM
QAM64 = qammod(0:63,64,'UnitAveragePower',true);

% ALL CONSTELLATIONS
consts = {PSK4, PSK8, QASK16, QAM16, QAM64};
names  = {'4-PSK','8-PSK','16-QASK','16-QAM','64-QAM'};

%% === Params ===
SNRdB = -2:2:30;
Nsamp = 1e5;                      % Monte Carlo samples
capacity = zeros(numel(consts), numel(SNRdB));

%% === MAIN ===
for m = 1:numel(consts)
    a = consts{m};
    N = numel(a);
    fprintf("Processing %s...\n", names{m});
    
    for s = 1:numel(SNRdB)
        rho = 10^(SNRdB(s)/10);   
        sigma2 = 1/(2*rho);       % average symbol power = 1
        
        % AWGN noise samples
        w = sqrt(sigma2) * (randn(Nsamp,1) + 1j*randn(Nsamp,1));
        
        Csum = 0;
        for k = 1:N
            ak = a(k);
            inner = zeros(Nsamp,1);
            
            % Summation over constellation
            for j = 1:N
                aj = a(j);
                inner = inner + exp(-(abs(ak + w - aj).^2 - abs(w).^2)/sigma2);
            end
            
            Csum = Csum + mean(log2(inner));
        end
        
        % Ungerboeck Eq. (5)
        capacity(m,s) = log2(N) - (1/N)*Csum;
    end
end

%% === Plot ===
figure; hold on; grid on;

for m = 1:numel(consts)
    plot(SNRdB, capacity(m,:), 'LineWidth', 2);
end

xlabel('SNR (dB)');
ylabel('Capacity C^* (bits/symbol)');
legend(names, 'Location','southeast');

%% ===========================================================
%  Compute BER crossing SNR values
% ============================================================

targetBER = 1e-3;
SNRtest = 0:0.1:40;
SNRlin  = 10.^(SNRtest/10);

% BER formulas
BER_QPSK  = qfunc(sqrt(2*SNRlin));
BER_16QAM = (3/4)*qfunc( sqrt((4/5)*SNRlin) );
BER_64QAM = (7/12)*qfunc( sqrt((12/63)*SNRlin*2) );

%% -------- Helper function to safely find SNR crossing -------------
findSNR = @(BERcurve) ...
    interp1( ...
        unique(flip(BERcurve)), ...            % remove duplicates
        flip(SNRtest(1:length(unique(BERcurve)))), ...
        targetBER, 'linear', 'extrap');

%% Compute thresholds
SNR_QPSK  = findSNR(BER_QPSK);
SNR_16QAM = findSNR(BER_16QAM);
SNR_64QAM = findSNR(BER_64QAM);

fprintf('Computed BER=1e-3 thresholds:\n');
fprintf('  QPSK   : %.2f dB\n', SNR_QPSK);
fprintf('  16-QAM : %.2f dB\n', SNR_16QAM);
fprintf('  64-QAM : %.2f dB\n', SNR_64QAM);

%% ============================================================
%  Add BER markers on the capacity plot
% ============================================================

Cap_QPSK  = interp1(SNRdB, capacity(strcmp(names,'4-PSK'),:), SNR_QPSK);
Cap_16QAM = interp1(SNRdB, capacity(strcmp(names,'16-QAM'),:), SNR_16QAM);
Cap_64QAM = interp1(SNRdB, capacity(strcmp(names,'64-QAM'),:), SNR_64QAM);

plot(SNR_QPSK,  Cap_QPSK,  'ko', 'MarkerFaceColor','k', 'MarkerSize',8);
plot(SNR_16QAM, Cap_16QAM, 'ko', 'MarkerFaceColor','k', 'MarkerSize',8);
plot(SNR_64QAM, Cap_64QAM, 'ko', 'MarkerFaceColor','k', 'MarkerSize',8);

text(SNR_QPSK+0.3,  Cap_QPSK,  'BER=10^{-3}', 'FontSize',8);
text(SNR_16QAM+0.3, Cap_16QAM, 'BER=10^{-3}', 'FontSize',8);
text(SNR_64QAM+0.3, Cap_64QAM, 'BER=10^{-3}', 'FontSize',8);

%% ============================================================
%  SNR where 16-QAM supports 2 bps ***
% ============================================================

targetRate = 2;  % bits per symbol
%% === Find SNR where 16-QAM capacity = 2 bps ===
Cap16 = capacity(strcmp(names,'16-QAM'),:);
uniqueCap = unique(Cap16);                          % strictly increasing
SNR_for_uniqueCap = SNRdB(ismember(Cap16, uniqueCap));

% ensure vectors have same length (sometimes duplicates remove mismatch)
L = min(length(uniqueCap), length(SNR_for_uniqueCap));
uniqueCap = uniqueCap(1:L);
SNR_for_uniqueCap = SNR_for_uniqueCap(1:L);

SNR_16QAM_rate2 = interp1(uniqueCap, SNR_for_uniqueCap, targetRate, 'linear');


fprintf('SNR(16-QAM achieving C=2 bps) = %.2f dB\n', SNR_16QAM_rate2);

%% Add marker for the "Rate = 2 bps" point
plot(SNR_16QAM_rate2, targetRate, 'rs', 'MarkerFaceColor','r', 'MarkerSize',8);
text(SNR_16QAM_rate2+0.3, targetRate, 'C=2 bps', 'Color','r', 'FontSize',8);

%% ============================================================
%  Compute Coding Gain (Required by HW)
% ============================================================

coding_gain = SNR_QPSK - SNR_16QAM_rate2;

fprintf('Coding Gain = SNR(QPSK @ BER=1e-3) - SNR(16-QAM @ C=2) = %.2f dB\n', coding_gain);


