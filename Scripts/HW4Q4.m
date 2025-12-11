clear; clc; close all;

% Parameters
SNR_dB = 0:0.1:30;
SNR = 10.^(SNR_dB/10);
K = 3;                             % Rice factor
AoA = pi/12;                    % theta
AoD   = pi/6;                     % phi

% Steering vector generator
steer = @(N,ang) exp(1j*(0:N-1)*ang).';

% --- Case 1: 1x1 MIMO (SISO Ricean) ---
N_r = 1; N_t = 1;
s_theta = steer(N_r,AoA);
s_phi   = steer(N_t,AoD);

C1 = zeros(size(SNR));
for k = 1:length(SNR)
    H_los  = sqrt(K/(K+1)) * (s_theta * s_phi');
    H_nlos = sqrt(1/(K+1)) * (randn(N_r,N_t)+1j*randn(N_r,N_t))/sqrt(2);
    H = H_los + H_nlos;
    C1(k) = log2(det( eye(N_r) + SNR(k)/N_t * (H*H') ));
end

% --- Case 2: 3x5 MIMO ---
N_r = 3; N_t = 5;
s_theta = steer(N_r,AoA);
s_phi   = steer(N_t,AoD);

C2 = zeros(size(SNR));
for k = 1:length(SNR)
    H_los  = sqrt(K/(K+1)) * (s_theta * s_phi');
    H_nlos = sqrt(1/(K+1)) * (randn(N_r,N_t)+1j*randn(N_r,N_t))/sqrt(2);
    H = H_los + H_nlos;
    C2(k) = log2(det( eye(N_r) + SNR(k)/N_t * (H*H') ));
end

% Plot
figure; hold on; grid on;
plot(SNR_dB, C1,'LineWidth',2);
plot(SNR_dB, C2,'LineWidth',2);
xlabel('SNR (dB)'); ylabel('Capacity (bits/s/Hz)');
%title('Capacity vs SNR for Ricean Channels (\theta=\pi/12, \phi=\pi/6)');
legend('1x1','3x5');

%% part b

angles = [pi/15 pi/7;   % θ1, φ1
          pi/12 pi/6;   % θ2, φ2
          pi/6  pi/3];  % θ3, φ3

figure; hold on; grid on;
for a = 1:size(angles,1)
    AoA = angles(a,1);
    AoD   = angles(a,2);

    s_theta = steer(3,AoA);
    s_phi   = steer(5,AoD);

    C = zeros(size(SNR));
    for k = 1:length(SNR)
        H_los  = sqrt(K/(K+1)) * (s_theta * s_phi');
        H_nlos = sqrt(1/(K+1)) * (randn(3,5)+1j*randn(3,5))/sqrt(2);
        H = H_los + H_nlos;
        C(k) = log2(det( eye(3) + SNR(k)/5 * (H*H') ));
    end

    plot(SNR_dB,C,'LineWidth',2);
end

xlabel('SNR (dB)'), ylabel('Capacity (bits/s/Hz)')
%title('3x5 Capacity for Different AoA/AoD')
legend('\theta=\pi/15,\phi=\pi/7', '\theta=\pi/12,\phi=\pi/6', '\theta=\pi/6,\phi=\pi/3');
