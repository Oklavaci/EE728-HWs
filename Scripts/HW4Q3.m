%% Ergodic Capacity vs SNR for 3x5 and 5x3 MIMO
% With and without CSIT (spatial water-filling)

clear; clc;

SNRdB = -10:2:30;                 % SNR range in dB
SNRlin = 10.^(SNRdB/10);          % linear SNR (rho = P/N0)
Nmc    = 3000;                    % Monte Carlo realizations

% Two systems: [Mt, Mr]
systems = [3 5;   % 3x5
           5 3];  % 5x3

C_noCSIT  = zeros(length(systems), length(SNRdB));
C_CSIT    = zeros(length(systems), length(SNRdB));

for s = 1:size(systems,1)
    Mt = systems(s,1);
    Mr = systems(s,2);
    M  = min(Mt, Mr);            % number of spatial modes

    fprintf('Simulating %dx%d MIMO...\n', Mt, Mr);

    for iSNR = 1:length(SNRlin)
        rho = SNRlin(iSNR);      % total transmit SNR (P/N0)

        C_no  = zeros(Nmc,1);
        C_yes = zeros(Nmc,1);

        for k = 1:Nmc
            % i.i.d. Rayleigh channel
            H = (randn(Mr,Mt) + 1j*randn(Mr,Mt))/sqrt(2);

            % ----- No CSIT: equal power across Tx antennas -----
            % Rx knows H, Tx uses Rx = (rho/Mt) * I
            C_no(k) = log2( det( eye(Mr) + (rho/Mt) * (H*H') ) );

            % ----- With CSIT: water-filling in space -----
            % Eigenvalues of H H^H
            HH = H*H';
            lambda = eig(HH);              % Mx1 vector (includes zero if rank-deficient)
            lambda = sort(real(lambda),'descend'); % sort descending
            lambda = lambda(1:M);          % only largest min(Mt,Mr) modes

            % Water-filling over spatial modes
            p = waterfill(lambda, rho);    % power allocation over M modes

            % Capacity for this realization
            C_yes(k) = sum( log2( 1 + p(:)'.*lambda(:)' ) );
        end

        % Ergodic capacities (bits/s/Hz)
        C_noCSIT(s,iSNR) = mean(C_no);
        C_CSIT(s,iSNR)   = mean(C_yes);
    end
end

%% Plot results
figure; hold on; grid on;
colors = lines(4);

% 3x5 curves
plot(SNRdB, C_noCSIT(1,:), '-',  'LineWidth', 1.8, 'Color', colors(1,:));
plot(SNRdB, C_CSIT(1,:),   '--', 'LineWidth', 1.8, 'Color', colors(1,:));

% 5x3 curves
plot(SNRdB, C_noCSIT(2,:), '-',  'LineWidth', 1.8, 'Color', colors(2,:));
plot(SNRdB, C_CSIT(2,:),   '--', 'LineWidth', 1.8, 'Color', colors(2,:));

xlabel('SNR (dB)');
ylabel('Ergodic Capacity (bits/s/Hz)');
legend('3x5 no CSIT', '3x5 with CSIT', ...
       '5x3 no CSIT', '5x3 with CSIT', ...
       'Location', 'NorthWest');
title('Ergodic Capacity of 3x5 and 5x3 MIMO (Rayleigh, B = 1 Hz)');

%% ---- Water-filling helper function ----
function p = waterfill(lambda, P)
% lambda : vector of eigenvalues (descending)
% P      : total power
    lambda = lambda(:);
    M = length(lambda);

    % Sort descending just in case
    [lambda_sorted, idx] = sort(lambda, 'descend');
    inv_lam = 1 ./ lambda_sorted;

    % Find active set size k and water level mu
    for k = 1:M
        mu = (P + sum(inv_lam(1:k))) / k;
        if k == M
            active = k;
            break;
        end
        % Stop when water level is not above next 1/lambda
        if mu <= inv_lam(k+1)
            active = k;
            break;
        end
    end

    % Final water level using active modes
    mu = (P + sum(inv_lam(1:active))) / active;
    p_sorted = max(mu - inv_lam, 0);      % allocated powers (sorted order)

    % Map back to original eigenvalue order
    p = zeros(M,1);
    p(idx) = p_sorted;
end
