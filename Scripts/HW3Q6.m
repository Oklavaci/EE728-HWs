%% PARAMETERS
N = 8;                         % OFDM size
h = [0.7 0.5 0 0.3];           % FIR channel taps (length 4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (a) NON-CIRCULANT LINEAR CONVOLUTION MATRIX (DMT matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = length(h);
H_linear = zeros(N);

for row = 1:N
    for col = 1:N
        if (row-col+1 >= 1) && (row-col+1 <= L)
            H_linear(row, col) = h(row-col+1);
        end
    end
end

disp('--- (a) Linear Convolution Matrix H (non-circulant) ---')
H_linear


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (b) CIRCULANT CONVOLUTION MATRIX + EIGEN-DECOMPOSITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Build circulant matrix by wrapping past the edges
h_circ = [h zeros(1, N-L)];
H_circ = gallery('circul', h_circ);

disp('--- (b) Circulant Convolution Matrix H ---')
H_circ

% DFT matrix
M  = fft(eye(N)) / sqrt(N);
Mh = M';

% Eigenvalues = DFT{channel taps padded to N}
lambda = fft(h_circ).';

disp('--- Eigenvalues of H (diagonal of Λ = MHM^H) ---')
lambda


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (c) FLAT-FADING GAINS FOR EACH SUBCARRIER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--- (c) Flat-fading channel gain per subcarrier |H[k]| ---')
abs(lambda)
