function f0 = f0_calc(filt_input_z, T, Fs)
% f0_calc  Berekent de fundamentele frequentie (f0) uit een gefilterd 
%          Z-richting signaal, op basis van FFT.
%
% INPUTS:
%   filt_input_z : gefilterd signaal van Z-richting (vector)
%   T            : indexen van het segment dat geanalyseerd wordt
%   Fs           : samplefrequentie in Hz
%
% OUTPUT:
%   f0           : fundamentele frequentie in Hz

    % Segment en lengte
    seg = filt_input_z(T);
    N = length(seg);

    % Preproc: mean removal & window
    seg = seg(:) - mean(seg);
    w = hann(N);
    seg_w = seg .* w;

    % FFT
    Y = fft(seg_w);
    A = abs(Y) / N;
    A = A / mean(w); % Compensatie voor window

    % Single-sided spectrum
    halfIdx = floor(N/2) + 1;
    A_single = A(1:halfIdx);
    if mod(N,2) == 0
        A_single(2:end-1) = 2 * A_single(2:end-1);
    else
        A_single(2:end) = 2 * A_single(2:end);
    end

    % Frequentievector
    f = (0:halfIdx-1) * (Fs / N);

    % Zoek fundament binnen verwacht bereik
    f_search_min = 0.5; 
    f_search_max = 3.0;
    search_idx = find(f >= f_search_min & f <= f_search_max);
    if isempty(search_idx)
        error('Geen fundamentele frequentie gevonden binnen zoekbereik.');
    end
    [~, relmax] = max(A_single(search_idx));
    idxFund = search_idx(relmax);

    % Output
    f0 = f(idxFund);
end
