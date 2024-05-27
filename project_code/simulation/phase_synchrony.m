function ps = phase_synchrony(tc1, tc2)
    % phase_synchrony calculates the cosine of the relative phase between two time series
    %
    % Inputs:
    %   tc1 - first input time series
    %   tc2 - second input time series
    %
    % Output:
    %   crp - cosine of the relative phase between tc1 and tc2
    
    % Calculate the analytic signal using the Hilbert transform
    hilbert_tc1 = hilbert(tc1);
    hilbert_tc2 = hilbert(tc2);
    
    % Compute the phase angles of the analytic signals
    phase_tc1 = angle(hilbert_tc1);
    phase_tc2 = angle(hilbert_tc2);
    
    % Calculate the relative phase
    relative_phase = phase_tc1 - phase_tc2;
    
    % Compute the cosine of the relative phase
    ps = cos(relative_phase);
end