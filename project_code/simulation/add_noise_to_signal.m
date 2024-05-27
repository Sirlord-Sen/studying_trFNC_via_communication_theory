function [noisy_tc, noise, snr_actual] = add_noise_to_signal(tc, dbSNR, noise_type)
% ADD_NOISE_TO_SIGNAL Add noise to a signal.
%
% Syntax:
%   noisy_tc = add_noise_to_signal(tc, dbSNR)
%
% Inputs:
%   - tc: Input signal.
%   - dbSNR: Desired signal-to-noise ratio in decibels (dB) (default: 100).
%   - noise_type: Type of noise to add ('gaussian', 'rayleigh', 'uniform', 'rician').
%
% Outputs:
%   - noisy_tc: Signal with added noise.
%   - noise: The noise signal.
%   - actual SNR after adding noise.
%
% Example:
%   % Generate a signal
%   dbSNR = 100;
%   t = 0:0.01:1;
%   tc = sin(2*pi*5*t);
%
%   % Add noise with SNR of 100 dB
%   noisy_tc = add_noise_to_signal(tc, dbSNR, 'rician');

    % Check if timecourse is given
    if nargin < 1
        error('Not enough input arguments. Please provide the input signal tc and optional desired SNR dbSNR.');
    end
    % Check if both input arguments are provided
    if nargin < 2
        dbSNR = 100; % Default value for dbSNR
    end

    % Check if noise_type is provided
    if nargin < 3
        noise_type = 'gaussian'; % Default to Gaussian noise
    end

     % Check if resamp_method is a valid interpolation method
    valid_types = {'gaussian', 'rayleigh', 'uniform', 'rician'};
    if ~ischar(char(noise_type)) || ~ismember(char(noise_type), valid_types)
        error('Invalid noise type. resamp_method must be one of: ''gaussian'', ''rayleigh'', ''uniform'', ''rician''.');
    end

    % Check if the SNR is non-negative
    if dbSNR < 0
        error('SNR value must be non-negative.');
    end
    
    % Check if the input signal is a vector
    if ~isvector(tc)
        error('Input signal must be a vector.');
    end
    
    % Calculate the power of the input signal
    tc_power = mean(tc.^2);
    
    % Calculate the noise power based on the SNR in dB
    noise_power = tc_power / (10^(dbSNR / 10));
    
    % Generate noise based on the specified type
    if strcmpi(noise_type, 'gaussian')
        % Generate Gaussian noise with the same length as the signal
        noise = sqrt(noise_power) * randn(size(tc));
    elseif strcmpi(noise_type, 'rayleigh')
        % Calculate the scale parameter for the Rayleigh distribution
        s= sqrt(noise_power/2);
        % Generate Rayleigh-distributed noise
        noise = raylrnd(s, size(tc));
    elseif strcmpi(noise_type, 'uniform')
        % Calculate the amplitude of uniform noise based on the SNR
        noise_amplitude = sqrt(12 * noise_power);
        % Generate uniform noise within [-amplitude/2, amplitude/2]
        noise = noise_amplitude * (rand(size(tc)) - 0.5);
    elseif strcmpi(noise_type, 'rician')
        % Calculate the scale parameter for the Rician distribution
        s = sqrt(noise_power/2);
        % Generate Rician-distributed noise
        noise = s * sqrt(randn(size(tc)).^2 + randn(size(tc)).^2);
    end

    % Add the noise to the signal
    noisy_tc = tc + noise;

    % Calculate and display the actual SNR
    snr_actual = 10*log10(mean(tc.^2)/mean(noise.^2));
end