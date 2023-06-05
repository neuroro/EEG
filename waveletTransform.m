function ...
[ spectralPower, phaseDirection, coefficients, frequencies, timeIndices ] ...
  = waveletTransform( signal, samplingRate, frequencyLimits, varargin )
%
% •.° Continuous Wavelet Transform °.•
% _________________________________________________________________________
%
% Spectrally decompose a signal using the continuous wavelet transform with
% a complex Morlet wavelet normalised to unit energy and convert to power
% in signal units squared and phase direction on the unit circle
%
% Edges are removed from each side of the time-spectra in proportion to the
% lowest frequency to ensure the absence of edge effects or distortions
%
% • Function •
% -------------------------------------------------------------------------
% >> [ spectralPower, phaseDirection, coefficients, frequencies, timeIndices ] ...
%      = waveletTransform( signal, samplingRate, frequencyLimits, ...
%                          'OptionName1', optionValue1, 'OptionName2', optionValue2, ... )
%
% Examples:
% >> spectralPower = waveletTransform( signal );
% >> [ ~, ~, coefficients, frequencies, timeIndices ] = waveletTransform( signal );
% >> [ pw, ph, c, f, t ] = waveletTransform( signal, 1000, [2 30] );
% >> [ pw, ph, c, f, t ] = waveletTransform( signal, 0, [], 'FrequencyResolution', 30 );
%
% • Inputs •
% -------------------------------------------------------------------------
% signal:               Time-series to spectrally decompose
% samplingRate:         Sampling rate of the signal in Hz
%                         (default 1000 Hz)
% frequencyLimits:      Frequencies in Hz to decompose as [minimum maximum]
%                         (default 2-60 Hz)
%
% Optional inputs:      Specified as 'OptionName1', optionValue1, ...
%                         names are case insensitive and can be incomplete
% 'FrequencyResolution' Number of frequencies to decompose per octave, 
%                         log-spaced in powers of 2
%                         (default 40)
% 'EdgeSize'            Edge size in number of wavelengths of the lowest
%                         frequency; edges are removed from each side of
%                         the time-spectra to ensure the absence of edge
%                         effects or distortions
%                         (default 1.5 wavelengths per side)
% 'CyclesMinimum'       Number of wavelet cycles at the minimum frequency
%                         and frequencies less than or equal to a cutoff
%                         (default pi * sqrt( 2/log(2) ) = 5.3364)
% 'CyclesCutoff'        Cutoff frequency above which the number of wavelet
%                         cycles increases logarithmically in octaves
%                         (default 12 Hz)
% 'CyclesSlope'         Rate of the logarithmic increase in the number of
%                         wavelet cycles above the cuffoff frequency in
%                         multiples of the minimum number of cycles
%                         (default 1)
%
% • Outputs •
% -------------------------------------------------------------------------
% spectralPower:       Power in signal units squared as frequencies x times
% phaseDirection:      Phase direction on the complex unit circle as 
%                        frequencies x times, used to calculate inter-trial
%                        phase coherence as abs( mean( phase directions ) )
% coefficients:        Wavelet complex coefficients as frequencies x times
% frequencies:         Frequencies in the decomposition
% timeIndices:         Time point indices of the input signal corresponding
%                        to the output time-spectra after excision of edges
%
% • Author •
% -------------------------------------------------------------------------
% Rohan O. C. King
% @neuroro
%
% Copyright 2023 Rohan King
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details (https://www.gnu.org/licenses/).
%
% Please cite this code if you use it


%% Inputs
% -------------------------------------------------------------------------

% Sampling rate in Hz
if nargin < 2 || ~samplingRate
    samplingRate        = 1000;
end

% Frequency limits in Hz
if nargin < 3 || isempty( frequencyLimits )
    frequencyLimits     = [2 60];
end

% Optional input name-value pairs
optionNames             = varargin(1:2:end-1);
optionValues            = varargin(2:2:end);
if length( optionNames ) ~= length( optionValues )
    error( 'Optional inputs are incomplete' )
end

% Frequencies per octave
optionExistence = contains( optionNames, { 'Frequency' 'Resolution' 'r' 'q'  }, 'IgnoreCase', true );
if any( optionExistence )
    iOption             = find( optionExistence, 1 );
    frequencyResolution = optionValues{iOption};
else
    frequencyResolution = 40;
end

% Edge size in cycles or wavelengths of the minimum frequency
optionExistence = contains( optionNames, { 'Edge' 'Size' 'd' 'g' 'z' }, 'IgnoreCase', true );
if any( optionExistence )
    iOption             = find( optionExistence, 1 );
    edgeSize            = optionValues{iOption};
else
    edgeSize            = 1.5;
end

% Minimum number of wavelet cycles
optionExistence = contains( optionNames, { 'Minimum' 'm' }, 'IgnoreCase', true );
if any( optionExistence )
    iOption             = find( optionExistence, 1 );
    cycles              = optionValues{iOption};
else
    cycles              = pi * sqrt( 2/log(2) );
end

% Cutoff frequency for the minimum number of wavelet cycles
optionExistence = contains( optionNames, { 'Cut' 'off' 'cu' 'of' 'ff' }, 'IgnoreCase', true );
if any( optionExistence )
    iOption             = find( optionExistence, 1 );
    frequencyCutoff     = optionValues{iOption};
else
    frequencyCutoff     = 12; % 12 Hz is an upper bound of the neural alpha band
end

% Rate of the logarithmic increase in the number of wavelet cycles
optionExistence = contains( optionNames, { 'Slope' 'sl' 'p' 'rate' }, 'IgnoreCase', true );
if any( optionExistence )
    iOption             = find( optionExistence, 1 );
    cyclesSlope         = optionValues{iOption};
else
    cyclesSlope         = 1; % Multiplied by the minimum number of cycles
end


%% Time-frequency decomposition parameters
% -------------------------------------------------------------------------

% Frequencies
octaves      = log2( frequencyLimits(2) / frequencyLimits(1) );
nFrequencies = ceil( frequencyResolution * octaves );
log2fLimits  = log2( frequencyLimits );
frequencies  = 2 .^ linspace( log2fLimits(1), log2fLimits(2), nFrequencies );

% Decomposition times
waveletTimes = -2:1/samplingRate:2;    % Time domain of the wavelet
nWavelet     = length( waveletTimes ); % Number of wavelet time points
nSignal      = length( signal );       % Number of signal time points
nConvolution = nSignal + nWavelet - 1; % Number of convolution time points
halfWavelet  = (nWavelet - 1)/2;       % Half-width of the wavelet

% Edges
longestCycle = samplingRate/frequencyLimits(1);
edge         = ceil( edgeSize * longestCycle );
if nSignal <= 2*edge
    error( 'The signal is too short for the specified lower frequency limit' )
end


%% Wavelet cycles
% -------------------------------------------------------------------------

% The number of cycles equals a minimum number for all frequencies below a
% cutoff then increases logarithmically in octaves with a sigmoid curve
% transition between constant and logarithmic to approximate continuity

% The width of the Gaussian envelope of the Morlet wavelet is given by
% sigma = n / ( 2*pi*f ), where n is number of cycles and f is frequency
%
% n = pi / sqrt( 2*log(2) ) = 2.6682 cycles gives a wavelet full width at
% half maximum of 1/f, and is recommended as the minimum (Cohen, 2019)
%
% n = 3 cycles at the lowest frequency is the default in EEGLAB 2019-2023
% (Delorme & Makeig, 2004)
%
% w0 = pi * sqrt( 2/log(2) ) = 5.3364 centre frequency gives a wavelet with
% a 2nd peak that is half of the main peak (Smith, 2022)
%
% Tallon-Baudry, Bertrand, Delpuech, & Permier (1997) argue that f/fSigma
% should be approximately 5 or greater (citing Grosmann et al., 1989),
% where fSigma = 1 / ( 2*pi*sigma ) is the Gaussian frequency shape and
% sigma is the width of the Gaussian
% If sigma = n / ( 2*pi*f ) then f/fSigma = n, ergo n should be >= 5
% 
% n = pi * sqrt( 2/log(2) ) cycles corresponds to a wavelet full width at
% half maximum of 2/f, and gives a 1 Hz sigma = 0.8493

% Minimum number of wavelet cycles
% Defined above as the optional input value or pi * sqrt( 2/log(2) )

% Frequencies above which the number of cycles increases
nAboveCutoff       = sum( frequencies > frequencyCutoff );

% Floating point error correction
if nAboveCutoff == 1 && round( frequencies(end) ) == round( frequencyCutoff )
    nAboveCutoff   = 0;
end

% Cycles constant at the minimum number for frequencies below the cutoff
waveletCycles      = cycles * ones( 1, nFrequencies - nAboveCutoff );

% Cycles increasing for frequencies above the cutoff
if nAboveCutoff

    % Magnitude of increase
    octavesAbove   = log2( frequencyLimits(2) / frequencyCutoff );
    maxCycles      = cycles * ( 1 + cyclesSlope * octavesAbove );

    % Sigmoid transition to approximate continuity
    sigmoidSlope   = 2;
    waveletCycles2 = sigmoidSpace( cycles, maxCycles, nAboveCutoff, sigmoidSlope );
    interceptRatio = exp(1)/10;
    nCycles2       = floor( nAboveCutoff * interceptRatio );
    waveletCycles2 = waveletCycles2(1:nCycles2);

    % Cycles increasing logarithmically in octaves
    nCycles3       = ceil( nAboveCutoff * (1 - interceptRatio) ) + 1;
    waveletCycles3 = 2 .^ linspace( log2( waveletCycles2(end) ), log2( maxCycles ), nCycles3 );
    waveletCycles3 = waveletCycles3(2:end); % Remove duplicate value

    % Number of cycles as a piece-wise function
    waveletCycles  = [ waveletCycles waveletCycles2 waveletCycles3 ];

end


%% Time-frequency decomposition
% -------------------------------------------------------------------------

% Pre-allocate
coefficients   = zeros( nFrequencies, nSignal );

% Signal in double precision
signal         = double( signal );

% Signal Fourier transform
signalSpectrum = fft( signal, nConvolution );

% Decompose each frequency
for f = 1:nFrequencies


    % Morlet wavelet at this frequency
    % ---------------------------------------------------------------------

    % Gaussian envelope
    gaussianWidth     = waveletCycles(f) / ( 2*pi*frequencies(f) );       % (Cohen, 2019)
    gaussianEnvelope  = exp( -waveletTimes.^2 ./ ( 2*gaussianWidth^2 ) ); % (Cohen, 2019; Tallon-Baudry, Bertrand, Delpuech, & Permier, 1997)

    % Complex sinusoid
    complexSinusoid   = exp( 1i * 2*pi*frequencies(f) .* waveletTimes );  % (Cohen, 2019; Tallon-Baudry, Bertrand, Delpuech, & Permier, 1997)

    % Normalisation
    %   Divide by the square root of the integral of the power of the
    %   wavelet, to give the wavelet unit energy (Ashmead, 2012; Torrence &
    %   Compo, 1998)
    %   Integral from -Inf to Inf of Morlet * conjugate( Morlet ) dt
    %   This simplifies to the Gaussian integral = sqrt( pi / a ), with 
    %   a = 2*pi^2*f^2/n^2, which gives the normalisation defined below
    normalisation     = ( 2*pi*frequencies(f)^2 / waveletCycles(f)^2 )^0.25; % For constant n = pi * sqrt( 2/log(2) ) this simplifies to ( 2*log(2) / pi )^0.25 * sqrt( frequencies(f) );

    % Morlet wavelet
    morletWavelet     = normalisation * gaussianEnvelope .* complexSinusoid;


    % Convolution in the frequency domain
    % ---------------------------------------------------------------------

    % Morlet wavelet Fourier transform
    morletSpectrum    = fft( morletWavelet, nConvolution );

    % Convolution
    convolution       = morletSpectrum .* signalSpectrum;
    convolution       = ifft( convolution, nConvolution );


    % Signal at this frequency
    % ---------------------------------------------------------------------

    % Complex coefficients of the signal at this frequency
    %   Excise half the wavelet width from each side of the convolution
    coefficients(f,:) = convolution(1+halfWavelet:end-halfWavelet);


end

% Excise a number of wavelengths from each side of the time-frequency
% coefficients to limit edge artifacts
timeIndices    = 1+edge:nSignal-edge;
coefficients   = coefficients(:,timeIndices);


%% Power and phase direction
% -------------------------------------------------------------------------

% Spectral power in signal units squared
spectralPower  = coefficients .* conj( coefficients );                     % Slightly cleaner than abs(z).^2 for simple cases

% Phase direction on the unit circle
%   Complex value equivalent to a unit vector
%   Can be used to calculate inter-trial phase coherence as
%   itpc = abs( mean( phase directions ) );
phaseDirection = coefficients ./ abs( coefficients );                      % Equivalent to exp( 1i * angle( coefficients ) ) different by about 1e-16 in R2022a


% _________________________________________________________________________
end


% • References •
% -------------------------------------------------------------------------
% Ashmead, J. (2012). Morlet Wavelets in Quantum Mechanics. Quanta, 1(1), 
%   58-70. https://doi.org/10.12743/quanta.v1i1.5
% Cohen, M. X. (2019). A better way to define and describe Morlet wavelets 
%   for time-frequency analysis. NeuroImage, 199, 81-86. 
%   https://doi.org/10.1016/j.neuroimage.2019.05.048
% Delorme, A. & Makeig, S. (2004). EEGLAB: an open-source toolbox for
%   analysis of single-trial EEG dynamics. Journal of Neuroscience Methods,
%   134, 9-21. https://doi.org/10.1016/j.jneumeth.2003.10.009
% Smith, J. O. (2022). Continuous Wavelet Transform. In Spectral Audio
%   Signal Processing (2011 ed.). Stanford.
%   https://ccrma.stanford.edu/~jos/sasp/Continuous_Wavelet_Transform.html
% Tallon-Baudry, C., Bertrand, O., Delpuech, C., & Permier, J. (1997).
%   Oscillatory gamma-band (30-70 Hz) activity induced by a visual search
%   task in humans. The Journal of neuroscience : the official journal of
%   the Society for Neuroscience, 17(2), 722–734.
%   https://doi.org/10.1523/JNEUROSCI.17-02-00722.1997
% Torrence, C., & Compo, G. P. (1998). A Practical Guide to Wavelet 
%   Analysis. Bulletin of the American Meteorological Society, 79, 61-78.



function y = sigmoidSpace( x1, x2, N, m )
%
% •.° Sigmoid Space °.•
% _________________________________________________________________________
% 
% Sigmoidally spaced vector between x1 and x2 over N points with slope m


% Default slope
if nargin < 4
    m = exp(1);
end

% Domain and range
sigmoidDomain  = linspace( -m, m, N );
sigmoidRange   = x2 - x1;

% Sigmoid curve about the origin 
sigmoidOrigin  = erf( sigmoidDomain );

% Sigmoid curve between 0 and 1
normalisation  = erf( m ) - erf( -m );
sigmoidZeroOne = ( sigmoidOrigin + erf( m ) ) / normalisation;

% Sigmoid space 
% 0-s-S-1 scaled by the range and shifted by the start point
y = sigmoidRange * sigmoidZeroOne + x1;


% _________________________________________________________________________
end
