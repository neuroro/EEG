function ...
[ spectralPower, phaseDirection, coefficients ] ...
  = spectralBlender( spectralPower, phaseDirection, coefficients, ...
                     frequencies, timePoints, blendTimes, blendDuration )
%
% •.° Time-Frequency Blender °.•
% _________________________________________________________________________
%
% Blend or limit time-frequency decompositions of EEG data around a time
% point
%
% Blend spectral power and phase direction up to full value or down to 0
% (baseline dB or coherence) around a time point over a finite duration
% and set coefficients outside-the-window to NaN + NaNi
%
% Alternatively, limit all spectra after or before a time point to NaN or
% NaN + NaNi (as appropriate) if the duration is + or - Inf
%
% Finite blending:
% 1. Sigmoid blend power from or to 0 dB around the blend time
% 3. Sigmoid blend phase direction from or to 0 + 0i around the blend time
% 4. Set coefficients up to or from the blend time to NaN + NaNi
% 
% • Function •
% -------------------------------------------------------------------------
% >> [ spectralPower, phaseDirection, coefficients ] ...
%    = spectralBlender( spectralPower, phaseDirection, coefficients, ...
%                       frequencies, timePoints, blendTimes, blendDuration )
% 
% spectralPower  as frequencies x time points
% phaseDirection as frequencies x time points
% coefficients   as frequencies x time points
%
% • Scalar blend time •
% -------------------------------------------------------------------------
%
% Positive blend duration     blendTime
%                                 |          
% coefficients      NaN..........NaN 1..........1
% phaseDirection    0.........0./////'1.........1
% spectralPower     0.........0./////'1.........1
%                             |_______|
%                           blendDuration
%
% Negative blend duration     blendTime
%                                 |          
% coefficients      1..........1 NaN..........NaN
% phaseDirection    1.........1`\\\\\.0.........0
% spectralPower     1.........1`\\\\\.0.........0
%                             |_______|
%                           blendDuration
%
% • Vector of blend times •
% -------------------------------------------------------------------------
%
% Positive blend duration   blendTime(1)        blendTime(2)
%                                |                   |
%                         ...0./////'1...........1`\\\\\.0...
%                            |_______|           |_______|
%                          blendDuration       blendDuration
%
% Negative blend duration   blendTime(1)        blendTime(2)
%                                |                   |
%                         ...1`\\\\\.0...........0./////'1...
%                            |_______|           |_______|
%                          blendDuration       blendDuration
%
% • Requires •
% -------------------------------------------------------------------------
% sspace.m (King, 2023)
%
% • Author •
% -------------------------------------------------------------------------
% Rohan O. C. King, 2023
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


% Set data to this value after or before the blend
exBlend = 0;

% Default blend duration
if nargin < 7
    blendDuration = 100;
end

% No blending
if isinf( blendDuration )
    noBlend = true;
else
    noBlend = false;
end

% Number of blend times
nBlends = length( blendTimes );

% Middle cut
if length( blendTimes ) > 1 && blendDuration < 0
    separate        = true;
    separateBlender = cell(1,nBlends);
    separateLimiter = cell(1,nBlends);
else
    separate        = false;
end

% Loop through: Blend times
for b = 1:nBlends

    blendTime = blendTimes(b);

    % Sanity check: Blend time is a number
    if ~isnan( blendTime )


        % Blend direction
        % -----------------------------------------------------------------

        % Reverse blender direction every other time
        if ~mod(b,2)
            blendDuration = -blendDuration;
        end

        % Blend in
        if blendDuration > 0
            proportionAtStart = 0;
            proportionAtEnd   = 1;

        % Blend out
        elseif blendDuration < 0
            proportionAtStart = 1;
            proportionAtEnd   = 0;

        end


        % Sigmoid blend times
        % -----------------------------------------------------------------

        % Centered on blend time
        halfTime           = round( abs( blendDuration )/2 );
        blendStart         = blendTime - halfTime;
        blendEnd           = blendTime + halfTime;

        % Find time point indices
        % Absolute differences closest to 0
        blendStart         = timePoints - blendStart;
        blendEnd           = timePoints - blendEnd;
        limitStart         = timePoints - blendTime;
        blendStart         = abs( blendStart );
        blendEnd           = abs( blendEnd   );
        limitStart         = abs( limitStart );
        [ ~, iBlendStart ] = min( blendStart );
        [ ~, iBlendEnd   ] = min( blendEnd   );
        [ ~, iLimitStart ] = min( limitStart );

        % Sigmoid blend time points
        blendTimePoints    = timePoints(iBlendStart:iBlendEnd);
        blendTimePoints    = length( blendTimePoints );


        % Build the blender and limiter
        % -----------------------------------------------------------------

        % Pre-allocate the blender and limiter
        nTimePoints        = length( timePoints  );
        nFrequencies       = length( frequencies );
        blender            = ones( nFrequencies, nTimePoints );
        limiter            = ones( nFrequencies, nTimePoints );

        % Loop through: Frequencies
        for f = 1:length( frequencies )

            % Build the sigmoid blend
            if abs( blendDuration ) < 200
                sigmoidSlope = exp(1);
            else
                sigmoidSlope = 1 + sqrt(5); % Twice the golden ratio
            end
            sigmoidBlend = sigmoidSpace( proportionAtStart, proportionAtEnd, blendTimePoints, sigmoidSlope );

            % Sigmoid blend the blender around the blend time
            blender(f,iBlendStart:iBlendEnd) = sigmoidBlend;

        end % for Frequencies

        % Build the limiter and zero the blender before or after the trial window
        switch proportionAtStart
            case 0
                blender(:,1:iBlendStart-1) = exBlend;      % Before sigmoid blend
                limiter(:,1:iLimitStart)   = NaN + 1i*NaN; % Up to blend time
            case 1
                blender(:,iBlendEnd+1:end) = exBlend;      % After sigmoid blend
                limiter(:,iLimitStart:end) = NaN + 1i*NaN; % From blend time
        end

        % Unplug the blender (chopping is cleaner)
        if noBlend
            blender = real( limiter );
        end


        % Put time-frequency in the blender (and limiter)
        % -----------------------------------------------------------------
        if ~separate
            spectralPower  = spectralPower  .* blender;
            phaseDirection = phaseDirection .* blender;
            coefficients   = coefficients   .* limiter;


        % Middle cut: Put the blender (and limiter) in the cupboard
        else
            separateBlender{b} = blender;
            separateLimiter{b} = limiter;

        end


    % Don't even build the blender
    else
        disp( 'Blend time does not make sense' )

    end % if Sanity check

end % for Blend times


% Middle cut
% -------------------------------------------------------------------------
if separate

    % Pre-allocate
    blender          = ones( nFrequencies, nTimePoints );
    limiter          = ones( nFrequencies, nTimePoints );
    blenderSharedNaN = ones( nFrequencies, nTimePoints );
    limiterSharedNaN = ones( nFrequencies, nTimePoints );

    % Shared NaN zones will be erased
    for b = 1:nBlends
        blenderSharedNaN = blenderSharedNaN .* isnan( separateBlender{b} );
        limiterSharedNaN = limiterSharedNaN .* isnan( separateLimiter{b} );
    end

    % Unshared NaN zones will be kept
    for b = 1:nBlends
        blenderAntiNaN = isnan( separateBlender{b} ) .* ~blenderSharedNaN;
        limiterAntiNaN = isnan( separateLimiter{b} ) .* ~limiterSharedNaN;
        separateBlender{b}(logical( blenderAntiNaN )) = 1;
        separateLimiter{b}(logical( limiterAntiNaN )) = 1 + 1i;
    end

    % Combine separate blenders and limiters
    for b = 1:nBlends
        blender = blender .* separateBlender{b};
        limiter = limiter .* separateLimiter{b};
    end

    % Put time-frequency in the blender (and limiter)
    spectralPower  = spectralPower  .* blender;
    phaseDirection = phaseDirection .* blender;
    coefficients   = coefficients   .* limiter;

end


% _________________________________________________________________________
end



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
