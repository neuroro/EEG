function eegData = linearDetrendZeroing( eegData )
%
% •.° Linear Detrend Zeroing °.•
% _________________________________________________________________________
%
% Zero the start and end points of a (long) continuous section of
% electroencephalographic (EEG) data (a block) using linear detrend of the
% entire block for each electrode separately. Doing so enables multiple
% blocks of data to be subsequently concatenated into unified continuous
% EEG signals per electrode without inter-block discontinuities, thereby
% annealing the EEG recording.
%
% • Function •
% -------------------------------------------------------------------------
% >> eegData = linearDetrendZeroing( eegData );
% 
% Input:
%   eegData: a continuous section of EEG data as electrodes x time points
%
% • Authors •
% -------------------------------------------------------------------------
% Rohan O. C. King & Chris C. King, 2023
%
% Copyright 2023 Rohan King & Chris King
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


% Number of electrodes
nElectrodes = size( eegData, 1 );

% Number of time points
nTimePoints = size( eegData, 2 );

% Voltages at the start and end per electrode
startPoints = eegData(:,1);
endPoints   = eegData(:,end);

% Pre-allocate
linearTrend = zeros( nElectrodes, nTimePoints );

% Start-to-end linear trend of the entire data per electrode
for e = 1:nElectrodes

    % Current electrode start and end voltages
    startVoltage = startPoints(e);
    endVoltage   = endPoints(e);

    % Current electrode linear trend
    linearTrend(e,:) = linspace( startVoltage, endVoltage, nTimePoints );

end

% Subtract the linear trend from the data
eegData = eegData - linearTrend;

% Near-zero the start and end points
%   Replace exact zeros with sign-matched near-zero random numbers to avoid
%   numerical issues arising in subsequent processing steps and analyses
startPolarity  = sign( startPoints );
endPolarity    = sign( endPoints   );
eegData(:,1)   = startPolarity .* rand( size( startPolarity ) ) .* 1e-15;
eegData(:,end) = endPolarity   .* rand( size( endPolarity   ) ) .* 1e-15;


% _________________________________________________________________________
end