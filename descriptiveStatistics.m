function Stats = descriptiveStatistics( data, dim )
%
% •.° Descriptive Statistics °.•
% _________________________________________________________________________
%
% Describe the distribution of numeric data in terms of central tendency,
% dispersion, standardised moments, and related measures
%
% • Usage •
% -------------------------------------------------------------------------
% >> Stats = descriptiveStatistics( data, dim )
%
% Examples:
% >> descriptiveStatistics( data )
% >> Stats = descriptiveStatistics( data, 1 )
%
% Inputs:
%   data: Numeric data in a vector, matrix, or array
%   dim:  Dimension to calculate statistics on
%         (optional input, default 1st dimension for matrices and arrays)
%
% Output: Struct of descriptive statistics in named fields
%
% • Statistics •
% -------------------------------------------------------------------------
% Central tendency:
% Mean
% Median
% Geometric mean or quasi-arithmetic log mean for positive ratio data or 
%   logarithmic ratio-level data (undefined for negative data and not valid 
%   for interval-level data)
%   = nth root of the product of Xn = e ^ E[ln(X)] for X > 0
%
% Dispersion:
% Standard deviation (SD)
% Normalised median absolute deviation, which is a robust estimate of the
%   SD for normally-distributed data (Akinshin, 2022)
%   = MAD / inverse cumulative distribution function (quantile function) of
%   the standard normal distribution at 0.75 with sample size bias
%   correction (Park et al., 2020)
% Normalised inter-quartile range, which is also a robust estimate of the
%   SD for normally-distributed data
%   = IQR / ( 2 root 2 * inverse error function at 0.5 ) (Ivezić, 2019)
% Geometric standard deviation for positive ratio data or logarithmic 
%   ratio-level data (undefined for negative data and not valid for 
%   interval-level data)
%   = e ^ square root of variance of ln(X)
% Standard error of the mean
% Mean of the absolute deviations from the mean
% Median of the absolute deviations from the median (MAD)
% Inter-quartile range (IQR)
% Range (minimum and maximum of the data)
% Coefficient of variation for normally-distributed ratio-level data (not 
%   valid for interval-level data)
%   = SD / absolute mean with sample size bias correction (%)
% Robust coefficient of variation for normally-distributed ratio-level data
%   (not valid for interval-level data)
%   = robust SD / absolute median with sample size bias correction (%)
% Geometric coefficient of variation for skewed or log-normally distributed
%   ratio-level data (not valid for interval-level data)
%   = square root of ( e ^ ( variance of ln(X) ) - 1 ) (%)
%
% Shape:
% Skewness
%   = E[(X - E[X])^3] / s^3 with sample size bias correction
% Non-parametric skew
%   = mean - median / SD
% Kurtosis, which is a measure of tailedness (the extent of outliers in the
%   tails) not of peakedness (Westfall, 2014)
%   = E[(X - E[X])^4] / s^4 with sample size bias correction
% Excess kurtosis, measuring tailedness relative to a normal distribution
%   = kurtosis - 3
%
% • Author •
% -------------------------------------------------------------------------
% Rohan O. C. King, 2023
% @neuroro
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


%% Default dimension
% -------------------------------------------------------------------------

if nargin < 2 || ~dim
    if isvector( data )
        [ ~, dim ] = max( size( data ) );
    else
        dim = 1;
    end
end


%% Descriptive statistic names
% -------------------------------------------------------------------------

centralTendency     = { 'Mean'          ...
                        'Median'        ...
                        'GeometricMean' };

dispersion          = { 'StandardDeviation'                 ...
                        'NormalisedMedianAbsoluteDeviation' ...
                        'NormalisedInterQuartileRange'      ...
                        'GeometricStandardDeviation'        ...
                        'StandardErrorOfTheMean'            ...
                        'MeanAbsoluteDeviation'             ...
                        'MedianAbsoluteDeviation'           ...
                        'InterQuartileRange'                ...
                        'Range'                             ...
                        'CoefficientOfVariation'            ...
                        'RobustCoefficientOfVariation'      ...
                        'GeometricCoefficientOfVariation'   };

standardisedMoment  = { 'Skewness'          ...
                        'NonParametricSkew' ...
                        'Kurtosis'          ...
                        'ExcessKurtosis'    };

statisticNames      = [ centralTendency    ...
                        dispersion         ...
                        standardisedMoment ];


%% Descriptive statistic function handles
% -------------------------------------------------------------------------

centralTendencyF    = { @( X, dim ) mean( X, dim, 'omitnan' ),   ...
                        @( X, dim ) median( X, dim, 'omitnan' ), ...
                        @dGeometricMean                          };

dispersionF         = { @( X, dim ) std( X, 0, dim, 'omitnan' ), ...        % Standard deviation
                        @dNormalisedMedianAbsoluteDeviation,     ...
                        @dNormalisedInterQuartileRange,          ...
                        @dGeometricStandardDeviation,            ...
                        @dStandardErrorOfTheMean,                ...
                        @( X, dim ) mad( X, 0, dim ),            ...        % Mean of the absolute deviations from the mean
                        @( X, dim ) mad( X, 1, dim ),            ...        % Median of the absolute deviations from the median
                        @( X, dim ) iqr( X, dim ),               ...        % Inter-quartile range
                        @dRange,                                 ...
                        @dCoefficientOfVariation,                ...
                        @dRobustCoefficientOfVariation,          ...
                        @dGeometricCoefficientOfVariation,       };

standardisedMomentF = { @( X, dim ) skewness( X, 0, dim ),       ...        % Skewness with sample size bias correction
                        @dNonParametricSkew,                     ...
                        @( X, dim ) kurtosis( X, 0, dim ),       ...        % Kurtosis with sample size bias correction
                        @( X, dim ) kurtosis( X, 0, dim ) - 3    };         % Excess kurtosis with sample size bias correction

statisticFunctions  = [ centralTendencyF    ...
                        dispersionF         ...
                        standardisedMomentF ];


%% Calculate descriptive statistics
% -------------------------------------------------------------------------

% Number of statistics
nStatistics = length( statisticNames );

% Loop through: Descriptive statistics
for s = 1:nStatistics

    % Current name and function
    statisticName         = statisticNames{s};
    statisticFunction     = statisticFunctions{s};

    % Calculate descriptive statistic
    descriptiveStatistic  = statisticFunction( data, dim );

    % Store in struct
    Stats.(statisticName) = squeeze( descriptiveStatistic );

end


% _________________________________________________________________________
end



%%
% •.° Descriptive Statistic Functions °.•
% _________________________________________________________________________


% • Geometric Mean or Quasi-Arithmetic Log Mean  •
% -------------------------------------------------------------------------
function m = dGeometricMean( X, dim )

% Check for negative values
if any( X < 0 )
    m = NaN;
    disp( [ '[' 8 'Warning: Geometric mean is undefined for negative values' ']' 8 ] )

% e ^ E[ln(X)]
else
    m = exp( mean( log(X), dim, 'omitnan' ) );

end

end


% • Geometric Standard Deviation •
% -------------------------------------------------------------------------
function s = dGeometricStandardDeviation( X, dim )

% Check for negative values
if any( X < 0 )
    s = NaN;
    disp( [ '[' 8 'Warning: Geometric standard deviation is undefined for negative values' ']' 8 ] )

% e ^ standard deviation of ln(X)
else
    s = exp( std( log(X), 0, dim, 'omitnan' ) );

end

end


% • Standard Error of the Mean •
% -------------------------------------------------------------------------
function s = dStandardErrorOfTheMean( X, dim )

% Sample size
N = size( X, dim );

% Standard deviation / square root of the sample size
s = std( X, 0, dim, 'omitnan' ) / sqrt(N);

end


% • Normalised Median Absolute Deviation •
% -------------------------------------------------------------------------
function normalisedMAD = dNormalisedMedianAbsoluteDeviation( X, dim )

% Median absolute deviation from the median
deviation     = mad( X, 1, dim ); 

% Normalisation
% MAD -> SD for normally distributed data
% Phi inverse = inverse cumulative distribution function (quantile 
% function) of the standard normal distribution at 0.75
normalisation = 1 / icdf( 'Normal', 0.75, 0, 1 );

% Sample size bias correction (Park et al., 2020)
N             = size( X, dim );
correction    = 1 / ( 1 - 0.804168866 * N ^ -1.008922 );

% NMAD = MAD * normalisation * bias correction
normalisedMAD = deviation * normalisation * correction;

end


% • Normalised Inter-Quartile Range •
% -------------------------------------------------------------------------
function normalisedIQR = dNormalisedInterQuartileRange( X, dim )

% NIQR = IQR / 2 root 2 inverse error function at 0.5 (Ivezić, 2019)
normalisedIQR = iqr( X, dim ) / ( 2 * sqrt(2) * erfinv(0.5) );

end


% • Range •
% -------------------------------------------------------------------------
function mM = dRange( X, dim )

m  = min( X, [], dim, 'omitnan' );
M  = max( X, [], dim, 'omitnan' );
mM = [ m M ];

end


% • Coefficient of Variation •
% -------------------------------------------------------------------------
function coefficient = dCoefficientOfVariation( X, dim )

% Coefficient valid for ratio-level data not for interval-level data
% = SD / mean
coefficient = std( X, 0, dim, 'omitnan' ) ./ mean( X, dim, 'omitnan' );

% Sample size bias correction for normally distributed data
N           = size( X, dim );
correction  = 1 + 1 / (4 * N);
coefficient = coefficient * correction;

% Percentage value
coefficient = coefficient * 100;

end


% • Robust Coefficient of Variation •
% -------------------------------------------------------------------------
function coefficient = dRobustCoefficientOfVariation( X, dim )

% Normalised IQR is a robust estimate of the standard deviation for
% normally distributed data
robustSD    = iqr( X, dim ) / (2 * sqrt(2) * erfinv(0.5));

% Coefficient valid for ratio-level data not for interval-level data
% = robust SD / median
coefficient = robustSD ./ median( X, dim, 'omitnan' );

% Sample size bias correction for normally distributed data
N           = size( X, dim );
correction  = 1 + 1 / (4 * N);
coefficient = coefficient * correction;

% Percentage value
coefficient = coefficient * 100;

end


% • Geometric Coefficient of Variation •
% -------------------------------------------------------------------------
function coefficient = dGeometricCoefficientOfVariation( X, dim )

% Coefficient for log-normally distrbiuted data 
% = square root of ( e ^ ( variance of ln(X) ) - 1 )
coefficient = sqrt( exp( var( log(X), 0, dim, 'omitnan' ) ) - 1 );

% Percentage value
coefficient = coefficient * 100;

end


% • Non-Parametric Skew •
% -------------------------------------------------------------------------
function skew = dNonParametricSkew( X, dim )

% Parameters
m    = mean( X, dim, 'omitnan' );
v    = median( X, dim, 'omitnan' );
SD   = std( X, 0, dim, 'omitnan' );

% Skew = mean - median / SD
skew = ( m - v ) ./ SD;

end


% _________________________________________________________________________



%% • References •
% -------------------------------------------------------------------------
%
% Akinshin, A. (2022). Finite-sample bias-correction factors for the median
%   absolute deviation based on the Harrell-Davis quantile estimator and
%   its trimmed modification [pre-print]. arXiv.
%   https://doi.org/10.48550/arXiv.2207.12005
%
% Ivezić, Ž., Connolly, A. J., VanderPlas, J. T., & Gray, A. (2019).
%   Statistics, data mining, and machine learning in astronomy: A practical
%   python guide for the analysis of survey data. Princeton University 
%   Press.
%
% Park, C., Kim, H., & Wang, M. (2022). Investigation of finite-sample
%   properties of robust location and scale estimators. Communications in
%   Statistics-Simulation and Computation, 51(5), 2619-2645.
%   https://doi.org/10.1080/03610918.2019.1699114
%
% Westfall, P. H. (2014). Kurtosis as peakedness, 1905–2014. RIP. The
%   American Statistician, 68(3), 191-195.
%   https://doi.org/10.1080/00031305.2014.917055


