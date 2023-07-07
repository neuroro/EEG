# EEG

EEG analysis in terms of event-related spectra and multi-scale entropy

Calculate time-frequency decompositions per EEGLAB dataset using eegTimeFrequency.m for tasks with stimulus and response, or eeg3TimeFrequency.m for tasks with three trial events; decompositions are calculated by waveletTransform.m and blended by spectralBlender.m

Collate time-frequency files and average across an electrode cluster using eegLocalSpectra.m

Plot event-related spectra using plotTimeFrequency.m

Find oscillatory peaks using spectralPeaks.m

Calculate composite multi-scale entropy per EEGLAB dataset using eegMultiScaleEntropy.m or eegMultiScaleEntropyData.m followed by eegMultiScaleEntropyAverage.m

Calculate descriptive statistics on any numeric array using descriptiveStatistics.m


# Coming soon...

Monte Carlo permutation test of paired data

Monte Carlo permutation test of unpaired data

Monte Carlo custer permutation tests with cluster formation and a distance function

Plot line graphs with shading representing the standard error of the mean and optional shading for a sub-domain window (an ERP component or frequency band, for example) using plotShadedLine.m


# Me
Rohan O. C. King, @neuroro

I hope my code finds use; if you use it please cite the function or the repository, for example in APA as, respectively:

King, R. (2023). functionName [MATLAB code]. GitHub. https://github.com/neuroro/EEG/functionName.m

King, R. (2023). EEG [GitHub repository]. GitHub. https://github.com/neuroro/EEG

Updated Sat 8th July 2023
