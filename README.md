# EEG

EEG analysis in terms of event-related spectra and multi-scale entropy

Calculate time-frequency decompositions per EEGLAB dataset using eegTimeFrequency.m for tasks with stimulus and response, or eeg3TimeFrequency.m for tasks with three trial events; decompositions are calculated by waveletTransform.m and blended by spectralBlender.m

Collate time-frequency files and average across an electrode cluster using eegLocalSpectra.m

Plot event-related spectra using plotTimeFrequency.m

Find oscillatory peaks using spectralPeaks.m

Calculate composite multi-scale entropy per EEGLAB dataset using eegMultiScaleEntropy.m or eegMultiScaleEntropyData.m followed by eegMultiScaleEntropyAverage.m

Calculate descriptive statistics on any numeric array using descriptiveStatistics.m

#
Rohan O. C. King, @neuroro

If you use this repository of code please cite it as:

King, R. (2023). EEG [GitHub repository]. GitHub. https://github.com/neuroro/EEG
