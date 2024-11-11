function temporalCrossBandEntropy = computeTemporalCrossBandEntropy(eegSignal)
% ------------------------------- Overview --------------------------------
% computeTemporalCrossBandEntropy: Computes the Temporal Cross-Band Entropy
%                                 of a given multi-channel EEG signal.

% ---------------------------- Function Inputs ----------------------------
    % eegSignal: a matrix of (numChannels x numSamples), where numChannels
    %            is the number of EEG leads, and numSamples is the number of samples.
    
% ---------------------------- Function Output ----------------------------
    % temporalCrossBandEntropy: returns the Temporal Cross-Band Entropy of the EEG signal.
    
% ========================================================================= 
% Number of EEG Channels and Samples
numChannels = size(eegSignal, 1);
numSamples = size(eegSignal, 2);

% Allow equal values to share the same index
indexedMatrix = zeros(size(eegSignal));

% Find the unique values of leads in each column over time
for sampleIdx = 1:numSamples
    [~, ~, indexedMatrix(:, sampleIdx)] = unique(eegSignal(:, sampleIdx), 'first');
end

% Assign unique numbers to each pattern using a base-n number system
patternIndices = numChannels.^(0:numChannels-1) * (indexedMatrix - 1);

% Sort the unique pattern numbers
[~, uniquePatternIndices, ~] = unique(sort(patternIndices), 'first');

% Compute Permutation entropy
patternCounts = diff([uniquePatternIndices; (length(patternIndices) + 1)]);
patternProbabilities = patternCounts / sum(patternCounts);
temporalCrossBandEntropy = -sum(patternProbabilities .* log2(patternProbabilities));
temporalCrossBandEntropy = temporalCrossBandEntropy / log2(factorial(numChannels));

end
