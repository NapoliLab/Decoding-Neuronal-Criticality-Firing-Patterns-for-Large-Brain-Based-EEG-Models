function [totalEngagement, dynamicEngagementIndex, sumDynamicEngagement, meanDynamicEngagement, stdDynamicEngagement] = calculateEngagementIndex(alphaBand, betaBand, thetaBand)
% ------------------------------- Overview --------------------------------
% calculateEngagementIndex: Calculates the EEG engagement index, a metric
%                           defined as the ratio of beta to (alpha + theta).
%                           The engagement index is used to identify and
%                           classify mental engagement.

% ---------------------------- Function Inputs ----------------------------
    % alphaBand: Output of the Activation Complexity Filter Bank in the
    %            Alpha Band
    % betaBand:  Output of the Activation Complexity Filter Bank in the
    %            Beta Band
    % thetaBand: Output of the Activation Complexity Filter Bank in the
    %            Theta Band

% --------------------------- Function Outputs ----------------------------
    % totalEngagement:   Total Engagement Index (non-dynamic)
    % dynamicEngagementIndex: Dynamic Engagement Index at each time
    % sumDynamicEngagement: Sum of all Dynamic Engagement Values
    % meanDynamicEngagement: Mean of all Dynamic Engagement Values
    % stdDynamicEngagement: Standard Deviation of all Dynamic Engagement Values

% ========================================================================= 
totalEngagement = sum(betaBand) / (sum(alphaBand) + sum(thetaBand));
dynamicEngagementIndex = betaBand ./ (alphaBand + thetaBand);
sumDynamicEngagement = sum(dynamicEngagementIndex);
meanDynamicEngagement = mean(dynamicEngagementIndex);
stdDynamicEngagement = std(dynamicEngagementIndex);

end
