function [Power] = optimizedEEGFilterBank(EEG_Signal,F_s, ...
    Sigma,PlotFilters,FiltCoef)
% ------------------------------- Overview --------------------------------
%   Activation_Complexity_Filter_Bank: This approach designs the wavelets
%   filters using a flat gaussian function in the frequency domain and
%   transforms the the EEG time series to the frequency domain to filter
%   the EEG time series. A gaussian filter is then implemented for 
%   smoothing. 


% ---------------------------- Function Inputs ----------------------------
    % EEG_Signal: EEG Time Series (N x 1) row vector 
    % F_s: Sampling Frequency
    % Sigma:  The gaussian Mask implementation is designed to smooth the
    %            the power out. The varible is scalar value that is
    %            described by sigma. 
    %            (Note: With a sampling fs=256 we have typically used a
    %            sigma around 60)
    % PlotFilters: A binary response, where the value of 1 would plot the
    %              designed filters and a value of 0 or leaving the input
    %              empty will not plot the filters.
    % FiltCoef: This is a binary response leaving the input empty or
    %           providing new design coefficients. If you want the option
    %           to use the default filter coefficients previously designed
    %           or a different set of coefficients. The design coefficient
    %           input should be set by f x 3 dimensions, where f is the
    %           number of filters and the columns represent parameters a, 
    %           b, and Cf. 

    
% ---------------------------- Function Output ----------------------------
    % Power: The output is the power of each designed filter. The
    %        variable dimensions are N x f. 
    
% =========================================================================  
% NOTES AND ASSUMPTIONS: The filter design was based on minimizing the
% standard deviation of the plateaut value of the filters, while
% maintaining the below EEG frequency cutoffs. The cutoff values were aimed
% to achieve a (1/e) where every other adjacent filter is orthogonal. 


% -----------------------------------
%         EEG Frequency Bands 
% -----------------------------------
     % Frequency Bands 
       % Delta       0.0 -   3.5 Hz
       % Theta       4.0 -   7.0 Hz
       % Alpha       8.0 -  12.0 Hz 
       % Low Beta   13.0 -  15.0 Hz
       % Mid Beta   15.0 -  18.0 Hz
       % High Beta  18.0 -  40.0 Hz
       % Gamma      32.0 - 100.0 Hz

%==========================================================================
% Contributions: This Code was designed and written by Matthew Demas and
%                Nicholas J. Napoli
%==========================================================================

%===========================================
% ----- Initialize Parameters for EEG ------
%===========================================
% Non-Default Filter Parameters
if nargin > 4 
            
    % ------ Row Vectors ------
    if size(FiltCoef,1) < size(FiltCoef,2)
        FiltCoef = FiltCoef.';
    end
    
    A  = FiltCoef(:,1);
    B  = FiltCoef(:,2);
    Cf = FiltCoef(:,3);
    
% Default Parameters 
else  
    A =  [0.150, 0.003, 0.180, 0.145, 0.120, 0.020, ...
          0.001, 0.001, 0.001, 0.001, 0.001, 0.001];
    B =  [0.075, 0.048, 0.110, 0.150, 0.190, 0.150, ...
          0.110, 0.075, 0.060, 0.050, 0.040, 0.030];
    Cf = [2.250, 5.630, 8.900, 11.54, 14.09, 16.78, ...
          19.79, 23.10, 26.70, 30.42, 34.37, 38.58];
end 

%==========================================
% ------------ Gaussian Mask --------------
%==========================================
radius = floor(3*Sigma);
sideSize = 1+2*radius;

Mask = zeros(1,sideSize);

den = 1 / (sqrt(2*pi*Sigma*Sigma));
for i = 1:sideSize
    x = i-radius-1;
    Mask(i) = den * exp( -0.5 * (x/Sigma) * (x/Sigma) );
end
%==============================================================
% ------ Initialize Time Series EEG Variable Assignment -------
%==============================================================
    %==============================================
    %======== Examine Data: Make Row Vector =======
    %==============================================
    nm = size(EEG_Signal);
    if nm(1) > 1
        EEG_Signal = EEG_Signal.';
    end

    % Correction Value for Convolution Padding
    GL_H = floor(length(Mask)/2);


%===============================================
%-------- Fourier Domain to Time Domain --------
%===============================================
N = round(length(EEG_Signal));
Fx = fft(EEG_Signal,N);
f = linspace(0, N-1,N).*F_s/N;
J = length(Cf);


    %=================================================
    %=============== Design Filter Bank ==============
    %=================================================
    WaveletsFreq = zeros(J,N);
    for j = 1:J
        WaveletsFreq(j,:) = exp(-A(j)*(f-Cf(j)).^2-B(j)*(f-Cf(j)).^4);
    end 

        %-------- Plot Filter Design in Frequency Domain --------
        if PlotFilters == 1 
            figure; 
            for j = 1:J
            plot(f,WaveletsFreq(j,:));
            hold on 
            end     
            plot(f,sum(WaveletsFreq))
            grid on 
            title('Wavelet Filter Bank Design')
            hold off
        end 

    %================================================
    %==== Filtering and Intensity Calculations ======
    %================================================
    
            % ------ Retrospective Analysis (Fourier Implementation) ------
            FreqFiltX = zeros(J, N);
            Y = zeros(J, N);
            Power = zeros(J, N);
            WL_H=floor(length(WaveletsFreq)/2);
            for j = 1:J
                TempWTA = fftshift(ifft(WaveletsFreq(j, :)));
                PowerTemp = conv(TempWTA, EEG_Signal);
                PowerTempB = 2 * abs(PowerTemp(WL_H + 1:end - WL_H));
                PowerTempC = conv(PowerTempB, Mask);
                Power(j, :) = PowerTempC(GL_H + 1:end - GL_H);
            end

end