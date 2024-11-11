function [hierarchicalArrangements] = createHierarchicalArrangements(X,Fs,sigma,PlotFilters,FiltCoef)
%-----Overview-------------------------
%   createHierarchicalArrangements: This approach designs the wavelets filters using a flat
%   gaussian function in the frequency domain and transforms the the EEG time series to the
%   frequency domain to filter the EEG time series. A gaussian filter is
%   then implemented for smoothing. It then finds the Hierarchical
%   Arrangements at each time sample based on the intensity of the wavelet
%   filter. 


% ------- Variable Inputs -------------
    % X: EEG Time Series (N x 1) row vector 
    % FS: Sampling Frequency
    % sigma:  The gaussian Mask implementation is designed to smooth the
    %            the power out. The varible is scalar value that is
    %            described by sigma. (Note: With a sampling fs=256 we have typically used a sigma around 60)
    % FiltCoef: This is a binary response leaving the input empty or providing new design
    %           coefficients . If you want the option to use the default filter coefficients
    %           previously designed or a different set of coefficients. The
    %           design coefficient input should be set by f x 3 dimensions,
    %           where f is the number of filters and the columns represent
    %           parameters a,b,Cf. 
    % PlotFilters: A binary response, where the value of 1 would plot the
    %              designed filters and a value of 0 or leaving the input empty will not
    %              plot the filters. 
    
% ------- Variable Output -------------
    % hierarchicalArrangements: The output is the specfic hierarchical arrangements for that sample in a f-based 
    %                       number system. The variable dimensions are N x 1. 
    
%=====================================================================================  
% NOTES  ASSUMPTIONS: The filter design was based on minimizing the
% standard deviation of the plateaut value of the filters, while
% maintaining the below EEG frequency cutoffs. The cutoff values were aimed
% to achieve a (1/e) where every other adjacent filter is orthogonal. 


%------------------------------------
% EGG Frequency Bands 
%------------------------------------
     % Frequency Bands 
       % Delta 0-3.5 Hz
       % Theta 4-7 Hz
       % Alpha 8-12 Hz 
       % Low Beta 13-15
       % Mid Beta 15-18 Hz
       % High Beta 18-40
       % Gamma 32-100 Hz    

%======================================================================================================================
% Contributions: The Respiration Complexity Code was designed and written by Matthew Demas and Nicholas J. Napoli at njn5fg@virginia.edu.
%=======================================================================================================================

%==========================================
% ------Initialize Parameters for EEG---------
%==========================================


        if nargin>4  % Non-Default Filter Parameters
                    
                         %------Column Vectors-----------
                    if size(FiltCoef,1)>size(FiltCoef,2)
                        A=FiltCoef(:,1);
                        B=FiltCoef(:,2);
                        Cf=FiltCoef(:,3);
                    else %-------Row Vectors ------------
                        FiltCoef=FiltCoef';
                        A=FiltCoef(:,1);
                        B=FiltCoef(:,2);
                        Cf=FiltCoef(:,3);
                    end 
            
        else  % Default Parameters 
            A=[.15, .003,   .18,   .145,  .12,   .02,  .001,  .001,  .001, .001,  .001,  .001];
            B=[.075, .048,   .11,   .15,   .19,   .15,   .11,  .075,  .06,  .05,   .04,   .03];
            Cf=[2.25, 5.63, 8.9 , 11.54, 14.09, 16.78, 19.79,23.1,  26.7, 30.42, 34.37, 38.58];
        end 

%==========================================
%------------- Gaussian Mask---------------
%==========================================
radius = floor(3*sigma);
sideSize = 1+2*radius;

Mask = zeros(1,sideSize);

den = 1 / (sqrt(2*pi*sigma*sigma));
for i=1:sideSize
  x = i-radius-1;
  Mask(i) = den * exp( -0.5 * (x/sigma) * (x/sigma) );
end
%==============================================================
%----------Initialize Time Series EEG Variable Assignment------
%==============================================================
                %==============================================
                % Examine Data: Make Row Vector
                %==============================================
                nm = size(X);
                if nm(1) > 1
                    X = X.';
                end
                GL_H=floor(length(Mask)/2); % Correction Value for Convolution Padding


%===============================================
%-------- Fourier Domain to Time Domain---------
%===============================================
    N = round(length(X));
    Fx = fft(X,N);
    f = linspace(0, N-1,N).*Fs/N;
    J=length(Cf);


%=================================================
%=============== Design Filter Bank ==============
%=================================================
WaveletsFreq=zeros(J,N);
for j=1:J
WaveletsFreq(j,:)=exp(-A(j)*(f-Cf(j)).^2-B(j)*(f-Cf(j)).^4);
end 
  

    %-------- Plot Filter Design in Frequency Domain--------
    if PlotFilters == 1 
        figure; 
        for j=1:J
        plot(f,WaveletsFreq(j,:));
        hold on 
        end     
        plot(f,sum(WaveletsFreq))
        grid on 
        title('Wavelet Filter Bank Design')
        hold off
    end 

%================================================
%====Filtering and Intensity Calculations =======
%================================================

        % Retrospective Analysis (Fourier Implementation)
                FreqFiltX = zeros(J,N);
                Y = zeros(J,N);
                Power = zeros(J,N);
                for j = 1:J
                    FreqFiltX(j,:) = Fx.*WaveletsFreq(j,:);
                    Y(j,:) = ifft(FreqFiltX(j,:));   
                    TempA=abs(Y(j,:)); 
                    TempB=abs(ifft(FreqFiltX(j,:).*f*2*pi*1i)/(2*pi*Cf(j)));
                    TempC=conv((TempA+TempB),Mask);
                    Power(j,:)=TempC(GL_H+1:end-GL_H);
                end
%================================================
%========= Hierarchical Arragements =============
%================================================

    n=length(Power(:,1)); 
    ind_mat = zeros(size(Power));
    for ii=1:size(ind_mat,2)
        [~,~,ind_mat(:,ii)]=unique(Power(:,ii),'first');
    end
    hierarchicalArragments = n.^(0:n-1) * (ind_mat-1);

end

