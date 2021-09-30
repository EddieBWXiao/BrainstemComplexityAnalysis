function brainstem_measures(nucleus,brainstemBOLDs,highpass_brainstemBOLDs)
%{
%based on BOLD signal from different brainstem ROIs, compute complexity
measures, neuronal variability, and fALFF 

%pre-requisits: 
%function sampen, by Víctor Martínez-Cagigal 
%link: https://uk.mathworks.com/matlabcentral/fileexchange/69381-sample-entropy
%function calc_lz_complexity, by Quang Thai
%link: https://uk.mathworks.com/matlabcentral/fileexchange/38211-calc_lz_complexity


%input: 
%nucleus: string, such as DRN, MRN, VTA or All5HT; must match table
variable names from the brainstemBOLDs file
%brainstemBOLDs: string for name of mat file (without .mat suffix), which
contains table of n_subject x n_nuclei, with each cell containing the BOLD
time series for each subject; 
%highpass_brainstemBOLDs: same as brainstemBOLDs, but for time series
without highpass filtering
%IMPORTANT: the two table above should be named as t

%output:
%n_subject x n_measure table (as file NucleiMetrics_%s.csv)
%columns: SD, sample entropy, Lempel-Ziv complexity, fALFF


%}


    %% prepare
    load(brainstemBOLDs,'t');
    hp = load(highpass_brainstemBOLDs);
    theChosen = sprintf('%s_BOLD',nucleus);
    close all
    
    %% compute the intrinsic measures
    
    %preallocate
    sigm = nan(height(t),1);
    SEn = nan(height(t),1);
    mssd = nan(height(t),1);
    LZCp = nan(height(t),1);
    fALFF = nan(height(t),1);
    
    for i = 1:height(t)
        TS = t.(theChosen){i};
        
        sigm(i) = std(TS);%compute standard deviation
        
        m = 3;
        r = 0.6; %note: r is the percentage; sampen automatically multiplies it with SD 
        
        SEn(i) = sampen(TS,m,r,'chebychev');%compute sample entropy
        mssd(i) = MSSD(TS);%compute MSSD
        [~, LZCp(i)] = myLZC(TS);%compute Lempel-Ziv complexity
        
        %compute fALFF
        TSraw = hp.t.(theChosen){i};
        fALFF(i) = sqrt(bandpower(TSraw,0.5,[0.01,0.08]))/sqrt(bandpower(TSraw,0.5,[0,0.25]));
        
    end
    
    
    SD = sigm;
    d = table(SD, SEn, LZCp, fALFF);
    
    %% output
    writetable(d,sprintf('NucleiMetrics_%s.csv',nucleus))

end
function outdata = MSSD(x)
    %calculates the mean square successive difference
    %Von Neuman et al. (1941)
    n = length(x);
    outdata = 0;
    for i = 1:n-1
        outdata = outdata + (x(i+1)-x(i))^2;
    end
    outdata = outdata/(n-1);

end
function [exhau, prim] = myLZC(S)
    S = hilbert(S);%get hilbert
    meanS = mean(S);
    S = S >= meanS;%binarise
    exhau = calc_lz_complexity(S, 'exhaustive', 0);
    prim = calc_lz_complexity(S, 'primitive', 0);
end

