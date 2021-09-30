function BrainstemPipeline
  
    %creates all the data tables necessary for statistical analysis in R
    %input: pre-processed ROI files for 81 participants

    % set initial input
    foldname = 'LegacyROIdata';%folder containing pre-processed time series
    hpname = 'Highpass';
    ROItable = 'ParcellationROINames_Schaefer.csv';
    behavtable = 'Legacy_clean_for_Eddie.csv';
    
    %convert time series to FC connectivity matrix files
    FCfold = ROItoFC(foldname,ROItable);
    %FCfold contains subject-specific FC matrices files, e.g., FC_Subject01.mat
    
    %extract brainstem nuclei specific info; do BOLD for highpass data
    getBrainstem_FCandBOLD(FCfold,foldname)%second arg is BOLD folder
    highpassBrainstemBOLD(hpname)
    %output:
    %1) folder: BrainstemFC (include n_region by n_nuclei tables, all Z values)
    %2) brainstemBOLDs.mat: contains bold signal, struct (n_nucleix81x300)
    %also: brainstemBOLDs_highpass.mat (same as above)
    
    %create MRN/DRNfcTable.csv (81 by number of nodes)
    %input from BrainstemFC folder
    genRapheFCtable
    %this is only done for 454 nodes
    
    %use...
    %brainstemBOLD.mat and Highpass data to calculate complexity & fALFF
    %also takes in the nuclei fc tables from previous fun for pos & neg
    %strength
    brainstemmets('MRN',behavtable)
    brainstemmets('DRN',behavtable)
    %output: NucleiMetrics_insertnucleusname.csv

    %combine receptor density data, behavioural data, and FC data
    medianFCgroup(behavtable)
    %output:
    %Nuclei_FCwithReceptor.csv etc.
    

end