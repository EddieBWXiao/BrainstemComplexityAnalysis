function t = getBrainstem_FCandBOLD(FCfold,BOLDfold)
    %takes in FC and raw BOLD files
    %extracts specific FC values for Raphe (DRN and MRN), ROI-to-ROI
    
    %{
    use:
    1) prepare for FC analysis (each region as rows, nuclei as columns; one file per subj)
    2) prepare for spectral analysis (subject as rows, each nucleus's time series as
    columns) => brainstemBOLDs.mat
    
    input:
    -FCfold: folder containing all FC files (output of ROItoFC)
    -BOLDfold: folder containing all original time series files, e.g., LegacyROIdata
    
    output:
    -t: table, n_subject x n_nuclei
    -also creates BrainstemBOLDs.mat file
    %}
    
    
    %% load data
    
    %raw bold
    nos = what(BOLDfold);
    fnames = nos.mat;
    %extract that number 
    snum = extractBetween(fnames,'Subject0','_Condition');
    %find those with no pair
    [num_occur, snum_id] = groupcounts(snum);
    nopair = snum_id(num_occur<2);
    disp(nopair)
    snum = unique(snum);%avoid that paired repetition!!
    
    %FC by region for each subject
    mkdir('BrainstemFC')
    
    %% iterate across subjects; append below
    %% preallocate
    %for the time series
    DRN_BOLD = cell(length(snum),1);
    MRN_BOLD = cell(length(snum),1);
    All5HT_BOLD = cell(length(snum),1);
    VTA_BOLD = cell(length(snum),1);
    LC_BOLD = cell(length(snum),1);
    Subject = nan(length(snum),1);
    
    iteration = 0;
    
    %add to big table
    t = table(Subject, DRN_BOLD,MRN_BOLD,All5HT_BOLD, VTA_BOLD, LC_BOLD);

    for subj = 1:length(snum)
        iteration = iteration + 1;
        
        %% skip faulty data
        %{
        exclu = [4,6,8,10,12,18,28,29,30,37,38,63,64,66,69];
        if any(exclu == subj)
            continue
        end
        %}
        %% get file & know indices
        load(replace('YYY/FC_SubjectXX.mat',{'YYY','XX'},{FCfold,snum{subj}}), 'fc')
        al = load(strrep('YYY/ROI_Subject0XX_Condition001.mat',{'YYY','XX'},{BOLDfold,snum{subj}}));
        %separate the name columns
        fc.shortn = fc.names(:,2);
        fc.longnames = fc.names(:,1);
        
        %% construct table for the FC file; extract nuclei ROI-to-ROI connectivity
        %DRN = nan(length(fc.Z),1);
        %MRN = nan(length(fc.Z),1);
        %All5HT = nan(length(fc.Z),1);
        %VTA = nan(length(fc.Z),1);
        %LC = nan(length(fc.Z),1);
        
        %use columns
        DRN = fc.Z(:,contains(fc.names,'AAN_DR'));
        MRN = fc.Z(:,contains(fc.names,'AAN_MR_'));
        All5HT = fc.Z(:,contains(fc.names,'Serotonergic_Whole'));
        VTA = fc.Z(:,contains(fc.names,'AAN_VTA'));
        LC = fc.Z(:,contains(fc.names,'AAN_LC'));
        
        %get region names
        RegionName = [fc.shortn(1:500);fc.longnames(501:end)];
        
        fct = table(RegionName,DRN,MRN,All5HT,VTA,LC);
        
        
        writetable(fct, replace('BrainstemFC/NucleiFC_SubjectXX.csv','XX',snum{subj}));
        
        %% extract BOLD
        t.DRN_BOLD{subj} = al.data{contains(al.names,'AAN_DR')};
        t.MRN_BOLD{subj} = al.data{contains(al.names,'AAN_MR_')}; %do NOT confuse with mRF       
        t.All5HT_BOLD{subj} = al.data{contains(al.names,'Serotonergic_Whole')};
        t.VTA_BOLD{subj} = al.data{contains(al.names,'AAN_VTA')};
        t.LC_BOLD{subj} = al.data{contains(al.names,'AAN_LC')};        
        t.Subject(subj) = str2num(snum{subj});
        
        fprintf('completed Subject%s \n',snum{subj})


        
    end
    save('brainstemBOLDs.mat','t')

end