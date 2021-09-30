function FCfoldname = ROItoFC(foldname,ROItable)
    %{    
    %compute Pearson's correlation between BOLD time series of different ROIs
    %converts to Z-transformed values, and provides names for each ROI
    %sort out the different parcellations to be examined separately
    
    inputs:
    1)foldname: string; name of folder containing ROI_...mat files, with pre-processed data
    2)ROItable: string; name of file with one column containing names of each ROI used, with .csv suffix
    
    requirements for contents of foldname:
    -Name: ROI_Subject0XX_Condition001, XX being subject number n=1,2,3...
    -ROIs must be identical to contents of ROItable
   

    output:
    -provides name of folder containing functional connectivity matrices
    -FCfolder: folder with .mat files (one for each participant),
    containing struct for different FC matrices

    
    %}
    %% load names for files and regions
    tname = readtable(ROItable);
    nos = what(foldname);
    fnames = nos.mat;

    %extract subject number
    snum = extractBetween(fnames,'Subject0','_Condition');

    %check for subjects with no pair in pre-processed files
    [num_occur, snum_id] = groupcounts(snum);
    nopair = snum_id(num_occur<2);
    disp(nopair)
    
    snum = unique(snum);%avoid repetition in loading files; get unique subject numbers
    
    %% create new folder with current date
    datestamp = erase(datestr(now, 'dd/mmm'),{'/'});
    FCfoldname = strrep('FCmatricesXX', 'XX',datestamp);
    mkdir(FCfoldname); 
    
    %% loop over files
    for subj = 1:length(snum)
        
        %% detect empty files
        fi1 = dir(strrep('FNAME/ROI_Subject0XX_Condition001.mat',{'FNAME','XX'},{foldname,snum{subj}}));
        if fi1.bytes == 0
            fprintf('\n No data: ROI_Subject0%s_Condition001.mat \n',snum{subj})
            continue
        end
        
        %% load matrix
        al = load(strrep('FNAME/ROI_Subject0XX_Condition001.mat',{'FNAME','XX'},{foldname,snum{subj}}));

        if length(al.names) == length(tname.ROIName)
            al.names = [al.names',tname.ROIName];
        elseif length(al.names) == length(tname.ROIName) - 4
            al.names = [al.names',tname.ROIName(1:end-4)];
        end
        
        al.xyz = al.xyz';%transpose in case of errors for merging with tables
        
        %% compute FC for all except the GM/WM/CSF and final few regressors

        %convert the cell-stored data into matrix of time series
        timeseries = cat(2,al.data{4:end-4}); %contatonate horizontally the series
        [fc.rPearson,~,~,~] = corrcoef(timeseries);

        %remove self correlations
        for i = 1:size(fc.rPearson,1)
            fc.rPearson(i,i) = NaN;
        end
        
        fc.Z = atanh(fc.rPearson);%compute fisher's Z
        fc.names = al.names(4:end-4,:);%add region names
        fc.xyz = al.xyz(4:end-4);%add region coordinates; note it is one dimensional        

        %% get matrix for each parcellation

        fc.S100.rPearson = fc.rPearson(1:100,1:100);
        fc.S100.Z = fc.rPearson(1:100,1:100);
        fc.S100.names = fc.names(1:100,:);
        fc.S100.xyz = fc.xyz(1:100);       

        fc.S400.rPearson = fc.rPearson(101:500,101:500);
        fc.S400.Z = fc.rPearson(101:500,101:500);
        fc.S400.names = fc.names(101:500,:);
        fc.S400.xyz = fc.xyz(101:500); 

        fc.T50.rPearson = fc.rPearson(501:550,501:550);
        fc.T50.Z = fc.rPearson(501:550,501:550);
        fc.T50.names = fc.names(501:550,:);
        fc.T50.xyz = fc.xyz(501:550);

        fc.T54.rPearson = fc.rPearson(551:604,551:604);
        fc.T54.Z = fc.rPearson(551:604, 551:604);
        fc.T54.names = fc.names(551:604,:);
        fc.T54.xyz = fc.xyz(551:604);

        %{
        %sanity check for if boundaries are correct
        fprintf('expect TSchaefer2018_100Parcels...cluster001, get: %s \n',r.S100.names{1})%insert name and create new line
        fprintf('expect TSchaefer2018_100Parcels...cluster100, get: %s \n',r.S100.names{end})
        fprintf('expect TSchaefer2018_400Parcels...cluster001, get: %s \n',r.S400.names{1})
        fprintf('expect TSchaefer2018_400Parcels...cluster400, get: %s \n',r.S400.names{end})

        fprintf('expect Tian_Subcortex_S3_3T.cluster001, get: %s \n',r.T50.names{1})
        fprintf('expect Tian_Subcortex_S3_3T.cluster050, get: %s \n',r.T50.names{end})    
        fprintf('expect Tian_Subcortex_S4_3T.cluster001, get: %s \n',r.T54.names{1})
        fprintf('expect Tian_Subcortex_S4_3T.cluster054, get: %s \n',r.T54.names{end})
        %}

        %% add a version for left and right next to each other
        fc.S100RL = structfun(@hemisphere_merge, fc.S100,'UniformOutput',false);
        fc.S400RL = structfun(@hemisphere_merge, fc.S400,'UniformOutput',false);
        fc.T50RL = structfun(@hemisphere_merge, fc.T50, 'UniformOutput',false);
        fc.T54RL = structfun(@hemisphere_merge, fc.T54, 'UniformOutput',false);

        %% combine cortical and subcortical
        %Schaefer 400 with Tian scale IV for 54 ROIs   
        fc.S400T54.rPearson = fc.rPearson([101:500,551:604],[101:500,551:604]);
        fc.S400T54.Z = fc.rPearson([101:500,551:604],[101:500,551:604]);
        fc.S400T54.names = fc.names([101:500,551:604],:);
        fc.S400T54.xyz = fc.xyz([101:500,551:604]);

        %% take out DMN-FPCN
        fc.DNCN100 = structfun(@takeDNCN, fc.S100,'UniformOutput',false);
        fc.DNCN400 = structfun(@takeDNCN, fc.S400,'UniformOutput',false);

        %% store in file

        %note: use replace, not strrep
        save(replace('FNAME/FC_SubjectXX.mat',{'FNAME','XX'},{FCfoldname,snum{subj}}),'fc')
        
        fprintf('completed Subject%s \n',snum{subj})
    end
end

function D = hemisphere_merge(x)
    xhalf = length(x)/2; 
    %disp(size(x,2))
    %A is left, B is right
    %treat name and connectivity matrix differently
    if size(x, 1) == size(x, 2)
        %shift all rows
        A = x(1:xhalf,1:end);%take all left
        B = x(xhalf+1:end,1:end);%take all right
        C = nan(size([A;B]));%create empty to put 
        C(1:2:end) = A;%assign to odd rows
        C(2:2:end) = B;%assign to even rows
        
        %shift all the columns
        D = nan(size(C));
        D(:,1:2:end) = C(:,1:xhalf);
        D(:,2:2:end) = C(:,xhalf+1:end);
        
    else
        %used to shuffle the name of regions
        A = x(1:xhalf,:);
        B = x(xhalf+1:end,:);
        D = cell(size([A;B]));
        %disp(size(A))
        %disp(size(B))
        %disp(size(D))
        D(1:2:end) = A;
        D(2:2:end) = B;
    end
end

function aus = takeDNCN(x)
    %gets the DNCN out of the bigger network
    %operate on the S100 and S400 bits of structure
    if length(x) == 100
        if size(x, 1) == size(x, 2)
            aus = x([34:50,81:100],[34:50,81:100]);
        elseif size(x,1) == 1 || size(x,2) == 1%if it is the xyz thing
            aus = x([34:50,81:100]);
        else
            aus = x([34:50,81:100],:);
        end
    elseif length(x) == 400
        if size(x, 1) == size(x, 2)
            aus = x([127:200,332:400],[127:200,332:400]);
        elseif size(x,1) == 1 || size(x,2) == 1%if it is the xyz thing
            aus = x([127:200,332:400]);
        else
            aus = x([127:200,332:400],:);
        end
    else
        disp('input not S100 or S400')
    end
end