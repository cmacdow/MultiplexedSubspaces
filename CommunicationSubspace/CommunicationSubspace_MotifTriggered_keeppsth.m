function CommunicationSubspace_MotifTriggered_keeppsth(motif,cur_rec,muaflag)
%Camden MacDowell - timeless
%Runs through a pipeline of CCA analyses for a given motif
%INPUTs
%win = [-5, 15] (default). The window around each motif spike to use.
%negative values are before onset. positive after. 
%EphysPath; the path of the ap_opts.mat file
%motif_fits; paths to the BasisMotifFits for a given mouse. 

if ispc
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Widefield_Imaging_Analysis'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\GenUtils'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Ephys'));
    savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\CommunicationSubspace_keeppsth\';
else
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/GenUtils'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Ephys'));
    savedir = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/CommunicationSubspace_keeppsth/';
end
tic
win=[-2 11]; %hardcoded write now. 
%starting
fprintf('Working on rec %d motif %d',cur_rec,motif);
%% Gathering Data
[rec_name,~,~,EphysPath,motif_fits] = LoadDataDirectories(cur_rec);

%load ephys data
[st_mat,~,st_depth] = LoadSpikes(EphysPath,'bindata',1,'offset',15,'mua',muaflag,'depth_type','probe'); 
% st_norm = cellfun(@(x) x./std(x),st_mat,'UniformOutput',0);
st_norm = st_mat;
st_norm = cellfun(@(x) x(1:2:end,:)+x(2:2:end,:),st_norm,'UniformOutput',0);

%get the anatomical locations
neu_area = LoadEphysAnatomicalLocations(EphysPath,st_depth);
neu_area = cat(2,neu_area{:});

%get the motif onsets for all recordings
[motif_onset,~] = CompileMotifOnsets(motif_fits); %return 'chunk_names' if you want to confirm ordering is correct
motif_onset = cellfun(@(x) floor(x/2), motif_onset,'UniformOutput',0);

if motif > numel(motif_onset) %get null periods that of window length    
    fprintf('\n\t Running on NULL')
    [~,trig_st] = ParseNullPeriods([],st_norm,motif_onset,win,10);
else %parse motif onsets
    fprintf('\n\t Running on Motifs')    
    [~,trig_st] = ParseByOnset([],st_norm,motif_onset,win,motif);
end

%parse activity per parent region 
[area_val, area_label] = ParseByArea(cat(1,trig_st{:}),neu_area,'general');

%remove neurons that fire less than 0.5spike/sec on average across trials
% [~, inactive_idx] = RemoveInactiveNeurons(area_val, 0.25/7.5);

%remove the rare edge case where a motif begins at the start (no baseline)
area_val = RemoveEdgeTrials(area_val);

%clean up areas %third input is the min # of spikes to keep area
[area_val, area_label] = CleanUpAreas(area_val, area_label, 10); 

%% Main
fprintf('USING RRR')
paired_areas = 1:numel(area_label);
%preallocate variables
rrr_b = cell(size(paired_areas,1),1);
rrr_V = cell(size(paired_areas,1),1);
rrr_B = cell(size(paired_areas,1),1);
cvl_rrr = cell(size(paired_areas,1),1);
rrr_B_xval = cell(size(paired_areas,1),1);
cvl_ridge = cell(size(paired_areas,1),1);
grouping = cell(size(paired_areas,1),1);
cvl_rrr_unique = cell(size(paired_areas,1),1);
contribution = cell(size(paired_areas,1),1);    
for i = 1:size(paired_areas,2)
    fprintf('\n\tWorking on subspace pairing %d of %d',i,size(paired_areas,2));
    idx = strcmp(area_label,area_label{paired_areas(i)});
    x = cat(1,area_val{idx==0,:});
    y = area_val{strcmp(area_label,area_label{paired_areas(i)}),:};

    %normalize to baseline
    x = normalizeToBaseline(x,1:2,'meansubtract');
    y = normalizeToBaseline(y,1:2,'meansubtract');

    %use post stimulus
    x = x(:,3:end,:);
    y = y(:,3:end,:);

    %while here, go ahead and parse beta by area
    grp = arrayfun(@(n) ones(size(area_val{n},1),1)*n,find(idx==0),'UniformOutput',0);
    grp = cat(1,grp{:});
    grouping{i} = grp;
    [rrr_b{i},rrr_B{i},rrr_V{i},cvl_rrr{i},cvl_ridge{i},cvl_rrr_unique{i},contribution{i},rrr_B_xval{i}] = RRR_full(x,y,0);
end %subspace identification loop
%% save off data
save([savedir,sprintf('%sregRRR_muaflag%d_GROUPEDmotif%d',rec_name,muaflag,motif)])
toc
fprintf('\ndone')
 
    
end %function end



