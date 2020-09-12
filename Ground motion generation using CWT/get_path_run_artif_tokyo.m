p = parpool('local',1);
tic;
SubFolder = GetFolders('.\Target Spectrum');
eqname={'EW.txt','EW.txt'};%
parfor ii=1:length(SubFolder)
    disp(ii);
    temp_string=strsplit(SubFolder{1,ii},'\');
    cur_foldername=temp_string{length(temp_string)};  %CHB001_0
    % temp_this_station=strsplit(cur_foldername,'_');
    % this_station=char(temp_this_station(1));
    this_station=char(cur_foldername);
    dsname   = 'PredictedGM.txt';                 % name of spectrum file (2 columns
    eqfolder = '.\groundmotion';   % directory with the accelerograms 输出地震动目录
%     eqname  = 'gm_322.txt';             % name of earthquake file gm_322
    ArtifQuakeLetII(SubFolder{1,ii},dsname,eqfolder,eqname{1,1},this_station);
end
toc;
delete(p);



