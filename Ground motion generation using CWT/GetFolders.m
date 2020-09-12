function [SubFolders] = GetFolders(ParentFolder)
%GetFolders
% 函数功能为获取父文件夹下所有子文件夹的路径
% 函数的输入为ParentFolder：父文件夹路径。eg: 'D:\Program Files'
% 函数的输出为SubFolders：子文件夹路径。为一个元胞数组，eg: {'D:\Program Files\FileZilla FTP Client\docs'}

SubFolderNames = dir(ParentFolder);
SubFolder(1).SubFolderName='';
for i=1:length(SubFolderNames)
        if( isequal( SubFolderNames( i ).name, '.' )||...
        isequal( SubFolderNames( i ).name, '..')||...
        ~SubFolderNames( i ).isdir) % 如果不是目录则跳过
            continue;
        end
        if strcmp(SubFolderNames( i ).name,'计算程序')~=true && strcmp(SubFolderNames( i ).name,'workspace')~=true...
            && strcmp(SubFolderNames( i ).name,'my_groundmotion')~=true
            SubFolder(i).SubFolderName = fullfile( ParentFolder, SubFolderNames( i ).name );
        end
end

temp = {SubFolder.SubFolderName};
idx = cellfun(@(x)~isempty(x),temp,'UniformOutput',true); % 利用cellfun函数得到元胞数组中所有非空元素的下标
SubFolders = temp(idx);
 
end