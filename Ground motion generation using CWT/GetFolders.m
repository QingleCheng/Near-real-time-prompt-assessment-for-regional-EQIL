function [SubFolders] = GetFolders(ParentFolder)
%GetFolders
% ��������Ϊ��ȡ���ļ������������ļ��е�·��
% ����������ΪParentFolder�����ļ���·����eg: 'D:\Program Files'
% ���������ΪSubFolders�����ļ���·����Ϊһ��Ԫ�����飬eg: {'D:\Program Files\FileZilla FTP Client\docs'}

SubFolderNames = dir(ParentFolder);
SubFolder(1).SubFolderName='';
for i=1:length(SubFolderNames)
        if( isequal( SubFolderNames( i ).name, '.' )||...
        isequal( SubFolderNames( i ).name, '..')||...
        ~SubFolderNames( i ).isdir) % �������Ŀ¼������
            continue;
        end
        if strcmp(SubFolderNames( i ).name,'�������')~=true && strcmp(SubFolderNames( i ).name,'workspace')~=true...
            && strcmp(SubFolderNames( i ).name,'my_groundmotion')~=true
            SubFolder(i).SubFolderName = fullfile( ParentFolder, SubFolderNames( i ).name );
        end
end

temp = {SubFolder.SubFolderName};
idx = cellfun(@(x)~isempty(x),temp,'UniformOutput',true); % ����cellfun�����õ�Ԫ�����������зǿ�Ԫ�ص��±�
SubFolders = temp(idx);
 
end