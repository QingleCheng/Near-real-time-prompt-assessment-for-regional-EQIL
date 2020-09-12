function [FileNames] = GetFileNames(Path,Format)
% GetFileNames
% �����Ĺ���Ϊ���ĳһ·���£�ĳ�ָ�ʽ�����ļ���
% ����������1ΪPath,Ҫ��ȡ��·����eg: 'D:\Program Files\FileZilla FTP Client\docs'
% ����������2ΪFormat��Ҫ��ȡ·�����ļ���ʽ��eg: '*.txt','*.docx','*.png'
 
fileFolder=fullfile(Path);
dirOutput=dir(fullfile(fileFolder,Format));
FileNames1={dirOutput.name};
j=1;
for i=1:length(FileNames1)   
    if strcmp(FileNames1{1,i}(1:5),'Comp-')~=true
        FileNames{1,j}=FileNames1{1,i};
        j=j+1;
    end
end

 
end