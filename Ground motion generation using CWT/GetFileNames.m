function [FileNames] = GetFileNames(Path,Format)
% GetFileNames
% 函数的功能为获得某一路径下，某种格式所有文件名
% 函数的输入1为Path,要获取的路径。eg: 'D:\Program Files\FileZilla FTP Client\docs'
% 函数的输入2为Format，要获取路径的文件格式。eg: '*.txt','*.docx','*.png'
 
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