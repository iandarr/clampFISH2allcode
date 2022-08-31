function listOfFiles=listFilesWithExtension(extension,varargin)
% listFilesWithExtension outputs the file names in a directory with a given extension, where listOfFiles is a cell array of character vectors
% 
% Usage
% listOfFiles=listFilesWithExtension(extension)
% listOfFiles=listFilesWithExtension(extension,'directoryPath',directoryPath)
% listOfFiles=listFilesWithExtension(extension,'outputWithExt',false)
% listOfFiles=listFilesWithExtension(extension,'ignoreFilesBeginningWithPeriod',false)
%
% 
% Example:
% 
% listFilesWithExtension('.nd2') looks in the current directory and returns
%
%     {'20200922_170759_414__ChannelDAPI,YFP,CY3,A594,CY5_Seq0000.nd2'}
%     {'20200922_222102_814__ChannelDAPI,YFP,CY3,A594,CY5_Seq0000.nd2'}
% 
% listFilesWithExtension('.nd2','directoryPath','data/20X','outputWithExt',false) % lists files in directoryPath where the extensions are not included
% 
%     {'20200922_170759_414__ChannelDAPI,YFP,CY3,A594,CY5_Seq0000'}
%     {'20200922_222102_814__ChannelDAPI,YFP,CY3,A594,CY5_Seq0000'}
% 
% listFilesWithExtension(extension,'outputWithExt',false)
%
%
% ignoreFilesBeginningWithPeriod
%   defaults to true (if starts with period it's a hidden file and isn't returned)


p=inputParser;
p.addParameter('outputWithExt', true, @islogical);
p.addParameter('directoryPath', pwd, @isfolder)
p.addParameter('ignoreFilesBeginningWithPeriod', true, @islogical)
p.parse(varargin{:})

outputWithExt=p.Results.outputWithExt;
directoryPath=p.Results.directoryPath;
ignoreFilesBeginningWithPeriod=p.Results.ignoreFilesBeginningWithPeriod;

if ~isfolder(directoryPath)
    error("could not find directoryPath provided")
end

temp=dir(directoryPath);
dirContents={temp.name}';

assert(ischar(extension))
if ~startsWith(extension,'.')
    extension=['.',extension];
end

listOfFiles=cell('');
number_of_files_with_ext=0;
for ii=1:length(dirContents)
    [~,fileName,fileExt]=fileparts(dirContents{ii});
    if strcmp(fileExt,extension)
        
        number_of_files_with_ext=number_of_files_with_ext+1;
        
    if outputWithExt
        listOfFiles(number_of_files_with_ext,1)={[fileName,extension]};
    else
        listOfFiles(number_of_files_with_ext,1)={fileName};
    end
    end
    
end

if ignoreFilesBeginningWithPeriod
    listOfFiles=listOfFiles(~startsWith(listOfFiles,'.'));
end
end