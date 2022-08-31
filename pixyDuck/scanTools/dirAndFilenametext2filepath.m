function [filepath,filename]=dirAndFilenametext2filepath(directory,filenametext)
        fileList=listFilesWithExtension('.nd2','directoryPath',directory,'ignoreFilesBeginningWithPeriod',true);
        
        filesInDirThatContainFilenametext=fileList(contains(fileList,filenametext));
        if length(filesInDirThatContainFilenametext)~=1
            error('more than one file was found in directory %s that contains filenametext %s',directory,filenametext)
        end
        filename=filesInDirThatContainFilenametext{1};
        filepath=fullfile(directory,filesep,filename);

end