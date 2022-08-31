function changeFilePathInMatFile(dataMatFilepath,newFilepathDir)
% changeFilePath(dataMatFilepath,newFilepathDir)
%
% Example:
%
%   changeFilePath('session1Folder/data001.mat','session7folder')
%
matObj = matfile(dataMatFilepath,'Writable',true);

objects=matObj.objects;

numObjects=length(objects);

for i=1:numObjects
   nodeLabels=objects(i).graph.labels;
   
   isNodeWithFilepath=~any([contains(nodeLabels,'imageObject');endsWith(nodeLabels,':Spots');endsWith(nodeLabels,':threshQC');endsWith(nodeLabels,'Proc');endsWith(nodeLabels,'Fitted')],1);
   %isNodeWithFilepath=ismember(nodeLabels,{'dapi','gfp','tmr','alexa','cy','nir','trans'});
   
   for ii=1:length(nodeLabels)
      if isNodeWithFilepath(ii)
          
            objects(i).graph.nodes{ii}.data.dirPath=newFilepathDir;
            
      elseif ~strcmp(nodeLabels{ii},'imageObject') % then you are a dependent node and need needsUpdate=1
          
            objects(i).graph.nodes{ii}.data.needsUpdate=true;
          
      else % then this is the objectHandle
      end
       
   end
   
   
end

matObj.objects=objects;

end