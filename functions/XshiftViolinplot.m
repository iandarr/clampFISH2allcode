function vs=XshiftViolinplot(vs,XshiftArray)
% shift the object from violinplot function
% Ian Dardani
% 
% violinplot is from:
%       Bechtold, Bastian, 2016. Violin Plots for Matlab, Github Project
%       https://github.com/bastibe/Violinplot-Matlab, DOI: 10.5281/zenodo.4559847

numViolins=length(vs);

if length(XshiftArray)~=numViolins
    error('second input XshiftArray is length %i but should be same as the number of violins %i', length(XshiftArray),numViolins)
end
            for iViolin=1:numViolins
                
                XDistanceToAdd=XshiftArray(iViolin);
                
                % move the violin
                violinVertices=vs(iViolin).ViolinPlot.Vertices;
                vs(iViolin).ViolinPlot.Vertices=[violinVertices(:,1) + XDistanceToAdd, violinVertices(:,2)];
                
                % move the scattered data
                vs(iViolin).ScatterPlot.XData=vs(iViolin).ScatterPlot.XData + XDistanceToAdd;
                
                % move the box plot
                vs(iViolin).BoxPlot.XData=vs(iViolin).BoxPlot.XData + XDistanceToAdd;
                
                % move the whisker plot
                vs(iViolin).WhiskerPlot.XData=vs(iViolin).WhiskerPlot.XData + XDistanceToAdd;
                
                % move the median plot
                vs(iViolin).MedianPlot.XData=vs(iViolin).MedianPlot.XData + XDistanceToAdd;
                
                % move the NotchPlots
                for iNotch=1:length(vs(iViolin).NotchPlots)
                    vs(iViolin).NotchPlots(iNotch).XData=vs(iViolin).NotchPlots(iNotch).XData + XDistanceToAdd;
                end
                
                % move meanplots
                vs(iViolin).MeanPlot.XData=vs(iViolin).MeanPlot.XData + XDistanceToAdd;
                
            end


end