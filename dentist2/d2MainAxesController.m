classdef d2MainAxesController < handle
    properties (Access = public)
        
        %Model objects
        spotTable
        scanObj
        maskObj
        nucleiObj
        
        %View
        viewObj
        
        %Other controllers
        threshCntrlr
        thumbCntrlr
        
        viewRect %Image coords not axes coords
        panRect
        cellViewRadius = 500
        zoomROI
        zoomStart
        zoomRect

        fixedZoom = false
        panStart
        channelIdx
        zoomMode = true
        scatterH %Not sure if we need these handle. 
        imageH
        spotScatterH
        
        nucleiScatterH
        maskH
        imagesInView
        nondapiChannelsInView
        dapiInView
        
        
        resizedImg %Not sure if we want to save these or regenerate
        resizedDapi
        
        numCells
    end
    
    methods
        function p = d2MainAxesController(view, scanObj, spotTable, maskObj, nucleiObj)
            p.viewObj = view;
            p.spotTable = spotTable;
            p.scanObj = scanObj;
            p.maskObj = maskObj;
            p.nucleiObj = nucleiObj;
            p.viewRect = [1, 1, min(p.scanObj.stitchDim, [25000 25000])];
            p.startup()
        end
        
        function startup(p)
            p.channelIdx = p.viewObj.channelPopup.Value;
            p.getNondapiChannelsInView()
            p.plotScatterMain();
            p.updateImageInView();
            spotChannelIdx=find(ismember(p.spotTable.spotChannels,p.nondapiChannelsInView{p.channelIdx}));
            if ~isempty(spotChannelIdx)
                p.numCells = height(p.spotTable.centroidLists{spotChannelIdx});
            else
                p.numCells=0;
            end
            %p.updateMainAxes()
        end
        
        function p = changeChannel(p, ~, ~)
            p.channelIdx = p.viewObj.channelPopup.Value;
            %Update centroid listbox
            p.updateCentroidListView();
            p.updateMainAxes();
            p.threshCntrlr.plotIntensityHistogram();
            p.threshCntrlr.plotIntensityThreshold();
        end
        
        function changeColormap(p, ~, ~)
            p.spotTable.paletteIdx = p.viewObj.colormapPopup.Value;
            p.spotTable.updateExpressionColors();
            p.updateMainAxes();
        end
        
        function p = updateMainAxes(p, ~, ~)
            if isvalid(p.viewObj.mainAxes.Children)
                delete(get(p.viewObj.mainAxes, 'Children'));
            end
            
            if logical(p.viewObj.scatterCheckBox.Value)
                p.plotScatterMain();
            else
                p.showImage();
                p.overlaySpots();
                p.overlayNuclei();
                p.overlayMasks();
            end
        end
        
        function updateCentroidListView(p)
            spotChannelIdx=find(ismember(p.spotTable.spotChannels,p.nondapiChannelsInView{p.channelIdx}));
            if ~isempty(spotChannelIdx)
                p.viewObj.centroidList.String = string(p.spotTable.centroidLists{spotChannelIdx}.GroupCount);
                p.numCells = height(p.spotTable.centroidLists{spotChannelIdx}); %Although this value doesn't depend on the channel. Easy to put this here rather than everywhere where cell # can change. 
            else
                p.viewObj.centroidList.String={};
                p.numCells=0;
            end
        end
        
        function p = centroidSelected(p, ~, ~)
            if strcmp(get(p.viewObj.figHandle, 'SelectionType'), 'open') %Respond to double mouse click
                spotChannelIdx=find(ismember(p.spotTable.spotChannels,p.nondapiChannelsInView{p.channelIdx}));
                cellIdx = get(p.viewObj.centroidList, 'Value');
                cellPos = p.spotTable.centroidLists{spotChannelIdx}{cellIdx, {'x', 'y'}};
                p.viewRect = d2utils.getRectAroundPoint(cellPos, 2 * p.cellViewRadius, 2 * p.cellViewRadius, p.scanObj.stitchDim);
                set(p.viewObj.scatterCheckBox, 'Value', 0)
                p.updateImageInView();
                p.updateMainAxes();
                p.thumbCntrlr.overlayThumbnailRect();
            end
        end
        
        function plotScatterMain(p)
            spotChannelIdx=find(ismember(p.spotTable.spotChannels,p.nondapiChannelsInView{p.channelIdx}));
            if ~isempty(spotChannelIdx)
            %centroidsInView = p.spotTable.centroidTableInRect(p.channelIdx, p.viewRect);
            centroidsInView = p.spotTable.centroidTableInRect(spotChannelIdx, p.viewRect);
            centroidsInView = flipud(centroidsInView); %In order to plot cells with 0 expression in the back
            set(p.viewObj.mainAxes, 'XLim', [p.viewRect(2)-0.5,  p.viewRect(2)+p.viewRect(4)-0.5])
            set(p.viewObj.mainAxes, 'YLim', [p.viewRect(1)-0.5,  p.viewRect(1)+p.viewRect(3)-0.5])
            hold(p.viewObj.mainAxes, 'on')
            %Note that it is slightly slower to specify RGB colors with expression_color
            %instead of changing the colormap on the axes however, specifying RGB colors
            %ensures that the colors don't change with the view 
            p.scatterH = scatter(centroidsInView.y, centroidsInView.x,...
                30, centroidsInView.expression_color, 'filled',...
                'Parent', p.viewObj.mainAxes, 'HitTest','off');
            hold(p.viewObj.mainAxes, 'off')
            %colorbar(p.viewObj.mainAxes, 'Location', 'eastoutside','HitTest','off')
            else % this is not a spot channel
                delete(p.scatterH)
            end
        end
        
        function zoomInPressed(p, ~, ~)
            p.zoomMode = true;
            %iptPointerManager(p.viewObj.figHandle, 'disable');
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
            %Below are old functions that draw rectangle using matlab
            %ROI instead of mainAxesZoom fcn. Possibly more efficient but doesn't allow using
            %the console until roi completed. 
            %set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDownZoomMode}) %need a zoomMode that doesn't use 'normal' clicks
%             while p.zoomMode
%                 if p.fixedZoom
%                     p.zoomROI = drawrectangle(p.viewObj.mainAxes, 'Color', 'r', 'FaceAlpha', 0, 'FixedAspectRatio', true);
%                 else 
%                     p.zoomROI = drawrectangle(p.viewObj.mainAxes, 'Color', 'r', 'FaceAlpha', 0);
%                 end
%                 
%                 if isempty(p.zoomROI.Position) %if escape pressed
%                     break
%                 elseif all(p.zoomROI.Position(3:4) > 4) %minimum view size
%                     p.viewRect  = d2utils.coordToPixelRect(round(p.zoomROI.Position));
%                     delete(p.zoomROI) %Possibly unnecessary since Children deleted in p.updateMainAxes(); 
%                     p.updateImageInView();
%                     p.updateMainAxes();
%                     p.fixedZoom = false;
%                 end
%             end
        end
    
        function figWindowDown(p, ~, ~)
            if d2utils.pointInSideViewRect(p.viewRect, get(p.viewObj.mainAxes, 'CurrentPoint'))
                switch(p.getSelectionType)
                    %'normal' click was used for drawing rectangle before using
                    %rectangle roi
                    case 'normal'
                        if p.zoomMode
                            p.zoomStart = get(p.viewObj.mainAxes, 'CurrentPoint');
                            set(p.viewObj.figHandle, 'WindowButtonUpFcn', {@p.stopDragFcn});
                            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', {@p.mainAxesZoom});
                        else
                            p.panStart =  get(p.viewObj.mainAxes, 'CurrentPoint');
                            set(p.viewObj.figHandle, 'WindowButtonUpFcn', {@p.stopPan})
                            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', {@p.panView});
                        end
                    case 'open'
                        p.viewRect = [1, 1, min(p.scanObj.stitchDim, [25000 25000])];
                        p.updateImageInView();
                        p.updateMainAxes();
                        p.thumbCntrlr.overlayThumbnailRect();
                    case 'alt'
                        p.viewRect = d2utils.expandView2x(p.viewRect, p.scanObj.stitchDim);
                        p.updateImageInView();
                        p.updateMainAxes();
                        p.thumbCntrlr.overlayThumbnailRect();
                end
            end
        end
        
        function mainAxesZoom(p, ~, ~)
            if ~isempty(p.viewObj.zoomH)
                delete(p.viewObj.zoomH)
            end
            
            if p.fixedZoom
                p.getSelectedRectangleCoordsFixed()
            else
                p.getSelectedRectangleCoords()
            end
            
            if all(p.zoomRect(3:4) > 4) %min view size
                p.viewObj.zoomH = rectangle(...
                    'Position', p.zoomRect, ...
                    'EdgeColor', 'r', ...
                    'Parent', p.viewObj.mainAxes,...
                    'Hittest', 'off');
            end

        end
        
        function stopDragFcn(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', '');
            if ~isempty(p.zoomRect) && all(p.zoomRect > 0)
                p.viewRect  = d2utils.coordToPixelRect(p.zoomRect); % Should notify event "viewChange"
                %Update imagesInView
                p.updateImageInView();
                %Update main axes
                p.updateMainAxes();
                p.thumbCntrlr.overlayThumbnailRect();
                p.zoomRect = [];
                p.fixedZoom = false;
            end
            delete(p.viewObj.zoomH)
            set(p.viewObj.figHandle, 'WindowButtonUpFcn', '');
        end
        
        function getSelectedRectangleCoords(p)
            currentPoint = get(p.viewObj.mainAxes, 'CurrentPoint');
            p.zoomRect = round([min(p.zoomStart(1,1:2), currentPoint(1,1:2)),...
                abs(p.zoomStart(1,1:2) - currentPoint(1,1:2))]);
        end
        
        function getSelectedRectangleCoordsFixed(p)
            currentPoint = get(p.viewObj.mainAxes, 'CurrentPoint');
            sz = max(abs(p.zoomStart(1,1:2) - currentPoint(1,1:2)));
            p.zoomRect = round([min(p.zoomStart(1,1:2), currentPoint(1,1:2)), [sz sz]]);
        end
        
        function panViewPressed(p, ~, ~)
            if any(p.viewRect(3:4) < p.scanObj.stitchDim) %No need to pan if view is of entire scan
                p.zoomMode = false; 
                %iptPointerManager(p.viewObj.figHandle, 'enable');
                %iptSetPointerBehavior(p.viewObj.mainAxes, @(hfig, currentPoint) set(hfig, 'Pointer', 'hand'));
                set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
            else
                p.zoomMode = false; 
                set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            end
        end
                
        function stopPan(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', '');
            set(p.viewObj.figHandle, 'WindowButtonUpFcn', '')
        end
        
        function panView(p, ~, ~)
            currentPoint = get(p.viewObj.mainAxes, 'CurrentPoint');
            displacement = p.panStart(1,1:2) - currentPoint(1,1:2);
            p.viewRect = d2utils.updateViewPanning(p.viewRect, displacement, p.scanObj.stitchDim);
            p.updateImageInView;
            p.updateMainAxes;
            p.thumbCntrlr.overlayThumbnailRect();
        end
        
        
        function getNondapiChannelsInView(p)
            p.nondapiChannelsInView=p.scanObj.channels(~ismember(p.scanObj.channelTypes,'dapi'));
        end
        
        function p = updateImageInView(p)
            %p.imagesInView = cell(0, numel(p.spotTable.spotChannels));
            channelsInView=p.nondapiChannelsInView;
            p.imagesInView = cell(0, numel(channelsInView));
            if prod(p.viewRect(3:4)) < 4000001
                for i = 1:numel(channelsInView)
                    p.imagesInView{i} = p.scanObj.getImageRect(channelsInView{i}, p.viewRect);
                end
                p.dapiInView = p.scanObj.getDapiImage(p.viewRect);
            elseif prod(p.viewRect(3:4)) < 64000001
                for i = 1:numel(channelsInView)
                    p.imagesInView{i} = p.scanObj.getSmallImageRect(channelsInView{i}, p.viewRect);
                end
                p.dapiInView = p.scanObj.getSmallDapiImage(p.viewRect);
            else
                disp('View is too large to display image. Plotting scatter')
                set(p.viewObj.scatterCheckBox, 'Value', 1)
            end
        end
        
        function showImage(p)
           %Adjust contrast
            contrastIn = [get(p.viewObj.lowerContrastSlider, 'Value'), get(p.viewObj.upperContrastSlider, 'Value')];
            tmpIm = imadjust(p.imagesInView{p.channelIdx}, contrastIn, []);
            %Decide if overlay with DAPI?
            tmpRGB = cat(3, tmpIm, tmpIm, tmpIm+p.dapiInView);
            xlimits = [p.viewRect(2)-0.5,  p.viewRect(2)+p.viewRect(4)-0.5];
            ylimits = [p.viewRect(1)-0.5,  p.viewRect(1)+p.viewRect(3)-0.5];
            %axis(p.viewObj.mainAxes, [xlimits, ylimits], 'square')
            set(p.viewObj.mainAxes, 'XLim', xlimits)
            set(p.viewObj.mainAxes, 'YLim', ylimits)
            hold(p.viewObj.mainAxes, 'on')
            %p.imageH = imshow(tmpRGB, 'XData', xlimits, 'YData', ylimits, 'Parent', p.viewObj.mainAxes);
            p.imageH = imshow(tmpRGB, 'XData', [p.viewRect(2), p.viewRect(2)+p.viewRect(4)-1], 'YData', [p.viewRect(1),  p.viewRect(1)+p.viewRect(3)-1], 'Parent', p.viewObj.mainAxes);
            %p.imageH = imshow(tmpRGB,'Parent', p.viewObj.mainAxes);
            %axis fill
            %pbaspect auto
            set(p.viewObj.mainAxes, 'Visible', 'on')
            hold(p.viewObj.mainAxes, 'off')
        end
        
        function overlaySpots(p, ~, ~) %update this
            if logical(p.viewObj.spotsCheckBox.Value) && ~logical(p.viewObj.scatterCheckBox.Value)
                spotChannelIdx=find(ismember(p.spotTable.spotChannels,p.nondapiChannelsInView{p.channelIdx}));
                %mainAxesIdxsWithSpots=find(ismember(p.nondapiChannelsInView,p.spotTable.spotChannels));
                if ~isempty(spotChannelIdx)
                    channel=p.nondapiChannelsInView{p.channelIdx};
                    [spotsInView, ~] = p.spotTable.getValidSpotsInRect(channel, p.viewRect);
                    hold(p.viewObj.mainAxes, 'on')
                    scatSZ=max(10,round(30000/max(p.viewRect(3:4))));
                    p.spotScatterH = scatter(spotsInView.y, spotsInView.x, scatSZ, spotsInView.colors,...
                        'Parent', p.viewObj.mainAxes, 'HitTest','off');
                    hold(p.viewObj.mainAxes, 'off')
                else
                    delete(p.spotScatterH)
                end
            else
                delete(p.spotScatterH)
            end
        end
        
        function overlayNuclei(p, ~, ~)
            spotChannelIdx=find(ismember(p.spotTable.spotChannels,p.nondapiChannelsInView{p.channelIdx}));
            if logical(p.viewObj.centroidsCheckBox.Value) && ~logical(p.viewObj.scatterCheckBox.Value) && ~isempty(spotChannelIdx)
                outTableTmp = p.spotTable.centroidTableInRect(spotChannelIdx, p.viewRect);
                hold(p.viewObj.mainAxes, 'on')
                p.nucleiScatterH = scatter(outTableTmp.y, outTableTmp.x, 35, outTableTmp.colors, 'filled',...
                    'Parent', p.viewObj.mainAxes, 'HitTest','off');
                hold(p.viewObj.mainAxes, 'off')
            else
                delete(p.nucleiScatterH)
            end
        end
        
        function overlayMasks(p, ~, ~)
            masks = findobj(p.viewObj.mainAxes, 'type', 'images.roi.freehand');
            delete(masks);
            if logical(p.viewObj.masksCheckBox.Value) && ~logical(p.viewObj.scatterCheckBox.Value)
                maskTableTmp = p.maskObj.getChannelMasksInRect(p.viewRect, p.nondapiChannelsInView{p.channelIdx});
                maskIDs = unique(maskTableTmp.maskID);
                maskIDs(maskIDs == 0) = [];
                for i = 1:numel(maskIDs)
                    drawfreehand(p.viewObj.mainAxes, 'Position', maskTableTmp{maskTableTmp.maskID == maskIDs(i), {'y', 'x'}},...
                        'Color', 'red', 'InteractionsAllowed', 'none');
                end
                dapiChannel=p.scanObj.channels(ismember(p.scanObj.channelTypes,'dapi'));
                cellMasksTmp = p.maskObj.getChannelMasksInRect(p.viewRect, dapiChannel);
                maskIDs = unique(cellMasksTmp.maskID);
                maskIDs(maskIDs == 0) = [];
                for i = 1:numel(maskIDs)
                    drawfreehand(p.viewObj.mainAxes, 'Position', cellMasksTmp{cellMasksTmp.maskID == maskIDs(i), {'y', 'x'}},...
                        'Color', 'blue', 'InteractionsAllowed', 'none');
                end
            end
        end
        
        function scatterCallback(p, ~, ~)
            if p.viewObj.scatterCheckBox.Value == 0
                p.updateImageInView();
            end
            p.updateMainAxes();
        end
        
        function shuffleColorsInView(p, ~, ~)
            if ~logical(p.viewObj.scatterCheckBox.Value)
                spotChannelIdx=find(ismember(p.spotTable.spotChannels,p.nondapiChannelsInView{p.channelIdx}));
                
                outTableTmp = p.spotTable.centroidTableInRect(spotChannelIdx, p.viewRect);
                randomColors  = single(d2utils.distinguishable_colors(50));
                outTableTmp.colors = randomColors(randi(50, height(outTableTmp), 1),:);
                
                %Update centroid list with new colors
                idx = ismember(p.spotTable.centroidLists{spotChannelIdx}.nucID, outTableTmp.nucID);
                p.spotTable.centroidLists{spotChannelIdx}.colors(idx, :) = outTableTmp.colors;
                %Update nuclei table wtih new colors
                [idxA, idxB] = ismember(p.nucleiObj.nuclei.nucID, outTableTmp.nucID);
                idxB(idxB == 0) = [];
                p.nucleiObj.nuclei.colors(idxA, :) = outTableTmp.colors(idxB,:);
                %Update spots table with new colors
                [idxA, idxB] = ismember(p.spotTable.spots.nearestNucID, outTableTmp.nucID);
                idxB(idxB == 0) = [];
                p.spotTable.spots.colors(idxA, :) = outTableTmp.colors(idxB,:);
                %Update plots
                p.overlayNuclei;
                p.overlaySpots;
            end
        end
        %Note that when masking/adding/deleting spots and cells using the
        %functions below, the threshold histogram is not automatically
        %updated. It will update once the 'filter masked spots' button is
        %pushed. Could add an automatic update but I don't think it's
        %necessary at this point. 
        function addSpotMask(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.addCellButton, p.viewObj.deleteCellButton,...
                p.viewObj.maskCellButton, p.viewObj.deleteMaskButton], 'Enable', 'off')
            channel = p.nondapiChannelsInView{p.channelIdx};
            p.maskH = drawfreehand(p.viewObj.mainAxes, 'Parent', p.viewObj.mainAxes);
            %addlistener(p.maskH, 'DrawingFinished', @p.maskSpots);
            if ~isempty(p.maskH.Position) && isvalid(p.maskH) && polyarea(p.maskH.Position(:,1), p.maskH.Position(:,2)) > 3 %min mask area
                if p.viewObj.maskSpotsAllChannelsCheckBox.Value
                    % then mask for all channels
                    p.maskSpots(p.maskH,p.spotTable.spotChannels)
                else
                    p.maskSpots(p.maskH, channel)
                end
            else
                delete(p.maskH)
            end
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.addCellButton, p.viewObj.deleteCellButton,...
                p.viewObj.maskCellButton, p.viewObj.deleteMaskButton], 'Enable', 'on')
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
        end
        
        function maskSpots(p, roi, channelOrChannels)
            %tmpPoly = roi.Position;
            set(p.viewObj.masksCheckBox, 'Value', true)
            
            newMaskID = p.maskObj.addMaskLocalCoords(roi.Position, channelOrChannels); % can handle channel array
            delete(roi)
            channelOrChannels=cellstr(channelOrChannels); % force cell array
            for i=1:length(channelOrChannels)
                channel=channelOrChannels{i};
                p.spotTable.addNewMask(channel, newMaskID);
                p.spotTable.updateSpotStatus(channel);
                p.spotTable.updateCentroidList(channel);
            end
            p.updateCentroidListView();
            p.updateMainAxes();
        end
        
        function addCellMask(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.addCellButton, p.viewObj.deleteCellButton,...
                p.viewObj.maskSpotsButton, p.viewObj.deleteMaskButton], 'Enable', 'off')
            p.maskH = drawfreehand(p.viewObj.mainAxes, 'Parent', p.viewObj.mainAxes);
            if ~isempty(p.maskH.Position) && isvalid(p.maskH) && polyarea(p.maskH.Position(:,1), p.maskH.Position(:,2)) > 3 %min mask area
                p.maskCells(p.maskH)
            else
                delete(p.maskH)
            end
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.addCellButton, p.viewObj.deleteCellButton,...
                p.viewObj.maskSpotsButton, p.viewObj.deleteMaskButton], 'Enable', 'on')
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
        end
        
        function maskCells(p, roi)
            %tmpPoly = roi.Position;
            newMaskID = p.maskObj.addMaskLocalCoords(roi.Position, 'dapi');
            delete(roi)
            p.nucleiObj.addNewMask(newMaskID);
            p.spotTable.assignSpotsInRect(p.viewRect); 
            p.spotTable.updateAllSpotStatus();
            p.spotTable.makeCentroidList();
            p.updateCentroidListView();
            set(p.viewObj.masksCheckBox, 'Value', true)
            p.updateMainAxes();
        end
        
        function deleteMask(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.maskCellButton, p.viewObj.deleteCellButton,...
                p.viewObj.maskSpotsButton, p.viewObj.addCellButton], 'Enable', 'off')
            channel = p.nondapiChannelsInView{p.channelIdx};
            [x, y] = getpts(p.viewObj.mainAxes); %Simple but not interruptible. Can make WindowButtonDownFcn if want something interruptiblef.   
            if ~isempty(x)  
                ptsInView = d2utils.getPtsInsideView([x, y], p.viewRect);
                p.maskObj.removeMasksByLocalPoints(ptsInView, p.viewRect);
                p.nucleiObj.removeMasks();
                p.spotTable.removeMasks2(channel, p.viewRect);
                if p.nucleiObj.nucleiChanged %Kinda ugly. Should write something better. Could make separate buttons for cell masks and spot masks
                    p.spotTable.assignSpotsInRect(p.viewRect);
                    p.spotTable.updateAllSpotStatus();
                    p.spotTable.makeCentroidList();
                else
                    p.spotTable.updateSpotStatus(channel);
                    p.spotTable.updateCentroidList(channel);
                end
                p.nucleiObj.nucleiChanged = false;
                p.updateCentroidListView();
                p.updateMainAxes();
            end
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.maskCellButton, p.viewObj.deleteCellButton,...
                p.viewObj.maskSpotsButton, p.viewObj.addCellButton], 'Enable', 'on')
        end
        
        function addCells(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.maskCellButton, p.viewObj.deleteCellButton,...
                p.viewObj.maskSpotsButton, p.viewObj.deleteMaskButton], 'Enable', 'off')
            [x, y] = getpts(p.viewObj.mainAxes); %Simple but not interruptible. Can make WindowButtonDownFcn if want something interruptiblef.   
            if ~isempty(x)
                ptsInView = d2utils.getPtsInsideView([x, y], p.viewRect);
                p.nucleiObj.addCell(ptsInView(:,2), ptsInView(:,1));
                p.spotTable.assignSpotsInRect(p.viewRect);
                p.spotTable.updateAllSpotStatus();
                p.spotTable.makeCentroidList();
                p.updateCentroidListView();
                p.updateMainAxes();
            end
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.maskCellButton, p.viewObj.deleteCellButton,...
                p.viewObj.maskSpotsButton, p.viewObj.deleteMaskButton], 'Enable', 'on')
        end
        
        function deleteCells(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.maskCellButton, p.viewObj.addCellButton,...
                p.viewObj.maskSpotsButton, p.viewObj.deleteMaskButton], 'Enable', 'off')
            [x, y] = getpts(p.viewObj.mainAxes); %Simple but not interruptible. Can make WindowButtonDownFcn if want something interruptiblef.
            if ~isempty(x)
                ptsInView = d2utils.getPtsInsideView([x, y], p.viewRect);
                p.nucleiObj.removeCell(ptsInView(:,2), ptsInView(:,1));
                p.spotTable.assignSpotsInRect(p.viewRect);
                p.spotTable.updateAllSpotStatus();
                p.spotTable.makeCentroidList();
                p.updateCentroidListView();
                p.updateMainAxes();
            end
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.maskCellButton, p.viewObj.addCellButton,...
                p.viewObj.maskSpotsButton, p.viewObj.deleteMaskButton], 'Enable', 'on')
        end
        
        function p = changeMaxDistance(p, distance)
            mustBePositive(distance)
            p.spotTable.updateMaxDistance(distance);
            p.spotTable.makeCentroidList();
            p.updateCentroidListView();
            p.updateMainAxes();
            p.threshCntrlr.startup();
            disp('done')
        end
        
        function p = adjustNucleiParameters(p, minSize, varargin)
            n = inputParser;
            n.addRequired('minSize', @isnumeric)
            n.addParameter('sensitivity', 0.1, @(x)validateattributes(x,{'numeric'}, {'scalar', '>=',0,'<=',1.0}));
            n.parse(minSize, varargin{:});
            
            p.nucleiObj.minNucleusSize = n.Results.minSize;
            
            disp('Finding nuclei')
            if ~(n.Results.sensitivity == p.nucleiObj.nucMaskSens) %Re-mask if sensitivity value different than current value.  
                if all(p.scanObj.scanDim == 1)
                    p.nucleiObj.stitchDAPImask2('sensitivity', n.Results.sensitivity);
                else
                    p.nucleiObj.stitchDAPImask(n.Results.sensitivity);
                end
            end
            disp('Finding nuclei boundaries.')
            p.nucleiObj.findNuclei();
            
            p.nucleiObj.addColors();
            p.nucleiObj.updateAllMasks();
            
            p.spotTable.assignSpotsToNuclei();
            p.spotTable.updateAllSpotStatus();
            
            p.spotTable.makeCentroidList();
            p.updateCentroidListView();
            p.updateMainAxes();
            p.threshCntrlr.startup();
            disp('done')
        end
        
        function p = changeSpotFilterThreshold(p, channel, threshFactor)
            %Check that channel and filter are valid
            mustBeMember(channel, p.spotTable.spotChannels)
            mustBePositive(threshFactor)
            %reload non-contrast scan
            fprintf('Loading stitched scan.\nThis may take several minutes.\n')
            p.scanObj.reloadChannelStitch(channel);
            %refind spots
            fprintf('Refinding spots.\nThis may take several minutes.\n')
            p.spotTable.findSpotsChannel(channel, threshFactor);
            %assign new spots to cells
            p.spotTable.assignSpotsToNucleiChannel(channel);
            %mask new spots. This will also update status. 
            p.spotTable.updateChannelMasks(channel);
            %update centroid list
            p.spotTable.updateCentroidList(channel);
            %re-contrast scans
            disp('Auto-contrasting stitched scan. This may take several minutes.')
            p.scanObj.contrastStitchedScanChannel(channel, [1 99], [0.9 3]);
            disp('Resizing stitched scan')
            p.scanObj.resizeStitchedScanChannel(channel);
            p.updateCentroidListView();
            p.updateMainAxes();
            p.threshCntrlr.startup();
            disp('done')
        end
        
        function p = subtractBackground(p, channel, value)
            %Check that channel and value
            mustBeMember(channel, p.spotTable.spotChannels)
            mustBeA(value, 'logical')
            if value
                fprintf('Measuring %s background.\nThis may take several minutes.\n', channel)
                p.scanObj.measureChannelBackground(channel);
            end
            fprintf('Stitching channel %s.\nThis may take several minutes.\n', channel)
            %Restitch channel
            p.scanObj.restitchChannel(channel, value);
            %Save stitch
            fprintf('Saving stitched scans.\nThis may take several minutes.\n')
            fishScans = p.scanObj.stitchedScans.stitches;
            save('stitchedScans.mat', 'fishScans', '-append')
            %Find spots
            fprintf('Refinding spots.\nThis may take several minutes.\n')
            p.spotTable.findSpotsChannel(channel);  
            %assign new spots to cells
            p.spotTable.assignSpotsToNucleiChannel(channel);
            %mask new spots. This will also update status. 
            p.spotTable.updateChannelMasks(channel);
            %update centroid list
            p.spotTable.updateCentroidList(channel);
            %re-contrast scans
            disp('Auto-contrasting stitched scan.')
            p.scanObj.contrastStitchedScanChannel(channel, [1 99], [0.9 3]);
            disp('Resizing stitched scan')
            p.scanObj.resizeStitchedScanChannel(channel);
            p.updateCentroidListView();
            p.updateMainAxes();
            p.threshCntrlr.startup();
            disp('done')
        end
        
        function type = getSelectionType(p)
            type = get(p.viewObj.figHandle, 'SelectionType');
        end
        
        function keyPressFunctions(p, ~, evt)
            keyPressed = evt.Key;
            modifierPressed = evt.Modifier;
            spotChannelIdx=find(ismember(p.spotTable.spotChannels,p.nondapiChannelsInView{p.channelIdx}));
            if isempty(modifierPressed)
                switch(keyPressed)
                    case 'z'
                        p.zoomInPressed();
                    case 'p'
                        p.panViewPressed();
                    case 's'
                        p.overlaySpots();
                    case 'n'
                        p.overlayNuclei();
                    case 'c'
                        p.scatterCallback();
                    case 'm'
                        p.addSpotMask();
                    case 'd'
                        p.deleteMask();
                    case 'f' 
                        p.fixedZoom = true;
                    case 'x'
                        set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
                    case 'uparrow'
                        cellIdx = max(1, get(p.viewObj.centroidList, 'Value')-1);
                        cellPos = p.spotTable.centroidLists{spotChannelIdx}{cellIdx, {'x', 'y'}};
                        p.viewRect = d2utils.getRectAroundPoint(cellPos, 2 * p.cellViewRadius, 2 * p.cellViewRadius, p.scanObj.stitchDim);
                        set(p.viewObj.scatterCheckBox, 'Value', 0)
                        p.updateImageInView();
                        p.updateMainAxes();
                        p.thumbCntrlr.overlayThumbnailRect();
                    case 'downarrow'
                        cellIdx = min(get(p.viewObj.centroidList, 'Value')+1, p.numCells);
                        cellPos = p.spotTable.centroidLists{spotChannelIdx}{cellIdx, {'x', 'y'}};
                        p.viewRect = d2utils.getRectAroundPoint(cellPos, 2 * p.cellViewRadius, 2 * p.cellViewRadius, p.scanObj.stitchDim);
                        set(p.viewObj.scatterCheckBox, 'Value', 0)
                        p.updateImageInView();
                        p.updateMainAxes();
                        p.thumbCntrlr.overlayThumbnailRect();
                end
            elseif strcmp(modifierPressed{1}, 'shift')
                switch(keyPressed)
                    case 'm'
                        p.addCellMask();
                    case 's'
                        p.viewObj.saveButtonPressed;
                    case 'e'
                        p.viewObj.exportButtonPressed;
                end
            end
        end
    end
end