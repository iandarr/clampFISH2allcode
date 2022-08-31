function Hs = layOutThresholdGUICore()
    backgColor = [0.247 0.247 0.247];
    
    Hs.figH = figure('Position',[650 100 255 500],...
        'NumberTitle','off',...
        'Name','improc2.ThresholdInspectorGUI',...
        'Resize','on',...
        'Toolbar','none',...
        'MenuBar','none',...
        'Color',backgColor,...
        'Visible','off');
    
    set(Hs.figH, 'Visible', 'on');
    
    Hs.threshPanel  = uipanel('Parent',Hs.figH,...
        'Units','normalized',...
        'Position',[0 0.3889 1 0.3889],...
        'BorderType', 'none', ...
        'BackgroundColor',backgColor);
    
    
    Hs.thresholdAx = axes('Parent',Hs.threshPanel,...
        'Units','normalized',...
        'Position',[0.0282 0.0829 0.9436 0.9029],...
        'YTick',[],...
        'XColor','w');
    
    Hs.channelMenu = uicontrol('Parent', Hs.threshPanel,...
        'Style', 'popupmenu', 'String', {'empty'}, ...
        'Value', 1, ...
        'Units', 'normalized',...
        'FontSize', 16, ...
        'Position', [0.6 0.85 0.39 0.15]);
    
    Hs.upperFreeArea = uipanel('Parent',Hs.figH,...
        'Units','normalized',...
        'BorderType', 'none', ...
        'Position',[0.7000 0.7889 0.2900 0.2111],...
        'BackgroundColor',backgColor);
    
    Hs.plotLogYCheck = uicontrol('Parent',Hs.figH,...
        'Style','checkbox',...
        'Units','normalized',...
        'BackgroundColor',backgColor,...
        'ForegroundColor',[1 1 1], ...
        'Position',[0.0282 0.8889 0.6301 0.0444],...
        'String','log Num Spots','FontSize',12);
    Hs.autoXAxisCheck = uicontrol('Parent',Hs.figH,...
        'Style','checkbox',...
        'Units','normalized',...
        'BackgroundColor',backgColor,...
        'ForegroundColor',[1 1 1], ...
        'Position',[0.0282 0.8389 0.6301 0.0444],...
        'String','auto X axis',...
        'ToolTip', ['otherwise set max to image saturation', ...
        'value in ''Fixed'' contrast mode'], ...
        'FontSize',12);
    Hs.hasClearThresholdPopup = uicontrol('Parent',Hs.figH,...
        'Style','popupmenu',...
        'Units','normalized',...
        'Position',[0.0282 0.7889 0.301 0.0444],...
        'String','-','FontSize',12);
    Hs.hasClearThresholdLabel = uicontrol('Parent',Hs.figH,...
        'Style','text',...
        'Units','normalized',...
        'Position',[0.32 0.7889 0.39 0.0444],...
        'BackgroundColor',backgColor,...
        'ForegroundColor',[1 1 1], ...
        'String','thresh is clear','FontSize',12);
    
    
    Hs.lowerFreeArea = uipanel('Parent',Hs.figH,...
        'Units','normalized',...
        'BorderType', 'none', ...
        'Position',[0.0345 0 0.9404 0.12],...
        'BackgroundColor',backgColor);
    
    Hs.naviPanel  = uipanel('Parent',Hs.figH,...
        'Units','normalized',...
        'Position',[0.0345 0.1333 0.9404 0.2249],...
        'BackgroundColor',[0.929 0.929 0.929]);
    Hs.goodCheck = uicontrol('Parent',Hs.naviPanel,...
        'Style','checkbox',...
        'Units','normalized',...
        'Position',[0.05 0.7663 0.4030 0.2262],...
        'String','Good Object',...
        'Value',0);
    Hs.prevObj = uicontrol('Parent',Hs.naviPanel,...
        'Style','pushbutton',...
        'Units','normalized',...
        'Position',[0.5250 0.5275 0.2160 0.400],...
        'String','Prev',...
        'FontSize',12,'FontWeight','bold');
    Hs.nextObj = uicontrol('Parent',Hs.naviPanel,...
        'Style','pushbutton',...
        'Units','normalized',...
        'Position',[0.7580 0.5275 0.2160 0.400],...
        'String','Next',...
        'FontSize',12,'FontWeight','bold');
    
    Hs.countsLabel = uicontrol('Parent',Hs.naviPanel,'Style','text',...
        'String','#Spots:',...
        'Units','normalized',...
        'Position',[0.0500 0.5375 0.2000 0.1962],...
        'FontSize',12,'FontWeight','bold');
    Hs.countsDisplay = uicontrol('Parent',Hs.naviPanel,'Style','text',...
        'String','0',...
        'Units','normalized',...
        'Position',[0.2500 0.5375 0.2500 0.1962],...
        'FontSize',12,'FontWeight','bold');
    
    Hs.fileLabel = uicontrol('Parent', Hs.naviPanel, ...
        'Style','text',...
        'String',sprintf('File:\n/ '),...
        'Units', 'normalized', ...
        'Position', [0.0250 0.1250 0.2160 0.3125], ...
        'FontSize', 12, 'FontWeight', 'bold');
    Hs.goToArray = uicontrol('Parent',Hs.naviPanel,...
        'Style','edit',...
        'Units','normalized',...
        'Position',[0.2750 0.1250 0.2160 0.3125],...
        'FontSize',12);
    Hs.objectLabel = uicontrol('Parent', Hs.naviPanel, ...
        'Style','text',...
        'String',sprintf('Obj:\n/ '),...
        'Units', 'normalized', ...
        'Position', [0.5250 0.1250 0.2160 0.3125], ...
        'FontSize', 12, 'FontWeight', 'bold');
    Hs.goToObj = uicontrol('Parent',Hs.naviPanel,...
        'Style','edit',...
        'Units','normalized',...
        'Position',[0.7750 0.1250 0.2160 0.3125],...
        'FontSize',12);
    
end

