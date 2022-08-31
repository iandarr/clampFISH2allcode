dentist.tests.cleanupForTests;

figH = figure(1);

thresholdRadio = uicontrol('Parent',figH,...
    'Style','radiobutton',...
    'Units','normalized',...
    'Position',[0.009,0.179,0.346,0.857],...
    'HandleVisibility','off',...
    'String','Threshold');
refObjRadio = uicontrol('Parent',figH,...
    'Style','radiobutton',...
    'Units','normalized',...
    'Position',[0.338,0.179,0.376,0.857],...
    'HandleVisibility','off',...
    'String','Ref. Object');
defaultRadio = uicontrol('Parent',figH,...
    'Style','radiobutton',...
    'Units','normalized',...
    'Position',[0.709,0.179,0.286,0.857],...
    'HandleVisibility','off',...
    'String','Default');

buttonStruct = struct('threshold', thresholdRadio, ...
    'refObj', refObjRadio, 'default', defaultRadio);

funcsOnToggleIn = struct('threshold', @() fprintf('threshold\n'), ...
    'refObj', @() fprintf('refObj\n'), ...
    'default', @() fprintf('default\n'));

funcsOnToggleOut = struct('threshold', @() fprintf('out threshold\n'), ...
    'refObj', @() fprintf('out refObj\n'), ...
    'default', @() fprintf('out default\n'));

x = dentist.utils.FunctionExecutingToggleGroup(buttonStruct, ...
    funcsOnToggleIn, funcsOnToggleOut);

x.initialize('threshold')

fprintf(['click on the radio buttons. \n',...
    'They should mutually exclude, and print switch status to console']);
