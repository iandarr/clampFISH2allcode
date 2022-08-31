
%% universal settings
iSeries=3;
iZ=1;
iC=1;
iT=1;
Nd2_filepath1='/Volumes/IAND_05/E159part2_scan/rawdata/s2_6well_cycle1/20210602_213621_459__WellA1_ChannelDAPI,YFP,CY3,A594,CY5_Seq0001.nd2';
%% Reader 1, With memoization
display('reader1 with memoization')
tic

% Construct an empty Bio-Formats reader
reader1 = bfGetReader();
% Decorate the reader with the Memoizer wrapper
reader1 = loci.formats.Memoizer(reader1); % if this line is commented out, you loose the speedup. Also, if you delete the hidden .filename.nd2.bfmemo file, then it'll again take long to initialize on the first attempt
% Initialize the reader with an input file
% If the call is longer than a minimal time, the initialized reader will
% be cached in a file under the same directory as the initial file
% name .large_file.bfmemo
reader1.setId(Nd2_filepath1);

% Perform work using the reader
% Read plane from series iSeries at Z, C, T coordinates (iZ, iC, iT)
% All indices are expected to be 1-based
reader1.setSeries(iSeries - 1);
iPlane = reader1.getIndex(iZ - 1, iC -1, iT - 1) + 1;
I = bfGetPlane(reader1, iPlane);

% Close the reader
reader1.close()

toc

%% Reader 2, Without memoization
display('reader2 without memoization')
tic
reader2=bfGetReader();
reader2.setId(Nd2_filepath1);

reader2.setSeries(iSeries - 1);
iPlane = reader2.getIndex(iZ - 1, iC -1, iT - 1) + 1;
I = bfGetPlane(reader2, iPlane);

% Close the reader
reader2.close()

toc

%%
%% test getting XY positions from cached vs. uncached
condIDtoCacheReaderMemoFor=1:7;

for iCond=1:length(condIDtoCacheReaderMemoFor)
    condID=condIDtoCacheReaderMemoFor(iCond);
    rowInTcond=find(Tcond.condID==condID); assert(length(rowInTcond)==1)
    Nd2_filepath=dirAndFilenametext2filepath(Tcond.Nd2_dir{rowInTcond},Tcond.Nd2_filenametext{rowInTcond});
    
    condID
    tic
    [Xpos,Ypos]=nd2positions(Nd2_filepath);
    length(Xpos)
    toc
end
