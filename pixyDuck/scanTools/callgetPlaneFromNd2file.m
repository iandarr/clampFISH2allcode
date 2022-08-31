%%
display('without closing reader') %slow at 1s
for tileID=1
    tic
    outImg = getPlaneFromNd2file(Nd2_filepath_ref, tileID, 'DAPI','iZ',1,'closeReaderAfterGettingPlane',false);
    toc
    imagesc(outImg)
end

display('with closing reader') %faster, 0.25s
for tileID=1
    tic
    outImg = getPlaneFromNd2file(Nd2_filepath_ref, tileID, 'DAPI','iZ',1,'closeReaderAfterGettingPlane',true);
    toc
    imagesc(outImg)
end

%% passing reader directly
display('giving reader directly, without closing reader') % this is the fastest (0.15 seconds when I don't close reader)
reader = bfGetReader();
    reader = loci.formats.Memoizer(reader,0);
    reader.setId(Nd2_filepath_ref);
tic
    outImg=getPlaneFromNd2file(reader,tileID, 'DAPI','iZ',1,'closeReaderAfterGettingPlane',false);
toc