% test bfGetPlane speed with cropped images

clear

Nd2_filepath = '/Volumes/IAND_04/E157_validation/ReadoutStrip/20X/s6_BeforeStrip_5x5/c1 20210311_183328_436__Well2_ChannelDAPI,YFP,YFP_1,YFP_2,YFP_3,CY3,CY3_1,CY3_2,CY3_3,A594,A594_1,A594_2,A594_3,CY5,CY5_1,CY5_3,CY5_2,Brightfield_Seq0001.nd2';
reader = bfGetReader();

reader = loci.formats.Memoizer(reader,0);
reader.setId(Nd2_filepath);

%%
img=[];
for i=1:4
    i
    tic
    img=getPlaneFromNd2file(reader,1,'DAPI','closeReaderAfterGettingPlane',false);
    toc
end
reader.close()
%%
img=[];
for i=1:4
    i
    tic
    img=getPlaneFromNd2file(reader,1,'DAPI','closeReaderAfterGettingPlane',false,'bbox',[1 1 200 200]);
    toc
end
reader.close()