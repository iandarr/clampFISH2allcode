
Tchan2=Tchan;
Tchan2.myNum=[1:height(Tchan2)]'
Tchan2.round(20)=7
Tchan2.channel(20)={'YFP'}

%%

T=formatChannelAssociatedInputs(Scan,Tchan2)
%T.channelLabel{1}
T.channelLabel{4}