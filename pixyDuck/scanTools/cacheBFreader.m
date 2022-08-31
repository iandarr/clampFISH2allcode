function cacheBFreader(Nd2_file)
% Creates a hidden file in the same directory as Nd2_file which stores a
% bioformats-matlab reader for faster access to the Nd2 file
% 
% Example:
%   cacheBFreader('/Volumes/IAND_05/data/20210719_183029_315__WellA1_ChannelDAPI.nd2')
%   will create a hidden file in the directory 'Volumes/IAND_05/data/' that
%   is named: 
%       .20210719_183029_315__WellA1_ChannelDAPI.nd2.bfmemo
%
%  On a mac, you'll need to go to the folder and hit Command+Shift+. to
%  toggle viewing hidden files
%
% 
%
% Taken from the bioformats-matlab documentation https://docs.openmicroscopy.org/bio-formats/5.7.3/developers/matlab-dev.html

% Construct an empty Bio-Formats reader
reader= bfGetReader();
% Decorate the reader with the Memoizer wrapper
%reader = loci.formats.Memoizer(reader);
reader = loci.formats.Memoizer(reader,0); % not sure what the zero does. Is this the DEFAULT_MINIMUM_ELAPSED variable? unclear to me


% Initialize the reader with an input file
% If the call is longer than a minimal time, the initialized reader will
% be cached in a file under the same directory as the initial file
% name .large_file.bfmemo
reader.setId(Nd2_file);

reader.setMetadataStore()

% Close the reader
reader.close() % close() doesn't seem to be necessary and adds time. I suppose it may prevent overuse of memory though, I'm not sure

end