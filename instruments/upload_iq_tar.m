function upload_iq_tar(Instr,file_source,file_destintion)
%function to upload an iq.tar file stored on the remote PC to an FSW or FPS
%
%Input: Instr           - Instrument handle (Instrument Control Toolbox)
%       file_source     - iq.tar file on remote PC
%       file_destintion - iq.tar file on instrument

filedata = dir(file_source);
filesize = filedata.bytes;

BlockSize = Instr.OutputBufferSize - 1024;

Instr.EOIMode = 'off';

header = ['#' num2str(length(num2str(filesize))) num2str(filesize)];
fwrite (Instr, ['MMEM:DATA ''' file_destintion ''',' header], 'uchar');

iFile = fopen(file_source,'rb');

% display a progress bar
bar = waitbar( 0, 'Sending IQ data...' );
    
%the number of bytes remaining is initially set to the total number of bytes to send
Blocks = 0;
remaining = filesize;

while( remaining > 0 )
    nextblock = remaining;
    if(nextblock > BlockSize)
        nextblock = BlockSize;
    end
    
    % set the new number of bytes remaining
    remaining = remaining - nextblock;
    Blocks = Blocks + 1;
    
    rawdata = fread (iFile, nextblock, 'uchar');
    
    % last block to send
    if (remaining==0)
        Instr.EOIMode = 'on';
    end
    
    % calculate state
    Progr = (1.0-remaining/filesize);
    waitbar( Progr, bar );
    
    %send binary data
    fwrite (Instr, rawdata, 'uchar');
end
close( bar );
fclose(iFile);