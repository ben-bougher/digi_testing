function wigglePlot(data, outfile)

fileID = fopen('tmp.bin','w');

fwrite(fileID, data, 'float32');
fclose(fileID);
nSamps = num2str(size(data,1));
['bash wigglePlot.sh tmp.bin ', nSamps,' ',outfile]

system(['bash wigglePlot.sh tmp.bin ', nSamps,' ',outfile]);

end