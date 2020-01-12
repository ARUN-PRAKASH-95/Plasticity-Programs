function [ ] = printdata(x,y,fname)
% printdata.m prints data into text file



%x = 0:.1:1;
data = [x; y];

% open the file with write permission
fid = fopen(fname, 'w');
%fprintf(fid, '%6.2f %12.8f\n', data);
fprintf(fid, '%12.8f %12.8f\n', data);
fclose(fid);

end

