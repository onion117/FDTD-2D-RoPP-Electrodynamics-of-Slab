function [  ] = write_dat( fname,x,y )

delete(fname)
fileID = fopen(fname,'w+');
[nrows] = length(x);
 for row = 1:nrows
fprintf(fileID,' %d %d \n ' ,x(row),y(row));
 end
fclose(fileID)


end

