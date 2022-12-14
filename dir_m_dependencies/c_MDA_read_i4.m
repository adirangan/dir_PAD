function output = c_MDA_read_i4(fname);
if (~exist(fname,'file')); disp(sprintf(' %% Warning! could not find %s',fname)); end;
fp = fopen(fname);
n_d = fread(fp,1,'int32');
d_ = fread(fp,n_d,'int32');
n_A = 1;
for na=1:n_d; n_A = n_A*d_(na); end;
A = fread(fp,n_A,'int');
if (n_d>1); A = reshape(A,transpose(d_)); end;
output = A;
fclose(fp);
