rseed=cast(3,'uint32');
rseed=RSEED_adv8(rseed); disp(sprintf(' %% rseed %d',rseed));
rseed=RSEED_adv8(rseed); disp(sprintf(' %% rseed %d',rseed));
n_iteration = 3;
v_ = zeros(n_iteration,1);
for niteration=0:n_iteration-1;
[v_(1+niteration),rseed] = R01GET(rseed);
end;%for niteration=0:n_iteration-1;
disp(sprintf(' %% v_: %f %f %f ... %f %f %f',v_([1,2,3,end-2,end-1,end-0])));
disp(sprintf(' %% rseed %d',rseed));