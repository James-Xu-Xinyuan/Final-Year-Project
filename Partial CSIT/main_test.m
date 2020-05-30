cvx_setup;
LASTN=maxNumCompThreads(1);
array_index=getenv('PBS_ARRAY_INDEX');
i_theta=str2num(array_index);
i1 = str2num(array_index);
rng(i1); 

for i  = 1:1000
    H = randn(10,10);
    H_inv = inv(H);
end

test = 123456;
save('outputtest','test','H')
