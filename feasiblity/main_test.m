LASTN=maxNumCompThreads(1);
array_index=getenv('PBS_ARRAY_INDEX');
i_test=str2num(array_index);
i1 = str2num(array_index);
rng(i1); 

load Feasibility_Jun1.mat

Feasibility_RS(1,i_test) = i_test
Feasibility_MULP(2,i_test) = i_test^2
save('Feasibility_Jun1.mat','Feasibility_RS','Feasibility_MULP')
