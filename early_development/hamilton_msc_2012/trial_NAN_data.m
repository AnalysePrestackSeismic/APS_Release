clear 

total_per_slice = 300*650;
skip = 0
for i = 1:466
    
    data(1+skip:skip+total_per_slice) = 1e30;
    skip = skip + total_per_slice;
    
end