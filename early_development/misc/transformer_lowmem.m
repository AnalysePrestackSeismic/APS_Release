clear all
clc

compute_lookup = input('Do you want to compute a new time lookup table? (1 = yes, 0 = no): ');

if compute_lookup == 1
    time_in = input('Enter time volume filename: ');
    nil = input('Enter number of inlines: ');
    nxl = input('Enter number of xlines: ');
    nt = input('Enter number of time samples: ');

    ns = nil*nxl*nt;
    ntrace = ns/nt;
    fid1 = fopen(time_in);

    time_flat_in = input('Enter flattened time volume filename: ');
    nt_flat = input('Enter number of time samples: ');

    fid2 = fopen(time_flat_in);
    
    attribute_flat_in = input('Enter flattened attribute volume filename: ');

    ns_flat = nil*nxl*nt_flat;
    fid3 = fopen(attribute_flat_in);
else
    index_mat_in = input('Enter time lookup table filename: ');
    nil = input('Enter number of inlines: ');
    nxl = input('Enter number of xlines: ');
    nt = input('Enter number of time samples: ');

    ns = nil*nxl*nt;
    ntrace = ns/nt;
    fid = fopen(index_mat_in);
    index_mat = fread(fid,ns,'uint');
    fclose(fid);
    
    attribute_flat_in = input('Enter flattened attribute volume filename: ');
    nt_flat = input('Enter number of time samples: ');
    
    ns_flat = nil*nxl*nt_flat;
    fid3 = fopen(attribute_flat_in);
end

if compute_lookup == 1
    index_mat = zeros(nt,ntrace);
    totalcols = 0;
    nblock = 1;
    nilblock = floor(nil/nblock);
    leftovers = nil - (nblock*nilblock);
    fprintf('\nCalculating time lookup table...\n');
    for k=1:nblock
        if k == nblock
            time = fread(fid1,[nt,(nilblock+leftovers)*nxl],'float32');
            time_flat = fread(fid2,[nt_flat,(nilblock+leftovers)*nxl],'float32');
        else
            time = fread(fid1,[nt,nilblock*nxl],'float32');
            time_flat = fread(fid2,[nt_flat,nilblock*nxl],'float32');
        end
        [rows,cols] = size(time);
        totalcols = totalcols + cols;
        index_mat_traces = zeros(nt,cols);
        for m=1:cols
            time_flat_rep = repmat(time_flat(:,m),1,nt);
            time_rep = repmat(time(:,m)',nt_flat,1);
            difference = abs(time_flat_rep-time_rep);
            [B,index] = min(difference);
            index_mat_traces(:,m) = index';
            if (ntrace/50)*floor((m+totalcols-cols)/(ntrace/50)) == m+totalcols-cols;
                fprintf('%d / %d = %.0f%%\n',m+totalcols-cols,ntrace,((m+totalcols-cols)*100)/(ntrace));
            end
        end
        index_mat(:,1+(totalcols-cols):totalcols) = index_mat_traces;
    end
    for n=1:ntrace
        index_mat(:,n) = index_mat(:,n) + (n-1)*nt_flat;
    end
    index_mat = reshape(index_mat,[],1);
    fclose(fid1);
    fclose(fid2);
    clearvars -except index_mat attribute_flat_in ns_flat fid3 ns compute_lookup time_flat_in time_in
end

fprintf('\nUnflattening volume ''%s''...\n',attribute_flat_in);

attribute_flat = fread(fid3,ns_flat,'float32');
fclose(fid3);
attribute_unflat = zeros(ns,1);
for l=1:ns
    attribute_unflat(l,1) = attribute_flat(index_mat(l,1),1);
    if (ns/50)*floor(l/(ns/50)) == l;
        fprintf('%d / %d = %.0f%%\n',l,ns,(l*100)/(ns));
    end
end

if compute_lookup == 1
    fprintf('\nSaving time lookup table to ''time_lookup_table_%s_to_%s''...\n',time_flat_in,time_in);
    fid4 = fopen(sprintf('time_lookup_table_%s_to_%s',time_flat_in,time_in),'w');
    fwrite(fid4,index_mat,'uint');
    fclose(fid4);
end

fprintf('\nSaving unflattened volume to ''unflat_%s''...\n',attribute_flat_in);

fid5 = fopen(sprintf('unflat_%s',attribute_flat_in),'w');
fwrite(fid5,attribute_unflat,'float32');
fclose(fid5);

fprintf('\nComplete\n');
