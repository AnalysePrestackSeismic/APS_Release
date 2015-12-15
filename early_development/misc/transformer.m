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
    fid = fopen(time_in);
    time = fread(fid,[nt,ntrace],'float32');
    fclose(fid);

    time = time';

    time_flat_in = input('Enter flattened time volume filename: ');
    nt_flat = input('Enter number of time samples: ');

    fid = fopen(time_flat_in);
    time_flat = fread(fid,[nt_flat,ntrace],'float32');
    fclose(fid);

    time_flat = time_flat';
    index_mat = zeros(nt,ntrace);
    
    attribute_flat_in = input('Enter flattened attribute volume filename: ');

    ns_flat = nil*nxl*nt_flat;
    fid = fopen(attribute_flat_in);
    attribute_flat = fread(fid,ns_flat,'float32');
    fclose(fid);
else
    index_mat_in = input('Enter time lookup table filename: ');
    nil = input('Enter number of inlines: ');
    nxl = input('Enter number of xlines: ');
    nt = input('Enter number of time samples: ');

    ns = nil*nxl*nt;
    ntrace = ns/nt;
    fid = fopen(index_mat_in);
    index_mat = fread(fid,ns,'float32');
    fclose(fid);
    
    attribute_flat_in = input('Enter flattened attribute volume filename: ');
    nt_flat = input('Enter number of time samples: ');
    
    ns_flat = nil*nxl*nt_flat;
    fid = fopen(attribute_flat_in);
    attribute_flat = fread(fid,ns_flat,'float32');
    fclose(fid);
end

attribute_unflat = zeros(ns,1);

if compute_lookup == 1
    fprintf('\nCalculating time lookup table...\n');
    for k=1:ntrace        
        time_flat_rep = repmat(time_flat(k,:)',1,nt);
        time_rep = repmat(time(k,:),nt_flat,1);
        difference = abs(time_flat_rep-time_rep);
        [B,index] = min(difference);
        index_mat(:,k) = index' + (k-1)*nt_flat;
        if (ntrace/50)*floor(k/(ntrace/50)) == k;
            fprintf('%d / %d = %.0f%%\n',k,ntrace,(k*100)/(ntrace));
        end
    end
    index_mat = reshape(index_mat,[],1);
end

fprintf('\nUnflattening volume ''%s''...\n',attribute_flat_in);

for l=1:ns
    attribute_unflat(l,1) = attribute_flat(index_mat(l,1),1);
    if (ns/50)*floor(l/(ns/50)) == l;
        fprintf('%d / %d = %.0f%%\n',l,ns,(l*100)/(ns));
    end
end

if compute_lookup == 1
    fprintf('\nSaving time lookup table to ''time_lookup_table_%s_to_%s''...\n',time_flat_in,time_in);
    fid = fopen(sprintf('time_lookup_table_%s_to_%s',time_flat_in,time_in),'w');
    fwrite(fid,index_mat,'float32');
    fclose(fid);
end

fprintf('\nSaving unflattened volume to ''unflat_%s''...\n',attribute_flat_in);

fid = fopen(sprintf('unflat_%s_2',attribute_flat_in),'w');
fwrite(fid,attribute_unflat,'float32');
fclose(fid);

fprintf('\nComplete\n');
