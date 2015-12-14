function status = mcmc_deconvolution(traces_in,output_location,proc_pos,warm_in,cycle_in)

% number of blocks to process
n_blocks = length(proc_pos)
dt = traces_in.s_rate/1e6;

% number of position to invert at once
n_pos = (2*traces_in.aperture+1)^2;

    for i_block = 1:1:n_blocks
        % Read all traces into memory
        traces = segy_read_traces(traces_in,proc_pos(i_block,4),proc_pos(i_block,5));     
        ntraces = size(traces.data,2);

        start_index = 1;
        n_groups = ntraces/n_pos;
        for ij=1:1:n_groups
            t = traces.data(:,start_index:n_pos*ij);
            % fprintf('Group %d\n', ij);
            % Estimate wavelet and noise level
            % parameters needed by bliss_kpe in order to find wavelet from data
            tw = 0.1; %input('Wavelet time length (s) = ');
            nw_h_in = round(tw/dt/2);
            get_phase_in = 1;
            damp_in = 3; % input('damping constant in Wiener filtering (>1) = ');
            [w0,vn0,~,~,~]=bliss_kpe(t,nw_h_in,get_phase_in,damp_in);
            
            % noise sd
            w0 = w0*sqrt(vn0);
            t = t(:);
                      
            [~,rm.data(:,ij)] = ...
                run_decon(t,w0,sqrt(vn0),warm_in,cycle_in);  
            
            rm.pos(1,ij) = traces.pos(1,(start_index+(ij*n_pos))/2); 
            rm.pos(2,ij) = traces.pos(2,(start_index+(ij*n_pos))/2);             
            start_index = start_index+n_pos;
        end
        
        worker = num2str(proc_pos(i_block,2));
        block = num2str(proc_pos(i_block,3));
        
        % Write data
        data_out = strcat(output_location,'/','mcmc_decon_worker_',worker,'_block_',block);
        file_out = fopen(data_out,'a');
        fwrite(file_out,rm.data,'float32');
        fclose(file_out);
        
        % Write position information
        meta_out = strcat(output_location,'/','mcmc_decon_pos_worker_',worker,'_block_',block);
        dlmwrite(meta_out,rm.pos','delimiter','\t');

    end
    
    status = 1;

end

function [rhat,rm] = run_decon(x,w,sn,n1,n2)


% Syntaxe [rhat, rm, rm0] = mcmcdecon(x,w,sn,n1,n2,verbose);
% Perform deconvolution by Markov Chain Monte-Carlo (MCMC) method under
% the model of independent sources (reflectivity) with Bernoulli-Gaussian
% distribution of unit variance, and Gaussian white noise.
% Input:
% - x: the data, assumed to obey the model x = conv(w,r) + noise
% - w: the wavelet (filter), assumed known and of short length
% - sn: the noise standard deviation, assumed known
% - n1: number of warm-up cycles in MCMC
% - n2: number of cycles in MCMC to compute the deconvolution result
% - verbose: set this variable is to true to get some printout
% Output:
% - rhat: Wiener deconvolution in the time domain (a column vector with
%	  length = length(x) + length(w)-1)
% - rm: MCMC deconvolution, same size as rhat
% * if n1, hence n2, are not given, only rhat is computed.
%
% Copyright Pham Dinh-Tuan Antoine, LJK, April 17 2007

%     if nargin < 6
%       verbose = false;
%     end

    m = numel(w);
    n = numel(x);
    vn = sn^2;

    np = n+m-1;						% length of rhat
    c = filter(w(m:-1:1),1,w(:)');			% ensure that c is a row vector
    C = spdiags(c(ones(np,1),m:-1:1),0:m-1,np,np);
    R = toeplitz([w(m); zeros(m-2,1)], w(m:-1:2));
    C(1:m-1,1:m-1) = R'*R;
    R = toeplitz(w(1:m-1), [w(1) zeros(1,m-2)]);
    C(n+1:np,n+1:np) = R'*R;
    R = chol(C + sparse(1:np,1:np,vn));		% R'*R = C + vn*speye(np)
    xp = filter(w(m:-1:1),1, [x(:); zeros(m-1,1)]);		% column vector
    rhat = R\(R'\xp);

    if nargin > 3
      phat = 3*(rhat'*rhat)^2/(np*sum(rhat.^4));	% must be < 1
      if phat >= 1
      error('Data is not super Gaussian enough for the Bernoulli-Gaussian model\n')
      end
      r = rhat;
      c1 = c(m-1:-1:1);
      c0 = c(m) + phat*vn;
      vr1 = vn*c0;
      sre = sn/sqrt(c0);
%       if verbose
%         fprintf('Warm up\n');
%       end
      for nrep = 1:n1
        for k = 1:m-1				% first m-1 indexes
          r1 = xp(k) - [C(k,1:k-1) C(k,k+1:k+m-1)]*[r(1:k-1); r(k+1:k+m-1)];
          ck = C(k,k) + phat*vn;
          lam = exp(-0.5*r1^2/(vn*ck))/sqrt(phat*vn/ck);
          if rand < phat/(phat + lam*(1-phat))
        r(k) = sqrt(vn/ck)*randn + r1/ck;
          else
        r(k) = 0;
          end
        end
        for k = m:n					% middle indexes
          r1 = xp(k) - c1*(r(k-1:-1:k+1-m) + r(k+1:k+m-1));
          lam = exp(-0.5*r1^2/vr1)/(sqrt(phat)*sre);
          if rand < phat/(phat + lam*(1-phat))
        r(k) = sre*randn + r1/c0;
          else
        r(k) = 0;
          end
        end
        for k = n+1:np				% last m-1 indexes
          r1 = xp(k) - [C(k-m+1:k-1,k)' C(k,k+1:np)]*[r(k-m+1:k-1); r(k+1:np)];
          ck = C(k,k) + phat*vn;
          lam = exp(-0.5*r1^2/(vn*ck))/sqrt(phat*vn/ck);
          if rand < phat/(phat + lam*(1-phat))
        r(k) = sqrt(vn/ck)*randn + r1/ck;
          else
        r(k) = 0;
          end
        end
        phat = sum(r~=0)/np;
%         if verbose
%           fprintf('cycle %3d percentage of nonzero %7.4f\r', nrep, phat)
%         end
        c0 = c(m) + phat*vn;
        vr1 = vn*c0;
        sre = sn/sqrt(c0);
      end

      rm = zeros(size(r));
      %rm0 = zeros(size(r));
%       if verbose
%         fprintf('\nCalculation\n');
%       end
      for nrep = 1:n2
        for k = 1:m-1				% first m-1 indexes
          r1 = xp(k) - [C(k,1:k-1) C(k,k+1:k+m-1)]*[r(1:k-1); r(k+1:k+m-1)];
          ck = C(k,k) + phat*vn;
          lam = exp(-0.5*r1^2/(vn*ck))/sqrt(phat*vn/ck);
          if rand < phat/(phat + lam*(1-phat))
        r(k) = sqrt(vn/ck)*randn + r1/ck;
          else
        r(k) = 0;
          end
        end
        for k = m:n					% middle indexes
          r1 = xp(k) - c1*(r(k+1:k+m-1) +r(k-1:-1:k+1-m));
          lam = exp(-0.5*r1^2/vr1)/(sqrt(phat)*sre);
          if rand < phat/(phat + lam*(1-phat))
        r(k) = sre*randn + r1/c0;
          else
        r(k) = 0;
          end
        end
        for k = n+1:np				% last m-1 indexes
          r1 = xp(k) - [C(k-m+1:k-1,k)' C(k,k+1:np)]*[r(k-m+1:k-1); r(k+1:np)];
          ck = C(k,k) + phat*vn;
          lam = exp(-0.5*r1^2/(vn*ck))/sqrt(phat*vn/ck);
          if rand < phat/(phat + lam*(1-phat))
        r(k) = sqrt(vn/ck)*randn + r1/ck;
          else
        r(k) = 0;
          end
        end
        idx = find(r~=0);
        nz = length(idx);
        phat = nz/np;
%         if verbose
%           fprintf('cycle %3d percentage of nonzero %7.4f\r', nrep, phat)
%         end
        rm = rm + r;
        %rm0 = rm0 + r;
        %R = chol(C(idx,idx) + sparse(1:nz,1:nz,vn*phat));
        %rm(idx) = rm(idx) + R\(R'\xp(idx));
        %This may inprove a little bit the estimator, but not worth the computation
        %%phat = 3*sum(r(idx).^2)/(np*sum(r(idx).^4));
        %%phat = phat + (nphat - phat)/(nrep+1);
        c0 = c(m) + phat*vn;
        vr1 = vn*c0;
        sre = sn/sqrt(c0);
      end
%       if verbose
%         fprintf('\n')
%       end
      %rm0 = rm0/n2;
      rm = rm/n2;
    end

end
