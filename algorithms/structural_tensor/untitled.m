        dy1 = -1/2*ones(nt*nxl*nil,1);
        dy2 = 1/2*ones(nt*nxl*nil,1);
        
        dy = spdiags([dy1,dy2],[0 2*nxl*nt],nt*nxl*nil,nt*nxl*nil);
        
        dy1 = -1*ones(nt*nxl,1);
        dy2 = ones(nt*nxl,1);
        
        dy_top = spdiags([dy1,dy2],[0 nt*nxl],nt*nxl,nt*nxl*nil);
        dy_bot = spdiags([dy1,dy2],[nsamp-(2*nt*nxl) nsamp-(nt*nxl)],nt*nxl,nt*nxl*nil);
        
        dy2 = [dy_top; dy(1:nsamp-(2*nt*nxl),:);dy_bot];
        dy = dy*data;
        
        
                nsamp = nt*nxl*nil; 
        dy1 = -1/2*ones(nt*nxl*nil,1);
        dy2 = 1/2*ones(nt*nxl*nil,1);
        
        dy = spdiags([dy1,dy2],[0 2*nxl*nt],nt*nxl*nil,nt*nxl*nil);
        
        dy1 = -1*ones(nt*nxl,1);
        dy2 = ones(nt*nxl,1);
        
        dy_top = spdiags([dy1,dy2],[0 nt*nxl],nt*nxl,nt*nxl*nil);
        dy_bot = spdiags([dy1,dy2],[nsamp-(2*nt*nxl) nsamp-(nt*nxl)],nt*nxl,nt*nxl*nil);
        
        dy = [dy_top; dy(1:nsamp-(2*nt*nxl),:);dy_bot];
        
        
                nsamp = nt*nxl*nil; 
        dy1 = [-1*ones(nt*nxl,1); zeros(nsamp-2*nt*nxl,1); ones(nt*nxl,1)];
        dy2 = [ones(2*nt*nxl,1); 1/2*ones(nt*nxl*nil-2*nt*nxl,1);];
        dy3 = [-1/2*ones(nt*nxl*nil-2*nt*nxl,1); -1*ones(2*nt*nxl,1)];       
      
        dy = spdiags([dy3,dy1,dy2],[-nxl*nt 0 nxl*nt],nt*nxl*nil,nt*nxl*nil);             
        dy = dy*data;
        
        
        
        
        
        
                nsamp = nt*nxl*nil; 
        dy1 = [-1*ones(nt*nxl,1); zeros(nsamp-(2*nt*nxl),1); ones(nt*nxl,1)];
        dy2 = [ones(nt*nxl,1); 1/2*ones(nsamp-nt*nxl,1)];
        dy2 = [ones(nt*nxl,1); dy2(1:end-nt*nxl,1)];
        
        dy3 = [-1/2*ones(nsamp-nt*nxl,1); -1*ones(nt*nxl,1)];
        dy3 = [dy3(1:end-2*nt*nxl+1,1); -1*ones(2*nt*nxl-1,1)];
        
        dy = spdiags([dy3,dy1,dy2],[-(nxl*nt)+1 0 (nxl*nt)+1],nt*nxl*nil,nt*nxl*nil);
        
        %dy1 = -1*ones(nt*nxl,1);
        %dy2 = ones(nt*nxl,1);
        
        %dy_top = spdiags([dy1,dy2],[0 nt*nxl],nt*nxl,nt*nxl*nil);
        %dy_bot = spdiags([dy1,dy2],[nsamp-(2*nt*nxl) nsamp-(nt*nxl)],nt*nxl,nt*nxl*nil);
        
        %dy = [dy_top; dy(1:nsamp-(2*nt*nxl),:);dy_bot];
        dy = dy*data;