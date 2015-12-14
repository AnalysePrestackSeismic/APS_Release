cjj = 1;
prevout = [];
totalloop = 400;
opts = spgSetParms('verbosity',0);
opts.quitPareto = 1;
iters_performed = 0;
tottime = 0;
iter = 300;
tol = 1e-3;
figure;
while cjj <  totalloop
    if cjj == 1
        tau = 0.1;
        prevout = [];
        %opts.decTol = 5e-1;
        opts.decTol = 5e-2;
    elseif cjj < 200
        tau = norm(prevout,1)/2;
        opts.decTol = 5e-3;
    else
        tau = norm(prevout,1)/2;
        opts.decTol = 5e-4;        
    end
    [ava_tmpp,r,g,info] = spgl1(IGiter, data,tau, 5.5e2, prevout,opts );
    prevout = ava_tmpp;
    cjj = cjj + 1;
    iters_performed = info.iter + iters_performed;
    tottime = info.timeTotal + tottime;
end
info
plot(ava_tmpp,'-g');
hold all;
disp(['(current cumulated total L1 gradient-step count: ' int2str(iters_performed) ')'])
disp(['(current cumulated time: ' num2str(tottime) ')'])

opts.quitPareto = 0;
% opts.decTol = 1e-4;
% [ava_tmpbra,r,g,info] = spgl1(IGiter, data, 0, 5.6e2, [],opts );
% plot(ava_tmpbra);
% info
opts.decTol = 1;
[ava_tmpbrb,r,g,info] = spgl1(IGiter, data, 0, 5.65e2, [],opts );
info
plot(ava_tmpbrb,'-b');


opts.decTol = 1;
[ava_tmpbrc,r,g,info] = spgl1(IGiter, data, [], [], [],opts );
info
plot(ava_tmpbrc,'-c');


[ava_tmp,lsqflag,~,fiternotmp] = lsqr(IGiter,data,tol,iter,[],[]);
plot(ava_tmp,'-r');
title('red = lsqr, blue = splgl1, green = spgl1 with restarts');