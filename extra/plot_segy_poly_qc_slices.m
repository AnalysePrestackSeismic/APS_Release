function plot_segy_poly_qc_slices( filepathin )

%PLOT_SEGY_POLY_QC_SLICES 

xv = [646474.45
683005.12
688841.39
689597.94
690030.25
690030.25
700621.99
700730.07
707647.12
710673.33
713807.62
713699.54
713699.54
719427.72
722994.32
725480.14
725912.45
726236.69
726236.69
725928.61
710139.01
700937.13
700949.41
667750.10
670213.87
674587.06
674463.87
674094.31
672953.78
672521.47
671548.76
670792.21
670359.89
670143.73
670035.65
669279.10
669279.10
669279.10
667982.15
667333.68
664218.30
664326.52
662471.04
662470.13
659119.68
657390.42
657386.52
655682.87
655661.16
648095.63
645501.74
644204.79
643340.16
642151.29
641286.66
640854.35
640422.03
639773.56
639125.08
638584.69
637936.22
634910.01
632424.19
630694.93
629722.22
628425.27
626263.69
623994.04
622264.77
620319.35
619238.56
617941.62
616644.67
615780.04
614158.85
612861.91
612429.59
611564.96
611024.57
610159.93
609608.09
609628.59
606765.42
605728.70
603134.80
601729.78
597514.70
595353.12
592860.00
592543.07
592867.31
592867.31
593083.46
593191.54
604854.85
641049.14
646474.45];

yv = [8949376.57
8949268.49
8949268.49
8947539.23
8947539.23
8946566.51
8921059.89
8924626.49
8927760.78
8929057.73
8930246.60
8931111.23
8931759.70
8934353.60
8935758.62
8937055.57
8937163.65
8936839.41
8910792.39
8894027.26
8894006.35
8894027.26
8889451.46
8889389.86
8880828.26
8863273.89
8862411.57
8862349.98
8859238.75
8857941.80
8856969.09
8856320.62
8855888.30
8856752.93
8857185.25
8856752.93
8855347.91
8854915.59
8854267.12
8853078.25
8852648.54
8875780.42
8875788.02
8875882.90
8875882.90
8875882.90
8875808.86
8875815.84
8876315.22
8876315.22
8876207.14
8876315.22
8876639.46
8876639.46
8877612.17
8878044.48
8878368.72
8879233.35
8880746.45
8881935.32
8883124.19
8883448.43
8884096.90
8884529.22
8884745.37
8885177.69
8885826.16
8886906.95
8887555.43
8888203.90
8888744.29
8889284.69
8889825.08
8890041.24
8891230.11
8892310.90
8892743.21
8893715.92
8894256.32
8895769.42
8896873.11
8903636.64
8903644.48
8906036.92
8912521.66
8915872.10
8926247.68
8931975.86
8938068.08
8942459.52
8945593.80
8947971.54
8949700.80
8949700.80
8949700.80
8949544.36
8949376.57];

% xv = [635575.77
% 688718.67
% 692404.66
% 657863.24
% 662270.00
% 671110.00
% 673570.00
% 671920.00
% 671960.00
% 664298.55
% 664326.52
% 662380.00
% 662340.00
% 653240.00
% 653230.00
% 649006.90
% 642290.93
% 635575.77];
% 
% yv = [8949635.57
% 8949430.30
% 8940689.92
% 8940662.15
% 8924600.00
% 8907900.00
% 8895100.00
% 8892800.00
% 8869800.00
% 8869800.00
% 8875780.41
% 8875788.41
% 8876900.00
% 8877100.00
% 8911500.00
% 8936290.37
% 8933662.39
% 8949635.57];

% Tanzania Block 3 Polygon
% xv = [638079.01
% 616591.06
% 606476.20
% 641453.01
% 672954.15
% 683680.74
% 668103.64
% 649928.29
% 644000.99
% 639741.26
% 636100.63
% 638079.01];
% 
% yv = [9021703.90
% 9016716.11
% 9060287.78
% 9060171.21
% 9060066.22
% 9013682.83
% 9009988.54
% 9023908.23
% 9049263.22
% 9030813.33
% 9029965.06
% 9021703.90];

% xv = [
%  962647
%  923218
%  851819
%  905189
%  962647
%  ];
% 
% yv = [
%  5929028
%  5883670
%  5945736
%  6007131
%  5929028
%  ];
% 
xv = xv.*100;
yv = yv.*100;
xv = int64(xv);
yv = int64(yv);

pltd = 200;



load(filepathin);

red = [repmat(1,50,1);(1:-0.02:0)'];
blue = [(0:0.02:1)';repmat(1,50,1)];
green = [(0:0.02:1)';(0.98:-0.02:0)'];
redblue = [red green blue];

scrsz = get(0,'ScreenSize');
cjf = figure('Position',[50 scrsz(4) scrsz(3)/1.1 (scrsz(4)/1.1)],'PaperType','A1','Visible','on','PaperPositionMode','auto','PaperOrientation','landscape');
colormap(redblue);

% colrngcj = (0.5*max(qcslice(3,1:2:endcj)));
% %props = {'LineStyle','none','Marker','o','MarkerEdge','b','MarkerSize',6};
% %line([qcslice(1,1:10:endcj),qcslice(1,1:10:endcj)],[qcslice(2,1:10:endcj),qcslice(2,1:10:endcj)],props{:});
% subplot(1,2,1)
% scatter(qcslice(4,1:pltd:endcj),qcslice(2,1:pltd:endcj),20,qcslice(3,1:pltd:endcj),'filled');
% caxis([(-1*colrngcj) colrngcj]);
% colorbar;
% hold on
% plot(xv,yv)
% hold off


%out_file_name_slice_tmp = strcat(out_file_name_slice,'_whole.pdf');

%print(cjf,'-dpdf','-noui','-r200',out_file_name_slice_tmp);


%tr_use = inpolygon(int64(qcslice(1,1:pltd:endcj)),int64(qcslice(2,1:pltd:endcj)),xv,yv);
%     %use logical to set the third dimension of the trace length to zero
%     tracesout = bsxfun(@times,traces,tr_use');

figure(1)
colrngcj = (0.5*max(qcslice(4,1:2:endcj)));
subplot(2,1,1); scatter(qcslice(1,1:pltd:endcj),qcslice(2,1:pltd:endcj),20,qcslice(4,1:pltd:endcj),'filled'); hold on; plot(xv,yv,'--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10); 
hold off
caxis([(-1*colrngcj) colrngcj]);
colorbar;

colrngcj = (0.5*max(qcslice(6,1:2:endcj)));
subplot(2,1,2); scatter(qcslice(1,1:pltd:endcj),qcslice(2,1:pltd:endcj),20,qcslice(6,1:pltd:endcj),'filled'); hold on; plot(xv,yv,'--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10);
hold off
caxis([(-1*colrngcj) colrngcj]);
colorbar;

% subplot(1,2,2)
% scatter(qcslice(1,1:pltd:endcj),qcslice(2,1:pltd:endcj),20,qcslice(4,1:pltd:endcj),'filled');
% %scatter(qcslice(1,1:pltd:endcj),qcslice(2,1:pltd:endcj),20,qcslice(5,1:pltd:endcj).*tr_use,'filled');
% caxis([(-1*colrngcj) colrngcj]);
% colorbar;
% hold on
% plot(xv,yv)
% hold off
%colorbar('peer',axes1);
%out_file_name_slice_tmp = strcat(out_file_name_slice,'_sub_polygon.pdf');
%print(cjf,'-dpdf','-noui','-r200',out_file_name_slice_tmp);

%close cjf;
end

