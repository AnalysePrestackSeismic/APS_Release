clearvars -except cube Lcube_true

tic;

npeaks1=1;%100;     % no. of peaks in first hough space to pick up per y-pane
npeaks2=1;%100;     % no. of peaks in second hough space to pick up per theta1-pane
max_peaks1=0.1;%0.3;  % min fraction of max peak in hough space to accept
max_peaks2=0.8;%0.3;
neigh_peaks1=10;%50; % exclusion zone around peaks picked, bigger is smaller neighbourhood
neigh_peaks2=50;
fill_gap1=2;
min_length1=7;
fill_gap=5;%20;     % bridging length allowed between segments for houghlines
min_length=7;%30;%70;   % minimum length of lines found with houghlines

theta1_bins=200;    % range of theta1 values possible and binnage (-theta1 -> +theta1)
theta1_max=40;%20;

theta2_bins=300;    % same for theta2 values
theta2_max=89;

% Hough transform in (x,z), fixed y
% get Hough peaks for each y pane
Hloc=zeros(size(cube,2),npeaks1,2)+1;
for yj=1:size(cube,2)
    
    % count
    [1,size(cube,2),yj]
    
    [H(:,yj,:),T1,R1]=hough(squeeze(cube(:,yj,:)),...
        'Theta',-theta1_max:theta1_max*2/(theta1_bins-1):theta1_max);
    holdpeaks=houghpeaks(squeeze(H(:,yj,:)),npeaks1,...
        'threshold',ceil(max_peaks1*max(max(H(:,yj,:)))),...
        'nhoodsize',2*round(((size(squeeze(H(:,yj,:)))/neigh_peaks1)+1)/2)-1);
    for loc=1:size(holdpeaks,1)
        Hloc(yj,loc,:)=holdpeaks(loc,:);
    end
end

% check plot
% for i=1:100,imagesc(squeeze(H(:,i,:))),colorbar,hold on,...
% scatter(Hloc(i,:,2),Hloc(i,:,1),'r'),pause(0.1),hold off,end

% put Hough peaks into cube, removing invalid contributions where no. of
% peaks < no. of peaks requested (returns [1,1] location)
Hpeaks=zeros(size(H));
for yj=1:size(cube,2)
    for loc=1:size(Hloc,2)
        if sum(Hloc(yj,loc,:))~=2
            Hpeaks(Hloc(yj,loc,1),yj,Hloc(yj,loc,2))=10;
        end
    end
    if sum(sum(Hpeaks(yj,:,:)))==10,Hpeaks(yj,:,:)=0;end  
    %empty slices with only 1 peak
end

% check
% for i=1:100,subplot(1,2,1),imagesc(squeeze(H(:,i,:))),colorbar,...
% subplot(1,2,2),imagesc(squeeze(Hpeaks(:,i,:))),pause(0.1),end

% Hough transform in (y,rho1), fixed theta1
% get Hough peaks
H2loc=zeros(size(H,3),npeaks2,2)+1;
for theta1=1:size(H,3)
    
    % count
    [2,size(H,3),theta1]    

    [H2(:,:,theta1),T2,R2]=hough(squeeze(Hpeaks(:,:,theta1)),...
        'Theta',-theta2_max:theta2_max*2/(theta2_bins-1):theta2_max);
    holdpeaks=houghpeaks(squeeze(H2(:,:,theta1)),npeaks2,...
        'threshold',ceil(max_peaks2*max(max(H2(:,:,theta1)))),...
        'nhoodsize',2*round(((size(H2(:,:,theta1))/neigh_peaks2)+1)/2)-1);
    for loc=1:size(holdpeaks,1)
        H2loc(theta1,loc,:)=holdpeaks(loc,:);
    end
end

% check
% for i=1:300,imagesc(squeeze(H2(:,:,i))),colorbar,hold on,...
% scatter(H2loc(i,:,2),H2loc(i,:,1),'r'),pause(0.02),hold off,end

% same as before with Hough peaks
H2peaks=zeros(size(H2));
for theta1=1:size(H,3)
    for loc=1:size(H2loc,2)
        if sum(H2loc(theta1,loc,:))~=2
            H2peaks(H2loc(theta1,loc,1),H2loc(theta1,loc,2),theta1)=10;
        end
    end
%    if sum(sum(H2peaks(theta1,:,:)))==10,H2peaks(theta1,:,:)=0;end  
    %empty slices with only 1 peak
end

% check
% for i=1:300,subplot(1,2,1),imagesc(squeeze(H2(:,:,i))),colorbar,...
% subplot(1,2,2),imagesc(squeeze(H2peaks(:,:,i))),pause(0.02),end

% work through each theta1 pane and extract Hough lines
Hlines.point1=[];Hlines.point2=[];Hlines.theta=[];Hlines.rho=[];
Hlines.theta1=[];Hlines.H2loc=[];
for theta1=1:size(H,3)
%     holdpeaks=[];
%     on=0;
%     for loc=1:size(H2loc,2)
%         if sum(H2loc(theta1,loc,:))~=2
%             holdpeaks=[holdpeaks;squeeze(H2loc(theta1,loc,:))'];
%             on=1;
%         end
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for loc=1:size(H2loc,2)
        holdline=houghlines(squeeze(H(:,:,theta1)),T2,R2,...
            squeeze(H2loc(theta1,loc,:))','FillGap',fill_gap1,'MinLength',min_length1);
        if isempty(fieldnames(holdline))~=1
            holdline(1).theta1=theta1;
            holdline(1).H2loc=[theta1,loc];
            Hlines=[Hlines,holdline(1)];
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % get H line for each H2 peak and store in Hlines, getting rid of any
    % results that are empty in the process
%     if on==1
%         holdline=houghlines(squeeze(H(:,:,theta1)),T2,R2,holdpeaks,...
%             'FillGap',10,'MinLength',20);
%         for lin=1:size(holdline,2)
%             if isempty(fieldnames(holdline(lin)))~=1
%                 holdline(lin).theta1=theta1;
%                 Hlines=[Hlines,holdline(lin)];
%             end
%         end
%     end
end

% check
% imagesc(squeeze(H(:,:,90))),colorbar,hold on,plot([Hlines(2).point1(1),...
% Hlines(2).point2(1)],[Hlines(2).point1(2),Hlines(2).point2(2)],'g'),hold off


% Go through all H lines, get all Hloc values lying within dthres of 
% H line entry (fHloc - [theta1,yj,rho1])
fHloc=[];
dthres=1.5;
count=0;
for lin=1:size(Hlines,2)
    
    % counter
    [3,size(Hlines,2),lin,toc/60]
    
    if isempty(fieldnames(Hlines(lin)))~=1 && isempty(Hlines(lin).point1)~=1
        Q1=[Hlines(lin).theta1,Hlines(lin).point1];
        Q2=[Hlines(lin).theta1,Hlines(lin).point2];
        holdloc=[];
        for yj=1:size(Hloc,1)
            for loc=1:size(Hloc,2)
                d=DistBetween2Segment(Q1,Q2,[Hloc(yj,loc,2),yj,Hloc(yj,loc,1)],...
                    [Hloc(yj,loc,2),yj,Hloc(yj,loc,1)+0.001]);
                
                %scatter3(Hloc(yj,loc,2),yj,Hloc(yj,loc,1),'r'),hold on
                
                if d<dthres
                    holdloc=[holdloc;[Hloc(yj,loc,2),yj,Hloc(yj,loc,1)]];
                end
            end
        end
        
        %         plot3([Q1(1),Q2(1)],[Q1(2),Q2(2)],[Q1(3),Q2(3)])
        %         if isempty(holdloc)~=1
        %             scatter3(holdloc(:,1),holdloc(:,2),holdloc(:,3))
        %         end
        
        if isempty(holdloc)~=1
            count=count+1;
            fHloc(count).loc=holdloc;
            fHloc(count).lin=lin;
        end
    end
end

% Use each fHloc entry as houghpeaks locations for houghlines to get fault
% lines in xyz space
lines.point1=[];lines.point2=[];lines.theta=[];lines.rho=[];lines.yj=[];
lines.faultID=[];
for i=1:size(fHloc,2)
    
    %count
    [4,size(fHloc,2),i]
    
    for j=1:size(fHloc(i).loc,1)
        
        
        
%         check
%         [A,B,C]=hough(squeeze(cube(:,fHloc(i).loc(j,2),:)),...
%             'Theta',-theta1_max:theta1_max*2/(theta1_bins-1):theta1_max);
%         subplot(1,2,1),imagesc(squeeze(cube(:,:,fHloc(i).loc(j,2)))),...
%             subplot(1,2,2),imagesc(A),hold on,...
%             subplot(1,2,2),scatter(fHloc(i).loc(j,1),fHloc(i).loc(j,3),'r'),hold off
%         pause(0.1)
        
        holdline=houghlines(squeeze(cube(:,fHloc(i).loc(j,2),:)),T1,R1,...
            [fHloc(i).loc(j,3),fHloc(i).loc(j,1)],...
            'FillGap',fill_gap,'MinLength',min_length);
        for lin=1:size(holdline,2)
            if isempty(fieldnames(holdline(lin)))~=1
                holdline(lin).yj=fHloc(i).loc(j,2);
                holdline(lin).faultID=i;
                lines=[lines,holdline(lin)];
            end
        end
    end
end

% paint fault lines into xyz space, store individual fault locations into
% Floc
Floc.loc=[];
for lin=1:size(lines,2)
    
    % counter
    [5,size(lines,2),lin]
    
    if isempty(lines(lin).point1)~=1
        len=sqrt((lines(lin).point1(1)-lines(lin).point2(1))^2+...
            (lines(lin).point1(2)-lines(lin).point2(2))^2);
        y=lines(lin).yj;
        for d=0:0.5/len:1
            p=round(lines(lin).point1+d*(lines(lin).point2-lines(lin).point1));
            if lines(lin).faultID>size(Floc,2)
                Floc(lines(lin).faultID).loc=[p(2),y,p(1)];
            else
                Floc(lines(lin).faultID).loc=...
                    [Floc(lines(lin).faultID).loc;[p(2),y,p(1)]];
            end
            Floc(lines(lin).faultID).ID=lines(lin).faultID;
        end
    end
end

% check
% for i=1:100,subplot(1,2,1),imagesc(squeeze(cube(:,i,:))),...
% subplot(1,2,2),imagesc(squeeze(Fcube(1,:,i,:))),pause(0.1),end

% find centroid of every fault and average scalar distance from points to
% centroid, then discard faults outside certain multiple of average scalar
% % distance
% thres_s=1.5;
% for i=1:size(Floc,2)
%     if isempty(Floc(i).loc)~=1
%         centroid=sum(Floc(i).loc,1)/size(Floc(i).loc,1);
%         d=Floc(i).loc-repmat(centroid,size(Floc(i).loc,1),1);
%         s=(d(:,1).^2+d(:,2).^2+d(:,3).^2).^0.5;
%         avg_s=sum(s,1)/size(s,1);
%         ind=[];
%         for j=1:size(Floc(i).loc,1)
%             if s(j)>avg_s*thres_s,ind=[ind;j];end
%         end
%         Floc(i).loc(ind,:)=[];
%     end
% end

% get fault height from each fault
for i=1:size(Floc,2)
    if isempty(Floc(i).loc)~=1
        Floc(i).height=max(Floc(i).loc(:,1))-min(Floc(i).loc(:,1));
    else
        Floc(i).height=0;
    end
end

% sort fault volumes by increasing fault height
[tmp ind]=sort([Floc.height],'descend');
Floc=Floc(fliplr(ind));

% get number of points in each fault plane and sort
%for i=1:size(Floc,2),Floc(i).size=size(Floc(i).loc,1);end
%[tmp ind]=sort([Floc.size]);
%Floc=Floc(fliplr(ind));

% combine Fcubes into one cube, overlaying onto original cube
Lcube_hough=zeros(size(cube));
for i=1:size(Floc,2)
    for j=1:size(Floc(i).loc,1)
        Lcube_hough(Floc(i).loc(j,1),Floc(i).loc(j,2),Floc(i).loc(j,3))=...
            Floc(i).height;%*...
            %cube(Floc(i).loc(j,1),Floc(i).loc(j,2),Floc(i).loc(j,3))/10;
    end
end

toc

% check
% for i=1:100,subplot(1,2,1),imagesc(squeeze(Lcube_true(:,i,:))),colorbar,...
% subplot(1,2,2),imagesc(squeeze(Lcube_hough(:,i,:))),colorbar,pause(0.2),end


