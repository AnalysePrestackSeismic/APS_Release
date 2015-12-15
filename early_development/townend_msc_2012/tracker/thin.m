function thin_slice=thin(slice)

%slice=traces.data(1:500,1:500);

%create thinned fault slice
thin_slice=zeros(size(slice));
for i=size(slice,1):-1:1
    maxj=0;
    for j=1:size(slice,2)
        if slice(i,j) == 10 && j>maxj
            on=1;
            maxj=j;
            if maxj+1<size(slice,2)
                while slice(i,maxj+1)==10 && maxj+1<size(slice,2)
                    maxj=maxj+1;
                end
            end
            maxj=maxj;
            width=maxj-j+1;
            if mod(width,2)==0
                thin_slice(i,j+width/2-1)=10;
                thin_slice(i,j+width/2)=10;
            else
                thin_slice(i,j+(width-1)/2)=10;
            end
        end
    end
end

% imagesc(thin_slice)
% figure
% imagesc(slice)