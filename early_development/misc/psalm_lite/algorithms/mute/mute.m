function mute(data)
% Click the left mouse button to define a point
% Drag the mouse to draw a line to the next point and
% left click again
% Right click the mouse to stop drawing
clear global bob 
clear global x_var 
clear global y_var
% 
%figure('WindowButtonDownFcn',@wbdcb)

% 
%imageHandle = imagesc(data);

f = figure;

%get(f)

global bob 
global x_var 
global y_var
bob = 1;
x_axis = size(data,2);
y_axis = size(data,1);
imageHandle = imagesc([0:1:x_axis],[0:1:y_axis],data);
colormap(gray)
%axis([0 5 0 2000])

set(imageHandle,'HitTest','off')
counter = 1;


set(gca,'ButtonDownFcn',{@mytestcallback,counter,f}); %'disp(''axis callback'')')

%uicontrol('Style', 'pushbutton','Position', [20 340 100 50],'String','Apply mute','Callback', 'uiresume(gcbf)');



%set(imageHandle,'ButtonDownFcn','disp(''axis callback'')')
% [x,y] = ginput;

% H_ax = gca;

%set(f,'WindowButtonDownFcn',@mytestcallback)
%set(imageHandle,'ButtonDownFcn',@mytestcallback);


%ah = axes('DrawMode','fast');
%axis ([1 10 1 10])

% % clc;clear;
% imObj = data;
% figure;
% hAxes = axes();
% imageHandle = imagesc(imObj);
% set(imageHandle,'ButtonDownFcn',@ImageClickCallback);
% 
% function ImageClickCallback ( objectHandle , eventData )
% axesHandle  = get(objectHandle,'Parent');
% coordinates = get(axesHandle,'CurrentPoint'); 
% coordinates = coordinates(1,1:2);
% message     = sprintf('x: %.1f , y: %.1f',coordinates (1) ,coordinates (2));
% helpdlg(message);
% end
% 
% end

function mytestcallback(hObject,evnt,count,f)
    
        %if strcmp(get(hObject,'SelectionType'),'alt')
        
        
        click_type=get(f,'SelectionType');        
        
        if strcmp(click_type,'normal')
        pos=get(hObject,'CurrentPoint');    

        
        

        x_var(bob) = pos(1,1);
        y_var(bob) = pos(1,2);
        
        str = ['(',num2str(x_var,'%0.3g'),', ',num2str(y_var,'%0.3g'),')'];
        
        if bob == 1
            h = line('XData',x_var,'YData',y_var,'Color','r','LineWidth',1,'Marker','+','MarkerSize',10,'MarkerEdgeColor','b');
        else
            h = line('XData',[x_var(bob) x_var(bob-1)],'YData',[y_var(bob) y_var(bob-1)],'Color','r','LineWidth',1,'Marker','+','MarkerSize',10,'MarkerEdgeColor','b');
        end   
        
%         if bob == 5;
%            delete(h(4)); redraw 
%         else
%         end
               
         disp(['event ',num2str(bob)]);
    
    disp(['You clicked X:',num2str(pos(1,1)),', Y:',num2str(pos(1,2))]);
    assignin('base', 'var1', x_var)
    assignin('base', 'var2', y_var)
   
%         uicontrol('Style', 'pushbutton','Position', [20 340 100 50],'String','Apply mute','Callback', @applymute);
% disp('Step one - please press button');
% drawnow; % Don't know, if this is needed.
% %uiwait(figH);
% %disp('Step two - please press button');
% %drawnow; % Don't know, if this is needed.
% %uiwait(figH);
% %disp('ready')


        
       
        %set(gca,'ButtonDownFcn',{@mytestcallback,count});
%     else
%         x(count) = pos(1,1);
%         y(count) = pos(1,2);
%         str = ['(',num2str(x,'%0.3g'),', ',num2str(y,'%0.3g'),')'];
%         line('XData',[x(end-1) x(end)],'YData',[y(end-1) y(end)],'Color','r','LineWidth',4,'Marker','.');
%         count = count + 1;
%         set(gca,'ButtonDownFcn',{@mytestcallback,count});
%     end
    
        %else
        %end
    
    %disp(['You clicked X:',num2str(x),', Y:',num2str(y)]);
    
%    set(gca,'ButtonDownFcn',{@mytestcallback,count});

        elseif strcmp(click_type,'alt')
            if bob == 1
                bob = 1;
                disp('Nothing to delete');
            else
                %pos=get(hObject,'CurrentPoint');
                
                new_x_var = x_var(1:end-1);
                new_y_var = y_var(1:end-1);
                
                % find nearest point
                disp(['You deleted X:',num2str(x_var(1:end)),', Y:',num2str(y_var(1:end))]);
               
                x_var = new_x_var;
                y_var = new_y_var;
                
                bob = bob - 1;
                %get(f)
                get(gca);
                h = findobj(gca,'Type','line');
                delete(h)
                line('XData',x_var,'YData',y_var,'Color','r','LineWidth',1,'Marker','+','MarkerSize',10,'MarkerEdgeColor','b');
            end
            
        end
        
         bob = bob  + count;

end
%     function applymute(hObject,evnt,data)
%         
%         f2 = figure;
%         
%         
%         
%         imageHandle = imagesc([0:1:x_axis],[0:1:y_axis],data);
%         
%         disp('Applying the mute')
%         
%     end


% function wbmcb(src,evnt)
%     [xn,yn,str] = disp_point(ah);
%     xdat = [x,xn];
%     ydat = [y,yn];
%     set(hl,'XData',xdat,'YData',ydat);
% end  
% 
% function [x,y,str] = disp_point(ah)
%     cp = get(ah,'CurrentPoint');  
%     x = cp(1,1);y = cp(1,2);
%     str = ['(',num2str(x,'%0.3g'),', ',num2str(y,'%0.3g'),')']; 
% end
% 
%    function wbdcb(src,evnt)
%       if strcmp(get(src,'SelectionType'),'normal')       
%          [x,y,str] = disp_point(ah);
%          hl = line('XData',x,'YData',y,'Marker','.');
%          text(x,y,str,'VerticalAlignment','bottom');drawnow
%          set(src,'WindowButtonMotionFcn',@wbmcb)
%       elseif strcmp(get(src,'SelectionType'),'alt')
%          set(src,'WindowButtonMotionFcn','')
%          [x,y,str] = disp_point(ah);
%          text(x,y,str,'VerticalAlignment','bottom');drawnow
%       end
%       function wbmcb(src,evnt)
%          [xn,yn,str] = disp_point(ah);
%          xdat = [x,xn];
%          ydat = [y,yn];
%          set(hl,'XData',xdat,'YData',ydat);
%       end  
%    end

end
