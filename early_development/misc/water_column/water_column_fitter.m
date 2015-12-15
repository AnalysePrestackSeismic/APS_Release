function fit = water_column_fitter(filename,sheet)

close all
% Import excel file
%file_name = 'vel_data.xls';
%vel_data = importdata(filename);
[vel_data.data] = xlsread(filename,sheet)

% Work out which sheets to use
% Sheets in xls file
sheet_names = fieldnames(vel_data);
for i=1:1:length(sheet_names)
   fprintf('%d - Name: %s\n',i,sheet_names{i}); 
end

% which sheets to use?
sheet_index = input('Enter number of sheet to use in brackets []: ');
for i=1:1:length(sheet_index)
    sheet{i}.name = cell2mat(sheet_names(sheet_index(i)));
    sheet{i}.data = getfield(vel_data.data,sheet_names{sheet_index(i)});
    fprintf('For sheet %s\n',sheet{i}.name);
    sheet{i}.indomain = input('Is the domain time or depth? [d - depth, t - time]: ');% Dependent variable
    fprintf('First row of data in sheet %s:\n',sheet{i}.name);
    disp(sheet{i}.data(1,:))
    %fprintf('%d %d \n',sheet{i}.data(1,:));
    sheet{i}.yvariable = input('Which column is velocity in? ');
    if (sheet{i}.indomain == 'd')
        sheet{i}.xvariable = input('Which column is depth in? ');
    else
        sheet{i}.xvariable = input('Which column is time in? ');
    end
    sheet{i}.outdomain = input('Which domain to fit in: [d - depth, t - time]: ');% Dependent variable
end

% Plot imported data
s = length(sheet_index);
figure(1)
for i=1:1:length(sheet_index)
    subplot(1,s,i), plot(sheet{i}.data(:,sheet{i}.xvariable),sheet{i}.data(:,sheet{i}.yvariable))
    grid on
end

% Take average
if (s > 1)
    take_average = input('Take average of plots for fitting? 1 - Yes, 0 - No: ');
    if (take_average == 1)
        sheet_i = s+1;
        sheet{sheet_i} = sheet{1};
        sheet{sheet_i}.name = 'average';
        sheet{sheet_i}.data = mean(sheet{1:s}.data);
    end
else
    sheet_i = input('Enter number of sheet to use: ');
end

%Conversion
if(sheet{sheet_i}.outdomain ~= sheet{sheet_i}.indomain)
   n = length(sheet{sheet_i}.data);
   fprintf('Converting to common domain\n') 
   if (sheet{sheet_i}.indomain == 'd')
       starti = 0;
       for i=1:1:n
            tmp_dist = (sheet{sheet_i}.data(i,sheet{sheet_i}.xvariable) - starti)/(sheet{sheet_i}.data(i,sheet{sheet_i}.yvariable));
            tmp_dist = 2*tmp_dist;
            starti = sheet{sheet_i}.data(i,sheet{sheet_i}.xvariable);
            if (i == 1)
                sheet{sheet_i}.data(i,end+1) = tmp_dist;
            else
                sheet{sheet_i}.data(i,end) = sheet{sheet_i}.data(i-1,end) + tmp_dist;
            end
       end
   elseif (sheet{sheet_i}.indomain == 't')
       starti = 0;
       for i=1:1:n
            tmp_dist = (sheet{sheet_i}.data(i,sheet{sheet_i}.xvariable) - starti)*(sheet{sheet_i}.data(i,sheet{sheet_i}.yvariable));
            starti = sheet{sheet_i}.data(i,sheet{sheet_i}.xvariable);
            if (i == 1)
                sheet{sheet_i}.data(i,end+1) = tmp_dist;
            else
                sheet{sheet_i}.data(i,end) = sheet{sheet_i}.data(i-1,end) + tmp_dist;
            end
       end
   end
   % create column with required domain values
else
    sheet{sheet_i}.data(:,end+1) = sheet{sheet_i}.data(:,sheet{sheet_i}.xvariable);
end

% Try working out sections for user
% deri = diff(sheet{sheet_i}.data(:,sheet{sheet_i}.yvariable));
% figure
% plot(sheet{sheet_i}.data(1:end-1,end),deri)
% hold all
% plot(sheet{sheet_i}.data(:,end),sheet{sheet_i}.data(:,sheet{sheet_i}.yvariable));

% y = filter(b,a,x); could perform moving average filter on data, then take
% difference to find break points (might even need second difference)

% Ask for plotting information
n_fits = input('How many sections do you want to divide the graph into? ');

for i=1:1:n_fits
    fprintf('For section %d ',i);
    fit_vals = input('please enter start and end x-values for sub-division []: ');
    fit{i}.xfit = [fit_vals(1) fit_vals(2)];
    fit{i}.index = [find(sheet{sheet_i}.data(:,sheet{sheet_i}.xvariable)==fit_vals(1)) find(sheet{sheet_i}.data(:,sheet{sheet_i}.xvariable)==fit_vals(2))];
    fit{i}.yvals = [sheet{sheet_i}.data(fit{i}.index(1),sheet{sheet_i}.yvariable) sheet{sheet_i}.data(fit{i}.index(2),sheet{sheet_i}.yvariable)];
    fit{i}.xvals = [sheet{sheet_i}.data(fit{i}.index(1),end) sheet{sheet_i}.data(fit{i}.index(2),end)];
    fit{i}.data = [sheet{sheet_i}.data(fit{i}.index(1):fit{i}.index(2),end) sheet{sheet_i}.data(fit{i}.index(1):fit{i}.index(2),sheet{sheet_i}.yvariable)];    
end

%wb_q = input('Do you want the function to reach a certain velocity at the water bottom?');

% Plot sections
figure(2)
for i=1:1:n_fits
   subplot(1,n_fits,i), plot(fit{i}.data(:,1),fit{i}.data(:,2))
   grid on
   axis tight
   fit{i}.type = input('Based on plot enter desired fit type [1 - average, 2 - linear, 3 - polynomial]: ');
   fit{i}.typename = get_fit_type_name(fit{i}.type);
   if fit{i}.type == 3
      fit{i}.order = input('Enter order of polynomial to fit: '); 
   end
   title({['Section ',int2str(i)];['Fit: ',fit{i}.typename]});
end

% Perform initial fits as independent sections
fprintf('Performing independent fits of sections \n');
for i=1:1:n_fits
   fit{i}.m{1} = perform_inp_fit(fit{i},1e-08,5000,0);
end

figure(3)
for i=1:1:n_fits
   subplot(1,n_fits,i), plot(fit{i}.data(:,1),fit{i}.data(:,2))
   axis tight
   title({['Section ',int2str(i)];['Fit: ',fit{i}.typename]});
   hold all
   y = plot_model(fit{i},1);
   plot(fit{i}.data(:,1),y);
   hold off
   %legend('on')
   grid on
end

print -painters -dpdf -r600 inp_result.pdf

fprintf('Performing dependent fits of sections \n');
fit = perform_dep_fit(fit,1e-08,5000,0);

figure(4)
h(1) = plot(sheet{sheet_i}.data(:,end),sheet{sheet_i}.data(:,sheet{sheet_i}.yvariable));
s = cell(1, n_fits+1);
s{1} = sprintf('Excel sheet: %s',sheet_names{sheet_i});
hold all
for i=1:1:n_fits
    y = plot_model(fit{i},2);
    h(i+1) = plot(fit{i}.data(:,1),y);
    if (fit{i}.type == 3)
        s{i+1} = sprintf('Section %d, fit type: %s order %d, range: %s, model: %s',i,fit{i}.typename,fit{i}.order,mat2str(fit{i}.bound_xvalue,4),mat2str(fit{i}.m{2}(1:fit{i}.order+1),8));
    elseif (fit{i}.type == 2)
        s{i+1} = sprintf('Section %d, fit type: %s,  range: %s, model: %s',i,fit{i}.typename,mat2str(fit{i}.bound_xvalue,4),mat2str(fit{i}.m{2}(1:2),8));
    else
        s{i+1} = sprintf('Section %d, fit type: %s, range: %s, model: %s',i,fit{i}.typename,mat2str(fit{i}.bound_xvalue,4),mat2str(fit{i}.m{2},8));
    end 
end
hold off;
title('Final fit');
ylabel('y')
xlabel('x')
legend(h, s,'Location','SouthOutside');
grid on

% print dtect statement

print -painters -dpdf -r600 dep_result.pdf
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fit_m = perform_inp_fit(fit,tol,maxit,~)
   
        switch fit.typename
            case 'average'
                fit_m = mean(fit.data(:,2));
            case 'linear'
                G = horzcat(ones(length(fit.data(:,1)),1),fit.data(:,1));
                [fit_m] = lsqr(G,fit.data(:,2),tol,maxit);
            case 'polynomial'
                G = ones(length(fit.data(:,1)),fit.order+1);
                for i=2:fit.order+1
                    G(:,i) = fit.data(:,1).^(i-1);
                end
                [fit_m] = lsqr(G,fit.data(:,2),tol,maxit);               
        end
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fit = perform_dep_fit(fit,tol,maxit,~)

        n_sections = length(fit);
            for i=1:1:n_sections
                % Fit average sections first
                if (fit{i}.type == 1)
                    avg = mean(fit{i}.data(:,2));
                    fit{i}.bound_yvalue = [avg avg];
                    fit{i}.bound_xvalue = [fit{i}.xvals(1) fit{i}.xvals(2)];
                    fit{i}.m{2} = avg;
                    % need to set avg value to fit before and after
                    if (i ~= 1 && i~= n_sections)
                        fit{i-1}.bound_yvalue(2) = avg;
                        fit{i+1}.bound_yvalue(1) = avg;
                    elseif i == n_sections
                        fit{i-1}.bound_yvalue(2) = avg;
                    elseif i == 1
                        fit{i+1}.bound_yvalue(1) = avg;
                    end
                    
                else % there are no averages to fit bound values are yvals and xvals
                    if (isfield(fit{i},'bound_yvalue') == 0)
                        fit{i}.bound_yvalue = [fit{i}.yvals(1) fit{i}.yvals(2)];
                    elseif (length(fit{i}.bound_yvalue) == 1)
                        fit{i}.bound_yvalue(2) = fit{i}.yvals(2);
                    end
                    
                    fit{i}.bound_xvalue = [fit{i}.xvals(1) fit{i}.xvals(2)];
                end
            end
            
            for i=1:1:n_sections
                if (fit{i}.type == 2) % linear
                    G = horzcat(ones(length(fit{i}.data(:,1)),1),fit{i}.data(:,1));
                    d = G'*fit{i}.data(:,2);
                    d = vertcat(d,fit{i}.bound_yvalue(1),fit{i}.bound_yvalue(2));
                    F = [1 fit{i}.bound_xvalue(1); 1 fit{i}.bound_xvalue(2)]; % 1 xvalue, bounds
                    V = zeros(size(G,2)+size(F,1),size(G,2)+size(F,1));
                    V(1:size(G,2),1:size(G,2)) = G'*G;
                    V(end-size(F,1)+1:end,1:size(F,2)) = F;
                    V(1:size(F,2),end-size(F,1)+1:end) = F';
                    [fit{i}.m{2}] = lsqr(V,d,tol,maxit);
                end
                if (fit{i}.type == 3) % poly
                    G = ones(length(fit{i}.data(:,1)),fit{i}.order+1);
                    for j=2:fit{i}.order+1
                        G(:,j) = fit{i}.data(:,1).^(j-1);
                    end
                    d = G'*fit{i}.data(:,2);
                    d = vertcat(d,fit{i}.bound_yvalue(1),fit{i}.bound_yvalue(2));
                    
                    F = ones(2,fit{i}.order+1); % 2 is number of bounds
                    for j=2:fit{i}.order+1
                        F(:,j) = fit{i}.bound_xvalue'.^(j-1);
                    end
                    V = zeros(size(G,2)+size(F,1),size(G,2)+size(F,1));
                    V(1:size(G,2),1:size(G,2)) = G'*G;
                    V(end-size(F,1)+1:end,1:size(F,2)) = F;
                    V(1:size(F,2),end-size(F,1)+1:end) = F';
                    [fit{i}.m{2}] = lsqr(V,d,tol,maxit);
                end
            end
            
end     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fit_name = get_fit_type_name(type)
% [1 - average, 2 - linear, 3 - poly]
    switch type
        case 1
            fit_name = 'average';
        case 2
            fit_name = 'linear';
        case 3
            fit_name = 'polynomial';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = plot_model(fit,index)

params = fit.type;

y = repmat(fit.m{index}(1),length(fit.data(:,1)),1);
    for i=2:1:params
       tmp = fit.m{index}(i).*(fit.data(:,1).^(i-1));
       y = y + tmp;
    end
end



% 
% 
% 
% depth_index = input('Based on the plot please enter values for sub-division []: ');
% for i=1:1:length(depth_index)
%     row_index(i) = find(sheet{sheet_i}.data(:,sheet{sheet_i}.xvariable)==depth_index(i));
% end
% 
% % Plot divisions
% n_divisions = length(row_index)+1;
% starti = 1;
% figure(2)
% for i=1:1:length(row_index)
%     endi = row_index(i);
%     subplot(1,n_divisions,i), plot(sheet{sheet_i}.data(starti:endi,sheet{sheet_i}.yvariable),sheet{sheet_i}.data(starti:endi,end))
%     fit{i}.type = input('Based on plot enter desired fit type [1 - average, 2 - linear, 3 - 2nd order poly., 4 - 3rd order poly,, 5 - 4th order poly]: ');
%     % get start and end from each division
%     starti = row_index(i);
%     fit{i}.division(1) = sheet{sheet_i}.data(starti,end);
%     fit{i}.division(2) = sheet{sheet_i}.data(starti,sheet{sheet_i}.yvariable);
%     
%     %starti = row_index(i);
%     
%     if i == length(row_index)
%         subplot(1,n_divisions,n_divisions), plot(sheet{sheet_i}.data(starti:end,sheet{sheet_i}.yvariable),sheet{sheet_i}.data(starti:end,end))
% %         division(n_divisions,1) = sheet{sheet_i}.data(end,end);
% %         division(n_divisions,2) = sheet{sheet_i}.data(end,sheet{sheet_i}.xvariable);
%         fit{i+1}.type = input('Based on plot enter desired fit type [1 - average, 2 - linear, 3 - 2nd order poly., 4 - 3rd order poly,, 5 - 4th order poly]: ');
%         fit{i+1}.division(1) = sheet{sheet_i}.data(starti,end);
%         fit{i+1}.division(2) = sheet{sheet_i}.data(starti,sheet{sheet_i}.yvariable);
%     
%     end
%     
% end
% 
% % Perform initial fits as independent sections
% start_index = 1;
% for i=1:1:length(row_index)
%    % find average values 
%    if (fit{i}.type == 1) % average
%         end_index = find(sheet{sheet_i}.data(:,end)==division(i,1));
%         fit{i}.values = mean(sheet{sheet_i}.data(start_index:end_index,sheet{sheet_i}.yvariable));
%         start_index = find(sheet{sheet_i}.data(:,end)==division(i,1));
%    end
% end
% 
% % Perform initial fits as independent sections
% 
% % Linear fit
% start_index = 1;
% for i=1:1:length(row_index)
% if (fit_type(i) == 2) % linear
% end_index = find(sheet{sheet_i}.data(:,end)==division(1,1));
% d = sheet{sheet_i}.data(start_index:end_index,sheet{sheet_i}.yvariable);
% G = horzcat(ones(length(sheet{sheet_i}.data(start_index:end_index,sheet{sheet_i}.xvariable)),1),sheet{sheet_i}.data(start_index:end_index,end));
% m = lsqr(G,d);
% 
% d = G'*d;
% d = vertcat(d,division(1,2));
% 
% % must fit
% %F = [1 division(1,1)];
% F = [1 3000];
% 
% V = zeros(length(G'*G)+length(F)-1,length(G'*G)+length(F)-1);
% V(1:length(G'*G),1:length(G'*G)) = G'*G;
% V(end,1:length(F)) = F;
% V(1:length(F),end) = F;
% 
% m2 = lsqr(V,d);
% 
% 
% 
% figure; plot(sheet{sheet_i}.data(start_index:end_index,end)*m(2)+m(1))
% hold all
% plot(sheet{sheet_i}.data(start_index:end_index,end)*m2(2)+m2(1))
% 
% 
% end
% 
% % end



% Save function as text file