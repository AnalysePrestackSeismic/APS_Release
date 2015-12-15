clear all

% gamma_data=0.5;
% k_data=1;
% m_data=0.5;
% 
% Rp_data=1;
% Rs_data=1;
% Rd_data=1;
% 
% counter=1;
% for theta=(0:5*pi/180:45*pi/180)
%     c1_data(counter,1)=1+(tan(theta)*tan(theta));
%     c2_data(counter,1)=(-8)*(gamma_data^2)*(tan(theta)*tan(theta));
%     c3_data(counter,1)=((-0.5)*(tan(theta)*tan(theta)))+(2*(gamma_data^2)*(sin(theta)*sin(theta)));
%     counter=counter+1;    
% end
% 
% c1_data=(0.5*c1_data)+(0.5*k_data*c2_data)+(m_data*c3_data);
% c2_data=(0.5*c2_data);
% 
% R_data=(c1_data*Rp_data)+(c2_data*Rs_data)+(c3_data*Rd_data);

DELIMITER = '\t';
HEADERLINES = 1;

% Import the file
newData1 = importdata('gather_1626-5639', DELIMITER, HEADERLINES);

% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

length(data)

R_data_all = data(:,4:length(data))';

gamma_inv=1.9;
k_inv=0.1;
m_inv=0.1;

counter=1;
for theta=(10*pi/180:10*pi/180:40*pi/180)
    c1_inv(counter,1)=1+(tan(theta)*tan(theta));
    c2_inv(counter,1)=(-8)*(gamma_inv^2)*(tan(theta)*tan(theta));
    c3_inv(counter,1)=((-0.5)*(tan(theta)*tan(theta)))+(2*(gamma_inv^2)*(sin(theta)*sin(theta)));
    counter=counter+1;
end

c1_inv=(0.5*c1_inv)+(0.5*k_inv*c2_inv)+(m_inv*c3_inv);
c2_inv=(0.5*c2_inv);

G=[c1_inv,c2_inv,c3_inv];

Rp0=1;
Rs0=0;
Rd0=0;

m0=[Rp0,Rs0,Rd0]';
R_model=zeros(length(m0),100);
R_est=zeros(length(m0),100);

model = cell(length(data)-3,1);

% data loop
for n=1:length(data)-3
    fprintf('Number of iterations: %d\n',n);
    
    R_data = R_data_all(n,:)';

    R_model(:,1)=G'*G*m0;
    R_est(:,1) = (inv(G'*G))*G'*R_data;

    res=zeros(length(m0),100);

    res(:,1)=(G'*R_data)-G'*G*R_model(:,1);
    alpha=(res(:,1)'*res(:,1))/(res(:,1)'*G'*G*res(:,1));
    res(:,2)=res(:,1)-(alpha*G'*G*res(:,1));
    dir=res(:,1);
    R_model(:,2)=R_model(:,1)+(alpha*dir);

    for i=3:100
        beta=(res(:,i-1)'*res(:,i-1))/(res(:,i-2)'*res(:,i-2));
        dir=res(:,i-1)+(beta*dir);
        alpha=(res(:,i-1)'*res(:,i-1))/(dir'*G'*G*dir);
        R_model(:,i)=R_model(:,i-1)+(alpha*dir);
        res(:,i)=res(:,i-1)-(alpha*G'*G*dir);
        if res(:,i) < 1e-14;
            fprintf('Number of iterations 2: %d\n',i);
            model{n} = R_model(:,i);
                               
%             figure(1)
%             plot(R_data)
%             hold all
%             plot(G*R_model(:,1))
%             plot(G*R_model(:,i))
%             hold off
            
        end
    end 
    data_model(n,:) = G*model{n};    
end

figure(1)
imagesc(R_data_all);
figure(2)
imagesc(data_model);
figure(3)
imagesc(R_data_all-data_model);




