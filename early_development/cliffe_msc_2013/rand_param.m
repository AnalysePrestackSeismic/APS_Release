clear all
% random parameters to be saved and used for the rest of the parameter
% testing

anomaly_multiplier = 1.5;
background_dominace = 0.8;  
nt = 1000;                                         
sparsity = 0.7;                                  
max_amplitude = 0.15;
ref=synref(nt,sparsity,max_amplitude); % make reflection coefficients, take out of loop, so randomising happens only once
random = rand(nt,1);
R1=rand(nt,1);        % series of nt random numbers, taken out of loop
R2 =rand(nt,1);
R3 = rand(nt,1);
R4=rand(2*nt,1);
R5=rand(nt*6,1);

class = nan(nt,1);                                      % class starts as a series of nt nans

class(ref~=0) = random(ref~=0);                         % where ref isn't zero, insert corresponding random number
class(class>=(1-background_dominace)) = 1;              % replace values in class vector with 1 or 2 depending on background dominance, ensures the right proportion of anomalous points
class(class<(1-background_dominace)) = 2;               
class_2_flag = 0;
for ii = 1:nt                                           %this loop makes anomalies, those with class 2 are flagged wtih a number 1, reflectivity is made negative, anomalous and squared
    if class(ii) == 2
        class_2_flag = 1;
        ref(ii) = -anomaly_multiplier*abs(ref(ii));
    end
    if class(ii) == 1                                   % class 1, 
        if class_2_flag == 1                            % if the previous loop has just made an anomalous value, the next one has class made =3, +ve anomalous and squared
            class(ii) = 3;
            ref(ii) = anomaly_multiplier*abs(ref(ii));  % ensures equal number of 2s and 3s, and equal +ve and -ve - white earth reflectivity assumption, could bias this at this point, break assumption
            class_2_flag = 0;                           % flag reset to zero and looped
        end
    end
end