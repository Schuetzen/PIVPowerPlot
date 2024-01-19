clear all
close all
clc

load velocity_field_removing_bubble_2023_12_01_11_32.mat

N = length(vel);
[m,n] = size(vel(1).u);
c_threshold = 0.7;

count = 0;
for i=2:N
    % read data with good pairs
    if ~isempty(vel(i).u)
        count = count + 1;
        ui = vel(i).u;
        vi = vel(i).v;
        ci = vel(i).c;
    
        % remove all data with low correlation
        flag = ci<c_threshold;
        ui(flag) = NaN;
        vi(flag) = NaN;
        
        % save all data in a three-dimensional array
        u(:,:,count) = ui;
        v(:,:,count) = vi;        
    end
end

load velocity_field_removing_bubble_2023_12_01_11_39.mat

N = length(vel);

for i=1:N
    % read data with good pairs
    if ~isempty(vel(i).u)
        count = count + 1;
        ui = vel(i).u;
        vi = vel(i).v;
        ci = vel(i).c;
    
        % remove all data with low correlation
        flag = ci<c_threshold;
        ui(flag) = NaN;
        vi(flag) = NaN;
        
        % save all data in a three-dimensional array
        u(:,:,count) = ui;
        v(:,:,count) = vi;        
    end
end

load velocity_field_removing_bubble_2023_12_01_11_42.mat

N = length(vel);

for i=1:N
    % read data with good pairs
    if ~isempty(vel(i).u)
        count = count + 1;
        ui = vel(i).u;
        vi = vel(i).v;
        ci = vel(i).c;
    
        % remove all data with low correlation
        flag = ci<c_threshold;
        ui(flag) = NaN;
        vi(flag) = NaN;
        
        % save all data in a three-dimensional array
        u(:,:,count) = ui;
        v(:,:,count) = vi;        
    end
end

load velocity_field_removing_bubble_2023_12_01_11_45.mat

N = length(vel);

for i=1:N
    % read data with good pairs
    if ~isempty(vel(i).u)
        count = count + 1;
        ui = vel(i).u;
        vi = vel(i).v;
        ci = vel(i).c;
    
        % remove all data with low correlation
        flag = ci<c_threshold;
        ui(flag) = NaN;
        vi(flag) = NaN;
        
        % save all data in a three-dimensional array
        u(:,:,count) = ui;
        v(:,:,count) = vi;        
    end
end

save velocity_field_all_after_clean u v xphy yphy