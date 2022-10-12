%% Script to assign interseismic slip to a fault

% from 15-24km depth

% Calculating total length of the fault
for i=1:length(utm_x)-1
   seg_length(i)=sqrt((utm_y(i)-utm_y(i+1))^2+(utm_x(i)-utm_x(i+1))^2);
   sum_length(1)=0;
   sum_length(i+1)=sum_length(i)+seg_length(i);
end
total_length=sum(seg_length);

% Importing the structural data
if exist(data_file)==2
    field_data=dlmread(data_file);
    data_lon=field_data(:,1);
    data_lat=field_data(:,2);
    throw=field_data(:,5);
    dip=field_data(:,11);
elseif exist(data_file)==0
    disp('No data file available')
else
end
% Omitting data lines with nan throw values
c=find(~isnan(throw));
data_lon=data_lon(c);
data_lat=data_lat(c);
throw=throw(c);
dip=dip(c);

if isnan(data_lon)==0 & isnan(data_lat)==0
    % Calculating the slip 
    for i=1:length(throw)
        if isnan(dip(i))==0
            slip(i)=throw(i)/sind(dip(i));
        else
            slip(i)=throw(i)/sind(shear_zone_dip);
        end
    end
    if length(utm_x)==2
        interseismic_slip=slip/15000;
    else
    % Finding the elements which contain the throw data
    for i=1:length(data_lon)
        a=find(utm_x>=data_lon(i),1,'first');
        b=find(utm_x<=data_lon(i),1,'last');
        if isempty(a)==1
            a=find(utm_y<=data_lat(i),1,'first');
            b=find(utm_y>=data_lat(i),1,'last');
        else
        end
        first_seg=sqrt((utm_x(b)-data_lon(i))^2+(utm_y(b)-data_lat(i))^2);
        second_seg=sqrt((utm_x(a)-data_lon(i))^2+(utm_y(a)-data_lat(i))^2);
        first_proportion=first_seg/(first_seg+second_seg);
        distances(i)=b*grid_size+first_proportion*grid_size;
    end

    % Calculating the mid points of all the boxes
    L=length(utm_x(1,:));
    d2=grid_size/2;
    if length(x_points(:,1))>3
        for i=1:L-3
            d(i)=d2+grid_size*i;
        end
    else
    end
    length_last=sqrt((utm_x(L-1)-utm_x(L))^2+(utm_y(L-1)-utm_y(L))^2);
    if length(x_points(:,1))>3
        distances_calc=[d2,d,((L-2)*grid_size+length_last/2)];
    else
        distances_calc=[d2,((L-2)*grid_size+length_last/2)];
    end
    distances_calc=sort(distances_calc).';

    data_dist=[0,distances,total_length].';
    data_throw=[0,(slip./15000),0].';

    interseismic_slip=interp1q(data_dist,data_throw,distances_calc);
    interseismic_slip=interseismic_slip.';
    end
elseif isnan(data_lon)==1 & isnan(data_lat)==1 % for faults where the throw has been calculated and assigned to the center
    % calculating the slip
    if isnan(dip)==0
    	slip=throw/sind(dip);
    else
    	slip=throw/sind(constant_dip);
    end
    % Calculating the mid points of all the boxes
    L=length(utm_x(1,:));
    d2=grid_size/2;
    if length(x_points(:,1))>3
        for i=1:L-3
            d(i)=d2+grid_size*i;
        end
    else
    end
    length_last=sqrt((utm_x(L-1)-utm_x(L))^2+(utm_y(L-1)-utm_y(L))^2);
    if length(x_points(:,1))>3
        distances_calc=[d2,d,((L-2)*grid_size+length_last/2)];
    else
        distances_calc=[d2,((L-2)*grid_size+length_last/2)];
    end
    distances_calc=sort(distances_calc).';

    data_dist=[0,(total_length/2),total_length].';
    data_throw=[0,(slip./15000),0].';
    interseismic_slip=interp1q(data_dist,data_throw,distances_calc);
    interseismic_slip=interseismic_slip.';
else
    disp('Interseismic slip calculations have gone wrong...')
end

