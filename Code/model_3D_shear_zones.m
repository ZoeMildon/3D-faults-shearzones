format long


%% Writing out the beginning of the Coulomb output file
output_data_path='Output_files/filename_grid_sizekm.inr';
gkm=num2str(grid_size/1000);
output_data_file=strrep(output_data_path,'filename',filename);
output_data_file=strrep(output_data_file,'grid_size',gkm);

fid=fopen(output_data_file, 'wt');
% comments about the data file
fprintf (fid,'This is a file created by rectangularly gridding the faults.\n');
fprintf (fid,'Interseismic loading, the grid size of faults is %2.0f km\n',(grid_size/1000));
fprintf (fid,'#reg1=  0  #reg2=  0   #fixed= 17787  sym=  1\n');
fprintf (fid,' PR1=       .250      PR2=       .250    DEPTH=        5.0\n');
fprintf (fid,'  E1=   0.800000E+06   E2=   0.800000E+06\n');
fprintf (fid,'XSYM=       .000     YSYM=       .000\n');
fprintf (fid,'FRIC=       .400\n');
fprintf (fid,'S1DR=    19.0001     S1DP=     -0.0001    S1IN=    100.000     S1GD=   .000000\n');
fprintf (fid,'S3DR=    89.9999     S3DP=      89.999    S3IN=     30.000     S3GD=   .000000\n');
fprintf (fid,'S2DR=   109.0001     S2DP=     -0.0001    S2IN=      0.000     S2GD=   .000000\n');
fprintf (fid,'\n');
fprintf (fid,'  #   X-start    Y-start     X-fin      Y-fin   Kode  rake    net slip   dip angle     top      bot\n');
fprintf (fid,'xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx\n');

%% Import the text file containing all the fault names
disp('Select the text file with the list of faults')
[fault_filename,fault_filepath]=uigetfile('*.txt');
if fault_filename==0
    disp('No file selected')
else
% disp('Select the folder where the fault kml files are located')
% [kml_folder]=uigetdir('');
kml_folder='Fault_traces';
boxes_with_data=0;
sum_difference=0;
fault_names=importdata([fault_filepath fault_filename]);
for i=1:length(fault_names)
    fault_name=fault_names{i}
    kml_path='folder/fault_name.kml';
    kml_file=strrep(kml_path,'folder',kml_folder);
    kml_file=strrep(kml_file,'fault_name',fault_name);
    data_path=strcat(FAULT_DATA_PATH,'fault_name.txt');
    data_file=strrep(data_path,'fault_name',fault_name);
    
    % Import kml file and convert to UTM coordinates
    fileid=fopen(kml_file);
    if fileid < 0
        msgkml=sprintf('Cannot find the kml file for %s.\nFault not plotted.',fault_name);
        disp(msgkml)
    else
    [lat,lon]=read_kml(kml_file);
    [utm_lon,utm_lat] = wgs2utm(lat,lon,UTM_zone,UTM_letter);
    if utm_lon(1)<utm_lon(end)
    else
        %disp('flipped coordinates')
        utm_lon=flip(utm_lon);
        utm_lat=flip(utm_lat);
    end
    % Grid the fault, starting at the NW end, by grid_size (as the crow flies)
    last_point=0;
    n=1;
    utm_x(1)=utm_lon(1);
    utm_y(1)=utm_lat(1);
    b=0;
    % Finding the next grid point by hypotenuse method
    while last_point<1
        [utm_x(n+1),utm_y(n+1),last_point,a]=nextpoint_hyp(utm_x(n),utm_y(n),grid_size,utm_lon,utm_lat,b);
        n=n+1;
        b=a;
    end
    %% Extending the fault to depth
    utm_z(:,length(utm_x))=0; % Assuming all faults come to the surface (0m depth)
    
    % extracting the relevant dip value to use
     A=exist('shear_zone_dip');
    if A==0
        [fault_name_slip,dip_values]=textread(DIPS_FILE,'%s %f');
        a=strmatch(fault_name,fault_name_slip,'exact');
        constant_dip=dip_values(a);
        shear_zone_dip=dip_values(a);
    elseif A==1
    end
% shear_zone_dip=constant_dip;
    if isempty(constant_dip)==1
        msgdip=sprintf('Missing dip information for %s.\nFault not plotted.',fault_name);
        disp(msgdip)
    else
        % extracting the relevant rake value to use
        [fault_name_slip,rake_values]=textread(RAKES_FILE,'%s %f');
        b=strmatch(fault_name,fault_name_slip,'exact');
        rake=rake_values(b);
        if isempty(rake)==1
            msgrake=sprintf('Missing rake information %s.\nFault not plotted.',fault_name);
            disp(msgrake)
        else
            % extracting the relevant depth to use. If length<15km then aspect ratio=1.
            [fault_name_depth,vardepth]=textread(SHORT_FAULT_LENGTHS_FILE,'%s %f');
            c=strmatch(fault_name,fault_name_depth,'exact');
            if isempty(c)==1
                fault_down_dip_length=depth/sind(constant_dip);
            else
                fault_down_dip_length=vardepth(c)*-1000;
            end
            % calculating the grid size to use to depth (to ensure a whole number
            % of boxes, resulting elements will be rectangular rather than square
            if grid_size<=abs((depth/sind(constant_dip)))
                    n=abs(round(fault_down_dip_length/1000)); % whole number of boxes that will fit into the fault_down_dip_length
                    grid_size_to_depth=-fault_down_dip_length/n;
            else 
                    grid_size_to_depth=-fault_down_dip_length;
            end
            grid_size_surface=grid_size_to_depth*cosd(constant_dip);
            grid_size_depth=grid_size_to_depth*sind(constant_dip);
    
            % extracting the relevant projection direction to use
            [fault_name_proj,proj_dir]=textread(PROJECTION_DIRECTION_FILE,'%s %f');
            e=strmatch(fault_name,fault_name_proj,'exact');
            proj_dir=proj_dir(e);
            if constant_dip~=90 && isempty(proj_dir)==1
                msgproj=sprintf('Missing projection direction information for %s.\nFault not plotted.',fault_name);
                disp(msgproj)
            else
                delta_z=grid_size_depth;
                if isempty(proj_dir)==1
                    delta_x=0;
                    delta_y=0;
                else
                    delta_x=abs(grid_size_surface*sind(proj_dir));
                    delta_y=abs(grid_size_surface*cosd(proj_dir));
                end
                x_points=utm_x;
                y_points=utm_y;
                z_points=utm_z;
                for j=1:length(utm_x) % calculating the brittle portion of the fault
                    for i=1:n
                        if isempty(proj_dir)==1 % strike-slip vertical faults
                            dx=0;
                            dy=0;
                        elseif proj_dir>=0 && proj_dir<90 % north east dipping faults
                            dx=1;
                            dy=1;
                        elseif proj_dir>=90 && proj_dir<180 % south east dipping faults
                            dx=1;
                            dy=-1;
                        elseif proj_dir>=180 && proj_dir<270 % south west dipping faults
                            dx=-1;
                            dy=-1;
                        elseif proj_dir>=270 && proj_dir<=360 % north west dipping faults
                            dx=-1;
                            dy=1;
                        else
                            disp('Error: problem with the projection direction')
                        end
                        x_points(1+i,j)=x_points(i,j)+dx*delta_x;
                        y_points(1+i,j)=y_points(i,j)+dy*delta_y;
                        z_points(1+i,j)=z_points(i,j)-delta_z;
                    end
                end
                % calculating the ductile portion of the fault
                if z_points(end)<(depth+1)
                shear_down_dip_length=-(depth-depth_ductile)/sind(shear_zone_dip);
                if grid_size<=abs(((depth-depth_ductile)/sind(shear_zone_dip)))
                    n=abs(round(shear_down_dip_length/grid_size)); % whole number of boxes that will fit into the fault_down_dip_length
                    grid_size_shear=-shear_down_dip_length/n;
                else 
                    grid_size_shear=-shear_down_dip_length;
                    %n=1;
                end
                grid_size_surface_i=grid_size_shear*cosd(shear_zone_dip);
                grid_size_depth_i=grid_size_shear*sind(shear_zone_dip);

                delta_x_i=abs(grid_size_surface_i*sind(proj_dir-180));
                delta_y_i=abs(grid_size_surface_i*cosd(proj_dir-180));
                delta_z_i=grid_size_depth_i;
                c=length(x_points(:,1));
                d=length(x_points(:,1))+n-1;
                for k=1:length(x_points(1,:)) 
                    for l=c:d
                        if isempty(proj_dir)==1 % strike-slip vertical faults
                            dx=0;
                            dy=0;
                        elseif proj_dir>=0 && proj_dir<90 % north east dipping faults
                            dx=1;
                            dy=1;
                        elseif proj_dir>=90 && proj_dir<180 % south east dipping faults
                            dx=1;
                            dy=-1;
                        elseif proj_dir>=180 && proj_dir<270 % south west dipping faults
                            dx=-1;
                            dy=-1;
                        elseif proj_dir>=270 && proj_dir<=360 % north west dipping faults
                            dx=-1;
                            dy=1;
                        else
                            disp('Error: problem with the projection direction ductile zone')
                        end
                        x_points(1+l,k)=x_points(l,k)+dx*delta_x_i;
                        y_points(1+l,k)=y_points(l,k)+dy*delta_y_i;
                        z_points(1+l,k)=z_points(l,k)-delta_z_i;
                    end
                end
                    interseismic_slip_distribution;
                    a=find(z_points(:,1)>=depth+1);
                    slip_distribution=zeros((length(z_points(a,1))),(length(x_points(1,:))-1)); % 0 slip in the brittle zone
                    b=find(z_points(:,1)<depth);
                    slip_distribution=[slip_distribution;interseismic_slip];
                    slip_distribution=[slip_distribution;repmat(interseismic_slip,length(b),1)];
                else
                    slip_distribution=zeros((length(z_points(:,1))-1),(length(x_points(1,:))-1)); % creates a slip of 0 for faults that don't penetrate seismogenic layer
                end
                patch_plotting
%                 for i=1:length(utm_x) % drawing the down-dip points
%                     scatter3(x_points(:,i),y_points(:,i),z_points(:,i),'m');
%                     hold on
%                 end
                % 2D plotting of the gridded fault traces
                plot(utm_lon,utm_lat,'r')
                axis equal
                hold on
                %% Writing the data to the Coulomb output file
                for i=1:length(z_points(:,1))-1
                    for j=1:length(x_points(1,:))-1
                        if isempty(proj_dir)==1 %for faults which are vertical
                            fprintf (fid,'  1    %4.3f   %4.3f    %4.3f   %4.3f 100     %2.2f      %2.5f    %2.0f     %2.2f     %2.2f    %s\n', x_points(i,j)/1000,y_points(i,j)/1000,x_points(i,j+1)/1000,y_points(i,j+1)/1000,rake,slip_distribution(i,j),constant_dip,abs(z_points(i,j)/1000),abs(z_points(i+1,j)/1000),fault_name);
                        elseif proj_dir>=90 && proj_dir<=270 % for south dipping faults
                            if abs(z_points(i,j)/1000)<15
                                fprintf (fid,'  1    %4.3f   %4.3f    %4.3f   %4.3f 100     %2.2f      %2.5f    %2.0f     %2.2f     %2.2f    %s\n', x_points(i,j)/1000,y_points(i,j)/1000,x_points(i,j+1)/1000,y_points(i,j+1)/1000,rake,slip_distribution(i,j),constant_dip,abs(z_points(i,j)/1000),abs(z_points(i+1,j)/1000),fault_name);
                            else
                                fprintf (fid,'  1    %4.3f   %4.3f    %4.3f   %4.3f 100     %2.2f      %2.5f    %2.0f     %2.2f     %2.2f    %s\n', x_points(i,j)/1000,y_points(i,j)/1000,x_points(i,j+1)/1000,y_points(i,j+1)/1000,rake,slip_distribution(i,j),shear_zone_dip,abs(z_points(i,j)/1000),abs(z_points(i+1,j)/1000),fault_name);
                            end
                       else % for north dipping faults
                           if abs(z_points(i,j)/1000)<15
                               fprintf (fid,'  1    %4.3f   %4.3f    %4.3f   %4.3f 100     %2.2f      %2.5f    %2.0f     %2.2f     %2.2f    %s\n', x_points(i,j+1)/1000,y_points(i,j+1)/1000,x_points(i,j)/1000,y_points(i,j)/1000,rake,slip_distribution(i,j),constant_dip,abs(z_points(i,j)/1000),abs(z_points(i+1,j)/1000),fault_name);
                           else
                               fprintf (fid,'  1    %4.3f   %4.3f    %4.3f   %4.3f 100     %2.2f      %2.5f    %2.0f     %2.2f     %2.2f    %s\n', x_points(i,j+1)/1000,y_points(i,j+1)/1000,x_points(i,j)/1000,y_points(i,j)/1000,rake,slip_distribution(i,j),shear_zone_dip,abs(z_points(i,j)/1000),abs(z_points(i+1,j)/1000),fault_name);
                           end
                       end
                    end
                end
                clearvars -except grid_size depth maximum_slip fault_names fault_slip_name fid output_data_file filename min_x max_x min_y max_y UTM_zone UTM_letter kml_folder COUL_GRID_SIZE constant_dip rake PROJECTION_DIRECTION_FILE SHORT_FAULT_LENGTHS_FILE DIPS_FILE shear_zone_dip depth_ductile % clears all data except variables required for each loop
            end
         end
      end
    end
end
end
%% Finishing off writing the Coulomb input file
fprintf (fid,'\n');
fprintf (fid,'\n');
fprintf (fid,'    Grid Parameters\n');
fprintf (fid,'  1  ----------------------------  Start-x =    %3.5f\n',min_x);
fprintf (fid,'  2  ----------------------------  Start-y =   %5.4f\n',min_y);
fprintf (fid,'  3  --------------------------   Finish-x =    %3.5f\n',max_x);
fprintf (fid,'  4  --------------------------   Finish-y =   %5.4f\n',max_y);
fprintf (fid,'  5  ------------------------  x-increment =      %2.4f\n',COUL_GRID_SIZE);
fprintf (fid,'  6  ------------------------  y-increment =      %2.4f\n',COUL_GRID_SIZE);
fprintf (fid,'     Size Parameters\n');
fprintf (fid,'  1  --------------------------  Plot size =     3.000000\n');
fprintf (fid,'  2  --------------  Shade/Color increment =     1.000000\n');
fprintf (fid,'  3  ------  Exaggeration for disp.& dist. =     150000.0\n');
fprintf (fid,'\n');
fprintf (fid,'Cross section default\n');
fprintf (fid,'  1  ----------------------------  Start-x =   %3.5f\n',min_x);
fprintf (fid,'  2  ----------------------------  Start-y =   %4.4f\n',min_y);
fprintf (fid,'  3  --------------------------   Finish-x =   %3.5f\n',max_x);
fprintf (fid,'  4  --------------------------   Finish-y =   %4.4f\n',max_y);
fprintf (fid,'  5  ------------------  Distant-increment =      %2.4f\n',COUL_GRID_SIZE);
fprintf (fid,'  6  ----------------------------  Z-depth =     20.00000\n');
fprintf (fid,'  7  ------------------------  Z-increment =      %2.4f\n',COUL_GRID_SIZE);
fclose(fid);
fclose('all');