% Script to plot data from SPH slosh simulation

function [] = save_vtu(particles,n_save, dirname)

% specify n_save = 0 for boundary particles

%% saving directory
% check for existence of paraviewfiles/vtu directory. this is the directory where
% the .vtu files will be stored. if it does not exist create it
% dirname='VTU_Results';
dirstatus=exist(dirname,'dir');
if(dirstatus==0)
    mkdir(dirname)
end

	
%% test vtu output (ascii)
i = n_save;
%for i = 1:size(save_pos,3)
	% specify file name   
    if n_save ~= 0
        name = strcat('slosh',num2str(i,'%3.3d'),'.dat');
        fid=eval(['fopen(''',dirname,'/PART' num2str(i,'%3.3d') '.vtu' ''',''w'' );']);
    else % for boundary particles
        %name = strcat('slosh',num2str(i,'%3.3d'),'.dat');
        fid=eval(['fopen(''',dirname,'/BOUND' num2str(i,'%3.3d') '.vtu' ''',''w'' );']);
        i = 1;
    end
        
    % specify data to store/ output
    % particles(:,:) = save_pos(:,:,i);
    np = size(particles,1);
    xp=particles(:,1); % position
    zp=particles(:,2);
    up=particles(:,6); % velocity
    wp=particles(:,7);
    rhop=particles(:,5); % density
    P=particles(:,8); % pressure

    % numbering
    idiv=2;
    nbeg1=1;
    nb=0;
    nend1=nb;
    nbeg2=nb+1;
    nend2=np;

    %   output to file in vtu format
    fprintf(fid,'<?xml version="1.0"?>\r\n');
    fprintf(fid,'<VTKFile type= "UnstructuredGrid"  version= "0.1"  byte_order= "BigEndian">\r\n');
    fprintf(fid,' <UnstructuredGrid>\r\n');
    fprintf(fid,'  <Piece NumberOfPoints="%d" NumberOfCells="%d">\r\n',np,np);

    % write in pressure data
    fprintf(fid,'   <PointData Scalars="Pressure" Vectors="Velocity">\r\n');
    fprintf(fid,'    <DataArray type="Float32" Name="Pressures" format="ascii">\r\n');
    for ii=1:np
        fprintf(fid,'%f\t',P(ii));
        fprintf(fid,'\n');
    end
    fprintf(fid,'\r\n');
    fprintf(fid,'    </DataArray>\r\n');

    % write density data
    fprintf(fid,'    <DataArray type="Float32" Name="Density" format="ascii">\r\n');
    for ii=1:np
        fprintf(fid,'%f\t',rhop(ii));
        fprintf(fid,'\n');
    end
    fprintf(fid,'\r\n');
    fprintf(fid,'    </DataArray>\r\n');

    % this section is used to color different particles based the input idiv specified above.
    fprintf(fid,'    <DataArray type="Float32" Name="Scalarplot" format="ascii">\r\n');
    for ii=1:idiv
        eval(['nbeg=nbeg' int2str(ii) ';'])
        eval(['nend=nend' int2str(ii) ';'])
        for jj=nbeg:nend
            fprintf(fid,'%f\t',ii);
            fprintf(fid,'\n');
        end
    end
    fprintf(fid,'\r\n');
    fprintf(fid,'    </DataArray>\r\n');

    % write velocity data
    fprintf(fid,'    <DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="ascii">\r\n');
    for ii=1:np
        vel=[up(ii) 0 wp(ii)];
        fprintf(fid,'%f\t %f\t %f\t',vel);
        fprintf(fid,'\n');
    end
    fprintf(fid,'\r\n');
    fprintf(fid,'    </DataArray>\r\n');
    fprintf(fid,'   </PointData>\r\n');

    % write particle position data
    fprintf(fid,'   <Points>\r\n');
    fprintf(fid,'    <DataArray type="Float32" NumberOfComponents="3" format="ascii">\r\n');
    for ii=1:np
        pos=[xp(ii) 0 zp(ii)];
        fprintf(fid,'%f\t %f\t %f\t',pos);
        fprintf(fid,'\n');
    end
    fprintf(fid,'\r\n');
    fprintf(fid,'    </DataArray>\r\n');
    fprintf(fid,'   </Points>\r\n');

    % write cell data. cell is of type vertex.
    fprintf(fid,'   <Cells>\r\n');
    fprintf(fid,'    <DataArray type="Int32" Name="connectivity" format="ascii">\r\n');
    for ii=1:np
        fprintf(fid,'%d\t',ii-1);
        fprintf(fid,'\n');
    end
    fprintf(fid,'\r\n');
    fprintf(fid,'    </DataArray>\r\n');
    fprintf(fid,'\r\n');
    fprintf(fid,'    <DataArray type="Int32" Name="offsets" format="ascii">\r\n');
    for ii=1:np
        fprintf(fid,'%d\t',ii);
        fprintf(fid,'\n');
    end
    fprintf(fid,'\r\n');
    fprintf(fid,'    </DataArray>\r\n');
    fprintf(fid,'\r\n');
    fprintf(fid,'    <DataArray type="Int32" Name="types" format="ascii">\r\n');
    for ii=1:np
        fprintf(fid,'%d\t',1);
        fprintf(fid,'\n');
    end
    fprintf(fid,'\r\n');
    fprintf(fid,'    </DataArray>\r\n');
    fprintf(fid,'   </Cells>\r\n');
    fprintf(fid,'  </Piece>\r\n');
    fprintf(fid,' </UnstructuredGrid>\r\n');
    fprintf(fid,'</VTKFile>');
    fclose(fid);

%end

	