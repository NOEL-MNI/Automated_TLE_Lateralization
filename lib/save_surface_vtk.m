function save_surface_vtk(surf,fname,fileType,valuesPerVertex)
% function save_surface_vtk(surf,fname,fileType,valuesPerVertex)
%
% surf            : A structure loaded by surfstat's reader
% fname           : filename.vtk
% fileType        : String, must be either 'ASCII' or 'BINARY'.
% valuesPerVertex : Data in a structure. Each field contains a 
%                   vector of size [nVertices 1]. You can have as many
%                   fields as you want. The name of each field gets
%                   reflected in the newly created vtk file.
%                   For example:
%                   valuesPerVertex.FA = [nVertices 1]
%                   valuesPerVertex.E1 = [nVertices 1]
%
% Luis Concha. BIC. July 2008.
%
%   See also SAVE_VOLUME_VTK, SAVE_TRACT_VTK



if nargin < 4
  writeFieldData = false;
else
  writeFieldData = true; 
end

if nargin < 3
   fileType = 'ASCII'; 
end

if ~strcmp(fileType,'BINARY') & ~strcmp(fileType,'ASCII')
   error('Invalid file type (ASCII or BINARY only)');
   return
end


%% Check that surf.tri has the coordinates the way we need them. Different
%% versions of surfstat transpose these.
if size(surf.tri,1) > size(surf.tri,2)
   surf.tri = surf.tri'; 
end


fid = fopen(fname,'w');

fprintf(fid,'%s\n','# vtk DataFile Version 3.0');
fprintf(fid,'%s\n','surface');
fprintf(fid,'%s\n',fileType);
fprintf(fid,'%s\n','DATASET POLYDATA');

nPoints = length(surf.coord);

fprintf(fid,'%s %d %s\n','POINTS',nPoints,'float');

if strcmp(fileType,'BINARY')
   fwrite(fid,surf.coord,'float','ieee-be');

else
   fprintf(fid,'%f %f %f\n',surf.coord);
end

nTriangles = length(surf.tri);
fprintf(fid,'%s %d %d\n','POLYGONS',nTriangles,4*size(surf.tri,2));

if strcmp(fileType,'BINARY')
   triangles = [zeros(1,length(surf.tri))+3;surf.tri-1];
   fwrite(fid,triangles,'int','ieee-be');
else
   triangles = [zeros(1,length(surf.tri))+3;surf.tri-1];
   fprintf(fid,'%d %d %d %d\n',triangles); 
end

% Writing normals not implemented yet
% fprintf(fid,'%s %s %s\n','NORMALS','point_normals','float');
% if strcmp(fileType,'BINARY')
%    fwrite(fid,surf.normal,'int','ieee-be');
% else
%    fprintf(fid,'%f %f %f\n',surf.normal); 
% end


if writeFieldData
    vars = fieldnames(valuesPerVertex);
    fprintf(fid,'%s %d\n','POINT_DATA',nPoints);
    for v = 1 : length(vars)
        var  = vars{v};
        eval(['data = valuesPerVertex.' var '(:,1);']);
        fprintf(fid,'%s %s %s\n','SCALARS',var,'float');
        fprintf(fid,'%s\n','LOOKUP_TABLE default');
        if strcmp(fileType,'BINARY')
           fwrite(fid,data,'float','ieee-be');
        else
           fprintf(fid,'%f\n',data);
           fprintf(fid,'\n');
        end
    end
end



fclose(fid);