function writeOpenSimControlFile(OutputData)
% 
% WriteOpenSimControlFile: Writes muscle controls to an OpenSim STO file
%
% Input:
%   OutputData: a Matlab structure containing simulation results
%
%   The expected stucture fields are:
%       name: A char array identifier of the data
%       nRows: the number of rows of data in the data field
%       nColumns: the number of columns of data in the data field
%       labels: an array of char arrays of data names from the header file
%       data: a nRows by nColumnss matrix of data values
%
% Output:
%  An OpenSim controls file (.sto)
%
% Authors: Brian Umberger & Leng-feng Lee, UMass Amherst
%
% This function was modeled on ReadOpenSimData.m by Daniel Jacobs (Stanford)
% 

% Open file with write permission
fid = fopen([OutputData.name,'.sto'],'w');
    

% Write the name to the first row of the STO file
fprintf(fid,'%s\n',OutputData.name);

% Write the version number to the second row
fprintf(fid,'version=1\n');

% Write the number of data rows to the third row
fprintf(fid,'%s\n',['nRows=', num2str(OutputData.nRows)]);

% Write the number of data coulumns to the fourth row
fprintf(fid,'%s\n',['nColumns=', num2str(OutputData.nColumns)]);

% Write the degress/radians flag to the fifth row
if OutputData.inDegrees == 0
   fprintf(fid,'inDegrees=no\n');
elseif OutputData.inDegrees == 1
   fprintf(fid,'inDegrees=yes\n');
else
   disp(['Warning: "inDegrees" flag not specified correctly'])
end

% Write the "end-of-header" flag to the sixth row
fprintf(fid,'endheader\n');

% Write the data column labels in the seventh row
for i = 1:OutputData.nColumns
   if i == OutputData.nColumns
      fprintf(fid,'%s\n',OutputData.labels{1,i});
   else
      fprintf(fid,'%s\t',OutputData.labels{1,i});
   end
end

% Write the data to the next 'nRows' rows and 'nColumns' columns
for i = 1:OutputData.nRows
   for j = 1:OutputData.nColumns
      if j == OutputData.nColumns
         fprintf(fid,'%16.8f\n',OutputData.data(i,j));
      else
         fprintf(fid,'%16.8f\t',OutputData.data(i,j));
      end
   end
end

% Close the file
fclose(fid);


end

