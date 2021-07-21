function [Vc,Fc] = read_vertices_and_faces_from_WRL_file(filename,n)
  % Reads a .wrl file e.g. for 67P and outputs the vertex and face list
  
  % Right now, requires user to edit the line that has the number of
  % vertices and number of facets
  % Take out the spaces so it looks like: # xxxx xxxx
  % I'll fix this later
  
  % Input:
  %  filename:  string of obj file's path
  %  n: no. digits of facets (3k = 4, 12k = 5 e.g.)
  %
  % Output:
  %  V  -- number of vertices x 3 array of vertex positions
  %  F  -- number of faces x 3 array of face indices
  
  Vc = zeros(0,3);
  Fc = zeros(0,3);
  
  fid = fopen(filename,'r');
  
  % Ignore the first 5 lines of the file
  for i=(1:5)
    line_2_discard = fgetl(fid);
  end
  
  % line containing number of vertices and no facets
  line_2_start_with = fgetl(fid);
  
  % n gives no. chars in physical number of verticies
  % Reads n verticies as 1:n of line_2_start_with, conv to float
  % Then takes a space and reads the rest of the line as n facets
  % embarrasing programming I will come back and fix this
  n_vertices = str2double(line_2_start_with(2:n+2))
  n_facets = str2double(line_2_start_with(n+3:length(line_2_start_with)))
  
  % read the appropriate bits of the file for each of V and F
  for i=(1:n_vertices)
    line = fgetl(fid);
    vertex = sscanf(line, '%f %f %f');
    Vc(i,:) = vertex;
  end
  
  %Now skip the next 3 lines and read in the coords of the faces in the same
  %way
  for i=(1:3)
      line_2_discard =fgetl(fid);
  end
  
  for i=(1:n_facets)
    line = fgetl(fid);
    facet = sscanf(line, '%f %f %f -1');
    Fc(i,:) = facet;
  end
  
end