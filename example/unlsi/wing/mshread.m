function [ msh_file_contents ] = mshread( filelocation, disp_loading_progress )
%MSHREAD Loads Fluent msh file.
%
% USAGE: (Loading files)
% [ msh_file_contents ] = mshread( filelocation )
% [ msh_file_contents ] = mshread( '' ); % asks for *.MSH file.
% OPTIONAL - show loading progress
% [ msh_file_contents ] = mshread( filelocation , 1); % load file // show progress
% [ msh_file_contents ] = mshread( '' , 1);           
%
% % VISUALISATION EXAMPLE
% [ msh_file_contents ] = mshread( '' , 1); % load file // show progress
%
% f1 = figure(1); set(f1,'color','w'); axis off; hold on
% colors = {'r','g','b','m','y','k','r','g','b','m','y','k','r','g','b','m','y','k','r','g','b','m','y','k','r','g','b','m','y','k','r','g','b','m','y','k','r','g','b','m','y','k'};
% fnames = fieldnames(msh_file_contents);
% fnames = fnames(~strcmpi(fnames,'vertices')'); lg=1;
% for i = 1:numel(fnames)
%     if strcmp(msh_file_contents.(fnames{i}).facetype,'wall')
%        patch('Faces',msh_file_contents.(fnames{i}).faces,'Vertices',msh_file_contents.vertices,'FaceColor',colors{i},'edgecolor','k'); hold on
%        legend_names{lg} = fnames{i}; lg = lg+1;
%     end
% end
% legend(legend_names{:})
% view(3); axis equal vis3d tight
% 
% --- Author information
% Wouter Potters
% Academic Medical Center, Amsterdam, The Netherlands
% w.v.potters@amc.nl
% Date: 14-November-2013
%% CHECK INPUTS
narginchk(0,2); nargoutchk(1,1);
switch nargin
    case 0, validfile = 0;
    case 1, if exist(filelocation,'file') == 2, validfile = 1; else validfile = 0; disp_loading_progress = 0; end
    case 2, if exist(filelocation,'file') == 2, validfile = 1; else validfile = 0; try if disp_loading_progress ~= 1,disp_loading_progress = 0; end, catch, disp_loading_progress = 0; end; end
end
if ~validfile % if no (valid) filepath is provided; let the user select one
    [filename,pathtofile] = uigetfile({'*.msh;*.MSH','msh files'},'Select *.msh file to load...');
    filelocation = fullfile(pathtofile,filename);
end
%% OPEN FILE
t1 = tic;
fid = fopen(filelocation,'r');
% ignore some crap
cur_line = fgetl(fid);
while (length(cur_line) < 6) || ~strcmp(cur_line(end-1:end),')(')
    cur_line = fgetl(fid);
end
%% read the vertices (elements)  - only one set of vertices supported
results         = sscanf(cur_line,'(%g (%g %g %x %g %g)(')';
% vertices_number = results(2);               % unique vertices_id
vertices_count  = results(4)-results(3)+1;  % # of vertices
vertices        = fscanf(fid,'%f',[3 vertices_count])'; % read vertices
if disp_loading_progress == 1, disp(['  vertices read:  ' num2str(length(vertices)) ', time: ' num2str(toc(t1)) ' seconds']), end
cur_line = fgetl(fid); % reset cur_line
%% read the faces - multiple faces enabled
all_faces = {}; %#ok<NASGU>
while ~feof(fid)
    % ignore some crap in the file
    niter=1;
    while ((length(cur_line) < 6) || ~strcmp(cur_line(end-1:end),')(')) && ~feof(fid)
        cur_line = fgetl(fid);
        if feof(fid) || any(strfind(cur_line,'Cells')), break; end % abort if cells or eof is found
        niter = niter+1; if niter>10, error('abort'); end
    end
    if feof(fid) || any(strfind(cur_line,'Cells')), break; end % abort if cells or eof is found
    
    % read the faces (elements)
    results = sscanf(cur_line,'(%g (%g %x %x %g %g)(')';
    faces_number = results(2);              % unique faces id
    faces_count  = results(4)-results(3)+1; % # of faces
    faces_type   = results(5);              % type of faces (3 connected points, 4 connected points)
    
    switch faces_type
        case 3 % case of 3 connected points
            faces = fscanf(fid,'%*f %x %x %x %*x %*x',[3 faces_count])';
            eval(['faces_' num2str(faces_number) ' = faces; all_faces = {[''faces_'' num2str(faces_number)] all_faces{:}};'])
            
        case 2 % % case of 4 connected points - volumetric elements - triangles
            faces = fscanf(fid,'%*f %x %x %x %x %*x',[4 faces_count])';
            eval(['faces_' num2str(faces_number) ' = faces; all_faces = {[''faces_'' num2str(faces_number)] all_faces{:}};'])
   
        otherwise % throw error if no supported face type was found
            error(['unsupported face type: ' num2str(faces_type)]);
            
    end
    if disp_loading_progress == 1, disp(['     faces read:  ' num2str(length(faces)) ', time: ' num2str(toc(t1)) ' seconds']), end
    cur_line = fgetl(fid);
end
% CEll part
cur_line = fgetl(fid); 
while ~isempty(cur_line) && ~feof(fid)
    result = sscanf(cur_line,'(%g (%g %x %x %g %g))');
    if disp_loading_progress == 1, disp(['number of cells: ' num2str(result(4))]), end
    cur_line = fgetl(fid); 
end
%% ZONES definition
cur_line = fgetl(fid); 
region_id = 1;
while ~isempty(cur_line) && ~feof(fid)
    cur_line = fgetl(fid);
    result = strsplit(cur_line,{')','(',' '}); %#ok<NASGU>
    eval('ZONE(region_id).number = num2str(result{3});')
    eval('ZONE(region_id).name   = result{5};')
    eval('ZONE(region_id).type   = result{4};')
    region_id = region_id + 1;
end
%% close file
if feof(fid)
    fclose(fid);
end
%% gather vertices and faces to msh_file_contents output variable
msh_file_contents.vertices = vertices;
for i = 1:length(ZONE)
    if ~strcmp(ZONE(i).type,'fluid') % somehow fluid part is not present in the msh file - is called default-interior
        msh_file_contents.(strrep(ZONE(i).name,'-','_')).faces    = eval(['faces_' num2str(ZONE(i).number)]);
        msh_file_contents.(strrep(ZONE(i).name,'-','_')).facetype = ZONE(i).type;
        msh_file_contents.(strrep(ZONE(i).name,'-','_')).number   = ZONE(i).number;
    end
end