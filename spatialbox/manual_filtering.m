function manual_filtering
% manual filtering

% read the file list:
[gui_files,gui_path,handles.dt,handles.scale,handles.state3d] = cil_uigetfiles;

handles.N = length(gui_files); % number of files selected
if  handles.N > 0
    handles.files = gui_files;
    handles.path = gui_path;
end

currentdir = pwd;

cd(handles.path);
if ~isempty(findstr(lower(handles.files{1}),'vec'))            % process .vec files
    % read the first file, determine the size
    for i = 1:handles.N
        [handles.xUnits,handles.velUnits,d] = vecread(fullfile(handles.path,handles.files{i}));
        [rows,cols,k] = size(d);
        %    [handles.u,handles.v] = deal(zeros(rows,cols,handles.N+1)); % 11.04.04, Alex
        handles.x           = d(:,:,1);
        handles.y           = d(:,:,2);
        handles.u           = d(:,:,3);
        handles.v           = d(:,:,4);

        % THIS IS THE FILTER CHECK
        renameifbadvecfile(handles,i)
    end
end
cd(currentdir);
end


function renameifbadvecfile(handles,i)
THRESHOLD = 1E-5;
% manual filter function - redefine
if sum(handles.u.^2 + handles.v.^2) < THRESHOLD
    oldfile = handles.files{i};
    newfile = [oldfile(1:end-4),'.badvec'];
    dos(sprintf('ren %s %s &',oldfile,newfile));
end
end

