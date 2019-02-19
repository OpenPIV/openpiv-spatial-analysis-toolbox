[gui_files,gui_path] = getVECfiles;

handles.N = length(gui_files);
if  handles.N > 1
    handles.files = gui_files;
    handles.path = gui_path;
%     set(handles.fig,'pointer','watch');
else
    warndlg('More than 2 files are required','Error','modal');
    set(handles.fig,'pointer','arrow');
    return
end
if ~isempty(findstr(lower(handles.files{1}),'vec')) & handles.N > 1
    [handles.xUnits,handles.velUnits,d] = vecread(fullfile(handles.path,handles.files{1}));
    [rows,cols,k] = size(d);
    [handles.u,handles.v] = deal(zeros(rows,cols,handles.N));
    handles.x           = d(:,:,1);
    handles.y           = d(:,:,2);
    handles.u(:,:,1)    = d(:,:,3);
    handles.v(:,:,1)    = d(:,:,4);
    for i = 2:handles.N
        d = vecread([handles.path,filesep,handles.files{i}],1,5);
        handles.u(:,:,i) = d(:,:,3);
        handles.v(:,:,i) = d(:,:,4);
    end
    clear d
end

quiver(handles.x,handles.y,handles.u(:,:,1),handles.v(:,:,1))