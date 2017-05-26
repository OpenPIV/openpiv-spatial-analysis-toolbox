function varargout = distrib(varargin)
% DISTRIB Application M-file for distrib.fig
%    FIG = DISTRIB launch distrib GUI.
%    DISTRIB('callback_name', ...) invoke the named callback.

% Copyright (c) 1998-2014 OpenPIV group
% See the file license.txt for copying permission.


if nargin == 1  % LAUNCH GUI
    
    fig = openfig(mfilename,'reuse');
    movegui(fig,'center');
    % set(gcf,'toolbar','figure');
    
    % Print option default : don't print uicontrols
    pt = printtemplate;
    pt.PrintUI = 0;
    set(fig, 'PrintTemplate', pt);
    
    % ------   Initialization --------------
    % Generate a structure of handles to pass to callbacks, and store it.
    handles                 = guihandles(fig);
    % in GUIDE the figure Tag is defined as 'fig', so
    % handles.fig is enough.
    %     handles.distribHandle   = fig;
    handles.data            = varargin{1};
    
    % --------------  Process only selected region---
    MatrixLeftMove  = min(handles.data.i)-1;
    MatrixUpMove    = min(handles.data.j)-1;
    
    % ------------- Cut selection and assign to handles.property ---------
    
    if isfield(handles,'property')
        handles  =  rmfield(handles,'property');
    end
    for i  =  1:length(handles.data.i)
        handles.property(handles.data.i(i) - MatrixLeftMove, handles.data.j(i) - MatrixUpMove) = ...
            handles.data.property(handles.data.i(i),handles.data.j(i));
    end
    
    
    if handles.data.by_row == true && handles.data.by_column == false
        % should be empty rows
        index = all(handles.property == 0,2);
        handles.property(index,:) = NaN;
        set(handles.UpdateGraph,'Value',1); % by rows
    elseif handles.data.by_row == false && handles.data.by_column == true
        %should be empty columns, replace by NaN
        index = all(handles.property == 0,1);
        handles.property(:,index) = NaN;
        set(handles.UpdateGraph,'Value',2); % by cols
    end
        
    
    %------------------------------- Delete empty rows ------------------
%     index = 1;
%     emptyrows               = [];
%     emptycols               = [];
    
%     for i = 1:size(handles.property,1)
%         if all(handles.property(i,:) == 0)
%             emptyrows(index) = i;  %#ok<AGROW>
%             index = index + 1;
%         end;
%     end;
%     if ~isempty(emptyrows) % there were empty rows, there should not be also empty columns
%         handles.property(emptyrows,:) = NaN;
%     else
%         
%         %--------------------Delete  Empty Columns ----------
%         index = 1;
%         for i = 1:size(handles.property,2)
%             if all(handles.property(:,i) == 0)
%                 emptycols(index) = i;  %#ok<AGROW>
%                 index = index+1;
%             end
%         end
%         if ~isempty(emptycols)
%             handles.property(:,emptycols) = NaN;
%             %           set(handles.UpdateGraph,'Value',2);
%             %            set(handles.UpdateGraph,'Enable','off');
%         end
%     end
    % --------------------------------
    
    handles.leftX = min(handles.data.j);
    handles.rightX = max(handles.data.j) ; %;(length(handles.data.j));
    handles.bottomY = size(handles.data.y,1) - min(handles.data.i) + 1;
    handles.topY = size(handles.data.y,1) - max(handles.data.i) + 1; %handles.data.i(length(handles.data.i));
    
    handles.swap = 0;   % default - no swap;
    handles.legend = 0; % default - no legend;
    
    handles.hold = 0;
    
    % handles.displayname = ''; %cellstr(num2str(handles.data.y(handles.topY:handles.bottomY,1)));
    handles.ydir = 'normal';
    
%     if length(unique(handles.data.j)) <= length(unique(handles.data.i))
%         set(handles.UpdateGraph,'Value',2); % by cols
%     else
%         set(handles.UpdateGraph,'Value',1); % by rows
%     end
    UpdateGraph_Callback(fig, [], handles);
    
    guidata(fig, handles);
        
%         plot(handles.data.x(1,handles.leftX:handles.rightX),handles.property,...
%             'DisplayName',handles.displayname);
%         %     legend show
%         ylabel(['x',handles.data.xUnits]); % [m]');
%         grid on;
%         xlabel(strcat(handles.data.previous_quantity,'  ',handles.data.units));

    
    %     plotedit on;
    if nargout > 0
        varargout{1}  =  fig;
    end
    
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
    
    try
        if (nargout)
            [varargout{1:nargout}]  =  feval(varargin{:}); % FEVAL switchyard
        else
            feval(varargin{:}); % FEVAL switchyard
        end
    catch
        disp(lasterr);
    end
    
end


% --------------------------------------------------------------------
function varargout  =  UpdateGraph_Callback(h, eventdata, handles, varargin)

val  =  get(handles.UpdateGraph,'Value');

switch val
    case 1
        x = handles.data.x(1,handles.leftX:handles.rightX);
        x_label = (['x ',handles.data.xUnits]); % 'x [m]';
        % y = handles.property';
        y = handles.property(~isnan(handles.property(:,1)),:)';
        % y = y(~isnan(y(:,1)),:);
        y_label = strcat(handles.data.previous_quantity,'  ',handles.data.units);
        xa = x;
        % ya = no_nan_mean(handles.property);
        ya = mean(y,2);
        % handles.displayname = cellstr(num2str(handles.data.y(handles.topY:handles.bottomY,1)));
        handles.displayname = cellstr( num2str( handles.data.y(unique(handles.data.i),1) ) );
    case 2
        y = handles.data.y(handles.topY:handles.bottomY,1);
        y_label = (['y ',handles.data.xUnits]); % 'y [m]';
        % x = handles.property;
        x = handles.property(:,~isnan(handles.property(1,:)));
        % x = x(:,~isnan(x(1,:)));
        x_label = strcat(handles.data.previous_quantity,'  ',handles.data.units);
        % xa = no_nan_mean(handles.property')';
        xa = mean(x,2);
        ya = y;
        % handles.displayname = cellstr(num2str(handles.data.x(1,handles.leftX:handles.rightX)'));
        
        handles.displayname = cellstr(num2str(handles.data.x(1,unique(handles.data.j))'));
end

if handles.swap
    tmp = x;
    x = y;
    y = tmp;
    tmp = x_label;
    x_label = y_label;
    y_label = tmp;
    tmp = xa;
    xa = ya;
    ya = tmp;
end



cla reset;
% ---------------------- Single is checked------------------
if get(handles.SnglCheckbox,'Value') == 1
    hold on;
    switch val
        case 1 % by rows
            for jj = 1:size(y,2)
                plot(x,y(:,jj),'DisplayName',handles.displayname{jj});
            end
        case 2 % by cols
            for jj = 1:size(x,2)
                plot(x(:,jj),y,'DisplayName',handles.displayname{jj});
            end
    end
    xlabel(x_label);
    ylabel(y_label);
    hold off;
end
% --------------------- Avge is checked -------------------------
if get(handles.AvgCheckbox,'Value') == get(handles.AvgCheckbox,'Max')
    hold on;
    plot(xa,ya,'r-','LineWidth',2);
    hold off;
    xlabel(x_label);
    ylabel(y_label);
end;
grid on;
if handles.legend
    legend('show');
else
    legend('hide');
end
% plotedit on;
set(gca,'ydir',handles.ydir);




% --------------------------------------------------------------------
function SwapXY_Callback(h, eventdata, handles, varargin)
if (get(h,'Value')  ==  get(h,'Max'))
    handles.swap = 1;
    
else
    % checkbox is not checked-take approriate action
    handles.swap = 0;
    
end
guidata(gcbo,handles);
UpdateGraph_Callback(gcbo,[],guidata(gcbo));


% --------------------------------------------------------------------
% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
delete(handles.fig);

% --------------------------------------------------------------------
% --- Executes on button press in AvgCheckbox.
function AvgCheckbox_Callback(hObject, eventdata, handles)
UpdateGraph_Callback(gcbo,[],guidata(gcbo));



% --- Executes on button press in SnglCheckbox.
function SnglCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to SnglCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SnglCheckbox
UpdateGraph_Callback(gcbo,[],guidata(gcbo));


% --------------------------------------------------------------------
% --- Executes on button press in savecsv.
function savecsv_Callback(hObject, eventdata, handles)
% hObject    handle to savecsv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hl  =  findobj(gca,'type','line');
if ~isempty(hl)
    xd  =  get(hl,'xdata');
    yd  =  get(hl,'ydata');
    if ~iscell(xd)
        matr(:,1)  =  xd(:);
        matr(:,2)  =  yd(:);
    elseif isequal(xd{1},xd{2})
        matr(:,1)  =  [xd{1}]';
        for i  =  1:length(hl),
            matr(:,i+1)  =  [yd{i}]';
        end
    elseif isequal(yd{1},yd{2})
        matr(:,1)  =  [yd{1}]';
        for i  =  1:length(hl),
            matr(:,i+1)  =  [xd{i}]';
        end
    else
        for i  =  1:2:length(hl),
            matr(:,i)  =  [xd{i}]';
            matr(:,i+1)  =  [yd{i}]';
        end
    end
    
    file  =  [];
    file  =  inputdlg('File Name','Input Name for CSV File');
    if ~isempty (file)
        csvwrite(file{1},matr);
    else
        errordlg ('Choose a valid file name !!! ');
    end;
end;


% --------------------------------------------------------------------
function export_Callback(hObject, eventdata, handles)
% hObject    handle to export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function export2figure_Callback(hObject, eventdata, handles)
% hObject    handle to export2figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.export_figure = figure;
% handles.export_axes   = axes;
handles.axes1 = get(handles.fig,'CurrentAxes');
copyobj(handles.axes1,handles.export_figure);
% copyobj(get(handles.axes_main,'Children'),handles.export_axes);

set(handles.export_figure,'Units','normalized');
handles.export_axes = get(handles.export_figure,'CurrentAxes');
set(handles.export_axes,'Units','normalized');
set(handles.export_axes,'Position',[0.13 0.11 0.775 0.815]);
set(handles.export_axes,'Box','on');
hl =  findobj(handles.export_axes,'type','line');
xd = get(hl,'xdata');
yd = get(hl,'ydata');
if iscell(xd)
    for i = 1:length(xd),
        if all(isnan(xd{i})) | all(isnan(yd{i})),
            delete(hl(i));
        end
    end
end


guidata(handles.fig, handles);


% --------------------------------------------------------------------
function export2CSV_Callback(hObject, eventdata, handles)
% hObject    handle to export2CSV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hl  =  findobj(gca,'type','line');
if ~isempty(hl)
    xd  =  get(hl,'xdata');
    yd  =  get(hl,'ydata');
    if ~iscell(xd)
        matr(:,1)  =  xd(:);
        matr(:,2)  =  yd(:);
    elseif isequal(xd{1},xd{2})
        matr(:,1)  =  [xd{1}]';
        for i  =  1:length(hl),
            matr(:,i+1)  =  [yd{i}]';
        end
    elseif isequal(yd{1},yd{2})
        matr(:,1)  =  [yd{1}]';
        for i  =  1:length(hl),
            matr(:,i+1)  =  [xd{i}]';
        end
    else
        for i  =  1:2:length(hl),
            matr(:,i)  =  [xd{i}]';
            matr(:,i+1)  =  [yd{i}]';
        end
    end
    
    file  =  [];
    file  =  inputdlg('File Name','Input Name for CSV File');
    if ~isempty (file)
        csvwrite(file{1},matr);
    else
        errordlg ('Choose a valid file name !!! ');
    end;
end;



function b = no_nan_mean(a)
anan = isnan(a);
indx = find(anan);
a(indx) = zeros(size(indx));
notnans = sum(~anan);
indx = find(notnans == 0);
notnans(indx) = 1;
b = sum(a)./notnans;
b(indx) = nan;


% --- Executes on button press in checkbox_legend.
function checkbox_legend_Callback(h, eventdata, handles, varargin)
if (get(h,'Value')  ==  get(h,'Max'))
    handles.legend = 1;
    
else
    % checkbox is not checked-take approriate action
    handles.legend = 0;
    
end
guidata(gcbo,handles);
UpdateGraph_Callback(gcbo,[],guidata(gcbo));




% --- Executes on button press in checkbox_legend.
function reverse_y_Callback(h, eventdata, handles, varargin)
if (get(h,'Value')  ==  get(h,'Max'))
    handles.ydir = 'reverse';
    
else
    % checkbox is not checked-take approriate action
    handles.ydir = 'normal';
    
end

guidata(gcbo,handles);
UpdateGraph_Callback(gcbo,[],guidata(gcbo));
