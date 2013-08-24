% SPATIALBOX - main GUI file, created by GUIDE
%

% Copyright (c) 1998-2012 OpenPIV group
% See the file license.txt for copying permission.

% Changes log:
% Last modified: 13-12-2004
% by Alex Liberzon,
% ydir is now reversed in update_gui for the Insight data. Check for the
% OpenPIV data
% - abs() is added to gridX and gridY - both have to be removed since
% it just doubles the handles.dx and handles.dy

function varargout = spatialbox(varargin)
warning off %#ok<WNOFF>

global orighandles current_path;

if nargin == 0   % LAUNCH GUI
    current_path = cd;
    addpath(current_path);
    
    fig = openfig(mfilename,'reuse','invisible');
    movegui(fig,'center')
    
    % Generate a structure of handles to pass to callbacks, and store it.
    handles = guihandles(fig);
    
    % Global several long strings
    handles.inst_list = '-|u|v|(u^2+v^2)^(1/2)|vorticity|sxx=du/dx|du/dy|dv/dx|syy=dv/dy|du/dx+dv/dy|sxy';
    handles.mean_list = '-|U|V|(U^2+V^2)^(1/2)|Vorticity|Sxx=dU/dx|dU/dy|dV/dx|Syy=dV/dy|dU/dx+dV/dy|Sxy';
    handles.fluct_list = '-|u''|v''|(u''^2+v''^2)^(1/2)|Vorticity''|sxx''=du''/dx|du''/dy|dv''/dx|syy''=dv''/dy|du''/dx+dv''/dy|sxy''';
    handles.fluct_mean_list = '-|u rms|v rms|Reynolds stress|Turb. intensity (u) |Turb. intensity (v)|Dissipation|Turb. Energy Production|TKE|Enstropy';
    
    handles.fig = fig;
    handles.previous_quantity = '-';
    orighandles = handles;          % backup for the later re-opening of the data
    guidata(handles.fig, handles);
    
    % Use system color scheme for figure:
    set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
    
    
    
    guidata(handles.fig, handles);
    
    if nargout > 0
        varargout{1} = fig;
    end
    set(findobj(handles.fig,'type','uicontrol'),'Enable','Off');
    
    % Show the figure
    set(handles.fig,'Visible','on');    % ensures that user gets only a prepared window.
    
elseif ischar(varargin{1}) % Invoke CallBack functions
    
    try
        if (nargout)
            [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
        else
            feval(varargin{:}); % FEVAL switchyard
        end
    catch ME
        % disp(lasterr);
        disp(ME.message); 
    end
    
end
% warning off

% --------------------------------------------------------------------
function checkbox_arrow_Callback(~, ~, handles, varargin) %#ok<DEFNU>
% arrow on / off callback
% if arrow - off is checked, hide the arrows.
if get(handles.checkbox_arrow,'Value') == 0
    handles.color=0;
    set(handles.checkbox_arrow_color,'Enable','off');
    set(handles.checkbox_arrow_color,'Value',0);
else
    set(handles.checkbox_arrow_color,'Enable','on');
    delete(get(handles.axes_main,'children'));
end
guidata(handles.fig,handles);
update_gui(handles.fig,[],handles);


% --------------------------------------------------------------------
function checkbox_arrow_color_Callback(hObject, ~, handles, varargin)
% Color arrows depends on the handles.property, chosen from the list of the
% avialable quantities

% color / black callback
if (get(hObject,'Value') == 1)
    handles.color = 1;
else
    % checkbox is not checked-take approriate action
    handles.color = 0;
    %     handles.color_quiver = get(get(handles.fig,'currentaxes'),'children');
    set(handles.color_quiver,'Visible','off');
end

guidata(handles.fig,handles);
update_gui(handles.fig,[],handles);


% --------------------------------------------------------------------
function popupmenu_quantity_Callback(~, ~, handles, varargin)
% Property Selection callback
% Assign handles.property to present as color or contour later, according
% to the flags of ensemble, fluctuatinons and possible quantities


% --------------Ensemble + Fluct --------------
if get(handles.checkbox_ensemble,'Value') == 1
    if get(handles.checkbox_fluct,'Value') == 1
        
        switch get(handles.popupmenu_quantity,'Value')
            case 1
                handles.property = [];
                handles.colorbar_flag = 0;
                
            case 2 % uf^2
                if isfield(handles,'uf2')
                    handles.property = handles.uf2;
                    
                else
                    handles.uf2 = sqrt(mean(handles.uf.^2,3)); %sqrt(u'^2) % 16.06.08. Alex
                    handles.property = handles.uf2;
                end
                %                                handles.units='[m/s]^2';
                handles.units = handles.velUnits; % 16.06.08 Alex
                % [handles.velUnits,'^2']; % '[m/s]^2';
            case 3 % vf^2
                
                if isfield(handles,'vf2')
                    handles.property = handles.vf2;
                else
                    handles.vf2 = sqrt(mean(handles.vf.^2,3)); %u'^2
                    handles.property = handles.vf2;
                end
                % Alex, 16.06.08
                handles.units = [handles.velUnits]; %'[m/s]^2';
            case 4
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units=['[',handles.velUnits,'^2/s^2]'];
                else
                    handles.units = ['[',handles.velUnits,'^2/\Delta t^2]'];
                end
                
                
                if isfield(handles,'rs')
                    handles.property = handles.rs;
                else
                    handles.rs = -1 * mean(handles.uf.*handles.vf,3); % Reynolds stress
                    handles.property = handles.rs;
                end
                % --------------------------
            case 5
                handles.units='[]';
                if isfield(handles,'Tu'),
                    handles.property = handles.Tu;
                else
                    meanU = abs(mean(mean(handles.u(:,:,handles.N+1))));
                    meanV = abs(mean(mean(handles.v(:,:,handles.N+1))));
                    if  meanU > meanV   % check which one is streamwise
                        handles.Tu = sqrt(mean(handles.uf.^2,3))/meanU;
                    else
                        handles.Tu = sqrt(mean(handles.uf.^2,3))/meanV;
                    end
                    handles.property = handles.Tu;
                end
                % ----------------------------
            case 6
                handles.units='[]';
                if isfield(handles,'Tv'),
                    handles.property = handles.Tv;
                else
                    meanU = abs(mean(mean(handles.u(:,:,handles.N+1))));
                    meanV = abs(mean(mean(handles.v(:,:,handles.N+1))));
                    if  meanU > meanV   % check which one is streamwise
                        handles.Tv = sqrt(mean(handles.vf.^2,3))/meanU;
                    else
                        handles.Tv = sqrt(mean(handles.vf.^2,3))/meanV;
                    end
                    handles.property = handles.Tv;
                end
            case 7 % dissipation, Alex, 01.04.04
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s^2]'  ;
                else
                    handles.units = '[1/\Delta t^2]';
                end
                
                if isfield(handles,'diss')
                    handles.property = handles.diss;
                else
                    handles.diss = zeros(size(handles.x));
                    for i = 1:handles.N
                        dudx = handles.dudx(:,:,i) - handles.dudx(:,:,handles.N+1);
                        dudy = handles.dudy(:,:,i) - handles.dudy(:,:,handles.N+1);
                        dvdx = handles.dvdx(:,:,i) - handles.dvdx(:,:,handles.N+1);
                        dvdy = handles.dvdy(:,:,i) - handles.dvdy(:,:,handles.N+1);
                        handles.diss = handles.diss + (dudx).^2 + (dvdy).^2 + 2*(dudy + dvdx).^2;
                    end
                    handles.diss = handles.diss/handles.N;
                    handles.property = handles.diss;
                end
            case 8 % energy production - <u'v'>S_xy
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[m^2/s^3]' ;
                else
                    handles.units = '[pixel^2/\Delta t^3]';
                end
                
                if isfield(handles,'prod')
                    handles.property = handles.prod;
                else
                    handles.prod = -1*( mean(handles.uf.*handles.vf,3).*( handles.dudy(:,:,handles.N+1) + handles.dvdx(:,:,handles.N+1) ) + ...
                        mean(handles.uf.*handles.uf,3).*handles.dudx(:,:,handles.N+1)  +...
                        mean(handles.vf.*handles.vf,3).*handles.dvdy(:,:,handles.N+1));
                    handles.property = handles.prod;
                end
            case 9 % TKE, 11.03.10 - Alex
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[m^2/s^2]' ;
                else
                    handles.units = '[pixel^2/\Delta t^2]';
                end
                
                if isfield(handles,'TKE')
                    handles.property = handles.TKE;
                else
                    handles.TKE = mean(handles.vf.^2 + handles.uf.^2,3);
                    handles.property = handles.TKE;
                end
            case 10 % Enstrophy, Alex, 11.03.10
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s^2]'  ;
                else
                    handles.units = '[1/\Delta t^2]';
                end
                
                if isfield(handles,'enst')
                    handles.property = handles.enst;
                else
                    handles.enst = zeros(size(handles.x));
                    for i = 1:handles.N
                        dudy = handles.dudy(:,:,i) - handles.dudy(:,:,handles.N+1);
                        dvdx = handles.dvdx(:,:,i) - handles.dvdx(:,:,handles.N+1);
                        handles.enst = handles.enst + (dvdx - dudy).^2;
                    end
                    handles.enst = handles.enst/handles.N;
                    handles.property = handles.enst;
                end
                
        end % of switch
        % ------------------------------------------- Ensemble ------
    else % of fluct, means only ensemble is chosen
        switch get(handles.popupmenu_quantity,'Value')
            case 1
                handles.property = [];    % default variable, choose properly later, alex, 16.02
                handles.colorbar_flag = 0;
                
            case 2
                handles.units = handles.velUnits;
                
                handles.property = handles.u(:,:,handles.current);
            case 3
                handles.units = handles.velUnits;
                handles.property = handles.v(:,:,handles.current);
            case 4
                handles.units = handles.velUnits;
                handles.property = sqrt(handles.u(:,:,handles.current).^2+handles.v(:,:,handles.current).^2);%Velocity Magnitude
            case 5 % vorticity;
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s]'  ;
                else
                    handles.units = '[1/\Delta t]';
                end
                
                handles.property = handles.dvdx(:,:,handles.current) - handles.dudy(:,:,handles.current);
            case 6
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s]'  ;
                else
                    handles.units = '[1/\Delta t]';
                end
                handles.property = handles.dudx(:,:,handles.current);
            case 7
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s]'  ;
                else
                    handles.units = '[1/\Delta t]';
                end
                handles.property = handles.dudy(:,:,handles.current);
            case 8
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s]'  ;
                else
                    handles.units = '[1/\Delta t]';
                end
                handles.property = handles.dvdx(:,:,handles.current);
            case 9
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s]'  ;
                else
                    handles.units = '[1/\Delta t]';
                end
                handles.property = handles.dvdy(:,:,handles.current);
            case 10
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s]'  ;
                else
                    handles.units = '[1/\Delta t]';
                end
                handles.property = handles.dudx(:,:,handles.current) + handles.dvdy(:,:,handles.current);
            case 11 % s_xy
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s]'  ;
                else
                    handles.units = '[1/\Delta t]';
                end
                handles.property = 0.5*(handles.dvdx(:,:,handles.current) + handles.dudy(:,:,handles.current));
        end
    end % if fluct, else
else % if no ensemble, it might be u or uf
    if get(handles.checkbox_fluct,'Value') == 1  % uf case
        
        % Derivative of the velocity fluctuations field
        dudx = handles.dudx(:,:,handles.current) - handles.dudx(:,:,handles.N+1);
        dudy = handles.dudy(:,:,handles.current) - handles.dudy(:,:,handles.N+1);
        dvdx = handles.dvdx(:,:,handles.current) - handles.dvdx(:,:,handles.N+1);
        dvdy = handles.dvdy(:,:,handles.current) - handles.dvdy(:,:,handles.N+1);
        
        
        switch get(handles.popupmenu_quantity,'Value')
            case 1
                handles.property = [];    % default variable, choose properly later, alex, 16.02
                handles.colorbar_flag = 0;
                
            case 2
                handles.units = handles.velUnits;
                handles.property = handles.uf(:,:,handles.current);
                handles.cmin = min(handles.uf(:)); handles.cmax = max(handles.uf(:));
            case 3
                handles.units = handles.velUnits;
                handles.property = handles.vf(:,:,handles.current);
                handles.cmin = min(handles.vf(:)); handles.cmax = max(handles.vf(:));
            case 4
                handles.units = handles.velUnits;
                handles.property = sqrt(handles.uf(:,:,handles.current).^2+handles.vf(:,:,handles.current).^2);%Velocity Magnitude
                handles.cmin = min(handles.uf(:).^2+handles.vf(:).^2);
                handles.cmax = max(handles.uf(:).^2+handles.vf(:).^2);
            case 5
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s]'  ;
                else
                    handles.units = '[1/\Delta t]';
                end
                handles.property = dvdx - dudy;
            case 6
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s]'  ;
                else
                    handles.units = '[1/\Delta t]';
                end
                handles.property = dudx;
            case 7
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s]'  ;
                else
                    handles.units = '[1/\Delta t]';
                end
                handles.property = dudy;
            case 8
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s]'  ;
                else
                    handles.units = '[1/\Delta t]';
                end
                handles.property = dvdx;
            case 9
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s]'  ;
                else
                    handles.units = '[1/\Delta t]';
                end
                handles.property = dvdy;
            case 10
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s]'  ;
                else
                    handles.units = '[1/\Delta t]';
                end
                handles.property = dudx+dvdy;
            case 11 % s_xy
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s]'  ;
                else
                    handles.units = '[1/\Delta t]';
                end
                handles.property = 0.5*(dvdx+dudy);
        end
    else % fluct, means nothing is selected, u values
        switch get(handles.popupmenu_quantity,'Value')
            case 1
                handles.property = [];
                handles.colorbar_flag = 0;
                
            case 2
                handles.units = handles.velUnits;
                handles.property = handles.u(:,:,handles.current);
            case 3
                handles.units = handles.velUnits;
                handles.property = handles.v(:,:,handles.current);
            case 4
                handles.units = handles.velUnits;
                handles.property = sqrt(handles.u(:,:,handles.current).^2 + handles.v(:,:,handles.current).^2);%Velocity Magnitude
            case 5
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s]'  ;
                else
                    handles.units = '[1/\Delta t]';
                end
                handles.property = handles.dvdx(:,:,handles.current) - handles.dudy(:,:,handles.current); % omega;
            case 6 % dudx
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s]'  ;
                else
                    handles.units = '[1/\Delta t]';
                end
                handles.property = handles.dudx(:,:,handles.current);
            case 7
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s]'  ;
                else
                    handles.units = '[1/\Delta t]';
                end
                handles.property = handles.dudy(:,:,handles.current);
            case 8
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s]'  ;
                else
                    handles.units = '[1/\Delta t]';
                end
                handles.property = handles.dvdx(:,:,handles.current);
            case 9
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s]'  ;
                else
                    handles.units = '[1/\Delta t]';
                end
                handles.property = handles.dvdy(:,:,handles.current);
            case 10
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s]'  ;
                else
                    handles.units = '[1/\Delta t]';
                end
                handles.property = handles.dudx(:,:,handles.current) + handles.dvdy(:,:,handles.current);
            case 11 % s_xy
                if ~isempty(findstr(handles.velUnits,'s'))
                    handles.units='[1/s]'  ;
                else
                    handles.units = '[1/\Delta t]';
                end
                handles.property = 0.5*(handles.dvdx(:,:,handles.current) + handles.dudy(:,:,handles.current));
        end
    end
end

if ~isempty(handles.property) % something was selected
    if strcmp( get(handles.edit_arrow_size,'Enable'),'on')
        set(handles.checkbox_arrow_color,'Enable','on');
    end
else
    handles.colorbar_flag = 0;
end

% Record the string which is in the quantity popupmenu in something
tmp = cellstr(get(handles.popupmenu_quantity,'String'));
if strcmp(handles.previous_quantity,tmp{get(handles.popupmenu_quantity,'Value')}) == 0 && ...
        get(handles.popupmenu_eachfield,'Value') == 3 % all fields, update it
    handles.previous_quantity = tmp{get(handles.popupmenu_quantity,'Value')};
    popupmenu_eachfield_Callback(handles.fig, [], handles);
else
    handles.previous_quantity = tmp{get(handles.popupmenu_quantity,'Value')};
    guidata(handles.fig,handles);
    update_gui(handles.fig,[],handles);
end

% if strcmp(get(handles.ed_text,'Visible'),'on')
%     if ~isempty(handles.property)
%         handles.ed_mean_value=num2str(mean(handles.property(:)));
%         handles.ed_std_value=num2str(std(handles.property(:)));
%         handles.ed_min_value=num2str(min(handles.property(:)));
%         handles.ed_max_value=num2str(max(handles.property(:)));
%         set(handles.ed_text,'String',handles.previous_quantity);
%         set(handles.ed_mean,'String',handles.ed_mean_value);
%         set(handles.ed_std,'String',handles.ed_std_value);
%         set(handles.ed_min,'String',handles.ed_min_value);
%         set(handles.ed_max,'String',handles.ed_max_value);
%         guidata(handles.fig,handles);
%     end;
% end;

% --------------------------------------------------------------------
function edit_numcolors_Callback(h, ~, handles, varargin)
% change number of colors
handles.numcolors = str2double(get(h,'String'));
if isempty(handles.numcolors),
    set(h,'String',10);
    handles.numcolors = 10;
end
guidata(handles.fig,handles);
update_gui(handles.fig,[],handles);


% --------------------------------------------------------------------
function update_gui(~, ~, handles, varargin)
% update_gui is responsible for update of the screen with current property and contour type

axes(handles.axes_main); %#ok<MAXES>
delete(get(handles.axes_main,'children'));

if get(handles.popupmenu_quantity,'Value') == 1 ,
    set(handles.checkbox_label,'Enable','Off');
    set(handles.checkbox_colorbar,'Enable','Off');
    set(handles.text_numberofcolors,'Enable','Off');
    set(handles.edit_numcolors,'Enable','Off');
    set(handles.popupmenu_eachfield,'Enable','Off');
    set(handles.edit_max_clim,'Enable','Off');
    set(handles.edit_min_clim,'Enable','Off');
    set(handles.pushbutton_set_clim,'Enable','Off');
    set(handles.popupmenu_contour_type,'Enable','Off');
else
    set(handles.checkbox_label,'Enable','On');
    set(handles.checkbox_colorbar,'Enable','On');
    set(handles.text_numberofcolors,'Enable','On');
    set(handles.edit_numcolors,'Enable','On');
    set(handles.popupmenu_eachfield,'Enable','On');
    set(handles.edit_max_clim,'Enable','On');
    set(handles.edit_min_clim,'Enable','On');
    set(handles.pushbutton_set_clim,'Enable','On');
    set(handles.popupmenu_contour_type,'Enable','On');
end

%-------------- Contour -----------------
if ~isempty(handles.property)
    
    switch (get(handles.popupmenu_contour_type,'Value'))
        case 1
            cla reset;
        case 2
            % Flood contour - type
            [handles.C,handles.CH] = contourf(handles.x,handles.y,handles.property,handles.numcolors);
            set(handles.CH,'edgecolor','none');
            % set (handles.CH,'FaceAlpha',0.5);    % if you want transparency, uncomment this.... very slow!!!
        case 3
            %     Color Line
            [handles.C,handles.CH] = contour(handles.x,handles.y,handles.property,handles.numcolors);
        case 4
            %Flood + Line
            [handles.C,handles.CH] = contourf(handles.x,handles.y,handles.property,handles.numcolors);
        case 5
            %Black Line
            [handles.C,handles.CH] = contour(handles.x,handles.y,handles.property,handles.numcolors);
            set(handles.CH,'edgecolor','black');
    end
    
    % -----------------------------------------
    handles.climit_prev = get(gca,'clim');
    % ----------------------------------------------------------
    if handles.alltodisp == 1                  % all to display checkbox is processed here
        set(gca,'CLim',handles.climit);
    elseif handles.allfields == 1               % all_fields checkbox is processed here
        if get(handles.popupmenu_eachfield,'Value')==4;  % in case of manual
            handles.cmin = str2double(get(handles.edit_min_clim,'String'));
            handles.cmax = str2double(get(handles.edit_max_clim,'String'));
        end
        set(gca,'CLim',[handles.cmin handles.cmax]);
    end
    % ----------------------------------------------------------
    if handles.labelit == 1  % label if necessary
        if get(handles.popupmenu_contour_type,'Value')>1
            clabel(handles.C,handles.CH);
        end;
    end
    % -----------------------------------------------------------
    if handles.colorbar_flag                       % colorbar if necessary
        handles.colorbar = colorbar('peer',handles.axes_main,'East');
    end
    % -------------------------------------------------
else    % if property is empty
    cla reset;   %bugfix 0304 ----- if property is empty dont dispay map...
    handles.color = 0;      % disable color
    set(handles.checkbox_arrow_color,'Enable','off');
    set(handles.checkbox_arrow_color,'Value',0);
end
% ----------------------end of ~isempty branch ---------------
% create arrows
if get(handles.checkbox_arrow,'Value') == 1
    if handles.color == 1                 % arrows colored / b&w
        hold on;
        if get(handles.checkbox_fluct,'Value') == 1 % checkbox_fluct On
            if ~isfield(handles ,'uf')
                for i = 1:handles.N
                    handles.uf(:,:,i) = handles.u(:,:,i) - handles.u(:,:,handles.N+1);
                end
            end
            if ~isfield(handles ,'vf')
                for i = 1:handles.N
                    handles.vf(:,:,i) = handles.v(:,:,i) - handles.v(:,:,handles.N+1);
                end
            end
            
            if get(handles.checkbox_ensemble,'Value') == 0
                handles.color_quiver = quiverc(handles.x,handles.y,handles.uf(:,:,handles.current),handles.vf(:,:,handles.current),handles.arrow_scale,handles.property);
            else
                handles.color_quiver = quiverc(handles.x,handles.y,handles.u(:,:,handles.current),handles.v(:,:,handles.current),handles.arrow_scale,handles.property);
                % nothing to display, arrows of ensemble + fluctuation = 0
                % by definition
                set(handles.color_quiver,'Visible','Off');
            end
        else
            handles.color_quiver = quiverc(handles.x,handles.y,...
                handles.u(:,:,handles.current),...
                handles.v(:,:,handles.current),...
                handles.arrow_scale,handles.property);
        end
        hold off;
        if handles.colorbar_flag == 1 % colorbar if necessary
            delete(handles.colorbar);
            handles.colorbar = colorbar('peer',handles.axes_main,'East');
% handles.colorbar = mcolorbar('EastOutside','peer',handles.axes_main);
        end
    else
        hold on;
        % ----------------------checknox_fluct is checked -------------
        if get(handles.checkbox_fluct,'Value') == 1
            if ~isfield(handles ,'uf')
                for i = 1:handles.N
                    handles.uf(:,:,i) = handles.u(:,:,i) - ...
                        handles.u(:,:,handles.N+1);
                end
            end
            if ~isfield(handles ,'vf')
                for i = 1:handles.N
                    handles.vf(:,:,i) = handles.v(:,:,i) - ...
                        handles.v(:,:,handles.N+1);
                end
            end
            
            if get(handles.checkbox_ensemble,'Value') == 0
                handles.quiverH = quiver(handles.x,handles.y, ...
                    handles.uf(:,:,handles.current), ...
                    handles.vf(:,:,handles.current),...
                    handles.arrow_scale,'k');
            else
                handles.quiverH = quiver(handles.x,handles.y,handles.u(:,:,handles.current),handles.v(:,:,handles.current),handles.arrow_scale,'k');
                % nothing to display, arrows of ensemble + fluctuation = 0
                % by definition
                set(handles.quiverH,'Visible','Off');
            end
        else
            handles.quiverH = quiver(handles.x,handles.y, ...
                handles.u(:,:,handles.current),...
                handles.v(:,:,handles.current),...
                handles.arrow_scale,'k');
        end
        hold off;
    end
end

set(handles.axes_main,'XLim',[min(handles.x(:)),max(handles.x(:))]);
set(handles.axes_main,'YLim',[min(handles.y(:)),max(handles.y(:))]);
xlabel(['x ',handles.xUnits]); % [m]');
ylabel(['y ',handles.xUnits]); % 'y [m]');

if isfield(handles,'colorbar') && get(handles.checkbox_colorbar,'Value') == 1
    axpos = handles.axpos; % 10.04.06
    set(handles.axes_main,'Units','normalized','Position',...
        [axpos(1),axpos(2),axpos(3)-.025,axpos(4)]);
    if ishandle(handles.colorbar)
        set(handles.colorbar,'Units','normalized','Position',...
            [axpos(1)+axpos(3)+.02,axpos(2),.025,axpos(4)]);
    elseif handles.colorbar_flag == 1
        handles.colorbar = colorbar('peer',handles.axes_main,'East');
        set(handles.colorbar,'Units','normalized','Position',...
            [axpos(1)+axpos(3)+.02,axpos(2),.025,axpos(4)]);
 %    handles.colorbar = mcolorbar('EastOutside','peer',handles.axes_main);
    end
    
end

set(gca,'ydir','reverse'); % Alex, 27.12.10, Hadar computer.
guidata(handles.fig,handles);

% --------------------------------------------------------------------
function edit_min_clim_Callback(h, ~, handles, varargin)
% first editbox for 'manual' checkbox callback
handles.cmin = str2double(get(h,'String')); % update cmin
handles.alltodisp = 0;
handles.allfields = 1;
guidata(handles.fig,handles);
% update_gui(handles.fig,[],handles);

% --------------------------------------------------------------------
function edit_max_clim_Callback(h, ~, handles, varargin)
% second editbox for 'manual' checkbox callback
handles.cmax = str2double(get(h,'String')); % get new value for cmax
handles.alltodisp = 0;
handles.allfields = 1;
guidata(handles.fig,handles);
% update_gui(handles.fig,[],handles);

% --------------------------------------------------------------------
function pushbutton_set_clim_Callback(~, ~, handles, varargin)
% set button callback
handles.cmin = eval(get(handles.edit_min_clim,'String'));
handles.cmax = eval(get(handles.edit_max_clim,'String'));
if handles.cmin >= handles.cmax
    errordlg('Wrong limits, min should be less than max', 'Error','modal');
    current_clim = get(handles.axes_main,'clim');
    handles.cmin = current_clim(1);
    handles.cmax = current_clim(2);
    set(handles.edit_min_clim,'String',sprintf('%3.2f',handles.cmin));
    set(handles.edit_max_clim,'String',sprintf('%3.2f',handles.cmax));
end
guidata(handles.fig,handles);
update_gui(handles.fig,[],handles);

% --------------------------------------------------------------------
function pushbutton_previous_Callback(~, ~, handles, varargin)
if handles.current > 1
    handles.current = handles.current - 1;  % update handles.current
    % display num. of current file being processed
    set(handles.edit_current,'String',handles.current); 
    delete(get(handles.axes_main,'children'));
    guidata(handles.fig,handles);
    popupmenu_quantity_Callback(handles.fig, [], handles);
else
    beep;
end

% --------------------------------------------------------------------
function pushbutton_next_Callback(~, ~, handles, varargin)
if handles.current < handles.N
    handles.current = handles.current + 1;      % update handles.current
    set(handles.edit_current,'String',handles.current);
    guidata(handles.fig,handles);
    delete(get(handles.axes_main,'children'));
    popupmenu_quantity_Callback(handles.fig, [], handles);
else
    beep;
end

% --------------------------------------------------------------------
function edit_current_Callback(~, ~, handles, varargin)
tmp = eval(get(handles.edit_current,'String'));
if tmp > 0 && tmp <= handles.N % valid number of map
    handles.current = tmp;
    guidata(handles.fig,handles);
    popupmenu_quantity_Callback(handles.fig, [], handles);
else
    beep
    set(handles.edit_current,'String',handles.current);
end

% --------------------------------------------------------------------
function edit_arrow_size_Callback(h, ~, handles, varargin)
% Arrows 'scale' editbox callback
handles.arrow_scale = eval(get(h,'String'));     % update handles.scale
if handles.arrow_scale == 0 || isempty(handles.arrow_scale) % no scaling
    handles.arrow_scale = [];
end
guidata(handles.fig,handles);
update_gui(handles.fig,[],handles);

% --------------------------------------------------------------------
function pushbutton_animate_Callback(~, ~, handles, varargin)
% animate button callback

if get(handles.pushbutton_animate,'Value') == 1
    %     handles.tmpenableprop = strcmp(get(handles.spatial_controls,'Enable'),'on')
    %     setenableprop(handles.spatial_controls,'off');
    %     setenableprop(handles.pushbutton_animate,'on');
    startpoint = handles.current;
    for i = startpoint:handles.N
        if get(handles.pushbutton_animate,'Value') == 0
            break;
        end
        handles.current = i;
        set(handles.edit_current,'String',handles.current);
        delete(get(handles.axes_main,'children'));
        guidata(handles.fig,handles);
        popupmenu_quantity_Callback(handles.fig, [], handles);
        drawnow;
    end
else
    
end
set(handles.pushbutton_animate,'Value',0);
% setenableprop(handles.spatial_controls(handles.tmpenableprop),'on');
% popupmenu_quantity_Callback(handles.fig, [], handles);
% guidata(handles.fig,handles);

% --------------------------------------------------------------------
function pushbutton_save_movie_Callback(~, ~, handles, varargin)

if get(handles.pushbutton_save_movie,'Value') == 1
    file = inputdlg('File Name','Input File Name for the movie');
    if isempty(file) || exist(file{1},'file') || exist([file{1},'.avi'],'file')
        set(handles.pushbutton_save_movie,'Value',0);
        return
    end
    handles.mov = avifile(file{1},'compression','none','quality',100,'fps',15);
    % It is possible to change the compression and the video codec, and the
    % frame rate
    % mov = avifile(file{1},'compression','Indeo5','quality',100,'fps',5);
    
    
    startpoint = handles.current;
    for i = startpoint:handles.N
        
        handles.current = i;
        set(handles.edit_current,'String',handles.current);
        delete(get(handles.axes_main,'children'));
        guidata(handles.fig,handles);
        popupmenu_quantity_Callback(handles.fig, [], handles);
        F = getframe(handles.axes_main);
        if get(handles.pushbutton_save_movie,'Value') == 0
            break;
        end
        handles.mov = addframe(handles.mov,F);
    end
else
    if isfield(handles,'mov')
        handles.mov = close(handles.mov);
        handles = rmfield(handles,'mov');
    end
    set(handles.pushbutton_save_movie,'Value',0);
end
set(handles.pushbutton_save_movie,'Value',0);
% guidata(handles.fig,handles);


% --------------------------------------------------------------------
function checkbox_label_Callback(h, ~, handles, varargin)
if (get(h,'Value') == get(h,'Max'))
    if get(handles.popupmenu_contour_type,'Value')>1
        handles.labelit = 1;
        
        clabel(handles.C,handles.CH);
        guidata(handles.fig,handles);
    else
        set(handles.checkbox_label,'Value',0);
    end;
    %update_gui(handles.fig,[],handles);
else
    handles.labelit = 0;
    guidata(handles.fig,handles);
    update_gui(handles.fig,[],handles);
end

% --------------------------------------------------------------------
function checkbox_colorbar_Callback(h, ~, handles, varargin)
if get(h,'Value') == 1
    handles.colorbar_flag = 1;
    guidata(handles.fig,handles);
    update_gui(handles.fig,[],handles);
else
    delete(findobj(handles.fig,'Tag','Colorbar'));
    handles.colorbar_flag = 0;
    set(handles.axes_main,'position',handles.axpos);
    guidata(handles.fig,handles);
end
% update_gui(handles.fig,[],handles);

% --------------------------------------------------------------------
function popupmenu_eachfield_Callback(~, ~, handles)
% Each Field, All to display, all fields & manual processing

% get value and assign selected property to handles.property which is default to display
val = get(handles.popupmenu_eachfield,'Value');
switch val
    case 1 % Each Field
        set(handles.edit_min_clim,'Visible','Off');
        set(handles.edit_max_clim,'Visible','Off');
        set(handles.pushbutton_set_clim,'Visible','Off');
        handles.alltodisp = 0;                             % unset all other flags
        handles.allfields = 0;
    case 2 % All to Display
        set(handles.edit_min_clim,'Visible','Off');
        set(handles.edit_max_clim,'Visible','Off');
        set(handles.pushbutton_set_clim,'Visible','Off');
        handles.alltodisp = 1;                % set alltodisp flag which is processed in popupmenu2_callback
        handles.climit = handles.climit_prev; % backup current climit which will be changed in popupmenu2_callback
    case 3 % All fields
        set(handles.edit_min_clim,'Visible','Off');
        set(handles.edit_max_clim,'Visible','Off');
        set(handles.pushbutton_set_clim,'Visible','Off');
        handles.alltodisp = 0;
        handles.allfields = 1;            % all_fields is processed in popupmenu2_callback
        
        if get(handles.checkbox_fluct,'Value') == 1  % fluctuations
            switch get(handles.popupmenu_quantity,'Value')
                case 1
                    handles.cmin = 0; handles.cmax = 1;
                case 2
                    handles.cmin = min(handles.uf(:)); handles.cmax = max(handles.uf(:));
                case 3
                    handles.cmin = min(handles.vf(:)); handles.cmax = max(handles.vf(:));
                case 4
                    handles.cmin = min(handles.uf(:).^2+handles.vf(:).^2);
                    handles.cmax = max(handles.uf(:).^2+handles.vf(:).^2);
                case 5
                    handles.cmin = Inf;
                    handles.cmax = -Inf;
                    for i = 1:handles.N
                        tmp = (handles.dvdx(:,:,i) - handles.dvdx(:,:,handles.N+1))-...
                            (handles.dudy(:,:,i) - handles.dudy(:,:,handles.N+1));
                        handles.cmin = min(handles.cmin,min(tmp(:)));
                        handles.cmax = max(handles.cmax,max(tmp(:)));
                    end
                case 6
                    handles.cmin = Inf;
                    handles.cmax = -Inf;
                    for i = 1:handles.N
                        tmp = handles.dudx(:,:,i) - handles.dudx(:,:,handles.N+1);
                        handles.cmin = min(handles.cmin,min(tmp(:)));
                        handles.cmax = max(handles.cmax,max(tmp(:)));
                    end
                case 7
                    %                 handles.property = dudy;
                    handles.cmin = Inf;
                    handles.cmax = -Inf;
                    for i = 1:handles.N
                        tmp = handles.dudy(:,:,i) - handles.dudy(:,:,handles.N+1);
                        handles.cmin = min(handles.cmin,min(tmp(:)));
                        handles.cmax = max(handles.cmax,max(tmp(:)));
                    end
                    
                case 8
                    %                 handles.property = dvdx;
                    handles.cmin = Inf;
                    handles.cmax = -Inf;
                    for i = 1:handles.N
                        tmp = handles.dvdx(:,:,i) - handles.dvdx(:,:,handles.N+1);
                        handles.cmin = min(handles.cmin,min(tmp(:)));
                        handles.cmax = max(handles.cmax,max(tmp(:)));
                    end
                    %
                case 9
                    %                 handles.property = dvdy;
                    handles.cmin = Inf;
                    handles.cmax = -Inf;
                    for i = 1:handles.N
                        tmp = handles.dvdy(:,:,i) - handles.dvdy(:,:,handles.N+1);
                        handles.cmin = min(handles.cmin,min(tmp(:)));
                        handles.cmax = max(handles.cmax,max(tmp(:)));
                    end
                    
                case 10
                    %                 handles.property = dudx+dvdy;
                    handles.cmin = Inf;
                    handles.cmax = -Inf;
                    for i = 1:handles.N
                        tmp = (handles.dudx(:,:,i) - handles.dudx(:,:,handles.N+1))+...
                            (handles.dvdy(:,:,i) - handles.dvdy(:,:,handles.N+1));
                        handles.cmin = min(handles.cmin,min(tmp(:)));
                        handles.cmax = max(handles.cmax,max(tmp(:)));
                    end
                case 11 % s_xy = 0.5*(dvdx+dudy);
                    handles.cmin = Inf;
                    handles.cmax = -Inf;
                    for i = 1:handles.N
                        tmp = 0.5*(handles.dudx(:,:,i) - handles.dudx(:,:,handles.N+1))-...
                            (handles.dvdy(:,:,i) - handles.dvdy(:,:,handles.N+1));
                        handles.cmin = min(handles.cmin,min(tmp(:)));
                        handles.cmax = max(handles.cmax,max(tmp(:)));
                    end
            end
        else                                    % no if, instantaneous
            switch get(handles.popupmenu_quantity,'Value')
                case 1
                    handles.cmin = 0; handles.cmax = 1;
                case 2
                    handles.cmin = min(handles.u(:)); handles.cmax = max(handles.u(:));
                case 3
                    handles.cmin = min(handles.v(:)); handles.cmax = max(handles.v(:));
                case 4
                    handles.cmin = min(handles.u(:).^2+handles.v(:).^2);
                    handles.cmax = max(handles.u(:).^2+handles.v(:).^2);
                case 5
                    handles.cmin = min(handles.dvdx(:)- handles.dudy(:));
                    handles.cmax = max(handles.dvdx(:)- handles.dudy(:));
                case 6
                    handles.cmin = min(handles.dudx(:));
                    handles.cmax = max(handles.dudx(:));
                case 7
                    handles.cmin = min(handles.dudy(:));
                    handles.cmax = max(handles.dudy(:));
                case 8
                    handles.cmin = min(handles.dvdx(:));
                    handles.cmax = max(handles.dvdx(:));
                case 9
                    handles.cmin = min(handles.dvdy(:));
                    handles.cmax = max(handles.dvdy(:));
                case 10
                    handles.cmin = min(handles.dudx(:)+handles.dvdy(:));
                    handles.cmax = max(handles.dudx(:)+handles.dvdy(:));
                case 11 % s_xy = 0.5*(dvdx+dudy);
                    handles.cmin = 0.5*min(handles.dvdx(:)+handles.dudy(:));
                    handles.cmax = 0.5*max(handles.dvdx(:)+handles.dudy(:));
            end
        end
        
    case 4 % Manual
        set(handles.edit_min_clim,'Visible','On'); % show 2 editboxes to enter the values
        set(handles.edit_max_clim,'Visible','On');
        set(handles.pushbutton_set_clim,'Visible','On');
        handles.alltodisp = 0;
        % we can use use allfields flag to process 'Manual'  but with cmin & cmax that are entered manually
        handles.allfields = 1;
        current_clim=get(handles.axes_main,'clim');
        set(handles.edit_min_clim,'String',sprintf('%3.2f',current_clim(1)));
        set(handles.edit_max_clim,'String',sprintf('%3.2f',current_clim(2)));
end

guidata(handles.fig,handles);
update_gui(handles.fig,[],handles);


%--------------------------------------------------------------------------
function checkbox_ensemble_Callback(~, ~, handles)

% global inst_list mean_list fluct_list fluct_mean_list;
set (handles.arrow_ctrls,'Enable','on');
if get(handles.checkbox_ensemble,'Value') == 1
    
    % No meaning for All fields and All to Display for Ensemble cases
    % However, we do not modify the actual subroutines, in order not to
    % destroy the code there
    set(handles.popupmenu_eachfield,'String','Each Field|---|---|Manual');
    if get(handles.checkbox_fluct,'Value') == 1
        set (handles.popupmenu_quantity,'String',handles.fluct_mean_list);
        set (handles.arrow_ctrls,'Enable','off');
    else
        set (handles.popupmenu_quantity,'String',handles.mean_list);
    end
    handles.current_index = handles.current;
    handles.current = handles.N + 1;
    set(handles.pushbutton_previous,'Enable','off');
    set(handles.pushbutton_next,'Enable','off');
    set(handles.pushbutton_animate,'Enable','off');
    set(handles.pushbutton_save_movie,'Enable','off');
    set(handles.edit_current,'Enable','off');
else
    set(handles.popupmenu_eachfield,'String',...
    'Each Field|All to Display|All Fields|Manual');
    if (get(handles.checkbox_fluct,'Value') == ...
            get(handles.checkbox_fluct,'Max'));
        set (handles.popupmenu_quantity,'String',handles.fluct_list);
    else
        set (handles.popupmenu_quantity,'String',handles.inst_list);
    end
    if handles.current == handles.N + 1;
        handles.current = handles.current_index;
    end
    set(handles.pushbutton_previous,'Enable','on');
    set(handles.pushbutton_next,'Enable','on');
    set(handles.pushbutton_animate,'Enable','on');
    set(handles.pushbutton_save_movie,'Enable','on');
    set(handles.edit_current,'Enable','on');
end;

set(handles.popupmenu_quantity,'Value',1);         % 2104 by Denis
handles.property=[];  handles.C=[]; handles.CH=[];  %0304 by Denis
guidata(handles.fig,handles);
update_gui(handles.fig,[],handles);



%--------------------------------------------------------------
function checkbox_fluct_Callback(~, ~, handles)
set (handles.arrow_ctrls,'Enable','on');
if (get(handles.checkbox_fluct,'Value') == 1)
    if (get(handles.checkbox_ensemble,'Value') == 1)
        set (handles.popupmenu_quantity,'String',handles.fluct_mean_list);
        set (handles.arrow_ctrls,'Enable','off');
    else
        set (handles.popupmenu_quantity,'String',handles.fluct_list);
    end
else
    if (get(handles.checkbox_ensemble,'Value') == 1)
        set (handles.popupmenu_quantity,'String',handles.mean_list);
    else
        set (handles.popupmenu_quantity,'String',handles.inst_list);
    end
end;
set(handles.popupmenu_quantity,'Value',1);     % 2104 by Denis
handles.property    =   [];
handles.C           =   [];
handles.CH          =   [];                        % 0304 by Denis

guidata(handles.fig,handles);
update_gui(handles.fig,[],handles);

% ------------------------------------------------------------------------
function File_Callback(~, ~, ~)




% ------------------   Load and prepare data module --------------------------------------------------
function loadVec_Callback(~, ~, handles)
global orighandles;
if isfield(handles,'restoreorig')
    handles = orighandles;
end
handles.restoreorig = 1;


try
    % [gui_files,gui_path,handles.dt,handles.scale,handles.state3d] = cil_uigetfiles;
    gui_files = uipickfiles;
    handles.dt = 1;
    handles.state3d = 0;
    handles.scale = 1;
    [gui_path,~,~] = fileparts(gui_files{1});
    
    
    handles.N = length(gui_files); % number of files selected
    if  handles.N > 0
        handles.files = gui_files;
        handles.path = gui_path;
        set(handles.fig,'pointer','watch');
    else
        return
    end
    
    switch handles.state3d
        case 1
            if ~isempty(findstr(handles.files{1},'v3d'))
                [handles.xUnits,handles.velUnits,d] = ...
                    svecread(fullfile(handles.path,handles.files{1}));
                [rows,cols,~] = size(d);
                [handles.u,handles.v,handles.w] = ...
                    deal(zeros(rows,cols,handles.N+1)); % 11.04.04, Alex
                % Bug fixes June 26, 2004, for the first release.
                handles.x = d(:,:,1);
                handles.y = d(:,:,2);
                handles.z = d(:,:,3);
                handles.u(:,:,1) = d(:,:,4);
                handles.v(:,:,1) = d(:,:,5);
                handles.w(:,:,1) = d(:,:,6);
%                 handles.x = d(:,:,1)*handles.scale/1000;
%                 handles.y = d(:,:,2)*handles.scale/1000;
%                 handles.z = d(:,:,3)*handles.scale/1000;
%                 handles.u(:,:,1) = d(:,:,4)*handles.scale/1000/handles.dt;
%                 handles.v(:,:,1) = d(:,:,5)*handles.scale/1000/handles.dt;
%                 handles.w(:,:,1) = d(:,:,6)*handles.scale/1000/handles.dt;
                
                for i = 2:handles.N
                    d = svecread([handles.path,filesep,handles.files{i}],1,8);
% handles.u(:,:,i) = d(:,:,4)*handles.scale/1000/handles.dt;
% handles.v(:,:,i) = d(:,:,5)*handles.scale/1000/handles.dt;
% handles.w(:,:,i) = d(:,:,6)*handles.scale/1000/handles.dt;
                    handles.u(:,:,i) = d(:,:,4);
                    handles.v(:,:,i) = d(:,:,5);
                    handles.w(:,:,i) = d(:,:,6);
                end
                clear d
            end
        case 0
            if ~isempty(findstr(lower(handles.files{1}),'vec'))            % process .vec files
                % read the first file, determine the size
                % [handles.xUnits,handles.velUnits,d] = vecread(fullfile(handles.path,handles.files{1}));
                [handles.xUnits,handles.velUnits,d] = vecread(handles.files{1});
                [rows,cols,~] = size(d);
                [handles.u,handles.v] = deal(zeros(rows,cols,handles.N+1)); % 11.04.04, Alex
                handles.x           = d(:,:,1);
                handles.y           = d(:,:,2);
                handles.u(:,:,1)    = d(:,:,3);
                handles.v(:,:,1)    = d(:,:,4);
                % Bug fixes, June 26, 2004
                %                 handles.x           = d(:,:,1)*handles.scale/1000;
                %                 handles.y           = d(:,:,2)*handles.scale/1000;
                %                 handles.u(:,:,1)    = d(:,:,3)*handles.scale/1000/handles.dt;
                %                 handles.v(:,:,1)    = d(:,:,4)*handles.scale/1000/handles.dt;
                %
                for i = 2:handles.N
                    % d = vecread([handles.path,filesep,handles.files{i}],1,5);
                    d = vecread(handles.files{i},1,5);
                    %                     handles.u(:,:,i) = d(:,:,3)*handles.scale/1000/handles.dt;
                    %                     handles.v(:,:,i) = d(:,:,4)*handles.scale/1000/handles.dt;
                    handles.u(:,:,i) = d(:,:,3);
                    handles.v(:,:,i) = d(:,:,4);
                end
                clear d
            elseif ~isempty(findstr(lower(handles.files{1}),'txt')) % new files, created for stratified
                % project, probably by Zach Taylor version of OpenPIV C++
                % the format is different from our ".txt" files which have
                % no headers, and different from VEC format of Insight 3G,
                % but has a header, single line that one can get out using:
%                 fid = fopen(handles.files{1},'r');
%                 header = fgetl(fid);
%                 fclose(fid);
                % get units - TODO. use findstr(header, '[') and ']'
                handles.xUnits = 'pixels';
                handles.velUnits = 'pixels';
                
                
                % and the read the data in the file using:
                d = dlmread(handles.files{1},'',1,0);
                
                % we need to know the reshape size:
                rows = find(diff(d(:,1))<0,1);
                cols = length(d(:,1))/rows;
                
                [handles.u,handles.v] = deal(zeros(rows,cols,handles.N+1)); % 11.04.04, Alex
                handles.x           = reshape(d(:,1),rows,cols);
                handles.y           = reshape(d(:,2),rows,cols);
                handles.u(:,:,1)    = reshape(d(:,3),rows,cols);
                handles.v(:,:,1)    = reshape(d(:,4),rows,cols);
                
                for i = 2:handles.N
                    d = dlmread(handles.files{i},'',1,0);
                    handles.u(:,:,i)    = reshape(d(:,3),rows,cols);
                    handles.v(:,:,i)    = reshape(d(:,4),rows,cols);
                end
                clear d
                
            end
            
            
            
            
    end
    
catch ME
    % errordlg('Something wrong with vector files');
    errordlg(ME.message)
    set(handles.fig,'pointer','arrow');
    return
end

handles.current = 1;                      % current file beeing displayed
% Display first file number, total number of files
set(handles.edit_current,'String',handles.current);
set(handles.edit_numfields,'String',handles.N);

% Initialize color, number of colors for contours
handles.numcolors = 10;                   % default number of colors
set(handles.edit_numcolors,'String', handles.numcolors);


% handles.dx = abs(handles.x(1,1) - handles.x(1,2));
% handles.dy = abs(handles.y(1,1) - handles.y(2,1));
% Bug, fixed at 16.06 after the email of Hai
handles.dx = handles.x(1,2) - handles.x(1,1);
handles.dy = handles.y(2,1) - handles.y(1,1);

% 19.03.08 - Alex found one case in which the VEC file
% is transposed (x is from top to bottom, etc.)
% therefore dx gave 0 and vorticity was NaN and Inf
% just for this case:

if handles.dx == 0 && handles.dy == 0
    handles.x = permute(handles.x,[2 1]);
    handles.y = permute(handles.y,[2 1]);
    handles.u = permute(handles.u,[2 1 3]);
    handles.v = permute(handles.v,[2 1 3]);
    handles.dx = handles.x(1,2) - handles.x(1,1);
    handles.dy = handles.y(2,1) - handles.y(1,1);
    tmprows = cols;
    cols = rows;
    rows = tmprows;
end



% Cat the mean values at the end
handles.u(:,:,handles.N+1) = mean(handles.u(:,:,1:handles.N),3);
handles.v(:,:,handles.N+1) = mean(handles.v(:,:,1:handles.N),3);
if handles.state3d
    handles.w(:,:,handles.N+1) = mean(handles.w(:,:,1:handles.N),3);
    handles.wf = zeros(rows,cols,handles.N);
    for i = 1:handles.N
        handles.wf(:,:,i) = handles.w(:,:,i) - handles.w(:,:,handles.N+1);
    end
end






% Preallocate memory and calculate all the necessary quantities
% Fluctuations (last one is zero, we do not need it)
handles.uf = zeros(rows,cols,handles.N);
handles.vf = zeros(rows,cols,handles.N);
%
for i = 1:handles.N
    handles.vf(:,:,i) = handles.v(:,:,i) - handles.v(:,:,handles.N+1);
end
%
for i = 1:handles.N
    handles.uf(:,:,i) = handles.u(:,:,i) - handles.u(:,:,handles.N+1);
end

% Derivatives
handles.dudx = zeros(rows,cols,handles.N+1);
handles.dudy = zeros(rows,cols,handles.N+1);
handles.dvdx = zeros(rows,cols,handles.N+1);
handles.dvdy = zeros(rows,cols,handles.N+1);
%
for i = 1:handles.N+1
    [handles.dudx(:,:,i),handles.dudy(:,:,i)] = ...
        lsgradient(handles.u(:,:,i),handles.dx, handles.dy);
    [handles.dvdx(:,:,i),handles.dvdy(:,:,i)] = ...
        lsgradient(handles.v(:,:,i),handles.dx, handles.dy);
end

% Possible future development, eliminating strong gradients
% on the borders
% handles.dudx(1:2,:,:)       = NaN;
% handles.dudx(end-1:end,:,:) = NaN;
% handles.dudx(:,1:2,:)       = NaN;
% handles.dudx(:,end-1:end,:) = NaN;

%
% More defaults
handles.arrow_scale = 1;                % default scale
set(handles.edit_arrow_size,'String',handles.arrow_scale);

% Default situation, instantaneous, not Ensemble, not fluctuations
set(handles.checkbox_ensemble,'Value',0);
set(handles.checkbox_fluct,'Value',0);

% No colors, no labels
handles.color = 0;                % display colored / black figure
handles.alltodisp = 0;            % default all_to display is unset
handles.allfields = 0;
handles.labelit = 0;              % label for contour
handles.colorbar_flag = 0;
handles.current_index = handles.current;
handles.distribOn=0;
handles.rowlock=0; handles.columnlock=0;
handles.previousSel=[];

% These are future values, for the TimeBox
handles.i=[];
handles.j=[];
handles.PointsH=[];
handles.Allselected=0;

% handles.dx = handles.x(1,2) - handles.x(1,1);
% handles.dy = handles.y(2,1) - handles.y(1,1);

% abs() is added on Dec.13,2004. by Alex on IHW_Liberzon computer
%
handles.gridX = abs(handles.x(1,2) - handles.x(1,1));
handles.gridY = abs(handles.y(1,1) - handles.y(2,1));

% Show some and hide some controls
set(handles.popupmenu_quantity,'Visible','on');
set(handles.popupmenu_quantity,'Value',1);
set(handles.popupmenu_contour_type,'Visible','on');
set(handles.popupmenu_contour_type,'Value',1);
set(handles.popupmenu_eachfield,'Visible','on');
set(handles.popupmenu_eachfield,'Value',1);

% Some only part of them in Enable mode
set(handles.checkbox_arrow,'Enable','On');
set(handles.edit_arrow_size,'Enable','On');
set(handles.checkbox_fluct,'Enable','On');
set(handles.checkbox_ensemble,'Enable','On');
set(handles.popupmenu_quantity,'Enable','On');
set(handles.pushbutton_previous,'Enable','On');
set(handles.pushbutton_next,'Enable','On');
set(handles.edit_current,'Enable','On');
set(handles.pushbutton_animate,'Enable','On');
set(handles.pushbutton_save_movie,'Enable','on');
% set(handles.pushbutton_stats,'Enable','on');


% Properties
set(handles.popupmenu_quantity,'String',handles.inst_list);
handles.property = [];

% ---------------- store all handles of spatial controls ------------------
% handles.spatial_controls=[handles.checkbox_ensemble,handles.checkbox_fluct,handles.checkbox_arrow,handles.checkbox_arrow_color,...
%         handles.checkbox_label,handles.checkbox_colorbar,handles.edit_arrow_size,handles.edit_numcolors,...
%         handles.edit_current,handles.edit_numfields,...
%         handles.text_contour_quantity, handles.text_contourtype, handles.text_numberofcolors, handles.text7,handles.text2,...
%         handles.text_arrow_size, handles.pushbutton_previous, handles.pushbutton_next, ...
%         handles.pushbutton_animate, handles.pushbutton_save_movie, ...
%         handles.frame_controls,handles.frame8,handles.frame_contour_quantity,handles.frame_contour_type,handles.frame7,...
%         handles.frame_arrow, handles.text5, handles.popupmenu_quantity,handles.popupmenu_contour_type,...
%         handles.popupmenu_eachfield];% ,handles.pushbutton_stats];

% -------------- store all handles of select controls ----------------
% handles.select_controls=[handles.pushbutton_selectpoints,handles.pushbutton_selectreg,handles.pushbutton_selectall,...
%          handles.pushbutton_profile1,handles.pushbutton_time,...
%          handles.pushbutton_reset, handles.rowpushbutton, handles.colpushbutton];

handles.arrow_ctrls=[handles.edit_arrow_size,handles.checkbox_arrow,handles.checkbox_arrow_color];
% handles.stat_controls=[handles.ed_mean,handles.ed_std,handles.ed_min,handles.ed_max,handles.ed_text,...
%         handles.ed_pushsavestl,handles.ed_frame,handles.ed_pushclose,handles.ed_textmean,...
%         handles.ed_textstd,handles.ed_textmax,handles.ed_textmin];


set(handles.pushbutton_spatial,'Enable','on');
set(handles.pushbutton_select,'Enable','on');

set(handles.fig,'pointer','arrow');

% added on 10.04.06 for R12SP3 version
handles.axpos = get(handles.axes_main,'Position');

% Update all handles structure
guidata(handles.fig,handles);

% Make default plot
update_gui(handles.fig,[],handles);



% --------------------------------------------------------------------
function exit_Callback(hObject, ~, ~)
% Find the highest parent - figure, and close it.
while ~strcmpi(get(hObject,'Type'),'figure'),
    hObject = get(hObject,'Parent');
end
delete(hObject);

% --- Executes on button press in pushbutton_spatial.
function pushbutton_spatial_Callback(~, ~, handles)

set(handles.spatial_controls,'Visible','on');
set(handles.select_controls,'Visible','off');
set(handles.pushbutton_spatial,'FontWeight','bold');
set(handles.pushbutton_select,'FontWeight','normal');
set(handles.pushbutton_spatial,'String','> Spatial <');
set(handles.pushbutton_select,'String','Select');
val = get(handles.popupmenu_eachfield,'Value');
if val==4   % in case if manual is selected
    set(handles.edit_min_clim,'Visible','On'); % show 2 editboxes to enter the values
    set(handles.edit_max_clim,'Visible','On');
    set(handles.pushbutton_set_clim,'Visible','On');
end

% save selection
handles.i = []; handles.j = [];
handles.rowlock=0; handles.columnlock=0;
handles.previousSel=[];

update_gui(gcbo,[],guidata(gcbo));
guidata(handles.fig,handles);
% update_gui(handles.fig,[],handles);


% -------------------------------------------------------------
%    SELECTION
% -------------------------------------------------------------

% --- Executes on button press in pushbutton_select.
function pushbutton_select_Callback(~, ~, handles)
% hObject    handle to pushbutton_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.spatial_controls,'Visible','off');
% set(handles.pushbutton_set_clim,'Visible','off');
set(handles.select_controls,'Visible','on');
% set(findobj(handles.spatial_controls,'type','uicontrol'),'Enable','Off');
set(findobj(handles.select_controls,'type','uicontrol'),'Enable','On');
set(handles.pushbutton_select,'String','> Select <');
set(handles.pushbutton_spatial,'String','Spatial');


set(handles.pushbutton_spatial,'FontWeight','normal');
set(handles.pushbutton_select,'FontWeight','bold');
% val = get(handles.popupmenu_eachfield,'Value');
handles.i = []; handles.j = [];
handles.rowlock=0; handles.columnlock=0;
handles.previousSel=[];

update_gui(gcbo,[],guidata(gcbo));
guidata(handles.fig,handles);




% --- Executes on button press in pushbutton_selectpoints.
function pushbutton_selectpoints_Callback(~, ~, handles)


set(handles.rowpushbutton,'Enable','off');
set(handles.pushbutton_selectreg,'Enable','off');
set(handles.colpushbutton,'Enable','off');
set(handles.pushbutton_selectall,'Enable','off');



set(handles.axes_main,'NextPlot','Add');
limX = xlim;
limY = ylim;
% leftcolX   = 1;
% bottomrowY = 1;


while 1
    
    [x1,y1, buttonNumber] = ginput(1);
    
    % When the right button is pressed, stop the loop
    if (buttonNumber == 2) || (buttonNumber==3)
        break
    end
    col = fix (( x1 - limX(1,1) )/ handles.gridX+0.5  )+1;
    row = fix(( y1 - limY(1,1) )/ handles.gridY+0.5  )+1;
    
    % check for errors ----------------
    if col<1 || col > (fix((limX(1,2)-limX(1,1))/handles.gridX)+1) ...
            || row < 1 ...
            || row >(fix((limY(1,2)-limY(1,1))/handles.gridY)+1);
        guidata(handles.fig,handles);
        return
    end
    % ---------------------------------
    sizeI = size(handles.i,1);
%     rightcolX = fix(( limX(1,2)-limX(1,1) )/  handles.gridX )+1;
%     uprowY    = fix(( limY(1,2)-limY(1,1) )/  handles.gridY )+1;
    sizeJ = size(handles.j,1);
%     numofcols = rightcolX - leftcolX + 1;
%     numofrows = uprowY - bottomrowY + 1;
    
    handles.i(sizeI+1,1) = row;
    handles.j(sizeJ+1,1) = col;
    
    
    
    line(limX(1,1)+(col-1)*...
        handles.gridX,limY(1,1)+(row-1)*handles.gridY,...
        'Marker','o','Color','k','MarkerSize',8);
    
end
disp('Points selected')
disp([handles.i,handles.j])
disp(handles.x(handles.i,handles.j))
guidata(handles.fig,handles);
% update_gui(handles.fig,[],handles);



% --- Executes on button press in pushbutton_selectreg.
function pushbutton_selectreg_Callback(~, ~, handles)
% check if all selected is off
set(handles.pushbutton_selectpoints,'Enable','off');
set(handles.rowpushbutton,'Enable','off');
set(handles.colpushbutton,'Enable','off');
set(handles.pushbutton_selectall,'Enable','off');

waitforbuttonpress;
point1  =   get(gca,'CurrentPoint');    % button down detected
% point1 is 2x3 matrix, first 2 elements are x,y
% finalRect = rbbox;                   % return figure units
rbbox;
point2  =    get(gca,'CurrentPoint');    % button up detected
point1  =    point1(1,1:2);              % extract x and y
point2  =    point2(1,1:2);
p1      =    min(point1,point2);             % calculate locations
offset  =    abs(point1-point2);         % and dimensions


limX = xlim; limY = ylim;                     % get Axis Limits
% -------------- calculate columns & rows ------------------
leftcolX = fix(( p1(1)-limX(1,1) )/ handles.gridX +1) + 1;
rightcolX = fix(( p1(1)+offset(1)-limX(1,1) )/ handles.gridX )+1;
bottomrowY =  fix(( p1(2)-limY(1,1) )/ handles.gridY+1 )+1;
uprowY = fix(( p1(2) + offset(2) - limY(1,1) )/ handles.gridY)+1;


% -------boundary check ----------------
plotstateX=0;plotstateY=0;
if leftcolX < 1
    leftcolX = 1;
    plotstateY = 1;
end
rightLimit = fix((limX(1,2)-limX(1,1))/handles.gridX)+1;

if rightcolX > rightLimit
    rightcolX = rightLimit; 
end
uprowLimit = fix((limY(1,2)-limY(1,1))/handles.gridY)+1;

if bottomrowY < 1
    bottomrowY = 1;
    plotstateX = 1;   % we need it to plot in right way
end
if uprowY>uprowLimit
    uprowY=uprowLimit; 
end

% --------- selection checking ---------------
sizeI = size(handles.i,1);
sizeJ = size(handles.j,1);
numofcols = rightcolX-leftcolX+1;
% numofrows = uprowY-bottomrowY+1;


% -------------- errorchecking --------
if ~isempty(handles.previousSel)
    a = handles.previousSel;
    if ((rightcolX-leftcolX) == a(2)-a(1) && a(2) == rightcolX && handles.rowlock~=1)
        handles.columnlock=1;
    elseif    ((uprowY-bottomrowY)==a(4)-a(3) && a(4)==uprowY && handles.columnlock~=1)
        handles.rowlock=1;
    else
        errordlg('Your Selection is Invalid...');
        return;
    end;
end;

if ismember([bottomrowY leftcolX],[handles.i handles.j],'rows') ...
        || ismember([bottomrowY rightcolX],[handles.i handles.j],'rows') ...
        || ismember([uprowY leftcolX],[handles.i handles.j],'rows') ...
        || ismember([uprowY rightcolX],[handles.i handles.j],'rows')
    errordlg('Your Selection is Invalid...');
    return;
end




% --------------------- Fill loop ----------------------------

for i1 = bottomrowY:uprowY
    handles.i(sizeI+1:sizeI+numofcols,1)    =   i1;
    handles.j(sizeJ+1:sizeJ+numofcols,1)    =   leftcolX:rightcolX;
    sizeI   =   sizeI+numofcols;
    sizeJ   =   sizeJ+numofcols;
end
% ------------- plot ----------------------

lx_box=limX(1,1)+(leftcolX-1)*handles.gridX*~plotstateY;
rx_box=limX(1,1)+(rightcolX-1)*handles.gridX;
uy_box=limY(1,1)+(uprowY-1)*handles.gridY;
by_box=limY(1,1)+(bottomrowY-1)*handles.gridY*~plotstateX;


x1 = [lx_box rx_box rx_box lx_box lx_box];
y1 = [by_box by_box uy_box uy_box by_box];
hold on
handles.selectionbox = plot(x1,y1,'--b','LineWidth',1.5);
hold off
handles.previousSel=[leftcolX rightcolX bottomrowY uprowY];

guidata(handles.fig,handles);


% --- Executes on button press in pushbutton_selectall.
function pushbutton_selectall_Callback(hObject, ~, handles)
set(handles.pushbutton_selectpoints,'Enable','off');
set(handles.rowpushbutton,'Enable','off');
set(handles.colpushbutton,'Enable','off');
set(handles.pushbutton_selectreg,'Enable','off');
set(handles.pushbutton_selectall,'Enable','off');
handles.Allselected=1;
update_gui(hObject,[],guidata(hObject));
handles.i=[]; handles.j=[]; handles.previousSel=[];
limX=xlim; limY=ylim;
x1 = [limX(1,1) limX(1,2) limX(1,2) limX(1,1) limX(1,1)];
y1 = [limY(1,1) limY(1,1) limY(1,2) limY(1,2) limY(1,1)];
hold on
handles.selectionbox = plot(x1,y1,':b','LineWidth',4);
hold off
leftcolX    =   1;
rightcolX   =   fix((limX(1,2)-limX(1,1))/handles.gridX)+1;
bottomrowY  =   1;
uprowY      =   fix((limY(1,2)-limY(1,1))/handles.gridY)+1;
% uprowY      =   fix((limY(1,2)-limY(1,1))/handles.gridY+1)+1;
numofcols=rightcolX-leftcolX+1;
sizeI=size(handles.i,1);
sizeJ=size(handles.j,1);
for i1 = bottomrowY:uprowY
    handles.i(sizeI+1:sizeI+numofcols,1)=i1;
    handles.j(sizeJ+1:sizeJ+numofcols,1)=leftcolX:rightcolX;
    sizeI=sizeI+numofcols; sizeJ=sizeJ+numofcols;
end
guidata(handles.fig,handles);


% --- Executes on button press in pushbutton_time.
function pushbutton_time_Callback(~, ~, handles)
if ~isempty(handles.i)
    timebox(handles); % timebox includes both versions with if ... else
    guidata(handles.fig,handles);
else
    errordlg('Select region of interest');
end



% --- Executes on button press in pushbutton_profile1.
function pushbutton_profile1_Callback(~, ~, handles)
if isempty(handles.property) || isempty(handles.i)
    errordlg('First, pick the quantity and region of interest !!!');
else
    distrib(handles);  % call to spatialbox
    handles.distribOn = 1;
end;

guidata(handles.fig,handles);


% --- Executes on button press in pushbutton_reset.
function pushbutton_reset_Callback(~, ~, handles)

handles.i = []; handles.j = [];
handles.rowlock=0; handles.columnlock=0;
handles.previousSel = [];
set(handles.pushbutton_selectpoints,'Enable','on');
set(handles.pushbutton_selectreg,'Enable','on');
set(handles.colpushbutton,'Enable','on');
set(handles.pushbutton_selectall,'Enable','on');
set(handles.rowpushbutton,'Enable','on');
guidata(handles.fig,handles);
update_gui(gcbo,[],guidata(gcbo));


% --- Executes during object creation, after setting all properties.
function figure_gradpiv_CreateFcn(hObject, ~, ~)

load cil_logo
image(im,'Parent',findobj(hObject,'type','axes')); %handles.axes_main);
axis off


% --- Executes during object creation, after setting all properties.
function axes_main_CreateFcn(hObject, ~, handles)

handles.axes_main = hObject;
guidata(hObject, handles);




function hh = quiverc(varargin)
% Modified QUIVER() that uses Patches instead of lines and allows
% to color arrows and have updated colorbar, according to the
% correct color data mapping. Arrows are of the same type as in quiver plot.
%
%   Modified version: (c) Alex Liberzon
%   $Revision: 1.02 $  $ Date: 05/10/2002 $
%   See also: HELP QUIVER (skipped here)


% Arrow head parameters
alpha = 0.33; % Size of arrow head relative to the length of the vector
beta = 0.33;  % Width of the base of the arrow head relative to the length
autoscale = 1; % Autoscale if ~= 0 then scale by this.
plotarrows = 1; % Plot arrows
% sym = '';

filled = 0;
% ls = '-';
ms = '';
col = '';

nin = nargin;
% Parse the string inputs
while ischar(varargin{nin}),
    vv = varargin{nin};
    if ~isempty(vv) && strcmpi(vv(1),'f')
        filled = 1;
        nin = nin-1;
    else
        [~,c,m,msg] = colstyle(vv);
        if ~isempty(msg),
            error('Unknown option "%s".',vv);
        end
        % if ~isempty(l), ls = l; end
        if ~isempty(c), col = c; end
        if ~isempty(m), ms = m; plotarrows = 0; end
        if isequal(m,'.'), ms = ''; end % Don't plot '.'
        nin = nin-1;
    end
end

error(nargchk(2,6,nin));

% Check numeric input arguments
if nin<4, % quiver(u,v) or quiver(u,v,s)
    [msg,x,y,u,v] = xyzchk(varargin{1:2});
else
    [msg,x,y,u,v] = xyzchk(varargin{1:4});
end
if ~isempty(msg), error(msg); end

if nin==4 || nin==6, % quiver(u,v,z,s) or quiver(x,y,u,v,z,s)
    autoscale = varargin{nin-1};
    z = varargin{nin};
end

% Scalar expand u,v
if numel(u)==1, u = u(ones(size(x))); end
if numel(v)==1, v = v(ones(size(u))); end

if autoscale,
    if min(size(x))==1
        n = sqrt(numel(x));
        m = n;
    else
        [m,n]=size(x);
    end
    delx = diff([min(x(:)) max(x(:))])/n;
    dely = diff([min(y(:)) max(y(:))])/m;
    del = delx.^2 + dely.^2;
    if del>0
        len = sqrt((u.^2 + v.^2)/del);
        maxlen = max(len(:));
    else
        maxlen = 0;
    end
    
    if maxlen>0
        autoscale = autoscale*0.9 / maxlen;
    else
        autoscale = autoscale*0.9;
    end
    u = u*autoscale; v = v*autoscale;
end


ax = newplot;
next = lower(get(ax,'NextPlot'));
hold_state = ishold;

% Make velocity vectors
x = x(:).'; y = y(:).';
u = u(:).'; v = v(:).';
uu = [x; x+u; NaN(size(u))];
vv = [y; y+v; NaN(size(u))];

% Prepare color matrix
z = [z(:)';z(:)';NaN*z(:)'];

h1 = patch([uu(:),uu(:)],[vv(:),vv(:)], [z(:),z(:)],'Parent',ax,'EdgeColor','Flat','FaceColor','None');

if plotarrows,
    % Make arrow heads and plot them
    hu = [x+u-alpha*(u+beta*(v+eps));x+u; ...
        x+u-alpha*(u-beta*(v+eps)); NaN(size(u))];
    hv = [y+v-alpha*(v-beta*(u+eps)); y+v; ...
        y+v-alpha*(v+beta*(u+eps)); NaN(size(v))];
    hold on
    % Modify color matrix
    z = [z(1,:); z];
    h2 = patch([hu(:),hu(:)],[hv(:),hv(:)], [z(:),z(:)],'Parent',ax,'EdgeColor','Flat','FaceColor','None');
    
else
    h2 = [];
end

if ~isempty(ms), % Plot marker on base
    hu = x; hv = y;
    hold on
    h3 = plot(hu(:),hv(:),[col ms]);
    if filled, set(h3,'markerfacecolor',get(h1,'color')); end
else
    h3 = [];
end

if ~hold_state, hold off, view(2); set(ax,'NextPlot',next); end

if nargout > 0, hh = [h1;h2;h3]; end


% --- Executes on button press in rowpushbutton.
function rowpushbutton_Callback(~, ~, handles)

set(handles.pushbutton_selectpoints,'Enable','off');
set(handles.pushbutton_selectreg,'Enable','off');
set(handles.colpushbutton,'Enable','off');
set(handles.pushbutton_selectall,'Enable','off');

set(handles.axes_main,'NextPlot','Add');
limX = xlim;
limY = ylim;
% leftcolX   = 1;
bottomrowY = 1;


while 1
    
    [x1,y1, buttonNumber] = ginput(1);
    
    % When the right button is pressed, stop the loop
    if (buttonNumber == 2) || (buttonNumber==3)
        break;
    end;
    col = fix (( x1 - limX(1,1) )/ handles.gridX+0.5  )+1;
    row = fix(( y1 - limY(1,1) )/ handles.gridY+0.5  )+1;
    % check for errors ----------------
    if col<1 || col>(fix((limX(1,2)-limX(1,1))/handles.gridX)+1) || row<1 || row...
            >(fix((limY(1,2)-limY(1,1))/handles.gridY)+1);
        guidata(handles.fig,handles);
        return;
    end;
    % ---------------------------------
    sizeI = size(handles.i,1);
    % rightcolX = fix(( limX(1,2)-limX(1,1) )/  handles.gridX )+1;
    uprowY    = fix(( limY(1,2)-limY(1,1) )/  handles.gridY )+1;
    % uprowY    = fix(( limY(1,2)-limY(1,1) )/  handles.gridY +1)+1; % Alex,
    % 26.10
    sizeJ = size(handles.j,1);
    % numofcols = rightcolX - leftcolX + 1;
    numofrows = uprowY - bottomrowY + 1;
    
    handles.i(sizeI+1:sizeI+numofrows,1) = 1:uprowY ;
    handles.j(sizeJ+1:sizeJ+numofrows,1) = col;
    row = 1:uprowY;
    % topLeft(1) = uprowY; 
    % bottomRight(1)=1;
    
    line(limX(1,1)+(col-1)*...
        handles.gridX,limY(1,1)+(row-1)*handles.gridY,...
        'Marker','o','Color','k','MarkerSize',8);
    
    
end;
guidata(handles.fig,handles);

% --- Executes on button press in colpushbutton.
function colpushbutton_Callback(~, ~, handles)
% hObject    handle to colpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbutton_selectpoints,'Enable','off');
set(handles.pushbutton_selectreg,'Enable','off');
set(handles.rowpushbutton,'Enable','off');
set(handles.pushbutton_selectall,'Enable','off');

set(handles.axes_main,'NextPlot','Add');
limX = xlim;
limY = ylim;
leftcolX   = 1;
% bottomrowY = 1;

while 1
    
    [x1,y1, buttonNumber] = ginput(1);
    
    % When the right button is pressed, stop the loop
    if (buttonNumber == 2) || (buttonNumber==3)
        break;
    end;
    col = fix (( x1 - limX(1,1) )/ handles.gridX+0.5  )+1;
    row = fix(( y1 - limY(1,1) )/ handles.gridY+0.5  )+1;
    % -------- find the corners of rectangle ----
    
    
    % check for errors ----------------
    if col<1 || col>(fix((limX(1,2)-limX(1,1))/handles.gridX)+1) ...
            || row < 1 ...
            || row >(fix((limY(1,2)-limY(1,1))/handles.gridY)+1);
        guidata(handles.fig,handles);
        return
    end
    % ---------------------------------
    sizeI = size(handles.i,1);
    rightcolX = fix(( limX(1,2)-limX(1,1) )/  handles.gridX )+1;
    % uprowY    = fix(( limY(1,2)-limY(1,1) ) /  handles.gridY )+1;
    sizeJ = size(handles.j,1);
    numofcols = rightcolX - leftcolX + 1;
    % numofrows = uprowY - bottomrowY + 1;
    %
    handles.i(sizeI+1:sizeI+numofcols,1) = row;
    handles.j(sizeJ+1:sizeJ+numofcols,1) = 1:rightcolX;
    col = 1:rightcolX;
    %
    line(limX(1,1)+(col-1)*...
        handles.gridX,limY(1,1)+(row-1)*handles.gridY,...
        'Marker','o','Color','k','MarkerSize',8);
end;


guidata(handles.fig,handles);


% --- Executes during object creation, after setting all properties.
function ed_max_CreateFcn(hObject, ~, ~)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor', ...
        get(0,'defaultUicontrolBackgroundColor'));
end



% --------------------------------------------------------------------
function export2figure_Callback(~, ~, handles)
handles.export_figure = figure;
copyobj(handles.axes_main,handles.export_figure);

set(handles.export_figure,'Units','normalized');
set(get(handles.export_figure,'children'),'Units','normalized');
set(get(handles.export_figure,'children'),'Position',[0.13 0.11 0.775 0.815]);
set(get(handles.export_figure,'children'),'Box','on');

if isfield(handles,'color_flag') && handles.colorbar_flag
    colorbar;
    % mcolorbar; % (get(handles.export_figure,'Children'));
end
guidata(handles.fig, handles);


% --- Executes on button press in pushbutton_stats.
function pushbutton_stats_Callback(~, ~, ~)
% hObject    handle to pushbutton_stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function loadMat_Callback(~, ~, handles)
% hObject    handle to loadMat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global orighandles;
if isfield(handles,'restoreorig')
    handles = orighandles;
end
handles.restoreorig = 1;
handles.dt = 1;
handles.state3d = 0;


try
    % Check the contents of the MAT file, if it's coordinates only or full file
    % curdir = cd;
    % [coordMatfile,coordMatpath] = uigetfile('*.mat','Choose Coordinates or EXPORTED MAT file');
    coordMatfile = uipickfiles('FilterSpec','*.mat');
    w = who('-file',coordMatfile{1});
    
    
    exportedMat = false;
    if sum(cellfun(@sum,strfind(w,'xUnits'))) > 0, exportedMat = true; end
    
    
    
    switch exportedMat
        case true
            % load "ready" dataset, no need for double loading of
            % coordinates and velocities
            tmp = load(coordMatfile{1});
            handles.x = tmp.x;
            handles.y = tmp.y;
            handles.u = tmp.u;
            handles.v = tmp.v;
            handles.uf = tmp.uf;
            handles.vf = tmp.vf;
            handles.files = tmp.files;
            handles.path = tmp.path;
            handles.dx = tmp.dx;
            handles.dy = tmp.dy;
            handles.dudx = tmp.dudx;
            handles.dvdx = tmp.dvdx;
            handles.dudy = tmp.dudy;
            handles.dvdy = tmp.dvdy;
            handles.gridX = tmp.gridX;
            handles.gridY = tmp.gridY;
            handles.N = tmp.N;
            handles.xUnits = tmp.xUnits;
            handles.velUnits = tmp.velUnits;
            clear tmp
            
        case false
            load(coordMatfile{1});
            
            
            if exist('x','var')
                handles.x = x;
                handles.y = y; % max(y(:)) - y;
                clear x y
            elseif exist('X','var')
                handles.x = X;
                handles.y = Y; % max(Y(:)) - Y;
                clear X Y
            else
                errordlg('Coordinates file does not include x or X variables');
                set(handles.fig,'pointer','arrow');
            end
            coordMatfile = uipickfiles('Type',{'*.mat','Velocity MAT file'});
            load(coordMatfile{1});
            if exist('fu','var')
                handles.u = fu;
                handles.v = fv;
                clear fu fv
            elseif exist('U','var')
                handles.u = U;
                handles.v = V;
                clear U V
            elseif exist('u','var')
                handles.u = u;
                handles.v = v;
                clear u v
            else
                errordlg('Wrong velocity names');
            end
            
            
            [rows,cols,handles.N] = size(handles.u);
            
            handles.xUnits = 'pix';
            handles.velUnits = 'pix/s';
        otherwise
            
    end
catch ME
    
    errordlg(ME.message); %'Something wrong with MAT files'); % vector -> MAT
    set(handles.fig,'pointer','arrow');
    return
end

try
    set(handles.data_info,'String',coordMatfile); % july 27, Alex's Laptop, SVN assembla
catch ME
    % setfield(handles,'data_info',coordMatfile); % 18.12.10, dropbox
    disp(ME.message);
    handles.('data_info') = coordMatfile; % Aug. 2013
end
handles.current = 1;                      % current file beeing displayed
% Display first file number, total number of files
set(handles.edit_current,'String',handles.current);
set(handles.edit_numfields,'String',handles.N);

% Initialize color, number of colors for contours
handles.numcolors = 10;                   % default number of colors
set(handles.edit_numcolors,'String', handles.numcolors);

% Cat the mean values at the end
if ~exportedMat
    handles.u(:,:,handles.N+1) = mean(handles.u(:,:,1:handles.N),3);
    handles.v(:,:,handles.N+1) = mean(handles.v(:,:,1:handles.N),3);
    if handles.state3d
        handles.w(:,:,handles.N+1) = mean(handles.w(:,:,1:handles.N),3);
        handles.wf = zeros(rows,cols,handles.N);
        for i = 1:handles.N
            handles.wf(:,:,i) = handles.w(:,:,i) - handles.w(:,:,handles.N+1);
        end
    end
    
    % handles.dx = abs(handles.x(1,1) - handles.x(1,2));
    % handles.dy = abs(handles.y(1,1) - handles.y(2,1));
    % Bug, fixed at 16.06 after the email of Hai
    handles.dx = handles.x(1,2) - handles.x(1,1);
    handles.dy = handles.y(2,1) - handles.y(1,1);
    
    % Preallocate memory and calculate all the necessary quantities
    % Fluctuations (last one is zero, we do not need it)
    handles.uf = zeros(rows,cols,handles.N);
    handles.vf = zeros(rows,cols,handles.N);
    %
    for i = 1:handles.N
        handles.vf(:,:,i) = handles.v(:,:,i) - handles.v(:,:,handles.N+1);
    end
    %
    for i = 1:handles.N
        handles.uf(:,:,i) = handles.u(:,:,i) - handles.u(:,:,handles.N+1);
    end
    
    % Derivatives
    handles.dudx = zeros(rows,cols,handles.N+1);
    handles.dudy = zeros(rows,cols,handles.N+1);
    handles.dvdx = zeros(rows,cols,handles.N+1);
    handles.dvdy = zeros(rows,cols,handles.N+1);
    %
    for i = 1:handles.N+1
        [handles.dudx(:,:,i),handles.dudy(:,:,i)] = ...
            lsgradient(handles.u(:,:,i),handles.dx, handles.dy);
        [handles.dvdx(:,:,i),handles.dvdy(:,:,i)] = ...
            lsgradient(handles.v(:,:,i),handles.dx, handles.dy);
    end
    
    % Possible future development, eliminating strong gradients
    % on the borders
    % handles.dudx(1:2,:,:)       = NaN;
    % handles.dudx(end-1:end,:,:) = NaN;
    % handles.dudx(:,1:2,:)       = NaN;
    % handles.dudx(:,end-1:end,:) = NaN;
    
end


%
% More defaults
handles.arrow_scale = 1;                % default scale
set(handles.edit_arrow_size,'String',handles.arrow_scale);

% Default situation, instantaneous, not Ensemble, not fluctuations
set(handles.checkbox_ensemble,'Value',0);
set(handles.checkbox_fluct,'Value',0);

% No colors, no labels
handles.color = 0;                % display colored / black figure
handles.alltodisp = 0;            % default all_to display is unset
handles.allfields = 0;
handles.labelit = 0;              % label for contour
handles.colorbar_flag = 0;
handles.current_index = handles.current;
handles.distribOn=0;
handles.rowlock=0; handles.columnlock=0;
handles.previousSel=[];

% These are future values, for the TimeBox
handles.i=[];
handles.j=[];
handles.PointsH=[];
handles.Allselected=0;

% handles.dx = handles.x(1,2) - handles.x(1,1);
% handles.dy = handles.y(2,1) - handles.y(1,1);

% abs() is added on Dec.13,2004. by Alex on IHW_Liberzon computer
%
handles.gridX = abs(handles.x(1,2) - handles.x(1,1));
handles.gridY = abs(handles.y(1,1) - handles.y(2,1));

% Show some and hide some controls
set(handles.popupmenu_quantity,'Visible','on');
set(handles.popupmenu_quantity,'Value',1);
set(handles.popupmenu_contour_type,'Visible','on');
set(handles.popupmenu_contour_type,'Value',1);
set(handles.popupmenu_eachfield,'Visible','on');
set(handles.popupmenu_eachfield,'Value',1);

% Some only part of them in Enable mode
set(handles.checkbox_arrow,'Enable','On');
set(handles.edit_arrow_size,'Enable','On');
set(handles.checkbox_fluct,'Enable','On');
set(handles.checkbox_ensemble,'Enable','On');
set(handles.popupmenu_quantity,'Enable','On');
set(handles.pushbutton_previous,'Enable','On');
set(handles.pushbutton_next,'Enable','On');
set(handles.edit_current,'Enable','On');
set(handles.pushbutton_animate,'Enable','On');
set(handles.pushbutton_save_movie,'Enable','on');
% set(handles.pushbutton_stats,'Enable','on');


% Properties
set(handles.popupmenu_quantity,'String',handles.inst_list);
handles.property = [];

% ---------------- store all handles of spatial controls ------------------
% handles.spatial_controls=[handles.checkbox_ensemble,handles.checkbox_fluct,handles.checkbox_arrow,handles.checkbox_arrow_color,...
%         handles.checkbox_label,handles.checkbox_colorbar,handles.edit_arrow_size,handles.edit_numcolors,...
%         handles.edit_current,handles.edit_numfields,...
%         handles.text_contour_quantity, handles.text_contourtype, handles.text_numberofcolors, handles.text7,handles.text2,...
%         handles.text_arrow_size, handles.pushbutton_previous, handles.pushbutton_next, ...
%         handles.pushbutton_animate, handles.pushbutton_save_movie, ...
%         handles.frame_controls,handles.frame8,handles.frame_contour_quantity,handles.frame_contour_type,handles.frame7,...
%         handles.frame_arrow, handles.text5, handles.popupmenu_quantity,handles.popupmenu_contour_type,...
%         handles.popupmenu_eachfield];% ,handles.pushbutton_stats];

% -------------- store all handles of select controls ----------------
% handles.select_controls=[handles.pushbutton_selectpoints,handles.pushbutton_selectreg,handles.pushbutton_selectall,...
%          handles.pushbutton_profile1,handles.pushbutton_time,...
%          handles.pushbutton_reset, handles.rowpushbutton, handles.colpushbutton];

handles.arrow_ctrls=[handles.edit_arrow_size,handles.checkbox_arrow,handles.checkbox_arrow_color];
% handles.stat_controls=[handles.ed_mean,handles.ed_std,handles.ed_min,handles.ed_max,handles.ed_text,...
%         handles.ed_pushsavestl,handles.ed_frame,handles.ed_pushclose,handles.ed_textmean,...
%         handles.ed_textstd,handles.ed_textmax,handles.ed_textmin];


set(handles.pushbutton_spatial,'Enable','on');
set(handles.pushbutton_select,'Enable','on');

set(handles.fig,'pointer','arrow');

% added on 10.04.06 for R12SP3 version
handles.axpos = get(handles.axes_main,'Position');

% Update all handles structure
guidata(handles.fig,handles);

% Make default plot
update_gui(handles.fig,[],handles);
return
% --------------------------------------------------------------------


% --------------------------------------------------------------------
function loadTXT_Callback(~, ~, handles)
% hObject    handle to loadTXT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global orighandles;
if isfield(handles,'restoreorig')
    handles = orighandles;
end

handles.restoreorig = 1;
handles.state3d = 0;
handles.dt  = 1;


try
    % [gui_files,gui_path] = uigetfile('*.txt','Select the location and the type of the TXT file');
    handles.files = uipickfiles; % April 2, 2010, AL
    
    
    handles.N = length(handles.files); % number of files selected
    if  handles.N > 0
        [handles.path,~,~] = fileparts(handles.files{1});
        set(handles.fig,'pointer','watch');
    else
        return
    end
    
    %%%%%%%%%%%%%%%%
    
    
    d = load(handles.files{1});
    
    
    x = d(:,1);
    x = x(x~=0);
    unX = unique(x);
    
    minX = min(unX);
    maxX = max(unX);
    dX = ceil((maxX-minX)/(length(unX)-1));
    
    y = d(:,2);
    y = y(y~=0);
    unY = unique(y);
    
    minY = min(unY);
    maxY = max(unY);
    dY = ceil((maxY-minY)/(length(unY)-1));
    
    [handles.x,handles.y] = meshgrid(minX:dX:maxX,minY:dY:maxY);
    [rows,cols] = size(handles.x);
    
    [handles.u,handles.v] = deal(zeros(rows,cols,handles.N+1)); % 11.04.04, Alex
    
    hwaitbar = waitbar(0,'Please wait...');
    
    
    
    tmp = d;
    tmp(tmp(:,1) == 0) = [];
    tmp(tmp(:,2) == 0) = [];
    % y = tmp(:,2);
    %  x = tmp(:,1);
    [m,n] = zeros(length(tmp(:,1)));
    for j = 1:length(tmp(:,1))
        [m(j),n(j)] = find(handles.x == tmp(j,1) & handles.y == tmp(j,2));
    end
    
    
    for i = 1:handles.N
        waitbar(i/handles.N,hwaitbar)
        % x = d(:,1,i);
        % tmp = d(x~=0,:,i);
        % presumably cleans bug values with zero in x and y - not clear
        % yet, 23.3.11, Alex
        
        tmp = load(handles.files{i});
        tmp(tmp(:,1) == 0) = [];
        tmp(tmp(:,2) == 0) = [];
        
        for j = 1:length(tmp(:,1))
            handles.u(m(j),n(j),i) = tmp(j,3);
            handles.v(m(j),n(j),i) = tmp(j,4);
        end
    end
    
    close(hwaitbar)
    handles.xUnits = 'pix';
    handles.velUnits = 'pix/dt';
    
    clear d tmp x y
    
catch ME
    errordlg(ME.message); %'Something wrong with TXT files');
    set(handles.fig,'pointer','arrow');
    return
end

handles.current = 1;                      % current file beeing displayed
% Display first file number, total number of files
set(handles.edit_current,'String',handles.current);
set(handles.edit_numfields,'String',handles.N);

% Initialize color, number of colors for contours
handles.numcolors = 10;                   % default number of colors
set(handles.edit_numcolors,'String', handles.numcolors);

% Cat the mean values at the end
handles.u(:,:,handles.N+1) = mean(handles.u(:,:,1:handles.N),3);
handles.v(:,:,handles.N+1) = mean(handles.v(:,:,1:handles.N),3);
if handles.state3d
    handles.w(:,:,handles.N+1) = mean(handles.w(:,:,1:handles.N),3);
    handles.wf = zeros(rows,cols,handles.N);
    for i = 1:handles.N
        handles.wf(:,:,i) = handles.w(:,:,i) - handles.w(:,:,handles.N+1);
    end
end

% handles.dx = abs(handles.x(1,1) - handles.x(1,2));
% handles.dy = abs(handles.y(1,1) - handles.y(2,1));
% Bug, fixed at 16.06 after the email of Hai
handles.dx = handles.x(1,2) - handles.x(1,1);
handles.dy = handles.y(2,1) - handles.y(1,1);

% Preallocate memory and calculate all the necessary quantities
% Fluctuations (last one is zero, we do not need it)
handles.uf = zeros(rows,cols,handles.N);
handles.vf = zeros(rows,cols,handles.N);
%
for i = 1:handles.N
    handles.vf(:,:,i) = handles.v(:,:,i) - handles.v(:,:,handles.N+1);
end
%
for i = 1:handles.N
    handles.uf(:,:,i) = handles.u(:,:,i) - handles.u(:,:,handles.N+1);
end

% Derivatives
handles.dudx = zeros(rows,cols,handles.N+1);
handles.dudy = zeros(rows,cols,handles.N+1);
handles.dvdx = zeros(rows,cols,handles.N+1);
handles.dvdy = zeros(rows,cols,handles.N+1);
%
for i = 1:handles.N+1
    [handles.dudx(:,:,i),handles.dudy(:,:,i)] = ...
        lsgradient(handles.u(:,:,i),handles.dx, handles.dy);
    [handles.dvdx(:,:,i),handles.dvdy(:,:,i)] = ...
        lsgradient(handles.v(:,:,i),handles.dx, handles.dy);
end

% Possible future development, eliminating strong gradients
% on the borders
% handles.dudx(1:2,:,:)       = NaN;
% handles.dudx(end-1:end,:,:) = NaN;
% handles.dudx(:,1:2,:)       = NaN;
% handles.dudx(:,end-1:end,:) = NaN;

%
% More defaults
handles.arrow_scale = 1;                % default scale
set(handles.edit_arrow_size,'String',handles.arrow_scale);

% Default situation, instantaneous, not Ensemble, not fluctuations
set(handles.checkbox_ensemble,'Value',0);
set(handles.checkbox_fluct,'Value',0);

% No colors, no labels
handles.color = 0;                % display colored / black figure
handles.alltodisp = 0;            % default all_to display is unset
handles.allfields = 0;
handles.labelit = 0;              % label for contour
handles.colorbar_flag = 0;
handles.current_index = handles.current;
handles.distribOn=0;
handles.rowlock=0; handles.columnlock=0;
handles.previousSel=[];

% These are future values, for the TimeBox
handles.i=[];
handles.j=[];
handles.PointsH=[];
handles.Allselected=0;

% handles.dx = handles.x(1,2) - handles.x(1,1);
% handles.dy = handles.y(2,1) - handles.y(1,1);

% abs() is added on Dec.13,2004. by Alex on IHW_Liberzon computer
%
handles.gridX = abs(handles.x(1,2) - handles.x(1,1));
handles.gridY = abs(handles.y(1,1) - handles.y(2,1));

% Show some and hide some controls
set(handles.popupmenu_quantity,'Visible','on');
set(handles.popupmenu_quantity,'Value',1);
set(handles.popupmenu_contour_type,'Visible','on');
set(handles.popupmenu_contour_type,'Value',1);
set(handles.popupmenu_eachfield,'Visible','on');
set(handles.popupmenu_eachfield,'Value',1);

% Some only part of them in Enable mode
set(handles.checkbox_arrow,'Enable','On');
set(handles.edit_arrow_size,'Enable','On');
set(handles.checkbox_fluct,'Enable','On');
set(handles.checkbox_ensemble,'Enable','On');
set(handles.popupmenu_quantity,'Enable','On');
set(handles.pushbutton_previous,'Enable','On');
set(handles.pushbutton_next,'Enable','On');
set(handles.edit_current,'Enable','On');
set(handles.pushbutton_animate,'Enable','On');
set(handles.pushbutton_save_movie,'Enable','on');
% set(handles.pushbutton_stats,'Enable','on');


% Properties
set(handles.popupmenu_quantity,'String',handles.inst_list);
handles.property = [];

handles.arrow_ctrls=[handles.edit_arrow_size,handles.checkbox_arrow,handles.checkbox_arrow_color];

set(handles.pushbutton_spatial,'Enable','on');
set(handles.pushbutton_select,'Enable','on');

set(handles.fig,'pointer','arrow');

% added on 10.04.06 for R12SP3 version
handles.axpos = get(handles.axes_main,'Position');

% Update all handles structure
guidata(handles.fig,handles);

% Make default plot
update_gui(handles.fig,[],handles);
return

% ---------------------------------------------------------------------
% function [filenames] = ReadTXTDir(dirname,data)
% if nargin < 2
%     data = 'txt';
% end
% 
% switch data
%     case{'_noflt.txt'} % a)
%         direc = dir([dirname,filesep,'*_noflt.txt']);
%     case{'_flt.txt'} %b)
%         direc = dir([dirname,filesep,'*_flt.txt']);
%     case{'txt'} % c)
%         direc = dir([dirname,filesep,'*.txt']);
%         tmp = struct('name',[]);
%         k = 0;
%         for i=1:length(direc)
%             if length(findstr(direc(i).name,'flt')) < 1
%                 k = k + 1;
%                 tmp(k).name = direc(i).name;
%             end
%         end
%         direc = tmp;
% end
% 
% if ~isempty(direc(1).name) && ~isempty(str2num(direc(1).name(1:length(direc(1).name)-4)))
%     for i = 1:length(direc)
%         n(i) = str2num(direc(i).name(1:length(direc(i).name)-4));
%     end
%     [junk,j] = sort(n);
%     direc = direc(j);
% end
% 
% filenames={};
% [filenames{1:length(direc),1}] = deal(direc.name);
% % filenames = sortrows(filenames);
% return


% --------------------------------------------------------------------
function Data_Callback(~, ~, ~)
% hObject    handle to Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Filter_Callback(~, ~, ~)
% hObject    handle to Filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Global_Callback(~, ~, handles)
% hObject    handle to Global (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%%%%%%%%%%%%%%%%%%%
% Reshape U and V matrices in two-dimensional grid and produce
% velocity vector in U + i*V form (real and imaginary parts):



[rows,cols] = size(handles.u(:,:,1));
%     % Remove outlayers - GLOBAL FILTERING
for i = 1:handles.N
    vector = sqrt(handles.u(:,:,i).^2 + handles.v(:,:,i).^2);
    % index = find(vector > (median(vector) + 6*std(vector)));
    [~,index,~] = deleteoutliers(vector(:),0.05);
    [k,m] = ind2sub([rows,cols],index);
    handles.u(k,m,i) = 0;
    handles.v(k,m,i) = 0;
end

% Update data
handles.u(:,:,handles.N+1) = mean(handles.u(:,:,1:handles.N),3);
handles.v(:,:,handles.N+1) = mean(handles.v(:,:,1:handles.N),3);
if handles.state3d
    handles.w(:,:,handles.N+1) = mean(handles.w(:,:,1:handles.N),3);
    handles.wf = zeros(rows,cols,handles.N);
    for i = 1:handles.N
        handles.wf(:,:,i) = handles.w(:,:,i) - handles.w(:,:,handles.N+1);
    end
end

handles.dx = handles.x(1,2) - handles.x(1,1);
handles.dy = handles.y(2,1) - handles.y(1,1);

%
for i = 1:handles.N
    handles.vf(:,:,i) = handles.v(:,:,i) - handles.v(:,:,handles.N+1);
    handles.uf(:,:,i) = handles.u(:,:,i) - handles.u(:,:,handles.N+1);
end
%
for i = 1:handles.N+1
    [handles.dudx(:,:,i),handles.dudy(:,:,i)] = ...
        lsgradient(handles.u(:,:,i),handles.dx, handles.dy);
    [handles.dvdx(:,:,i),handles.dvdy(:,:,i)] = ...
        lsgradient(handles.v(:,:,i),handles.dx, handles.dy);
end

update_gui(handles.fig,[],handles);
guidata(handles.fig,handles);


update_gui(handles.fig,[],handles);
guidata(handles.fig,handles)

% --------------------------------------------------------------------
function Median_Callback(~, ~, handles)
% hObject    handle to Median (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Adaptive Local Median filtering

kernel = [-1 -1 -1; -1 8 -1; -1 -1 -1];
[rows,cols] = size(handles.u(:,:,1));

for i = 1:handles.N
    tmpv = abs(conv2(handles.v(:,:,i),kernel,'same'));
    tmpu = abs(conv2(handles.u(:,:,i),kernel,'same'));
    
    ind = tmpv~=0 & tmpu~=0;
    lmtv = mean(tmpv(ind)) + 3*std(tmpv(ind));
    lmtu = mean(tmpu(ind)) + 3*std(tmpu(ind));
    
    out = find(tmpu > lmtu  & tmpv > lmtv);
    
    [k,m] = ind2sub([rows,cols],out);
    handles.u(k,m,i) = 0;
    handles.v(k,m,i) = 0;
end

% Update data

handles.u(:,:,handles.N+1) = mean(handles.u(:,:,1:handles.N),3);
handles.v(:,:,handles.N+1) = mean(handles.v(:,:,1:handles.N),3);
if handles.state3d
    handles.w(:,:,handles.N+1) = mean(handles.w(:,:,1:handles.N),3);
    handles.wf = zeros(rows,cols,handles.N);
    for i = 1:handles.N
        handles.wf(:,:,i) = handles.w(:,:,i) - handles.w(:,:,handles.N+1);
    end
end

handles.dx = handles.x(1,2) - handles.x(1,1);
handles.dy = handles.y(2,1) - handles.y(1,1);

%
for i = 1:handles.N
    handles.vf(:,:,i) = handles.v(:,:,i) - handles.v(:,:,handles.N+1);
    handles.uf(:,:,i) = handles.u(:,:,i) - handles.u(:,:,handles.N+1);
end
%
for i = 1:handles.N+1
    [handles.dudx(:,:,i),handles.dudy(:,:,i)] = ...
        lsgradient(handles.u(:,:,i),handles.dx, handles.dy);
    [handles.dvdx(:,:,i),handles.dvdy(:,:,i)] = ...
        lsgradient(handles.v(:,:,i),handles.dx, handles.dy);
end

update_gui(handles.fig,[],handles);
guidata(handles.fig,handles);


update_gui(handles.fig,[],handles);
guidata(handles.fig,handles);


% --------------------------------------------------------------------
function Crop_Callback(~, ~, handles)
% hObject    handle to Crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% axes(handles.axes_main);
handles.rect = getrect(handles.fig);
rectangle('position',handles.rect,'Curvature',[0 0],'Edgecolor','blue','Linestyle',':');
% handles.xind = find(handles.x(1,:) > handles.rect(1) & handles.x(1,:) < (handles.rect(1) + handles.rect(3)));
% handles.yind = find(handles.y(:,1) > handles.rect(2) & handles.y(:,1) < (handles.rect(2) + handles.rect(4)));
handles.xind = find(handles.x(1,:) < handles.rect(1) | handles.x(1,:) > (handles.rect(1) + handles.rect(3)));
handles.yind = find(handles.y(:,1) < handles.rect(2) | handles.y(:,1) > (handles.rect(2) + handles.rect(4)));

% w = whos('handles');
% if w.bytes > 60000000
%
%     errordlg('Matlab does not know how to release memory, call Mathworks');
%
%
% else
%     handles.x = handles.x(handles.yind,handles.xind);
%     handles.y = handles.y(handles.yind,handles.xind);
%     handles.u = handles.u(handles.yind,handles.xind,:);
%     handles.v = handles.v(handles.yind,handles.xind,:);
%     handles.dudx = handles.dudx(handles.yind,handles.xind,:);
%     handles.dudy = handles.dudy(handles.yind,handles.xind,:);
%     handles.dvdx = handles.dvdx(handles.yind,handles.xind,:);
%     handles.dvdy = handles.dvdy(handles.yind,handles.xind,:);
%
%     handles.uf = handles.uf(handles.yind,handles.xind,:);
%     handles.vf = handles.vf(handles.yind,handles.xind,:);
%
% end


try
    handles.x(handles.yind,:) = [];
    handles.y(handles.yind,:)= [];
    handles.u(handles.yind,:,:)= [];
    handles.v(handles.yind,:,:)= [];
    handles.dudx(handles.yind,:,:)= [];
    handles.dudy(handles.yind,:,:)= [];
    handles.dvdx(handles.yind,:,:)= [];
    handles.dvdy(handles.yind,:,:)= [];
    
    handles.uf(handles.yind,:,:)= [];
    handles.vf(handles.yind,:,:)= [];
    
    
    handles.x(:,handles.xind) = [];
    handles.y(:,handles.xind)= [];
    handles.u(:,handles.xind,:)= [];
    handles.v(:,handles.xind,:)= [];
    handles.dudx(:,handles.xind,:)= [];
    handles.dudy(:,handles.xind,:)= [];
    handles.dvdx(:,handles.xind,:)= [];
    handles.dvdy(:,handles.xind,:)= [];
    
    handles.uf(:,handles.xind,:)= [];
    handles.vf(:,handles.xind,:)= [];
    
    if ~isempty(handles.property)
        handles.property(handles.yind,:,:) = [];
        handles.property(:,handles.xind,:) = [];
    end
    
    if ~isempty(handles.uf2)
        handles.uf2(handles.yind,:,:) = [];
        handles.uf2(:,handles.xind,:) = [];
        handles.vf2(handles.yind,:,:) = [];
        handles.vf2(:,handles.xind,:) = [];
        
        handles.rs(handles.yind,:,:) = [];
        handles.rs(:,handles.xind,:) = [];
        handles.Tu(handles.yind,:,:) = [];
        handles.Tu(:,handles.xind,:) = [];
        handles.diss(handles.yind,:,:) = [];
        handles.diss(:,handles.xind,:) = [];
        handles.prod(handles.yind,:,:) = [];
        handles.prod(:,handles.xind,:) = [];
    end
catch ME
    errordlg(ME.message); % 'Something went wrong with crop');
end

update_gui(handles.fig,[],handles);
guidata(handles.fig,handles);


% --------------------------------------------------------------------
function interpolate_Callback(~, ~, handles)
% hObject    handle to interpolate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[rows,cols] = size(handles.u(:,:,1));





for i = 1:handles.N
    
    u = handles.u(:,:,i);
    u(u==0) = NaN;
    handles.u(:,:,i) = inpaint_nans(u,3);
    u = handles.v(:,:,i);
    u(u==0) = NaN;
    handles.v(:,:,i) = inpaint_nans(u,3);
    
    %     vector = handles.u(:,:,i) + sqrt(-1)*handles.v(:,:,i);
    %
    %     [indx,indy] = find(abs(vector) == 0);
    %
    %     while ~isempty(indx)
    %         for z=1:length(indx)
    %             k = [max(3,indx(z))-2:min(rows-2,indx(z))+2];
    %             m = [max(3,indy(z))-2:min(cols-2,indy(z))+2];
    %             tmpvec = vector(k,m);
    %             tmpvec = tmpvec(find(tmpvec));
    %             vector(indx(z),indy(z)) = mean(real(tmpvec))+ sqrt(-1)*mean(imag(tmpvec));
    %         end
    %         try
    %             [indx,indy] = find(abs(vector) == 0);
    %         catch
    %             indx = [];
    %         end
    %     end
end


% Update data

handles.u(:,:,handles.N+1) = mean(handles.u(:,:,1:handles.N),3);
handles.v(:,:,handles.N+1) = mean(handles.v(:,:,1:handles.N),3);
if handles.state3d
    handles.w(:,:,handles.N+1) = mean(handles.w(:,:,1:handles.N),3);
    handles.wf = zeros(rows,cols,handles.N);
    for i = 1:handles.N
        handles.wf(:,:,i) = handles.w(:,:,i) - handles.w(:,:,handles.N+1);
    end
end

handles.dx = handles.x(1,2) - handles.x(1,1);
handles.dy = handles.y(2,1) - handles.y(1,1);

%
for i = 1:handles.N
    handles.vf(:,:,i) = handles.v(:,:,i) - handles.v(:,:,handles.N+1);
    handles.uf(:,:,i) = handles.u(:,:,i) - handles.u(:,:,handles.N+1);
end
%
for i = 1:handles.N+1
    [handles.dudx(:,:,i),handles.dudy(:,:,i)] = lsgradient(handles.u(:,:,i),handles.dx, handles.dy);
    [handles.dvdx(:,:,i),handles.dvdy(:,:,i)] = lsgradient(handles.v(:,:,i),handles.dx, handles.dy);
end

update_gui(handles.fig,[],handles);
guidata(handles.fig,handles);


% --------------------------------------------------------------------
function export2MAT_Callback(~, ~, ~)
% hObject    handle to export2MAT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% file = [];
file = inputdlg('File Name','Input Name for CSV File');
if ~isempty (file)
    eval(['save ',file{1},' -struct handles']);
else
    errordlg ('Choose a valid file name !!! ');
end;
% save exported -struct handles


% --------------------------------------------------------------------
function delete_current_Callback(~, ~, handles)
% hObject    handle to delete_current (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% axes(handles.axes_main);


handles.u(:,:,handles.current) = [];
handles.v(:,:,handles.current) = [];
handles.dudx(:,:,handles.current) = [];
handles.dudy(:,:,handles.current) = [];
handles.dvdx(:,:,handles.current) = [];
handles.dvdy(:,:,handles.current) = [];
handles.uf(:,:,handles.current) = [];
handles.vf(:,:,handles.current) = [];

handles.N = handles.N - 1;

update_gui(handles.fig,[],handles);
guidata(handles.fig,handles);


% --- Executes on button press in pod_pushbutton.
function pod_pushbutton_Callback(~, ~, handles)
% hObject    handle to pod_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.property) || isempty(handles.i)
    errordlg('First, pick the quantity and region of interest !!!');
else
    podbox(handles);  % call to podbox
    handles.podOn = 1;
end;

guidata(handles.fig,handles);
