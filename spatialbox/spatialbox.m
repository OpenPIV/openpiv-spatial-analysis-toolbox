function varargout = spatialbox(varargin)
warning off
% SPATIALBOX - Spatial Toolbox for the PIV data analysis
%    SPATIALBOX launch SpatialBox GUI.
%    SPATIALBOX('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 16-May-2003 21:09:40
% modified by Alex 18-Apr-2004 22:12:00
% - the cil_uigetfiles is implemented
% modified: Alex, 19/04/04, 24:40
% - handles.rs is added for the Reynolds stress
% Last modified: Alex 20/04/04
% - animation could run from every point, not only from the first file
% - see startpoint variable in the animation_Callback routine
% - major bug: for long series it is impossible to save every animation
% that user might want only to see on the screen. The animation button
% is not longer connected to the movie button. Animation runs animation. 
% if you want movie, you put again the first file as you like, and then
% press movie, it works exactly as animation, but at the end asks about
% AVI file name. The only difference between two, is that in order to have
% very long movie - user might need more than 1Gb of memory. I don't know
% how to check. But we didn't promise 1000 vector maps movies. There are
% probably other ways to make short movies, and then combine them together.
% or write directly to disk? look at avifile, it might work.
% - AVIFILE is introduced, instead of MOVIE2AVI - see 
% http://www.mathworks.com/support/tech-notes/1200/1204.html
% Modified: 21/04/04
% - no STOP button, animation and movie are toggle buttons, press one to
% run, press second to stop. every time movie is pressed, it asks for a new
% file
% - All fields is not called each time, but only when the name of the
% property changes. new handles.previous_property is introduced
% 
% Last modified: May 5, 2004, 01:31AM
% - all units are converted to SI: seconds, meter per second, etc.
% - there is no need to convert later, in timebox
% - 
%
%

global orighandles current_path;

if nargin == 0   % LAUNCH GUI
    current_path = cd;
    addpath(current_path);
    
    fig = openfig(mfilename,'reuse','invisible');
    movegui(fig,'center')
    
    % Generate a structure of handles to pass to callbacks, and store it. 
    handles = guihandles(fig);
    
    % Global several long strings
%     handles.inst_list = '-|u|v|(u^2+v^2)^(1/2)|vorticity|s_xx=du/dx|du/dy|dv/dx|s_yy=dv/dy|du/dx+dv/dy|s_xy';
%     handles.mean_list = '-|U|V|(U^2+V^2)^(1/2)|Vorticity|S_xx=dU/dx|dU/dy|dV/dx|S_yy=dV/dy|dU/dx+dV/dy|S_xy';
%     handles.fluct_list = '-|u''|v''|(u''^2+v''^2)^(1/2)|Vorticity''|s_xx''=du''/dx|du''/dy|dv''/dx|s_yy''=dv''/dy|du''/dx+dv''/dy|s_xy''';
%     handles.fluct_mean_list = '-|<u''^2>|<v''^2>|<u''v''>|T_u|T_v|diss|<u''v''>S_xy';
    handles.inst_list = '-|u|v|(u^2+v^2)^(1/2)|vorticity|s_xx=du/dx|du/dy|dv/dx|s_yy=dv/dy|du/dx+dv/dy|s_xy';
    handles.mean_list = '-|U|V|(U^2+V^2)^(1/2)|Vorticity|S_xx=dU/dx|dU/dy|dV/dx|S_yy=dV/dy|dU/dx+dV/dy|S_xy';
    handles.fluct_list = '-|u''|v''|(u''^2+v''^2)^(1/2)|Vorticity''|s_xx''=du''/dx|du''/dy|dv''/dx|s_yy''=dv''/dy|du''/dx+dv''/dy|s_xy''';
    handles.fluct_mean_list = '-|u rms|v rms|Reynolds stress|Turb. intensity (u) |Turb. intensity (v)|Dissipation|Turb. Energy Production';

    handles.fig = fig;
    handles.previous_quantity = '-';
    orighandles = handles;          % backup for the later re-opening of the data
    guidata(handles.fig, handles);  
    
    % Use system color scheme for figure:
    set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
    
    %     % Add a new toolbar with standard Matlab icons
    if isempty(findall(fig,'type','uitoolbar'))
        hToolbar = uitoolbar('Parent',fig);
        
        load tsi_icons;
        
        uipushtool('parent',hToolbar,'Click','spatialbox(''load_Callback'',gcbo, [], guidata(gcbo))',...
            'cdata',openfile,'Tag','openfilebtn');
        
        handles.hzoomin=uitoggletool('Parent',hToolbar,'OnCallback',...
    'spatialbox(''toggle_zoomin'',gcbo, [], guidata(gcbo))',...
 'OffCallback','spatialbox(''toggle_off'',gcbo, [], guidata(gcbo))',...   
'cdata',zoomin,'Tag','zoominbtn','state','off','ToolTipString','Zoom In'); 

handles.hzoomout=uitoggletool('Parent',hToolbar,'OnCallback',...
    'spatialbox(''toggle_zoomout'',gcbo, [], guidata(gcbo))',...
'OffCallback','spatialbox(''toggle_off'',gcbo, [], guidata(gcbo))',...     
'cdata',zoomout,'Tag','zoomoutbtn','state','off','ToolTipString','Zoom Out'); 
   guidata(handles.fig, handles); 

uipushtool('parent',hToolbar,'ClickedCallback',...
'spatialbox(''push_zoomreset'',gcbo, [], guidata(gcbo))',...
'cdata',zoomall,'Tag','zoomallbtn', 'ToolTipString','View All');
  
   handles.hzoomx=uitoggletool('Parent',hToolbar,'OnCallback',...
    'spatialbox(''toggle_zoomx'',gcbo, [], guidata(gcbo))',...
'OffCallback','spatialbox(''toggle_off'',gcbo, [], guidata(gcbo))',... 
'cdata',zoomx,'Tag','zoomxbtn','state','off','ToolTipString','Zoom X'); 

handles.hzoomy=uitoggletool('Parent',hToolbar,'OnCallback',...
    'spatialbox(''toggle_zoomy'',gcbo, [], guidata(gcbo))',...
'OffCallback','spatialbox(''toggle_off'',gcbo, [], guidata(gcbo))',... 
'cdata',zoomy,'Tag','zoomxbtn','state','off','ToolTipString','Zoom Y'); 

   guidata(handles.fig, handles); 
        
            
        
     
    end
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
    catch
        disp(lasterr);
    end
    
end
warning on

% --------------------------------------------------------------------
function varargout = checkbox_arrow_Callback(h, eventdata, handles, varargin)
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
function varargout = checkbox_arrow_color_Callback(hObject, eventdata, handles, varargin)
% Color arrows depends on the handles.property, chosen from the list of the
% avialable quantities

% color / black callback
if (get(hObject,'Value') == 1) 
    handles.color = 1;
else
    % checkbox is not checked-take approriate action
    handles.color = 0;
    set(handles.color_quiver,'Visible','off');
end

guidata(handles.fig,handles);
update_gui(handles.fig,[],handles);


% --------------------------------------------------------------------
function varargout = popupmenu_quantity_Callback(h, eventdata, handles, varargin)
% Property Selection callback
% Assign handles.property to present as color or contour later, according
% to the flags of ensemble, fluctuatinons and possible quantities

% --------------Ensemble + Fluct --------------
if get(handles.checkbox_ensemble,'Value') == 1  
    if get(handles.checkbox_fluct,'Value') == 1 
        
        switch get(handles.popupmenu_quantity,'Value')
            case 1
                handles.property = [];    
                
                
            case 2 % uf^2 
                if isfield(handles,'uf2')
                    handles.property = handles.uf2;
                    
                else
                    handles.uf2 = mean(handles.uf.^2,3); %u'^2
                    handles.property = handles.uf2;
                end
                handles.units='[m/s]^2';
            case 3 % vf^2 
                if isfield(handles,'vf2')
                    handles.property = handles.vf2;
                else
                    handles.vf2 = mean(handles.vf.^2,3); %u'^2
                    handles.property = handles.vf2;
                end
                handles.units='[m/s]^2';
            case 4
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
                handles.units='[1/s^2]'  ;
                
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
                handles.units='[m^2/s^3]' ;
                
                if isfield(handles,'prod')
                    handles.property = handles.prod;
                else
                    handles.prod = -1*( mean(handles.uf.*handles.vf,3).*( handles.dudy(:,:,handles.N+1) + handles.dvdx(:,:,handles.N+1) ) + ...
                        mean(handles.uf.*handles.uf,3).*handles.dudx(:,:,handles.N+1)  +...
                        mean(handles.vf.*handles.vf,3).*handles.dvdy(:,:,handles.N+1));
                    handles.property = handles.prod;
                end
        end % of switch
        % ------------------------------------------- Ensemble ------
    else % of fluct, means only ensemble is chosen
        switch get(handles.popupmenu_quantity,'Value')
            case 1
                handles.property = [];    % default variable, choose properly later, alex, 16.02
                
            case 2
                handles.units='[m/s]' ;   
                
                handles.property = handles.u(:,:,handles.current);
            case 3
                handles.units='[m/s]' ;   
                handles.property = handles.v(:,:,handles.current);
            case 4
                handles.units='[m/s]' ;   
                handles.property = sqrt(handles.u(:,:,handles.current).^2+handles.v(:,:,handles.current).^2);%Velocity Magnitude
            case 5 % vorticity; 
                handles.units='[1/s]' ;   
                
                handles.property = handles.dvdx(:,:,handles.current) - handles.dudy(:,:,handles.current); 
            case 6
                handles.units='[1/s]' ;
                handles.property = handles.dudx(:,:,handles.current);
            case 7    
                handles.units='[1/s]' ;
                handles.property = handles.dudy(:,:,handles.current);
            case 8
                handles.units='[1/s]' ;
                handles.property = handles.dvdx(:,:,handles.current);
            case 9
                handles.units='[1/s]' ;
                handles.property = handles.dvdy(:,:,handles.current);
            case 10
                handles.units='[1/s]' ;
                handles.property = handles.dudx(:,:,handles.current) + handles.dvdy(:,:,handles.current);
            case 11 % s_xy
                handles.units='[1/s]' ;
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
                
            case 2
                handles.units='[m/s]' ;  
                handles.property = handles.uf(:,:,handles.current);
                handles.cmin = min(handles.uf(:)); handles.cmax = max(handles.uf(:));
            case 3
                handles.units='[m/s]' ;  
                handles.property = handles.vf(:,:,handles.current);
                handles.cmin = min(handles.vf(:)); handles.cmax = max(handles.vf(:));
            case 4
                handles.units='[m/s]' ;  
                handles.property = sqrt(handles.uf(:,:,handles.current).^2+handles.vf(:,:,handles.current).^2);%Velocity Magnitude
                handles.cmin = min(handles.uf(:).^2+handles.vf(:).^2); 
                handles.cmax = max(handles.uf(:).^2+handles.vf(:).^2);
            case 5
                handles.units='[1/s]' ;  
                handles.property = dvdx - dudy;
                %                 handles.cmin = Inf;
                %                 handles.cmax = -Inf;
                %                 for i = 1:handles.N
                %                     tmp = (handles.dvdx(:,:,i) - handles.dvdx(:,:,handles.N+1))-...
                %                         (handles.dudy(:,:,i) - handles.dudy(:,:,handles.N+1));
                %                     handles.cmin = min(handles.cmin,min(tmp,3));
                %                     handles.cmax = max(handles.cmax,max(tmp,3));
                %                 end
            case 6
                handles.units='[1/s]' ; 
                handles.property = dudx;
            case 7    
                handles.units='[1/s]' ; 
                handles.property = dudy;
            case 8
                handles.units='[1/s]' ; 
                handles.property = dvdx;
            case 9
                handles.units='[1/s]' ; 
                handles.property = dvdy;
            case 10
                handles.units='[1/s]' ; 
                handles.property = dudx+dvdy;
            case 11 % s_xy
                handles.units='[1/s]' ; 
                handles.property = 0.5*(dvdx+dudy);      
        end
    else % fluct, means nothing is selected, u values
        switch get(handles.popupmenu_quantity,'Value')
            case 1
                handles.property = [];    % default variable, choose properly later, alex, 16.02
                
            case 2
                handles.units='[m/s]' ; 
                handles.property = handles.u(:,:,handles.current);
            case 3
                handles.units='[m/s]' ; 
                handles.property = handles.v(:,:,handles.current);
            case 4
                handles.units='[m/s]' ; 
                handles.property = sqrt(handles.u(:,:,handles.current).^2 + handles.v(:,:,handles.current).^2);%Velocity Magnitude
            case 5
                handles.units='[1/s]' ; 
                handles.property = handles.dvdx(:,:,handles.current) - handles.dudy(:,:,handles.current); % omega;
            case 6
                handles.units='[1/s]' ;
                handles.property = handles.dudx(:,:,handles.current);
            case 7    
                handles.units='[1/s]' ;
                handles.property = handles.dudy(:,:,handles.current);
            case 8
                handles.units='[1/s]' ;
                handles.property = handles.dvdx(:,:,handles.current);
            case 9
                handles.units='[1/s]' ;
                handles.property = handles.dvdy(:,:,handles.current);
            case 10
                handles.units='[1/s]' ;
                handles.property = handles.dudx(:,:,handles.current) + handles.dvdy(:,:,handles.current);
            case 11 % s_xy
                handles.units='[1/s]' ;
                handles.property = 0.5*(handles.dvdx(:,:,handles.current) + handles.dudy(:,:,handles.current));      
        end
    end
end

if ~isempty(handles.property) % something was selected
    if strcmp( get(handles.edit_arrow_size,'Enable'),'on')
        set(handles.checkbox_arrow_color,'Enable','on');
    end
end

% % Modified by Alex - it is stupid, we update the whole dataset check each 
% % animation frame, or anything, if All Fields is selected, 21.04.04
% % We need more 'intelligent' check - it should update all fields, only 
% if the call was from the popupmenu_eachfield_Button, and not from the
% animation or something. Only if you change the quantity, then it should
% update the the color for the All fields.
% Record the string which is in the quantity popupmenu in something
tmp = cellstr(get(handles.popupmenu_quantity,'String'));
if strcmp(handles.previous_quantity,tmp{get(handles.popupmenu_quantity,'Value')}) == 0 & ...
        get(handles.popupmenu_eachfield,'Value') == 3 % all fields, update it
    handles.previous_quantity = tmp{get(handles.popupmenu_quantity,'Value')};
    popupmenu_eachfield_Callback(handles.fig, [], handles);
else
    handles.previous_quantity = tmp{get(handles.popupmenu_quantity,'Value')};
    guidata(handles.fig,handles);
    update_gui(handles.fig,[],handles);
end

if strcmp(get(handles.ed_text,'Visible'),'on')
  if ~isempty(handles.property)
      handles.ed_mean_value=num2str(mean(handles.property(:)));
      handles.ed_std_value=num2str(std(handles.property(:)));
      handles.ed_min_value=num2str(min(handles.property(:)));
      handles.ed_max_value=num2str(max(handles.property(:)));
    set(handles.ed_text,'String',handles.previous_quantity);
    set(handles.ed_mean,'String',handles.ed_mean_value);
    set(handles.ed_std,'String',handles.ed_std_value);
    set(handles.ed_min,'String',handles.ed_min_value);
    set(handles.ed_max,'String',handles.ed_max_value);
    guidata(handles.fig,handles);
   end;
end;

% --------------------------------------------------------------------
function varargout = edit_numcolors_Callback(h, eventdata, handles, varargin)
% change number of colors
handles.numcolors = str2num(get(h,'String'));
if isempty(handles.numcolors),
    set(h,'String',10);
    handles.numcolors = 10;
end
guidata(handles.fig,handles);
update_gui(handles.fig,[],handles);


% --------------------------------------------------------------------
function varargout = update_gui(h, eventdata, handles, varargin)
% update_gui is responsible for update of the screen with current property and contour type

axes(handles.axes_main);
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
            handles.cmin = str2num(get(handles.edit_min_clim,'String'));
            handles.cmax = str2num(get(handles.edit_max_clim,'String'));
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
        handles.colorbar = mcolorbar('vert','peer',handles.axes_main);
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
            handles.color_quiver = quiverc(handles.x,handles.y,handles.u(:,:,handles.current),handles.v(:,:,handles.current),handles.arrow_scale,handles.property);
        end
        hold off;
        if handles.colorbar_flag                       % colorbar if necessary
            handles.colorbar = mcolorbar('vert','peer',handles.axes_main);
        end
    else
        hold on;
        % ----------------------checknox_fluct is checked -------------
        if get(handles.checkbox_fluct,'Value') == 1
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
                handles.quiver = quiver(handles.x,handles.y,handles.uf(:,:,handles.current),handles.vf(:,:,handles.current),handles.arrow_scale,'k');
            else
                handles.quiver = quiver(handles.x,handles.y,handles.u(:,:,handles.current),handles.v(:,:,handles.current),handles.arrow_scale,'k');
                % nothing to display, arrows of ensemble + fluctuation = 0
                % by definition
                set(handles.quiver,'Visible','Off');
            end
        else
            handles.quiver = quiver(handles.x,handles.y,handles.u(:,:,handles.current),handles.v(:,:,handles.current),handles.arrow_scale,'k');
        end
        hold off;
    end
end

set(handles.axes_main,'XLim',[min(handles.x(:)),max(handles.x(:))]);
set(handles.axes_main,'YLim',[min(handles.y(:)),max(handles.y(:))]);
xlabel('x [m]');
ylabel('y [m]');
guidata(handles.fig,handles);

% --------------------------------------------------------------------
function varargout = edit_min_clim_Callback(h, eventdata, handles, varargin)
% first editbox for 'manual' checkbox callback
handles.cmin = str2num(get(h,'String')); % update cmin
handles.alltodisp = 0;
handles.allfields = 1;
guidata(handles.fig,handles);
% update_gui(handles.fig,[],handles);

% --------------------------------------------------------------------
function varargout = edit_max_clim_Callback(h, eventdata, handles, varargin)
% second editbox for 'manual' checkbox callback
handles.cmax=str2num(get(h,'String')); % get new value for cmax
handles.alltodisp = 0;
handles.allfields = 1;
guidata(handles.fig,handles);
% update_gui(handles.fig,[],handles);

% --------------------------------------------------------------------
function varargout = pushbutton_set_clim_Callback(h, eventdata, handles, varargin)
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
function varargout = pushbutton_previous_Callback(h, eventdata, handles, varargin)
if handles.current > 1
    handles.current = handles.current - 1;          % update handles.current
    set(handles.edit_current,'String',handles.current);        % display num. of current file being processed
    delete(get(handles.axes_main,'children'));
    guidata(handles.fig,handles);
    popupmenu_quantity_Callback(handles.fig, [], handles);
else
    beep;
end

% --------------------------------------------------------------------
function varargout = pushbutton_next_Callback(h, eventdata, handles, varargin)
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
function varargout = edit_current_Callback(h, eventdata, handles, varargin)
tmp = eval(get(handles.edit_current,'String'));
if tmp > 0 & tmp <= handles.N % valid number of map
    handles.current = tmp;
    guidata(handles.fig,handles);
    popupmenu_quantity_Callback(handles.fig, [], handles);
else
    beep
    set(handles.edit_current,'String',handles.current);
end

% --------------------------------------------------------------------
function varargout = edit_arrow_size_Callback(h, eventdata, handles, varargin)
% 'Scale' editbox callback
handles.arrow_scale = eval(get(h,'String'));     % update handles.scale
if handles.arrow_scale == 0 | isempty(handles.arrow_scale) % no scaling
    handles.arrow_scale = [];
end
guidata(handles.fig,handles);
update_gui(handles.fig,[],handles);

% --------------------------------------------------------------------
function varargout = pushbutton_animate_Callback(h, eventdata, handles, varargin)
% animate button callback
% set(handles.pushbutton_stop,'Enable','on');
if get(handles.pushbutton_animate,'Value') == 1
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
    ;
end
set(handles.pushbutton_animate,'Value',0);
% set(handles.pushbutton_save_movie,'Enable','on');
% set(handles.pushbutton_stop,'Enable','off');
guidata(handles.fig,handles);

% --------------------------------------------------------------------
function varargout = pushbutton_save_movie_Callback(h, eventdata, handles, varargin)

if get(handles.pushbutton_save_movie,'Value') == 1 
    file = [];
    file = inputdlg('File Name','Input File Name for the movie');   
    if isempty(file) | exist(file{1},'file') | exist([file{1},'.avi'],'file') 
        set(handles.pushbutton_save_movie,'Value',0);
        return
    end
    handles.mov = avifile(file{1},'compression','none','quality',100,'fps',15); 
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
guidata(handles.fig,handles);


% % --------------------------------------------------------------------
% function varargout = pushbutton_stop_Callback(hObject, eventdata, handles, varargin)
% % stop button callback
% % global stop1;       
% % stop1 = 1;        % set stop flag
% guidata(hObject,handles);


% --------------------------------------------------------------------
function varargout = checkbox_label_Callback(h, eventdata, handles, varargin)
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
function varargout = checkbox_colorbar_Callback(h, eventdata, handles, varargin)
if get(h,'Value') == 1
    handles.colorbar_flag = 1;
    guidata(handles.fig,handles);
    update_gui(handles.fig,[],handles);
else
    delete(handles.colorbar);
    handles.colorbar_flag = 0;
    guidata(handles.fig,handles);
end
% update_gui(handles.fig,[],handles);

% --------------------------------------------------------------------



% --------------------------------------------------------------------
function popupmenu_eachfield_Callback(hObject, eventdata, handles)
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
function checkbox_ensemble_Callback(hObject, eventdata, handles)

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
    set(handles.popupmenu_eachfield,'String','Each Field|All to Display|All Fields|Manual');
    if (get(handles.checkbox_fluct,'Value') == get(handles.checkbox_fluct,'Max'))
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
function checkbox_fluct_Callback(hObject, eventdata, handles)
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
function pushbutton_export_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.export_figure = figure;
% handles.export_axes   = axes;
copyobj(handles.axes_main,handles.export_figure);
% copyobj(get(handles.axes_main,'Children'),handles.export_axes);

set(handles.export_figure,'Units','normalized');
set(get(handles.export_figure,'children'),'Units','normalized');
set(get(handles.export_figure,'children'),'Position',[0.13 0.11 0.775 0.815]);
set(get(handles.export_figure,'children'),'Box','on');

if handles.colorbar_flag
    % copyobj(handles.colorbar,handles.export_figure);
    mcolorbar; % (get(handles.export_figure,'Children'));
end
guidata(handles.fig, handles);

% ------------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% ------------------   Load and prepare data module --------------------------------------------------
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global orighandles;
if isfield(handles,'restoreorig') 
    handles = orighandles;
end  
handles.restoreorig = 1;

[gui_files,gui_path,handles.dt,handles.scale,handles.state3d] = cil_uigetfiles;

% Alex, 04.05.04, cil_uigetfiles returns this parameter
% if ~isempty(findstr(lower(gui_files{1}),'.v3d')) 
%     handles.state3d = 1;
% end;

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
            [header,d] = svecread(fullfile(handles.path,handles.files{1}));
            [rows,cols,k] = size(d);
            [handles.u,handles.v,handles.w] = deal(zeros(rows,cols,handles.N+1)); % 11.04.04, Alex
            handles.x = d(:,:,1)*handles.scale/1000;
            handles.y = d(:,:,2)*handles.scale/1000;    
            handles.z = d(:,:,3)*handles.scale/1000;
            handles.u(:,:,1) = d(:,:,4)*handles.scale/1000/handles.dt;
            handles.v(:,:,1) = d(:,:,5)*handles.scale/1000/handles.dt;    
            handles.w(:,:,1) = d(:,:,6)*handles.scale/1000/handles.dt;
            for i = 2:handles.N
                d = svecread([handles.path,filesep,handles.files{i}],1,8);
                handles.u(:,:,i) = d(:,:,4)*handles.scale/1000/handles.dt;
                handles.v(:,:,i) = d(:,:,5)*handles.scale/1000/handles.dt;
                handles.w(:,:,i) = d(:,:,6)*handles.scale/1000/handles.dt;
            end
            clear d
        end
    case 0
        if ~isempty(findstr(handles.files{1},'vec'))            % process .vec files
            % short copy of readexpdir.m
            % read first file, determine the size
            [header,d] = vecread(fullfile(handles.path,handles.files{1}));
            [rows,cols,k] = size(d);
            [handles.u,handles.v] = deal(zeros(rows,cols,handles.N+1)); % 11.04.04, Alex
            handles.x           = d(:,:,1)*handles.scale/1000;
            handles.y           = d(:,:,2)*handles.scale/1000;
            handles.u(:,:,1)    = d(:,:,3)*handles.scale/1000/handles.dt;
            handles.v(:,:,1)    = d(:,:,4)*handles.scale/1000/handles.dt;
            
            for i = 2:handles.N
                d = vecread([handles.path,filesep,handles.files{i}],1,5);
                handles.u(:,:,i) = d(:,:,3)*handles.scale/1000/handles.dt;
                handles.v(:,:,i) = d(:,:,4)*handles.scale/1000/handles.dt;
            end
            clear d
        end
end

handles.current = 1;                      % current file beeing displayed
% Display first file number, total number of files
set(handles.edit_current,'String',handles.current);
set(handles.edit_numfields,'String',handles.N);
% Initialize color, number of colors for contours
% handles.stop = 0;
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

% Alex 21.02.04
handles.dx = abs(handles.x(1,1) - handles.x(1,2));
handles.dy = abs(handles.y(1,1) - handles.y(2,1));

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
    [handles.dudx(:,:,i),handles.dudy(:,:,i)] = lsgradient(handles.u(:,:,i),handles.dx, handles.dy);  
    [handles.dvdx(:,:,i),handles.dvdy(:,:,i)] = lsgradient(handles.v(:,:,i),handles.dx, handles.dy); 
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

% Default plot is quiver of the first map.
handles.quiver = quiver(handles.x,handles.y,handles.u(:,:,handles.current),handles.v(:,:,handles.current),handles.arrow_scale,'k');
xlabel('x [m]');
ylabel('y [m]');

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
handles.gridX = handles.x(1,2) - handles.x(1,1);    
handles.gridY = handles.y(1,1) - handles.y(2,1);

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
set(handles.pushbutton_export,'Enable','on');
set(handles.pushbutton_stats,'Enable','on');


% Properties
set(handles.popupmenu_quantity,'String',handles.inst_list);
handles.property = []; 

% Very nice implementation of the Tab Panel for later TimeBox usage
% ---------------- store all handles of spatial controls ------------------
handles.spatial_controls=[handles.checkbox_ensemble,handles.checkbox_fluct,handles.checkbox_arrow,handles.checkbox_arrow_color,...
        handles.checkbox_label,handles.checkbox_colorbar,handles.edit_arrow_size,handles.edit_numcolors,...
        handles.edit_current,handles.edit_numfields,...
        handles.text_contour_quantity, handles.text_contourtype, handles.text_numberofcolors, handles.text7,handles.text2,...
        handles.text_arrow_size, handles.pushbutton_previous, handles.pushbutton_next, ...
        handles.pushbutton_animate, handles.pushbutton_save_movie, handles.pushbutton_export,...
        handles.frame_controls,handles.frame8,handles.frame_contour_quantity,handles.frame_contour_type,handles.frame7,...
        handles.frame_arrow, handles.text5, handles.popupmenu_quantity,handles.popupmenu_contour_type,...
        handles.popupmenu_eachfield,handles.pushbutton_stats];

% -------------- store all handles of select controls ----------------
handles.select_controls=[handles.pushbutton_selectpoints,handles.pushbutton_selectreg,handles.pushbutton_selectall,...
        handles.pushbutton_profile1,handles.pushbutton_time,handles.frame_select,handles.frame_region,...
        handles.pushbutton_reset, handles.rowpushbutton, handles.colpushbutton,handles.frame_sel];

handles.arrow_ctrls=[handles.edit_arrow_size,handles.checkbox_arrow,handles.checkbox_arrow_color];
handles.stat_controls=[handles.ed_mean,handles.ed_std,handles.ed_min,handles.ed_max,handles.ed_text,...
        handles.ed_pushsavestl,handles.ed_frame,handles.ed_pushclose,handles.ed_textmean,...
        handles.ed_textstd,handles.ed_textmax,handles.ed_textmin];

set(handles.pushbutton_spatial,'Enable','on');
set(handles.pushbutton_select,'Enable','on');


% Update all handles structure
guidata(handles.fig,handles);

set(handles.fig,'pointer','arrow');


% --------------------------------------------------------------------
function exit_Callback(hObject, eventdata, handles)
% Find the highest parent - figure, and close it.
while ~strcmpi(get(hObject,'Type'),'figure'),
    hObject = get(hObject,'Parent');
end
delete(hObject);

% --- Executes on button press in pushbutton_spatial.
function pushbutton_spatial_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_spatial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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
function pushbutton_select_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.spatial_controls,'Visible','off');
set(handles.pushbutton_set_clim,'Visible','off');
set(handles.select_controls,'Visible','on');
set(handles.select_controls,'Enable','on');
set(handles.pushbutton_select,'String','> Select <');
set(handles.pushbutton_spatial,'String','Spatial');


set(handles.pushbutton_spatial,'FontWeight','normal');
set(handles.pushbutton_select,'FontWeight','bold');
val = get(handles.popupmenu_eachfield,'Value');

guidata(handles.fig,handles);


% --- Executes on button press in pushbutton_selectpoints.
function pushbutton_selectpoints_Callback(hObject, eventdata, handles)


set(handles.rowpushbutton,'Enable','off');
set(handles.pushbutton_selectreg,'Enable','off');
set(handles.colpushbutton,'Enable','off');
set(handles.pushbutton_selectall,'Enable','off');



set(handles.axes_main,'NextPlot','Add');
limX = xlim; 
limY = ylim;
leftcolX   = 1; 
bottomrowY = 1; 


while 1
    
    [x1,y1, buttonNumber] = ginput(1);
    
    % When the right button is pressed, stop the loop
    if (buttonNumber == 2) | (buttonNumber==3)
        break;
    end;
    col = fix (( x1 - limX(1,1) )/ handles.gridX+0.5  )+1; 
    row = fix(( y1 - limY(1,1) )/ handles.gridY+0.5  )+1; 
    
    % check for errors ----------------
    if col<1 | col>(fix((limX(1,2)-limX(1,1))/handles.gridX)+1) | row<1 | row...
            >(fix((limY(1,2)-limY(1,1))/handles.gridY)+1);
        guidata(handles.fig,handles);
        return;
    end;
    % ---------------------------------
    sizeI = size(handles.i,1);
    rightcolX = fix(( limX(1,2)-limX(1,1) )/  handles.gridX )+1;
    uprowY    = fix(( limY(1,2)-limY(1,1) )/  handles.gridY )+1;
    sizeJ = size(handles.j,1);
    numofcols = rightcolX - leftcolX + 1; 
    numofrows = uprowY - bottomrowY + 1;
    
    handles.i(sizeI+1,1) = row; 
    handles.j(sizeJ+1,1) = col;
    
    line(limX(1,1)+(col-1)*...
        handles.gridX,limY(1,1)+(row-1)*handles.gridY,'Marker','o','Color','k','MarkerSize',8);
    
end;
guidata(handles.fig,handles);
% update_gui(handles.fig,[],handles);



% --- Executes on button press in pushbutton_selectreg.
function pushbutton_selectreg_Callback(hObject, eventdata, handles)
% check if all selected is off
set(handles.pushbutton_selectpoints,'Enable','off');
set(handles.rowpushbutton,'Enable','off');
set(handles.colpushbutton,'Enable','off');
set(handles.pushbutton_selectall,'Enable','off');
% if handles.Allselected==1
%    handles.i = []; handles.j = [];
% handles.rowlock=0; handles.columnlock=0;
% handles.previousSel=[];
% update_gui(gcbo,[],guidata(gcbo));
%     handles.Allselected=0;
% end;  
% ------------------------------
k = waitforbuttonpress;
point1 = get(gca,'CurrentPoint');    % button down detected
% point1 is 2x3 matrix, first 2 elements are x,y
finalRect = rbbox;                   % return figure units
point2 = get(gca,'CurrentPoint');    % button up detected
point1 = point1(1,1:2);              % extract x and y
point2 = point2(1,1:2);
p1 = min(point1,point2);             % calculate locations
offset = abs(point1-point2);         % and dimensions


limX = xlim; limY = ylim;                     % get Axis Limits
% -------------- calculate columns & rows ------------------
leftcolX = fix(( p1(1)-limX(1,1) )/ handles.gridX +1) + 1;  
rightcolX = fix(( p1(1)+offset(1)-limX(1,1) )/ handles.gridX )+1;
bottomrowY =  fix(( p1(2)-limY(1,1) )/ handles.gridY+1 )+1;
uprowY = fix(( p1(2) + offset(2) - limY(1,1) )/ handles.gridY)+1;


% -------boundary check ----------------
plotstateX=0;plotstateY=0;
if leftcolX<1     leftcolX=1;
    plotstateY= 1;
end;
rightLimit=fix((limX(1,2)-limX(1,1))/handles.gridX)+1;
if rightcolX>rightLimit  rightcolX=rightLimit; end;
uprowLimit=fix((limY(1,2)-limY(1,1))/handles.gridY)+1;
if bottomrowY<1  bottomrowY=1; 
    plotstateX= 1;   % we need it to plot in right way
end
if uprowY>uprowLimit     uprowY=uprowLimit; end;
% --------- selection checking ---------------  

%---------------------------------------         
% leftcolX=1;rightcolX=fix((limX(1,2)-limX(1,1))/handles.gridX);
% bottomrowY=1; uprowY=fix((limY(1,2)-limY(1,1))/handles.gridY);

sizeI = size(handles.i,1);
sizeJ = size(handles.j,1);
% -------------------- check if constrainX or constrainY is checked
% if  get(handles.checkbox_constrainx,'Value') == 1
%     bottomrowY = 1; 
%     uprowY = fix(( limY(1,2)-limY(1,1) )/ handles.gridY )+1; % select all the rows
% elseif get(handles.checkbox_constrainy,'Value') == 1
%     leftcolX = 1;
%     rightcolX = fix(( limX(1,2)-limX(1,1) ) / handles.gridX )+1; %select all the columns
% end;
numofcols=rightcolX-leftcolX+1; 
numofrows=uprowY-bottomrowY+1;
% -------------- errorchecking --------
if ~isempty(handles.previousSel)
    a=handles.previousSel;
    if ((rightcolX-leftcolX)==a(2)-a(1) & a(2)==rightcolX & handles.rowlock~=1)
        handles.columnlock=1;
    elseif    ((uprowY-bottomrowY)==a(4)-a(3) & a(4)==uprowY & handles.columnlock~=1)
        handles.rowlock=1;
    else
        errordlg('Your Selection is Invalid...');
        return;
    end; 
end;
% --------------------- Fill loop ----------------------------

for i1 = bottomrowY:uprowY
    handles.i(sizeI+1:sizeI+numofcols,1)=i1;
    handles.j(sizeJ+1:sizeJ+numofcols,1)=leftcolX:rightcolX;
    sizeI=sizeI+numofcols; sizeJ=sizeJ+numofcols;
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

% update_gui(handles.fig,[],handles);



% --- Executes on button press in pushbutton_selectall.
function pushbutton_selectall_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_selectall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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
leftcolX=1;rightcolX=fix((limX(1,2)-limX(1,1))/handles.gridX)+1;
bottomrowY=1; uprowY=fix((limY(1,2)-limY(1,1))/handles.gridY)+1;
numofcols=rightcolX-leftcolX+1;
sizeI=size(handles.i,1);
sizeJ=size(handles.j,1);
for i1 = bottomrowY:uprowY
    handles.i(sizeI+1:sizeI+numofcols,1)=i1;
    handles.j(sizeJ+1:sizeJ+numofcols,1)=leftcolX:rightcolX;
    sizeI=sizeI+numofcols; sizeJ=sizeJ+numofcols;
end
guidata(handles.fig,handles);
% update_gui(handles.fig,[],handles);


% --- Executes on button press in pushbutton_time.
function pushbutton_time_Callback(hObject, eventdata, handles)
% if handles.state3d
%     timeboxHandle = timebox_v3d(handles);
% else
%     
% end
timeboxHandle = timebox(handles); % timebox includes both versions with if ... else
guidata(handles.fig,handles);



% --- Executes on button press in pushbutton_profile1.
function pushbutton_profile1_Callback(hObject, eventdata, handles)
if isempty(handles.property)
    errordlg('First, pick the quantity !!!');
else
    distribHandle = distrib(handles);  % call to spatialbox
    handles.distribOn = 1;
end; 

guidata(handles.fig,handles);


% --- Executes on button press in pushbutton_reset.
function pushbutton_reset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%if ~isempty(handles.i) & ~isempty(handles.j)
handles.i = []; handles.j = [];
handles.rowlock=0; handles.columnlock=0;
handles.previousSel=[];
set(handles.pushbutton_selectpoints,'Enable','on');
set(handles.pushbutton_selectreg,'Enable','on');
set(handles.colpushbutton,'Enable','on');
set(handles.pushbutton_selectall,'Enable','on');
set(handles.rowpushbutton,'Enable','on');
%guidata(gcbo,handles);        % even if some bug with deleting, first of all we have to save handles
%H=findobj('type','line','marker','o');

%     delete(H);

%if isfield(handles,'selectionbox')     %bugfix 280404
% delete(handles.selectionbox);
%  rmfield(handles,'selectionbox');
%end;
guidata(handles.fig,handles);
update_gui(gcbo,[],guidata(gcbo));
%end;


% --- Executes during object creation, after setting all properties.
function figure_gradpiv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure_gradpiv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
load cil_logo
image(im,'Parent',handles.axes_main);
axis off


% --- Executes during object creation, after setting all properties.
function axes_main_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.axes_main = hObject;
guidata(hObject, handles);


% % --- Executes on button press in checkbox_constrainx.
% function checkbox_constrainx_Callback(hObject, eventdata, handles)
% if get(handles.checkbox_constrainx,'Value') == 1
%     set(handles.checkbox_constrainy,'value',0);
% end;
% guidata(handles.fig,handles);
% % update_gui(handles.fig,[],handles);



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
sym = '';

filled = 0;
ls = '-';
ms = '';
col = '';

nin = nargin;
% Parse the string inputs
while isstr(varargin{nin}),
    vv = varargin{nin};
    if ~isempty(vv) & strcmp(lower(vv(1)),'f')
        filled = 1;
        nin = nin-1;
    else
        [l,c,m,msg] = colstyle(vv);
        if ~isempty(msg), 
            error(sprintf('Unknown option "%s".',vv));
        end
        if ~isempty(l), ls = l; end
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

% Modified, Alex Liberzon, 2002
if nin==4 | nin==6, % quiver(u,v,z,s) or quiver(x,y,u,v,z,s)
    autoscale = varargin{nin-1};
    z = varargin{nin};
end

% Scalar expand u,v
if prod(size(u))==1, u = u(ones(size(x))); end
if prod(size(v))==1, v = v(ones(size(u))); end

if autoscale,
    % Base autoscale value on average spacing in the x and y
    % directions.  Estimate number of points in each direction as
    % either the size of the input arrays or the effective square
    % spacing if x and y are vectors.
    if min(size(x))==1, n=sqrt(prod(size(x))); m=n; else [m,n]=size(x); end
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
uu = [x;x+u;repmat(NaN,size(u))];
vv = [y;y+v;repmat(NaN,size(u))];

% h1 = plot(uu(:),vv(:),[col ls]);

% Prepare color matrix
z = [z(:)';z(:)';NaN*z(:)'];

h1 = patch([uu(:),uu(:)],[vv(:),vv(:)], [z(:),z(:)],'Parent',ax,'EdgeColor','Flat','FaceColor','None');

if plotarrows,
    % Make arrow heads and plot them
    hu = [x+u-alpha*(u+beta*(v+eps));x+u; ...
            x+u-alpha*(u-beta*(v+eps));repmat(NaN,size(u))];
    hv = [y+v-alpha*(v-beta*(u+eps));y+v; ...
            y+v-alpha*(v+beta*(u+eps));repmat(NaN,size(v))];
    hold on
    %  h2 = plot(hu(:),hv(:),[col ls]);
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
function rowpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to rowpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.pushbutton_selectpoints,'Enable','off');
set(handles.pushbutton_selectreg,'Enable','off');
set(handles.colpushbutton,'Enable','off');
set(handles.pushbutton_selectall,'Enable','off');

set(handles.axes_main,'NextPlot','Add');
limX = xlim; 
limY = ylim;
leftcolX   = 1; 
bottomrowY = 1; 


while 1
    
    [x1,y1, buttonNumber] = ginput(1);
    
    % When the right button is pressed, stop the loop
    if (buttonNumber == 2) | (buttonNumber==3)
        break;
    end;
    col = fix (( x1 - limX(1,1) )/ handles.gridX+0.5  )+1; 
    row = fix(( y1 - limY(1,1) )/ handles.gridY+0.5  )+1; 
    % check for errors ----------------
    if col<1 | col>(fix((limX(1,2)-limX(1,1))/handles.gridX)+1) | row<1 | row...
            >(fix((limY(1,2)-limY(1,1))/handles.gridY)+1);
        guidata(handles.fig,handles);
        return;
    end;
    % ---------------------------------
    sizeI = size(handles.i,1);
    rightcolX = fix(( limX(1,2)-limX(1,1) )/  handles.gridX )+1;
    uprowY    = fix(( limY(1,2)-limY(1,1) )/  handles.gridY )+1;
    sizeJ = size(handles.j,1);
    numofcols = rightcolX - leftcolX + 1; 
    numofrows = uprowY - bottomrowY + 1;
    
    handles.i(sizeI+1:sizeI+numofrows,1) = 1:uprowY ;
    handles.j(sizeJ+1:sizeJ+numofrows,1) = col;
    row = 1:uprowY;
    topLeft(1)=uprowY; bottomRight(1)=1;
    
    line(limX(1,1)+(col-1)*...
        handles.gridX,limY(1,1)+(row-1)*handles.gridY,'Marker','o','Color','k','MarkerSize',8);
    
    % plot(handles.i*handles.gridX,handles.j*handles.gridY,'Marker','o','Color','k','MarkerSize',8);
    
end;
guidata(handles.fig,handles);
% update_gui(handles.fig,[],handles);





% --- Executes on button press in colpushbutton.
function colpushbutton_Callback(hObject, eventdata, handles)
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
bottomrowY = 1; 

while 1
    
    [x1,y1, buttonNumber] = ginput(1);
    
    % When the right button is pressed, stop the loop
    if (buttonNumber == 2) | (buttonNumber==3)
        break;
    end;
    col = fix (( x1 - limX(1,1) )/ handles.gridX+0.5  )+1; 
    row = fix(( y1 - limY(1,1) )/ handles.gridY+0.5  )+1; 
    % -------- find the corners of rectangle ----
    
    
    % check for errors ----------------
    if col<1 | col>(fix((limX(1,2)-limX(1,1))/handles.gridX)+1) | row<1 | row...
            >(fix((limY(1,2)-limY(1,1))/handles.gridY)+1);
        guidata(handles.fig,handles);
        return;
    end;
    % ---------------------------------
    sizeI = size(handles.i,1);
    rightcolX = fix(( limX(1,2)-limX(1,1) )/  handles.gridX )+1;
    uprowY    = fix(( limY(1,2)-limY(1,1) )/  handles.gridY )+1;
    sizeJ = size(handles.j,1);
    numofcols = rightcolX - leftcolX + 1; 
    numofrows = uprowY - bottomrowY + 1;
    %
    handles.i(sizeI+1:sizeI+numofcols,1) = row;
    handles.j(sizeJ+1:sizeJ+numofcols,1) = 1:rightcolX;
    col = 1:rightcolX;
    %
    line(limX(1,1)+(col-1)*...
        handles.gridX,limY(1,1)+(row-1)*handles.gridY,'Marker','o','Color','k','MarkerSize',8);
end;


guidata(handles.fig,handles);

%---- zoom handling functions -----------------
function toggle_zoomin(hObject, eventdata, handles)
set(handles.hzoomout,'state','off');
set(handles.hzoomx,'state','off');
set(handles.hzoomy,'state','off');
putdowntext('zoomin',gcbf);
guidata(handles.fig,handles);
 
function toggle_zoomout(hObject, eventdata, handles)
set(handles.hzoomin,'state','off');
set(handles.hzoomx,'state','off');
set(handles.hzoomy,'state','off');
putdowntext('zoomout',gcbf);
guidata(handles.fig,handles);

function toggle_zoomx(hObject, eventdata, handles)

set(handles.hzoomin,'state','off');
set(handles.hzoomy,'state','off');
set(handles.hzoomout,'state','off');
zoom(gcbf,'xon')
guidata(handles.fig,handles);

function toggle_zoomy(hObject, eventdata, handles)

set(handles.hzoomin,'state','off');
set(handles.hzoomx,'state','off');
set(handles.hzoomout,'state','off');
zoom(gcbf,'yon')
guidata(handles.fig,handles);

function push_zoomreset(hObject, eventdata, handles)
zoom(gcbf,'reset');


function toggle_off(hObject, eventdata, handles)
zoom(gcbf,'off')
% ----------------------------------------------------


% --- Executes during object creation, after setting all properties.
function ed_mean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ed_mean_Callback(hObject, eventdata, handles)
% hObject    handle to ed_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_mean as text
%        str2double(get(hObject,'String')) returns contents of ed_mean as a double


% --- Executes during object creation, after setting all properties.
function ed_std_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ed_std_Callback(hObject, eventdata, handles)
% hObject    handle to ed_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_std as text
%        str2double(get(hObject,'String')) returns contents of ed_std as a double


% --- Executes during object creation, after setting all properties.
function ed_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ed_min_Callback(hObject, eventdata, handles)
% hObject    handle to ed_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_min as text
%        str2double(get(hObject,'String')) returns contents of ed_min as a double


% --- Executes during object creation, after setting all properties.
function ed_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ed_max_Callback(hObject, eventdata, handles)
% hObject    handle to ed_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_max as text
%        str2double(get(hObject,'String')) returns contents of ed_max as a double


% --- Executes on button press in ed_pushsavestl.
function ed_pushsavestl_Callback(hObject, eventdata, handles)
% hObject    handle to ed_pushsavestl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = inputdlg('File Name','Save File Name ');  
fid = fopen(file{1},'w');
 if fid ~= -1, 
header_string='TSI_STL_VERSION 1.1';
a=handles.previous_quantity;  % name of the current property
names_string=strcat(['File, ' a '_Mean, ' a '_StdDev, ' a '_Min, ' a '_Max']);
 fprintf(fid,'%s\n', header_string);
 fprintf(fid,'%s\n', names_string);
  if get(handles.checkbox_ensemble,'Value') == 1  % if ensemble is checked - no file name
     filename='Ensemble'; 
  else
     filename=handles.files{handles.current};
  end;   
 fprintf(fid,'%s, %s, %s, %s, %s\n',filename,...
 handles.ed_mean_value,handles.ed_std_value,handles.ed_min_value,handles.ed_max_value);
 fclose(fid);
else
    errordlg('Please, select a valid file name');
end;  

% --- Executes on button press in pushbutton_stats.
function pushbutton_stats_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pushbutton_stats
if isempty(handles.property)
    errordlg('First, pick the quantity !!!');
 else  
    set(handles.stat_controls,'Visible','on');
    set(handles.stat_controls,'Enable','on');
    handles.ed_mean_value=num2str(mean(handles.property(:)));
      handles.ed_std_value=num2str(std(handles.property(:)));
      handles.ed_min_value=num2str(min(handles.property(:)));
      handles.ed_max_value=num2str(max(handles.property(:)));
    set(handles.ed_text,'String',handles.previous_quantity);
    set(handles.ed_mean,'String',handles.ed_mean_value);
    set(handles.ed_std,'String',handles.ed_std_value);
    set(handles.ed_min,'String',handles.ed_min_value);
    set(handles.ed_max,'String',handles.ed_max_value);
    guidata(handles.fig,handles);
end;
% --- Executes on button press in ed_pushclose.
function ed_pushclose_Callback(hObject, eventdata, handles)
% hObject    handle to ed_pushclose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


set(handles.stat_controls,'Visible','off');


% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


