function varargout = timebox_v3d(varargin)
%TIMEBOX_V3D First prototype of the time toolbox
% Last modified: 18-04-04, 22:28, Alex
% - see if nargin < 1 development part.
% Last modified: April 20, 01:50, Alex
% - majou changes, all in loops, see update_gui
% - at least time axes works now
% Last modified: April 24, 17:01
% - handles.property is for time series
% - handles.quantity is for frequency space or correlation
% Last modified: April 28, 2004, 12:00
% - real(fft) is actually the cosine transform
% according to Pope. p. 69 and p. 680 it should work
% Last modified: April 28, 2004, 12:00
% - bug with uf, vf
% 

if strcmp(inputname(1),'handles') 
    % GUI
    
    
    fig = openfig(mfilename,'reuse');
    movegui(fig, 'center');
    
    % Generate a structure of handles to pass to callbacks, and store it. 
    handles = guihandles(fig);
    
    
    set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
    
    % ---------------------------------------------------------------------
    if isempty(findall(fig,'type','uitoolbar'))
        
        % Add a new toolbar
        hToolbar = uitoolbar('Parent',fig);
        
        load tsi_icons;
        
        uipushtool('parent',hToolbar,'Click','process1(''load_Callback'',gcbo, [], guidata(gcbo))',...
            'cdata',openfile,'Tag','openfilebtn');
        uipushtool('parent',hToolbar,'Click','savefile_Callback',...
            'cdata',savefile,'Tag','savefilebtn');
        uitoggletool('parent',hToolbar,'ClickedCallback','putdowntext(''zoomin'',gcbf)',...
            'cdata',zoomin,'Tag','zoominbtn','state','off','ToolTipString','Zoom In');
        uitoggletool('parent',hToolbar,'Click','putdowntext(''zoomout'',gcbf)',...
            'cdata',zoomout,'Tag','zoomoutbtn','state','off', 'ToolTipString','Zoom Out');
        uipushtool('parent',hToolbar,'Click','zoom(gcbf,''reset'')',...
            'cdata',zoomall,'Tag','zoomallbtn', 'ToolTipString','View All');
        uitoggletool('parent',hToolbar,'Click','zoom(gcbf,''xon'')',...
            'cdata',zoomx,'Tag','zoomxbtn','state','off', 'ToolTipString','Zoom X');
        uitoggletool('parent',hToolbar,'Click','zoom(gcbf,''yon'')',...
            'cdata',zoomy,'Tag','zoomybtn','state','off', 'ToolTipString','Zoom Y');
        %         uitoggletool('parent',hToolbar,'Click','grid',...
        %             'Tag','gridbtn','state','on', 'ToolTipString','Grid On/Off');
    end
    
    
    % ---------------------------------------------------------------------
    
    % Out data coming from process1 in the initial stage
    handles.data = varargin{1};
    handles.t = (0:handles.data.N-1)*handles.data.dt;
    handles.Fs = 1/handles.data.dt;      
    
    handles.maxX = max(handles.data.x(:));
    handles.maxY = max(handles.data.y(:));
    [r,c,k] = size(handles.data.u);
    handles.linestyle = {'r-','b-','k-','m-','r--','b--','k--','m--','r-.','b-.','k-.','m-.','r:','b:','k:','m:'};
    
    handles.indm = zeros(r,c);
    handles.ind = sub2ind([r,c],handles.data.i,handles.data.j);     % long vector of indices
    handles.indm(handles.ind) = true;                               % 2D matrix, could visualize the selection
    handles.rows  = max(handles.data.i) - min(handles.data.i) + 1;
    handles.cols  = max(handles.data.j) - min(handles.data.j) + 1;
    
    % Calculate u and v average velocities over spatial regions
    handles.umean = squeeze( handles.data.u(handles.data.i(1),handles.data.j(1),1:end-1) );
    handles.vmean = squeeze( handles.data.v(handles.data.i(1),handles.data.j(1),1:end-1) );
    handles.wmean = squeeze( handles.data.w(handles.data.i(1),handles.data.j(1),1:end-1) );
    handles.ufmean = squeeze( handles.data.uf(handles.data.i(1),handles.data.j(1),1:end) );
    handles.vfmean = squeeze( handles.data.vf(handles.data.i(1),handles.data.j(1),1:end) );
    handles.wfmean = squeeze( handles.data.wf(handles.data.i(1),handles.data.j(1),1:end) );

    
    for i = 2:length(handles.data.i)
        handles.umean  = (1/i)*(handles.umean*(i-1) + squeeze( handles.data.u(handles.data.i(i),handles.data.j(i),1:end-1)));
        handles.vmean  = (1/i)*(handles.vmean*(i-1) + squeeze( handles.data.v(handles.data.i(i),handles.data.j(i),1:end-1)));
        handles.wmean  = (1/i)*(handles.wmean*(i-1) + squeeze( handles.data.w(handles.data.i(i),handles.data.j(i),1:end-1)));
        handles.ufmean = (1/i)*(handles.ufmean*(i-1) + squeeze( handles.data.uf(handles.data.i(i),handles.data.j(i),1:end)));
        handles.vfmean = (1/i)*(handles.vfmean*(i-1) + squeeze( handles.data.vf(handles.data.i(i),handles.data.j(i),1:end)));
        handles.wfmean = (1/i)*(handles.wfmean*(i-1) + squeeze( handles.data.wf(handles.data.i(i),handles.data.j(i),1:end)));
    end
    
    switch get(handles.velocity_component,'Value')
        case 1 % u
            handles.property  = handles.data.u(:,:,1:end-1); % we do not need the average
            handles.mean_property = handles.umean;
        case 2 % v
            handles.property = handles.data.v(:,:,1:end-1);
            handles.mean_property = handles.vmean;
        case 3 % w
            handles.property = handles.data.w(:,:,1:end-1);
            handles.mean_property = handles.wmean;
        case 4 % uf
            handles.property  = handles.data.uf;
            handles.mean_property = handles.ufmean;
        case 5 % vf
            handles.property  = handles.data.vf;
            handles.mean_property = handles.vfmean;
        case 6 % wf
            handles.property  = handles.data.wf;
            handles.mean_property = handles.wfmean;
    end
    
    guidata(handles.fig, handles);
    update_time_axis(handles.fig,[],handles);
    
    if nargout > 0
        varargout{1} = fig;
    end
    
    set(fig, 'Visible', 'on')
    
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
    
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


% -------------------------------------------------------------------
function varargout = update_time_axis(hObject, eventdata, handles)
% Time series plot
axes(handles.time_axes)
delete(get(handles.time_axes,'Children'));

if get(handles.spatial_average,'Value') == 1
    handles.time_average_plot = plot(handles.t, handles.mean_property,'Linewidth',2);
end 
if get(handles.checkbox_series,'Value') == 1
    set(handles.time_axes,'NextPlot','add');
    for i = 1:length(handles.data.i)
        handles.time_plots(i) = plot( handles.t,...
            squeeze( handles.property(handles.data.i(i),handles.data.j(i),:) ), handles.linestyle{mod(i,16)+1});
    end
    %     set(handles.time_axes,'NextPlot','replace');
end
set(handles.time_axes,'NextPlot','replace');
grid on
set(get(handles.time_axes,'xlabel'),'string','Time [s]');        
set(get(handles.time_axes,'ylabel'),'string','Velocity [m/s]');
guidata(hObject, handles);





% -------------------------------------------------------------------------
function varargout = update_spectrum_axis(hObject, eventdata, handles)

warning off

if get(handles.popupmenu_quantity,'Value') > 1
    % Frequency domain plot
    axes(handles.spectrum_axes);
    delete(get(handles.spectrum_axes,'Children'));
    
    switch get(handles.popupmenu_quantity,'Value')
        case {5}
            [cs,h] = contour(handles.spectrum_yaxis,handles.spectrum_xaxis,handles.quantity,'k-');
            set(h(find(cell2mat(get(h(:),'userdata'))< 0)),'linestyle',':');
            % clabel(cs,h,'fontsize',10,'color','r','rotation',0,'labelspacing',200)
        case {6,7,8} 
            handles.spectrum_plot = loglog(handles.spectrum_xaxis,2*handles.quantity);
        case {2,3,4}
            handles.spectrum_plot = plot(handles.spectrum_xaxis,2*handles.quantity);
        otherwise
            errodlg('Wrong selection');
    end 
    set(get(handles.spectrum_axes,'xlabel'),'string',handles.spectrum_xlabel);        
    set(get(handles.spectrum_axes,'ylabel'),'string',handles.spectrum_ylabel);
    grid on
    
end
guidata(hObject, handles);

% --------------------------------------------------------------------
function varargout = velocity_component_Callback(hObject, eventdata, handles)
% Change velocity_component
% handles.property = repmat(0,[size(handles.data.i),size(handles.data.u,3)]);
    switch get(handles.velocity_component,'Value')
        case 1 % u
            handles.property  = handles.data.u(:,:,1:end-1); % we do not need the average
            handles.mean_property = handles.umean;
        case 2 % v
            handles.property = handles.data.v(:,:,1:end-1);
            handles.mean_property = handles.vmean;
        case 3 % w
            handles.property = handles.data.w(:,:,1:end-1);
            handles.mean_property = handles.wmean;
        case 4 % uf
            handles.property  = handles.data.uf;
            handles.mean_property = handles.ufmean;
        case 5 % vf
            handles.property  = handles.data.vf;
            handles.mean_property = handles.vfmean;
        case 6 % wf
            handles.property  = handles.data.wf;
            handles.mean_property = handles.wfmean;
    end
tmp = cellstr(get(hObject,'String'));
handles.property_name = tmp{get(hObject,'Value')};
guidata(hObject,handles);
update_time_axis(hObject,[],handles);




% --------------------------------------------------------------------
function edit_nfft_Callback(hObject, eventdata, handles)
hanldes.nfft = str2num(get(hObject,'String'));
guidata(hObject,handles);
% update_gui(hObject,handles);

% --------------------------------------------------------------------
function varargout = quit_Callback(hObject, eventdata, handles)

if strcmpi(get(hObject,'Type'),'figure'),
    fig = hObject;
else
    fig = get(hObject,'Parent');
end

delete(fig);


% --------------------------------------------------------------------
function spatial_average_Callback(hObject, eventdata, handles)

if get(hObject,'Value') == 0
    set(handles.time_average_plot,'Visible','off')
end
guidata(hObject,handles);
update_time_axis(hObject,[],handles);


% --------------------------------------------------------------------
% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(get(hObject,'Parent'));


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function distribution_Callback(hObject, eventdata, handles)
% hObject    handle to distribution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function timebox_v3d_Callback(hObject, eventdata, handles)
% hObject    handle to timebox_v3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function time_axes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate time_axes


% --- Executes on button press in checkbox_series.
function checkbox_series_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_series (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_series

update_time_axis(hObject,[],handles);


% --- Executes on selection change in popupmenu_quantity.
function calculate_quantity(hObject, eventdata, handles)
% hObject    handle to popupmenu_quantity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.nfft        = fix(str2num(get(handles.edit_nfft,'String')));
handles.noverlap    = fix(str2num(get(handles.edit_overlap,'String')));
handles.method      = get(handles.popupmenu_method,'value');
handles.detrend     = get(handles.popupmenu_detrend,'value') - 1; % switch to 0,1,2. see psd_est()

switch get(handles.popupmenu_quantity,'Value')
    case 1 % none
        
    case 2 % autocorrelation in time
        [tmp,lags] = crosscorr(squeeze(handles.property(handles.data.i(1),handles.data.j(1),:)));    % time direction
        k = 1;
        for i = 2:length(handles.data.i)
            k = k + 1;
            tmp = 1/k*(tmp*(k-1) + crosscorr(squeeze(handles.property(handles.data.i(i),handles.data.j(i),:))));
        end
        i = find(lags > 0);
        handles.quantity = tmp(i);
        handles.spectrum_xaxis = lags(i)*handles.data.dt;
        handles.spectrum_xlabel = 'Time lag [sec]';
        handles.spectrum_ylabel = 'R_{ii}(t)';
    case 3 % autocorr in x
        uniqueI = unique(handles.data.i);           % unique rows
        indx = find(handles.data.i == uniqueI(1));
        [tmp,lags] = crosscorr(handles.property(uniqueI(1),handles.data.j(indx),1));
        k = 1;
        for i = 2:length(uniqueI)                     % for each unique row
            indx = find(handles.data.i == uniqueI(i));
            for j = 1:handles.data.N
                k = k + 1;
                tmp = 1/k*(tmp*(k-1)+ crosscorr(handles.property(uniqueI(i),handles.data.j(indx),j) ) );
            end
        end
        i = find(lags > 0);
        handles.quantity = tmp(i);
        handles.spectrum_xaxis = lags(i)*handles.data.scale*handles.data.gridX;
        handles.spectrum_xlabel = 'Lag [mm]';
        handles.spectrum_ylabel = 'R_{ii}(x)';
        
    case 4 % autocorrelation in y
        uniqueJ = unique(handles.data.j);           % unique rows
        indx = find(handles.data.j == uniqueJ(1));
        [tmp,lags] = crosscorr(handles.property(handles.data.i(indx),uniqueJ(1),1));
        k = 1;
        for i = 2:length(uniqueJ)                     % for each unique row
            indx = find(handles.data.j == uniqueJ(1));
            for j = 1:handles.data.N
                k = k + 1;
                tmp = 1/k*(tmp*(k-1) + crosscorr(handles.property(handles.data.i(indx),uniqueJ(i),j) ) );
            end
        end
        i = find(lags > 0);
        handles.quantity = tmp(i);
        handles.spectrum_xaxis = lags(i)*handles.data.scale*handles.data.gridY;
        handles.spectrum_xlabel = 'Lag [mm]';
        handles.spectrum_ylabel = 'R_{ii}(y)';
        
    case 5 % cross-correlation
        ind = sub2ind(size(handles.property),handles.data.i,handles.data.j,ones(size(handles.data.i)));
        %        Nice, but not correct.
        tmp = reshape(handles.property(ind),handles.rows,handles.cols);
        %        tmp = repmat(0,[length(handles.ind),1]
        handles.quantity = corrcoef(tmp);     % 2-point covariance map;
        for i = 1:handles.data.N
            ind = sub2ind(size(handles.property),handles.data.i,handles.data.j,i*ones(size(handles.data.i)));
            tmp = reshape(handles.property(ind),handles.rows,handles.cols);
            handles.quantity = 1/i*(handles.quantity*(i-1) + corrcoef(tmp));     % averaging
        end
        handles.spectrum_xaxis = (0:size(handles.quantity,1)-1)*handles.data.scale*handles.data.gridX;
        handles.spectrum_yaxis = (0:size(handles.quantity,2)-1)*handles.data.scale*handles.data.gridY;
        handles.spectrum_xlabel = 'Lag [mm]';
        handles.spectrum_ylabel = 'Lag [mm]';
        
    case 6 % spectrum in t
        switch handles.method
            case 1 % fft(corr)
                [tmp,lags] = crosscorr(squeeze(handles.property(handles.data.i(1),handles.data.j(1),:)));    % time direction
                for i = 2:length(handles.data.i)
                    tmp = 1/i*(tmp*(i-1) + crosscorr(squeeze(handles.property(handles.data.i(i),handles.data.j(i),:))));
                end
                i = find(lags > 0);
                w = windowvector(length(tmp(i)),get(handles.window,'Value'));
                tmp(i) = tmp(i).*w;
                handles.quantity = real(fft(tmp(i),handles.nfft))./(w'*w);
                handles.spectrum_xaxis = 2*pi*(0:ceil((handles.nfft+1)/2-1))*handles.Fs/handles.nfft;
                handles.quantity = handles.quantity(1:ceil((handles.nfft+1)/2));
            case 2 % psd
                for i = 1:length(handles.data.i)
                    tmp = handles.property(handles.data.i(i),handles.data.j(i),:);
                    w = windowvector(length(tmp(:)),get(handles.window,'Value'));
                    [psd,handles.spectrum_xaxis] = psd_est(tmp,handles.nfft,w,handles.Fs,handles.detrend);
                    handles.quantity = 1/i*(handles.quantity*(i-1) + psd);
                end
        end
        
        handles.spectrum_xlabel = 'Frequency [1/sec]';
        handles.spectrum_ylabel = 'E_{ii}(f)';
    case 7 % spectrum in k1
        switch handles.method
            case 1 % fft(corr)
                uniqueI = unique(handles.data.i);           % unique rows
                indx = find(handles.data.i == uniqueI(1));
                [tmp,lags] = crosscorr(handles.property(uniqueI(1),handles.data.j(indx),1));
                k = 1;
                for i = 2:length(uniqueI)                     % for each unique row
                    indx = find(handles.data.i == uniqueI(i));
                    for j = 1:handles.data.N
                        k = k + 1;
                        tmp = 1/k*(tmp*(k-1) + crosscorr(handles.property(uniqueI(i),handles.data.j(indx),j) ) );
                    end
                end
                i = find(lags > 0);
                w = windowvector(length(tmp(i)),get(handles.window,'Value'));
                tmp(i) = tmp(i).*w(:)';
                handles.quantity = real(fft(tmp(i),handles.nfft))./(w'*w);
                handles.spectrum_xaxis = 2*pi*(0:ceil((handles.nfft+1)/2-1))/(handles.data.scale*handles.data.gridX)/handles.nfft;
                handles.quantity = handles.quantity(1:ceil((handles.nfft+1)/2));
                
            case 2 % psd
                uniqueI = unique(handles.data.i);           % unique rows
                indx = find(handles.data.i == uniqueI(1));
                w = windowvector(length(indx),get(handles.window,'Value'));
                [handles.quantity,handles.spectrum_xaxis] = psd_est(handles.property(uniqueI(1),handles.data.j(indx),1),...
                    handles.nfft,w,handles.Fs,handles.detrend);
                k = 1;
                for i = 2:length(uniqueI)                     % for each unique row
                    indx = find(handles.data.i == uniqueI(i));
                    for j = 1:handles.data.N
                        k = k + 1;
                        handles.quantity = 1/k*(handles.quantity*(k-1) + psd_est(handles.property(uniqueI(i),handles.data.j(indx),j),...
                            handles.nfft,w,handles.Fs,handles.detrend));
                    end
                end
        end
        handles.spectrum_xlabel = 'k_1 [rad/s]';
        handles.spectrum_ylabel = 'E_{ii}(k_1)';
        
    case 8 % spectrum in k2
        switch handles.method
            case 1 % fft(corr)
                
                uniqueJ = unique(handles.data.j);           % unique rows
                indx = find(handles.data.j == uniqueJ(1));
                [tmp,lags] = crosscorr(handles.property(handles.data.i(indx),uniqueJ(1),1));
                k = 1;
                for i = 2:length(uniqueJ)                     % for each unique row
                    indx = find(handles.data.j == uniqueJ(1));
                    for j = 1:handles.data.N
                        k = k + 1;
                        tmp = 1/k*(tmp*(k-1) + crosscorr(handles.property(handles.data.i(indx),uniqueJ(i),j) ) );
                    end
                end
                i = find(lags > 0);
                w = windowvector(length(tmp(i)),get(handles.window,'Value'));
                tmp(i) = tmp(i).*w(:);
                handles.quantity = real(fft(tmp(i),handles.nfft))./(w'*w);
                handles.spectrum_xaxis = 2*pi*(0:ceil((handles.nfft+1)/2-1))/(handles.data.scale*handles.data.gridY)/handles.nfft;
                handles.quantity = handles.quantity(1:ceil((handles.nfft+1)/2));
                
            case 2 % psd
                
                uniqueJ = unique(handles.data.j);           % unique rows
                indx = find(handles.data.j == uniqueJ(1));
                w = windowvector(length(indx),get(handles.window,'Value'));
                
                [handles.quantity,handles.spectrum_xaxis] = psd_est(handles.property(handles.data.i(indx),uniqueJ(1),1),...
                    handles.nfft,w,handles.Fs,handles.detrend);
                k = 1;
                for i = 2:length(uniqueJ)                     % for each unique row
                    indx = find(handles.data.j == uniqueJ(1));
                    for j = 1:handles.data.N
                        k = k + 1;
                        handles.quantity = 1/k*(handles.quantity*(k-1) + psd_est(handles.property(handles.data.i(indx),uniqueJ(i),j),...
                            handles.nfft,w,handles.Fs,handles.detrend));
                    end
                end
        end
        
        handles.spectrum_xlabel = 'k_2 [rad/s]';
        handles.spectrum_ylabel = 'E_{ii}(k_2)';

end
guidata(hObject,handles);
% update_gui(hObject,[],handles);



% --- Executes on selection change in Window.
function window_Callback(hObject, eventdata, handles)
% guidata(hObject,handles);
% popupmenu_quantity_Callback(hObject,[],handles);

% ------------------------------------------------------------------------
function [w] = windowvector(n,window,varargin)
% WINDOWVECTOR returns the vector of the
% WINDOW of the proper length N, of one of
% the types: 'boxcar','hamming','hanning',
% 'blackman', numbered in Window Listbox:
% 1 - boxcar
% 2 - hamming
% 3 - hanning
% 4 - blackman
% 5 - bartlett
% 
% Based on Mathwork's HAMMING, HANNING, BOXCAR, ...

if window == 1 % boxcar
    w = ones(n,1);
elseif window  == 5 % bartlett
    w = 2*(0:(n-1)/2)/(n-1);
    if ~rem(n,2) % even
        w = [w w(n/2:-1:1)]';
    else % odd
        w = [w w((n-1)/2:-1:1)]';
    end
else
    switch window
        case 2 % 'hamming'
            % Hamming window: w = (54 - 46*cos(2*pi*(0:m-1)'/(n-1)))/100;
            a0 = 0.54;
            a1 = 0.46;
            a2 = 0;
            a3 = 0;
            a4 = 0;
        case 3 % 'hanning'
            % Hanning window: w = 0.5 * (1 - cos(2*pi*(0:m-1)'/(n-1))); 
            a0 = 0.5;
            a1 = 0.5;
            a2 = 0;
            a3 = 0;
            a4 = 0; 
        case 4 % 'blackman'
            % Blackman window: w = (42 - 50*cos(2*pi*(0:m-1)/(n-1)) + 8*cos(4*pi*(0:m-1)/(n-1)))'/100;
            a0 = 0.42;
            a1 = 0.5;
            a2 = 0.08;
            a3 = 0;
            a4 = 0;
    end
    
    if ~rem(n,2)
        x = (0:n/2-1)'/(n-1);
        w = a0 - a1*cos(2*pi*x) + a2*cos(4*pi*x) - a3*cos(6*pi*x) + a4*cos(8*pi*x);
        w = [w; w(end:-1:1)];
    else % Odd length window
        x = (0:(n+1)/2-1)'/(n-1);
        w = a0 - a1*cos(2*pi*x) + a2*cos(4*pi*x) - a3*cos(6*pi*x) + a4*cos(8*pi*x);
        w = [w; w(end-1:-1:1)];
    end
end
% [EOF] windowvector


% --- Executes on button press in pushbutton_calc.
function pushbutton_calc_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text_tip,'String','Calculating ...');
calculate_quantity(hObject,[],guidata(handles.fig));
set(handles.text_tip,'String','Plotting ...');
% handles = guihandles(hObject);
update_spectrum_axis(hObject,[],guidata(handles.fig));
set(handles.text_tip,'String','Ready ...');
guidata(handles.fig,handles);





%--------------------------------------------------------------------------
% --- Executes on selection change in popupmenu_quantity.
function popupmenu_quantity_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_quantity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%--------------------------------------------------------------------------
function [c,lags] = crosscorr(a,b)
% Calculates Cross-Correlation with
% the Mathwork's built-in filter() function
% Based on CONV
%
if nargin == 1
    b = a;
end
na = length(a);
a(2*na-1) = 0;
c = filter(b(end:-1:1), 1, a);
lags = [-(na-1):na-1];



%--------------------------------------------------------------------------
function [pxx,freq_vector] = psd_est(x,nfft,window,noverlap,Fs,detrendflag)
% PSD_EST Power Spectral Density estimate, based on Mathworks' PSD.M
% The method is modified periodogram, Welch method.
% 
%
%   Pxx = PSD(X,NFFT,Fs,WINDOW) estimates the Power Spectral Density of 
%   a discrete-time signal vector X using Welch's averaged, modified
%   periodogram method.  
%   
%   X is divided into overlapping sections, each of which is detrended 
%   (according to the detrending flag, if specified), then windowed by 
%   the WINDOW parameter, then zero-padded to length NFFT.  The magnitude 
%   squared of the length NFFT DFTs of the sections are averaged to form
%   Pxx.  Pxx is length NFFT/2+1 for NFFT even, (NFFT+1)/2 for NFFT odd,
%   or NFFT if the signal X is complex.  If you specify a scalar for 
%   WINDOW, a Hanning window of that length is used.  Fs is the sampling
%   frequency which doesn't affect the spectrum estimate but is used 
%   for scaling the X-axis of the plots.
%
%   [Pxx,F] = PSD(X,NFFT,Fs,WINDOW,NOVERLAP) returns a vector of frequen-
%   cies the same size as Pxx at which the PSD is estimated, and overlaps
%   the sections of X by NOVERLAP samples.
%
%
%   detrendflag = 0 -  'none', 1 - 'linear', 2 - 'mean'
%   
%   PSD with no output arguments plots the PSD in the current figure window,
%   with confidence intervals if you provide the P parameter.
%
%   The default values for the parameters are NFFT = 256 (or LENGTH(X),
%   whichever is smaller), NOVERLAP = 0, WINDOW = HANNING(NFFT), Fs = 2, 
%   P = .95, and DFLAG = 'none'.  You can obtain a default parameter by 
%   leaving it off or inserting an empty matrix [], e.g. PSD(X,[],10000).
%
%   NOTE: For Welch's method implementation which scales by the sampling
%         frequency, 1/Fs, see PWELCH.
%
%   See also PWELCH, CSD, COHERE, TFE.
%   ETFE, SPA, and ARX in the System Identification Toolbox.

%   Author(s): T. Krauss, 3-26-93
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.12 $  $Date: 2002/03/28 17:30:13 $

%   NOTE 1: To express the result of PSD, Pxx, in units of
%           Power per Hertz multiply Pxx by 1/Fs [1].
%
%   NOTE 2: The Power Spectral Density of a continuous-time signal,
%           Pss (watts/Hz), is proportional to the Power Spectral 
%           Density of the sampled discrete-time signal, Pxx, by Ts
%           (sampling period). [2] 
%       
%               Pss(w/Ts) = Pxx(w)*Ts,    |w| < pi; where w = 2*pi*f*Ts

%   References:
%     [1] Petre Stoica and Randolph Moses, Introduction To Spectral
%         Analysis, Prentice hall, 1997, pg, 15
%     [2] A.V. Oppenheim and R.W. Schafer, Discrete-Time Signal
%         Processing, Prentice-Hall, 1989, pg. 731
%     [3] A.V. Oppenheim and R.W. Schafer, Digital Signal
%         Processing, Prentice-Hall, 1975, pg. 556

% Default values
if nargin < 2 | isempty(nfft)
    nfft = min(length(x),256);
end
if nargin < 3 | isempty(window)
    window = hanning(nfft);
end 
if nargin < 4 | isempty(noverlap)
    noverlap = 0;
end
if nargin < 5 | isempty(Fs)
    Fs = 2;
end
if nargin < 6 | isempty(detrendflag)
    detrendflag = 0;
end

% compute PSD
window = window(:);
n = length(x);		    
nwind = length(window); 
% zero-pad x if it has length less than the window length
if n < nwind            
    x(nwind)=0;  n=nwind;
end
x = x(:);		

% Number of overlaps
k = fix((n-noverlap)/(nwind-noverlap));	

index = 1:nwind;

% Normalizing scale factor ==> asymptotically unbiased
scale = k*norm(window)^2;
% scale = k*sum(window)^2;

pxx = zeros(nfft,1);
switch detrendflag
    case 0 % no detrend
        for i=1:k
            xw = window.*(x(index));
            index = index + (nwind - noverlap);
            pxx = pxx + abs(fft(xw,nfft)).^2;
        end
    case 1 % detrend
        for i=1:k
            xw = window.*detrend(x(index));
            index = index + (nwind - noverlap);
            pxx = pxx + abs(fft(xw,nfft)).^2;
        end
    case 2 % remove mean
        for i=1:k
            xw = window.*detrend(x(index),'constant');
            index = index + (nwind - noverlap);
            pxx = pxx + abs(fft(xw,nfft)).^2;
        end
end


if rem(nfft,2),    % nfft odd
    indx = (1:(nfft+1)/2)';
else
    indx = (1:nfft/2+1)';
end
pxx = pxx(indx);
freq_vector = (indx-1)*Fs/nfft;
% [EOF] psd_est


% --- Executes during object creation, after setting all properties.
function edit_overlap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit_overlap_Callback(hObject, eventdata, handles)
% hObject    handle to edit_overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_overlap as text
%        str2double(get(hObject,'String')) returns contents of edit_overlap as a double


% --- Executes during object creation, after setting all properties.
function popupmenu_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu_method.
function popupmenu_method_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_method


% --- Executes during object creation, after setting all properties.
function popupmenu_detrend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_detrend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu_detrend.
function popupmenu_detrend_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_detrend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_detrend contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_detrend


% --------------------------------------------------------------------
function export_timeplot_Callback(hObject, eventdata, handles)
% hObject    handle to export_timeplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hl = findobj(handles.time_axes,'type','line');
if ~isempty(hl)
    xd = get(hl,'xdata');
    yd = get(hl,'ydata');
    if isequal(xd{1},xd{2})
        matr(:,1) = [xd{1}]';
        for i = 1:length(hl), 
            matr(:,i+1) = [yd{i}]';
        end
    elseif isequal(yd{1},yd{2})
        matr(:,1) = [yd{1}]';
        for i = 1:length(hl), 
            matr(:,i+1) = [xd{i}]';
        end
    else
        for i = 1:2:length(hl), 
            matr(:,i) = [xd{i}]';
            matr(:,i+1) = [yd{i}]';
        end
    end
    
    file = [];
    file = inputdlg('File Name','Input Name for CSV File');   
    if ~isempty (file)
        csvwrite(file{1},matr);  
    else 
        errordlg ('Choose a valid file name !!! ');
    end;
end;

% --------------------------------------------------------------------
function export_spectrum_Callback(hObject, eventdata, handles)
% hObject    handle to export_spectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hl = findobj(handles.spectrum_axes,'type','line');
if ~isempty(hl)
    xd = get(hl,'xdata');
    yd = get(hl,'ydata');
    if isequal(xd{1},xd{2})
        matr(:,1) = [xd{1}]';
        for i = 1:length(hl), 
            matr(:,i+1) = [yd{i}]';
        end
    elseif isequal(yd{1},yd{2})
        matr(:,1) = [yd{1}]';
        for i = 1:length(hl), 
            matr(:,i+1) = [xd{i}]';
        end
    else
        for i = 1:2:length(hl), 
            matr(:,i) = [xd{i}]';
            matr(:,i+1) = [yd{i}]';
        end
    end
    
    file = [];
    file = inputdlg('File Name','Input Name for CSV File');   
    if ~isempty (file)
        csvwrite(file{1},matr);  
    else 
        errordlg ('Choose a valid file name !!! ');
    end;
end;


% --- Executes during object creation, after setting all properties.
function edit_winlen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_winlen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit_winlen_Callback(hObject, eventdata, handles)
% hObject    handle to edit_winlen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_winlen as text
%        str2double(get(hObject,'String')) returns contents of edit_winlen as a double


