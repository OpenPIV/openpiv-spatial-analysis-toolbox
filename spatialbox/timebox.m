function varargout = timebox(varargin)
%TIMEBOX First prototype of the time toolbox
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
% Last modified: May 5, 2004, 01:31AM
% - all units are in m/s, therefore 2pi should work
% window length is introduced, as well as detrend and overlaping
% Last modified: June 21, 2004
% handles.data.dx and handles.data.dy are replaced with abs()
% due to the bug in SpatialBox
% 

if strcmp(inputname(1),'handles') 
    % GUI
    
    
    fig = openfig(mfilename,'new');
    movegui(fig, 'center');
    
    % Generate a structure of handles to pass to callbacks, and store it. 
    handles = guihandles(fig);
    
    
    set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
    
    
    
    
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
    
    % Prevent two point correlation if rows or cols are selected
    if length(handles.ind) ~= handles.rows*handles.cols
        %    set(handles.popupmenu_quantity,'string','---|Rii (t)|Rii (r 1)|Rii (r 2)|Rij (r 1, r 2)|E (f)|E (k 1)|E (k 2)');
        set(handles.popupmenu_quantity,'string','---|Rii (t)|Rii (r 1)|Rii (r 2)|---|E (f)|E (k 1)|E (k 2)');
    end
    
    % Calculate u and v average velocities over spatial regions
    handles.umean = squeeze( handles.data.u(handles.data.i(1),handles.data.j(1),1:end-1) );
    handles.vmean = squeeze( handles.data.v(handles.data.i(1),handles.data.j(1),1:end-1) );
    handles.ufmean = squeeze( handles.data.uf(handles.data.i(1),handles.data.j(1),1:end) );
    handles.vfmean = squeeze( handles.data.vf(handles.data.i(1),handles.data.j(1),1:end) );
    
    for i = 2:length(handles.data.i)
        handles.umean  = (1/i)*(handles.umean*(i-1) + squeeze( handles.data.u(handles.data.i(i),handles.data.j(i),1:end-1)));
        handles.vmean  = (1/i)*(handles.vmean*(i-1) + squeeze( handles.data.v(handles.data.i(i),handles.data.j(i),1:end-1)));
        handles.ufmean = (1/i)*(handles.ufmean*(i-1) + squeeze( handles.data.uf(handles.data.i(i),handles.data.j(i),1:end)));
        handles.vfmean = (1/i)*(handles.vfmean*(i-1) + squeeze( handles.data.vf(handles.data.i(i),handles.data.j(i),1:end)));
    end
    
    if handles.data.state3d
        handles.wmean = squeeze( handles.data.w(handles.data.i(1),handles.data.j(1),1:end-1) );
        handles.wfmean = squeeze( handles.data.wf(handles.data.i(1),handles.data.j(1),1:end) );
        for i = 2:length(handles.data.i)
            handles.wmean  = (1/i)*(handles.wmean*(i-1) + squeeze( handles.data.w(handles.data.i(i),handles.data.j(i),1:end-1)));
            handles.wfmean = (1/i)*(handles.wfmean*(i-1) + squeeze( handles.data.wf(handles.data.i(i),handles.data.j(i),1:end)));
        end
    end
    
    if handles.data.state3d
        set(handles.velocity_component,'String','u|v|w|u''|v''|w''');
    else
        set(handles.velocity_component,'String','u|v|u''|v''');
    end
    
    % Default values
    handles.property  = handles.data.u(:,:,1:end-1); % we do not need the average
    handles.mean_property = handles.umean;
    
    
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

set(handles.fig,'Pointer','watch');
set(handles.text_tip,'String','Plotting ...');
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

set(handles.text_tip,'String','Ready ...');
set(handles.fig,'Pointer','arrow');

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
            if ~isempty(h)
                set(h(find(cell2mat(get(h(:),'userdata'))< 0)),'linestyle',':');
            end
            % clabel(cs,h,'fontsize',10,'color','r','rotation',0,'labelspacing',200)
        case {6,7,8} 
            handles.spectrum_plot = loglog(handles.spectrum_xaxis,handles.quantity);
        case {2,3,4}
            handles.spectrum_plot = plot(handles.spectrum_xaxis,handles.quantity);
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

if handles.data.state3d
    switch get(handles.velocity_component,'Value')
        case 1 % u
            handles.property  = handles.data.u(:,:,1:end-1); % we do not need the average
            handles.mean_property = handles.umean;
            set(handles.pushbutton_calc,'Enable','On');
            
        case 2 % v
            handles.property = handles.data.v(:,:,1:end-1);
            handles.mean_property = handles.vmean;
            set(handles.pushbutton_calc,'Enable','On');
        case 3 % w
            handles.property = handles.data.w(:,:,1:end-1);
            handles.mean_property = handles.wmean;
            set(handles.pushbutton_calc,'Enable','Off');
        case 4 % uf
            handles.property  = handles.data.uf;
            handles.mean_property = handles.ufmean;
            set(handles.pushbutton_calc,'Enable','On');
        case 5 % vf
            handles.property  = handles.data.vf;
            handles.mean_property = handles.vfmean;
        case 6 % wf
            handles.property  = handles.data.wf;
            handles.mean_property = handles.wfmean;
            set(handles.pushbutton_calc,'Enable','Off');
    end
else
    set(handles.pushbutton_calc,'Enable','On');
    
    switch get(handles.velocity_component,'Value')
        case 1 % u
            handles.property  = handles.data.u(:,:,1:end-1); % we do not need the average
            handles.mean_property = handles.umean;
        case 2 % v
            handles.property = handles.data.v(:,:,1:end-1);
            handles.mean_property = handles.vmean;
        case 3 % uf
            handles.property  = handles.data.uf;
            handles.mean_property = handles.ufmean;
        case 4 % vf
            handles.property  = handles.data.vf;
            handles.mean_property = handles.vfmean;
    end
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
function timebox_Callback(hObject, eventdata, handles)
% hObject    handle to timebox (see GCBO)
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
        i = find(lags >= 0);
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
        i = find(lags >= 0);
        handles.quantity = tmp(i)./tmp(i(1));   % normalize to 1.
        handles.spectrum_xaxis = lags(i)*abs(handles.data.dx);
        handles.spectrum_xlabel = 'Lag [m]';
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
        i = find(lags >= 0);
        handles.quantity = tmp(i)./tmp(i(1));   % normalize to 1
        handles.spectrum_xaxis = lags(i)*abs(handles.data.dy);
        handles.spectrum_xlabel = 'Lag [m]';
        handles.spectrum_ylabel = 'R_{ii}(y)';
        
    case 5 % cross-correlation
        ind = sub2ind(size(handles.property),handles.data.i,handles.data.j,ones(size(handles.data.i)));
        if length(handles.property(ind))~= handles.rows*handles.cols
            set(handles.text_tip,'String','Wrong selection !!!');
            pause(2)
            handles.quantity = zeros(handles.rows,handles.cols);
            handles.spectrum_xaxis = (0:size(handles.quantity,1)-1)*abs(handles.data.dx);
            handles.spectrum_yaxis = (0:size(handles.quantity,2)-1)*abs(handles.data.dy);
            handles.spectrum_xlabel = 'Lag [m]';
            handles.spectrum_ylabel = 'Lag [m]';
        else
            tmp = reshape(handles.property(ind),handles.rows,handles.cols);
            handles.quantity = crosscorr2d(tmp);     % 2-point covariance map;
%                         handles.quantity = corrcoef(tmp);     % 2-point covariance map;

            for i = 2:handles.data.N
                ind = sub2ind(size(handles.property),handles.data.i,handles.data.j,i*ones(size(handles.data.i)));
                tmp = reshape(handles.property(ind),handles.rows,handles.cols);
                handles.quantity = 1/i*(handles.quantity*(i-1) + crosscorr2d(tmp));     % averaging
%                                 handles.quantity = 1/i*(handles.quantity*(i-1) + corrcoef(tmp));     % averaging

            end
            handles.quantity = handles.quantity./max(max(handles.quantity));
            handles.spectrum_xaxis = (0:size(handles.quantity,1)-1)*abs(handles.data.dx);
            handles.spectrum_yaxis = (0:size(handles.quantity,2)-1)*abs(handles.data.dy);
            handles.spectrum_xlabel = 'Lag [m]';
            handles.spectrum_ylabel = 'Lag [m]';
        end
    case 6 % spectrum in t
        switch handles.method
            case 1 % fft(corr)
                [tmp,lags] = crosscorr(squeeze(handles.property(handles.data.i(1),handles.data.j(1),:)));    % time direction
                for i = 2:length(handles.data.i)
                    tmp = 1/i*(tmp*(i-1) + crosscorr(squeeze(handles.property(handles.data.i(i),handles.data.j(i),:))));
                end
                i = find(lags >= 0);
                tmp = 2*real(tmp(i))/length(tmp);
                w = windowvector(length(tmp),get(handles.window,'Value'));
                tmp = tmp(:).*w(:);
                handles.quantity = (1/pi/handles.Fs)*real(fft(tmp,handles.nfft))./(w'*w);
                handles.spectrum_xaxis = 2*pi*(0:ceil((handles.nfft+1)/2-1))*handles.Fs/handles.nfft;
                handles.quantity = handles.quantity(1:ceil((handles.nfft+1)/2));
                handles.quantity = handles.quantity(2:end);
                handles.spectrum_xaxis = handles.spectrum_xaxis(2:end);
            case 2 % psd
                %                 handles.quantity = 0;
                %                 for i = 1:length(handles.data.i)
                %                     tmp = handles.property(handles.data.i(i),handles.data.j(i),:);
                %                     winlen = str2double(get(handles.edit_winlen,'String'));
                %                     if winlen <= 0 & winlen > length(tmp), winlen = length(tmp(:)); end 
                %                     w = windowvector(winlen,get(handles.window,'Value'));
                %                     [pxx,handles.spectrum_xaxis] = psd_est(tmp,handles.nfft,w,handles.noverlap,handles.Fs,handles.detrend);
                %                     handles.quantity = 1/i*(handles.quantity*(i-1) + pxx);
                %                 end
        end
        handles.quantity = handles.quantity(2:end);
        handles.spectrum_xaxis = handles.spectrum_xaxis(2:end);
        handles.spectrum_xlabel = 'Frequency [rad/sec]';
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
                i = find(lags >= 0);
                tmp = 2*real(tmp(i))/length(tmp);
                w = windowvector(length(tmp),get(handles.window,'Value'));
                tmp = tmp(:).*w(:);
                handles.quantity = (1/pi/handles.Fs)*real(fft(tmp,handles.nfft))./(w'*w);
                handles.spectrum_xaxis = 2*pi*(0:ceil((handles.nfft+1)/2-1))/(abs(handles.data.dx))/handles.nfft;
                handles.quantity = handles.quantity(1:ceil((handles.nfft+1)/2));
                handles.quantity = handles.quantity(2:end);
                handles.spectrum_xaxis = handles.spectrum_xaxis(2:end);
                
            case 2 % psd
                %                 uniqueI = unique(handles.data.i);           % unique rows
                %                 indx = find(handles.data.i == uniqueI(1));
                %                 winlen = str2double(get(handles.edit_winlen,'String'));
                %                     if winlen <= 0 & winlen > length(indx), winlen = length(indx); end 
                %                     w = windowvector(winlen,get(handles.window,'Value'));
                % %                w = windowvector(length(indx),get(handles.window,'Value'));
                %                 [handles.quantity,handles.spectrum_xaxis] = psd_est(handles.property(uniqueI(1),handles.data.j(indx),1),...
                %                     handles.nfft,w,handles.noverlap,handles.Fs,handles.detrend);
                %                 k = 1;
                %                 for i = 2:length(uniqueI)                     % for each unique row
                %                     indx = find(handles.data.i == uniqueI(i));
                %                     for j = 1:handles.data.N
                %                         k = k + 1;
                %                         handles.quantity = 1/k*(handles.quantity*(k-1) + psd_est(handles.property(uniqueI(i),handles.data.j(indx),j),...
                %                             handles.nfft,w,handles.noverlap,handles.Fs,handles.detrend));
                %                     end
                %                 end
        end
        handles.quantity = handles.quantity(2:end);
        handles.spectrum_xaxis = handles.spectrum_xaxis(2:end);
        
        handles.spectrum_xlabel = 'k_1 [rad/m]';
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
                        tmp = 1/k*(tmp*(k-1) + crosscorr(handles.property(handles.data.i(indx),uniqueJ(i),j)));
                    end
                end
                i = find(lags >= 0);
                tmp = 2*real(tmp(i))/length(tmp);
                w = windowvector(length(i),get(handles.window,'Value'));
                tmp = tmp(:).*w(:);
                handles.quantity = (1/pi/handles.Fs)*real(fft(tmp,handles.nfft))./(w'*w);
                handles.spectrum_xaxis = 2*pi*(0:ceil((handles.nfft+1)/2-1))/(abs(handles.data.dy))/handles.nfft;
                handles.quantity = handles.quantity(1:ceil((handles.nfft+1)/2));
                handles.quantity = handles.quantity(2:end);
                handles.spectrum_xaxis = handles.spectrum_xaxis(2:end);
                
            case 2 % psd
                
                %                 uniqueJ = unique(handles.data.j);           % unique rows
                %                 indx = find(handles.data.j == uniqueJ(1));
                % %                 w = windowvector(length(indx),get(handles.window,'Value'));
                %                 winlen = str2double(get(handles.edit_winlen,'String'));
                %                 if winlen <= 0 & winlen > length(indx), winlen = length(indx); end 
                %                 w = windowvector(winlen,get(handles.window,'Value'));
                %                 [handles.quantity,handles.spectrum_xaxis] = psd_est(handles.property(handles.data.i(indx),uniqueJ(1),1),...
                %                     handles.nfft,w,handles.noverlap,handles.Fs,handles.detrend);
                %                 k = 1;
                %                 for i = 2:length(uniqueJ)                     % for each unique row
                %                     indx = find(handles.data.j == uniqueJ(1));
                %                     for j = 1:handles.data.N
                %                         k = k + 1;
                %                         handles.quantity = 1/k*(handles.quantity*(k-1) + psd_est(handles.property(handles.data.i(indx),uniqueJ(i),j),...
                %                             handles.nfft,w,handles.noverlap,handles.Fs,handles.detrend));
                %                     end
                %                 end
        end
        handles.quantity = handles.quantity(2:end);
        handles.spectrum_xaxis = handles.spectrum_xaxis(2:end);
        
        handles.spectrum_xlabel = 'k_2 [rad/m]';
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
set(handles.fig,'Pointer','watch');
calculate_quantity(hObject,[],guidata(handles.fig));
set(handles.text_tip,'String','Plotting ...');
% handles = guihandles(hObject);
update_spectrum_axis(hObject,[],guidata(handles.fig));
set(handles.text_tip,'String','Ready ...');
set(handles.fig,'Pointer','arrow');
guidata(handles.fig,handles);





%--------------------------------------------------------------------------
% --- Executes on selection change in popupmenu_quantity.
function popupmenu_quantity_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_quantity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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
    if ~iscell(xd) 
        matr(:,1)  =  xd(:);
        matr(:,2)  =  yd(:);
    elseif isequal(xd{1},xd{2})
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
    if ~iscell(xd) 
        matr(:,1)  =  xd(:);
        matr(:,2)  =  yd(:);
    elseif isequal(xd{1},xd{2})
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

% --- Executes during object creation, after setting all properties.
function export_Callback(hObject, eventdata, handles)



function edit_winlen_Callback(hObject, eventdata, handles)
% hObject    handle to edit_winlen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_winlen as text
%        str2double(get(hObject,'String')) returns contents of edit_winlen as a double


function [c,lags] = crosscorr(a,b)
% XCORRF - replaced CROSSCORR due to 
% slow performance with long time series
%
% Based on XCORRF2.M 
%
if nargin == 1
    b = a;
end
a = a(:);
b = b(:);
na = length(a);
nb = length(b);
% make reverse conjugate of one array
b = conj(b(nb:-1:1));
% use power of 2 transform lengths
nf = 2^nextpow2(na+nb);
% take Fourier transform with zero padding
at =  fft(a,nf);
bt =  fft(b,nf);
% multiply transforms then inverse transform
c = ifft(at.*bt);
%  trim to standard size (2*N-1)
c(na+nb:nf) = [];
lags = [-(na-1):na-1];


% --------------------------------------------------------------------
function export2csv_Callback(hObject, eventdata, handles)
% hObject    handle to export2csv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function export2matlab_Callback(hObject, eventdata, handles)
% hObject    handle to export2matlab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function timeplot2figure_Callback(hObject, eventdata, handles)
% hObject    handle to timeplot2figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.export_figure = figure;
% handles.export_axes   = axes;
copyobj(handles.time_axes,handles.export_figure);
% copyobj(get(handles.axes_main,'Children'),handles.export_axes);

set(handles.export_figure,'Units','normalized');
set(get(handles.export_figure,'children'),'Units','normalized');
set(get(handles.export_figure,'children'),'Position',[0.13 0.11 0.775 0.815]);
set(get(handles.export_figure,'children'),'Box','on');

guidata(handles.fig, handles);


% --------------------------------------------------------------------
function spectrum2figure_Callback(hObject, eventdata, handles)
% hObject    handle to spectrum2figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.export_figure = figure;
% handles.export_axes   = axes;
copyobj(handles.spectrum_axes,handles.export_figure);
% copyobj(get(handles.axes_main,'Children'),handles.export_axes);

set(handles.export_figure,'Units','normalized');
set(get(handles.export_figure,'children'),'Units','normalized');
set(get(handles.export_figure,'children'),'Position',[0.13 0.11 0.775 0.815]);
set(get(handles.export_figure,'children'),'Box','on');

guidata(handles.fig, handles);

% --------------------------------------------------------------------
function c = crosscorr2d(x,y)
if nargin == 1
	y = x;
end
c = conv2(x, y(end:-1:1,end:-1:1));


