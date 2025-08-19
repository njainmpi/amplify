%% signal_change_map_matlab.m
% Axial percent-change viewer (8 slices/page), controls at bottom panel
% - Underlay: Mean(4D) or Baseline mean (popup)
% - Overlay: ((Signal - Baseline)/Baseline)*100
% - Paging: Prev 8 / Next 8
% - Scale slider: symmetric ±v (0.1 .. 15 %)
% - Colormap: blue -> black -> hot (zero = black)
clear; clc;

%% ---- SETTINGS ----
in_nii       = 'cleaned_mc_func.nii.gz';
baseline_idx = 350:550;
signal_idx   = 1300:1500;
eps_baseline = 1e-6;
alpha_overlay = 0.6;
slices_per_page = 8;            % 8 slices (2x4 grid)

%% ---- LOAD 4D NIFTI (.nii/.nii.gz) ----
is_gz = endsWith(in_nii,'.gz','IgnoreCase',true);
tmp_nii = '';
if is_gz
    [tmpdir,~,~] = fileparts(tempname); if ~exist(tmpdir,'dir'), mkdir(tmpdir); end
    gunzip(in_nii,tmpdir);
    [~,b,~] = fileparts(in_nii); [~,b2,e2] = fileparts(b);
    tmp_nii = fullfile(tmpdir,[b2,e2]);
    niiPathToRead = tmp_nii;
else
    niiPathToRead = in_nii;
end
info = niftiinfo(niiPathToRead);
Y = double(niftiread(info));                 % X x Y x Z x T
[X,Ydim,Z,T] = size(Y);
fprintf('Loaded: %d x %d x %d x %d\n', X, Ydim, Z, T);

%% ---- MAPS ----
baseline_mean = mean(Y(:,:,:,baseline_idx),4);
signal_mean   = mean(Y(:,:,:,signal_idx),4);

den = baseline_mean; den(abs(den)<eps_baseline) = eps_baseline;
pc = ((signal_mean - baseline_mean) ./ den) * 100;     % 3D percent change
underlay_mean4D = mean(Y,4);                            % 3D

% Scale defaults
pc_flat = pc(isfinite(pc)); if isempty(pc_flat), pc_flat = 0; end
pc_lo = prctile(pc_flat,2); pc_hi = prctile(pc_flat,98);
pc_mag_default = max(5, max(abs([pc_lo pc_hi])));

% Underlay windowing
[und_lo_mean, und_hi_mean] = robust_window(underlay_mean4D);
[und_lo_base, und_hi_base] = robust_window(baseline_mean);

% Mask from 4D mean
try
    uflat = underlay_mean4D(isfinite(underlay_mean4D));
    u_lo = prctile(uflat,2); u_hi = prctile(uflat,98);
    thr = graythresh(mat2gray(underlay_mean4D));
    brainMask = underlay_mean4D > (u_lo + thr*(u_hi - u_lo));
catch
    brainMask = underlay_mean4D > 0;
end

%% ---- FIGURE LAYOUT (top: images, bottom: controls) ----
fig = figure('Color','w','Name','Axial %Change Viewer (8 per page)', ...
             'Units','normalized','Position',[0.06 0.06 0.88 0.86]);

% Panels: slices (top) and controls (bottom)
panelSlices   = uipanel('Parent',fig,'Units','normalized','Position',[0.03 0.18 0.94 0.79], ...
                        'BorderType','none');
panelControls = uipanel('Parent',fig,'Units','normalized','Position',[0.03 0.03 0.94 0.12], ...
                        'BorderType','none','BackgroundColor','w');

% Tiled grid for 8 slices
tl = tiledlayout(panelSlices,2,4,'Padding','compact','TileSpacing','compact');
ax  = gobjects(slices_per_page,1);
imU = gobjects(slices_per_page,1);
imP = gobjects(slices_per_page,1);

for k = 1:slices_per_page
    ax(k) = nexttile(tl,k); hold(ax(k),'on'); axis(ax(k),'image'); axis(ax(k),'off');
end

% Colormap + shared colorbar (attach to layout for recent MATLAB;
% for older MATLAB, fallback to figure-level colorbar)
cmap = gbhot_diverging(256);
colormap(fig, cmap);
try
    cb = colorbar(tl,'Location','southoutside'); % newer versions
catch
    cb = colorbar('southoutside');               % fallback
end
cb.Label.String = '% change';

% ---------- Controls (bottom panel) ----------
% Paging
btn_prev = uicontrol('Parent',panelControls,'Style','pushbutton', ...
    'Units','normalized','Position',[0.01 0.12 0.10 0.76], ...
    'String','◀ Prev 8','FontWeight','bold','Callback',@prevPage);
btn_next = uicontrol('Parent',panelControls,'Style','pushbutton', ...
    'Units','normalized','Position',[0.12 0.12 0.10 0.76], ...
    'String','Next 8 ▶','FontWeight','bold','Callback',@nextPage);
lbl_page = uicontrol('Parent',panelControls,'Style','text', ...
    'Units','normalized','Position',[0.23 0.12 0.20 0.76], ...
    'String','','BackgroundColor','w','FontWeight','bold','HorizontalAlignment','center');

% Underlay chooser
uicontrol('Parent',panelControls,'Style','text','Units','normalized', ...
    'Position',[0.44 0.58 0.12 0.32], 'String','Underlay:', ...
    'BackgroundColor','w','HorizontalAlignment','right','FontWeight','bold');
pop_under = uicontrol('Parent',panelControls,'Style','popupmenu','Units','normalized', ...
    'Position',[0.57 0.58 0.18 0.32], ...
    'String',{'Mean(4D)','Baseline mean'}, 'Callback',@changeUnderlay);

% Scale slider
uicontrol('Parent',panelControls,'Style','text','Units','normalized', ...
    'Position',[0.44 0.12 0.12 0.32], 'String','Scale (±%)', ...
    'BackgroundColor','w','HorizontalAlignment','right','FontWeight','bold');
sld_scale = uicontrol('Parent',panelControls,'Style','slider','Units','normalized', ...
    'Position',[0.57 0.12 0.18 0.32], 'Min',0.1, 'Max',15, ...
    'Value',min(15, max(1, pc_mag_default)), ...
    'SliderStep',[1/150 1/15], 'Callback',@adjustScale);
lbl_scale = uicontrol('Parent',panelControls,'Style','text','Units','normalized', ...
    'Position',[0.76 0.12 0.06 0.32], 'String',sprintf('%.1f',get(sld_scale,'Value')), ...
    'BackgroundColor','w','HorizontalAlignment','left');

% ---- Pack state ----
S.pc = pc; S.mask = brainMask;
S.U_mean = underlay_mean4D; S.U_base = baseline_mean;
S.und_lo_mean=und_lo_mean; S.und_hi_mean=und_hi_mean;
S.und_lo_base=und_lo_base; S.und_hi_base=und_hi_base;

S.Z = Z; S.spp = slices_per_page; S.startSlice = 1;
S.alpha_overlay = alpha_overlay;

S.fig = fig; S.tl = tl; S.ax = ax; S.imU = imU; S.imP = imP; S.cb = cb;
S.pop_under = pop_under; S.lbl_page = lbl_page;
S.sld_scale = sld_scale; S.lbl_scale = lbl_scale;

guidata(fig,S);

% Initial render + scale
renderPage(fig);
adjustScale(sld_scale,[]);

% Cleanup temp
if ~isempty(tmp_nii) && exist(tmp_nii,'file')
    try, delete(tmp_nii); catch, end
end

%% =========================
%% === CALLBACKS / FUNCS ===
%% =========================

function prevPage(src,~)
    fig = ancestor(src,'figure'); if isempty(fig), return; end
    S = guidata(fig); if ~isstruct(S), return; end
    S.startSlice = max(S.startSlice - S.spp, 1);
    guidata(fig,S); renderPage(fig);
end

function nextPage(src,~)
    fig = ancestor(src,'figure'); if isempty(fig), return; end
    S = guidata(fig); if ~isstruct(S), return; end
    S.startSlice = min(S.startSlice + S.spp, max(S.Z - S.spp + 1, 1));
    guidata(fig,S); renderPage(fig);
end

function changeUnderlay(src,~)
    fig = ancestor(src,'figure'); if isempty(fig), return; end
    renderPage(fig);
end

function adjustScale(src,~)
    fig = ancestor(src,'figure'); if isempty(fig), return; end
    S = guidata(fig); if ~isstruct(S), return; end
    v = get(S.sld_scale,'Value'); if v<=0, v=0.1; set(S.sld_scale,'Value',v); end
    for k = 1:numel(S.ax)
        if isgraphics(S.ax(k)), caxis(S.ax(k),[-v v]); end
    end
    S.cb.Limits = [-v v];
    S.cb.Ticks  = [-v 0 v];
    S.cb.TickLabels = {num2str(-v,'%.0f'),'0',num2str(v,'%.0f')};
    set(S.lbl_scale,'String',sprintf('%.1f',v));
    guidata(fig,S);
end

function renderPage(fig)
    S = guidata(fig); if ~isstruct(S), return; end

    % Underlay source + contrast
    if get(S.pop_under,'Value')==1
        U = S.U_mean; lo = S.und_lo_mean; hi = S.und_hi_mean;
    else
        U = S.U_base; lo = S.und_lo_base; hi = S.und_hi_base;
    end

    first = S.startSlice; last = min(first + S.spp - 1, S.Z);
    slices = first:last;
    set(S.lbl_page,'String',sprintf('Slices %d–%d of %d', first, last, S.Z));

    % (Re)draw tiles
    for k = 1:S.spp
        axk = S.ax(k);
        if k <= numel(slices)
            idx = slices(k);
            Uslc  = rot90(U(:,:,idx),-1);
            PCslc = rot90(S.pc(:,:,idx),-1);
            Mslc  = rot90(S.mask(:,:,idx),-1);

            if ~isgraphics(S.imU(k))
                axes(axk); cla(axk); hold(axk,'on'); axis(axk,'image'); axis(axk,'off');
                S.imU(k) = image(axk, to_rgb_gray(Uslc,lo,hi));
                S.imP(k) = imagesc(axk, PCslc);
                set(S.imP(k),'AlphaData', S.alpha_overlay * double(Mslc));
                title(axk, sprintf('Z = %d', idx));
            else
                set(S.imU(k),'CData', to_rgb_gray(Uslc,lo,hi));
                set(S.imP(k),'CData', PCslc, 'AlphaData', S.alpha_overlay * double(Mslc));
                title(axk, sprintf('Z = %d', idx));
            end
            set(axk,'Visible','on');
        else
            cla(axk); set(axk,'Visible','off');
        end
    end

    guidata(fig,S);
    % keep current scale
    adjustScale(S.sld_scale,[]);
end

%% ---- HELPERS ----
function [lo,hi]=robust_window(U)
    U = U(isfinite(U));
    if isempty(U), lo=0; hi=1;
    else, lo=prctile(U,2); hi=prctile(U,98); if hi<=lo, hi=lo+eps; end
    end
end

function RGB=to_rgb_gray(slice2d,lo,hi)
    if hi<=lo, hi=lo+eps; end
    N=(double(slice2d)-lo)/(hi-lo); N=max(0,min(1,N));
    RGB=repmat(N,[1 1 3]); % truecolor grayscale
end

function cmap = gbhot_diverging(n)
    % Green/Cyan/Blue -> Black -> Hot (red->yellow->white), zero is black.
    if mod(n,2)==1, n=n+1; end
    half = n/2;
    % Negative side: green->cyan->blue, darkening to black at zero
    t = linspace(0,1,half)';     % from -v to 0
    g = (1-t);                   % green fades
    c = t.*(1-t);                % cyan mid
    b = t;                       % blue increases
    neg = [0*g + 0*c, 1*g + 1*c, 0*g + 1*b]; % blend G/C->B
    neg = neg .* (1 - t);        % darken toward zero (black at center)
    % Positive side: MATLAB 'hot' from black->white
    pos = hot(half);
    cmap = [neg; pos];
end
