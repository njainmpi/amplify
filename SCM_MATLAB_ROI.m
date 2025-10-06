function SCM_MATLAB_ROI_Interactive
%% SCM_MATLAB_ROI_Interactive.m — Axial %Change Viewer (3D ROI + interactive UI)
% PSC policy (updated):
% - Compute PSC ONLY inside brainMask. Outside => NaN (never shown).
% - Replace near-zero baseline inside mask by a small, dynamic epsilon:
%   eps_den = 1e-4 * max(1, median(U_mean4D(brainMask)))
%   This keeps math finite without inventing huge % values.

clear; clc;

%% ---- SETTINGS ----
in_nii        = 'cleaned_mc_func.nii.gz';  % path or will prompt
baseline_idx  = 100:300;
signal_idx    = 1100:1300;
alpha_base    = 0.6;
default_tiles = 8;                         % 4 or 8

% Plot settings
PSC_YLIM        = [-5 15];   % hardcoded PSC y-limits
SHOW_VOXEL_PSC  = false;     % don't plot clicked-voxel PSC (ROI-only)

% Global UI style
UI_FONT     = 'Calibri';
UI_FONTSIZE = 18;            % minimum; helper enforces >= this

%% ---- LOAD 4D NIfTI (.nii/.nii.gz) ----
if isstring(in_nii), in_nii = char(in_nii); end
in_nii = strtrim(in_nii);
if ~isempty(in_nii)
    if startsWith(in_nii, ''''), in_nii = in_nii(2:end); end
    if ~isempty(in_nii) && endsWith(in_nii, ''''), in_nii = in_nii(1:end-1); end
    if startsWith(in_nii, '"'),  in_nii = in_nii(2:end); end
    if ~isempty(in_nii) && endsWith(in_nii, '"'),  in_nii = in_nii(1:end-1); end
end
if ~exist(in_nii,'file')
    warning('Input NIfTI not found: %s', in_nii);
    [f,p] = uigetfile({'*.nii;*.nii.gz','NIfTI files (*.nii, *.nii.gz)'}, ...
                      'Select your 4D functional NIfTI');
    if isequal(f,0), error('No NIfTI selected.'); end
    in_nii = fullfile(p,f);
end
is_gz   = endsWith(in_nii,'.gz','IgnoreCase',true);
tmp_nii = '';
if is_gz
    [tmpdir,~,~] = fileparts(tempname);
    if ~exist(tmpdir,'dir'), mkdir(tmpdir); end
    gunzip(in_nii,tmpdir);
    [~,b,~] = fileparts(in_nii);
    [~,b2,e2] = fileparts(b);
    tmp_nii = fullfile(tmpdir,[b2, e2]);
    niiPathToRead = tmp_nii;
else
    niiPathToRead = in_nii;
end

info = niftiinfo(niiPathToRead);
Y4D  = double(niftiread(info));  % X x Y x Z x T
[X,Ydim,Z,T] = size(Y4D);
fprintf('Loaded: %d x %d x %d x %d\n', X, Ydim, Z, T);

%% ---- GLOBAL UNDERLAY + MASK ----
U_mean4D = mean(Y4D,4);

% robust underlay windows for mean and baseline (computed later)
[und_lo_mean, und_hi_mean] = robust_window(U_mean4D);

% quick Otsu+percentile gating for a brain-like mask
try
    uflat = U_mean4D(isfinite(U_mean4D)); u_lo = prctile(uflat,2); u_hi = prctile(uflat,98);
    thr   = graythresh(mat2gray(U_mean4D));
    brainMask = U_mean4D > (u_lo + thr*(u_hi - u_lo));
catch
    brainMask = U_mean4D > 0;
end

% global median signal level inside mask (used for dynamic epsilon)
global_med = median(U_mean4D(brainMask),'omitnan');
if ~isfinite(global_med) || global_med<=0, global_med = 1; end
eps_den_global = 1e-4 * max(1, global_med);  % dynamic epsilon seed

%% ---- BASELINE / SIGNAL WINDOWS & PSC MAP (MASK-AWARE) ----
baseline_mean = mean(Y4D(:,:,:,baseline_idx),4,'omitnan');
signal_mean   = mean(Y4D(:,:,:,signal_idx),4,'omitnan');

% denominator handling (inside mask only)
den = baseline_mean;
tiny = abs(den) < eps_den_global;
den(tiny & brainMask) = eps_den_global;

pc = ((signal_mean - baseline_mean) ./ den) * 100;
pc(~brainMask) = NaN;  % PSC undefined outside mask
pc_abs = abs(pc);

% robust ranges for alpha modulation from in-mask finite values
pc_flat         = pc(isfinite(pc) & brainMask);
if isempty(pc_flat), pc_flat = 0; end
pc_lo           = prctile(pc_flat,2); pc_hi = prctile(pc_flat,98);
pc_mag_default  = max(5, max(abs([pc_lo pc_hi])));

pc_abs_flat  = pc_abs(isfinite(pc_abs) & brainMask);
if isempty(pc_abs_flat), pc_abs_flat = 0; end
pc_abs_lo    = prctile(pc_abs_flat,2); pc_abs_hi = prctile(pc_abs_flat,98);
if pc_abs_hi <= pc_abs_lo, pc_abs_hi = pc_abs_lo + eps; end

% baseline underlay window (for background when "Baseline mean" is selected)
[und_lo_base, und_hi_base] = robust_window(baseline_mean);

%% ---- FIGURE & PANELS ----
fig = figure('Color','w','Name','Axial %Change Viewer (Interactive)', ...
             'Units','pixels','Position',[100 80 1400 900]);

panelLeft  = uipanel('Parent',fig,'Units','pixels','BorderType','none');
panelRight = uipanel('Parent',fig,'Units','pixels','BorderType','none','BackgroundColor','w');
panelCtrl  = uipanel('Parent',fig,'Units','pixels','BorderType','none','BackgroundColor','w');

% --- Right-side axes (time-series / PSC) ---
axTS = axes('Parent',panelRight);
box(axTS,'on'); grid(axTS,'on');
xlabel(axTS,'Time (in sec)');
ylabel(axTS,'Raw MRI Signal');
title(axTS,'Time series / PSC (click a slice)');

% --- PSC controls (right panel) ---
chk_psc = uicontrol('Parent',panelRight,'Style','checkbox','Units','pixels', ...
    'String','Plot PSC', 'Value',1, 'BackgroundColor','w', ...
    'Callback',@(~,~)refreshTSLabel(axTS), 'TooltipString','Toggle PSC (%) vs raw signal');

lbl_bmode = uicontrol('Parent',panelRight,'Style','text','Units','pixels', ...
    'String','Baseline mode:', 'BackgroundColor','w','HorizontalAlignment','right','FontWeight','bold');

pop_bmode = uicontrol('Parent',panelRight,'Style','popupmenu','Units','pixels', ...
    'String',{'Global baseline_idx','Manual [b1 b2]'}, 'Callback',@onBaselineCtrlChanged);

lbl_b1 = uicontrol('Parent',panelRight,'Style','text','Units','pixels', ...
    'String','b1:', 'BackgroundColor','w','HorizontalAlignment','right');

edt_b1 = uicontrol('Parent',panelRight,'Style','edit','Units','pixels', ...
    'String',num2str(min(baseline_idx)), 'Callback',@onBaselineCtrlChanged, ...
    'TooltipString','Baseline start (frame)');

lbl_b2 = uicontrol('Parent',panelRight,'Style','text','Units','pixels', ...
    'String','b2:', 'BackgroundColor','w','HorizontalAlignment','right');

edt_b2 = uicontrol('Parent',panelRight,'Style','edit','Units','pixels', ...
    'String',num2str(max(baseline_idx)), 'Callback',@onBaselineCtrlChanged, ...
    'TooltipString','Baseline end (frame)');

lbl_bhint = uicontrol('Parent',panelRight,'Style','text','Units','pixels', ...
    'String','(Used if Manual)', 'BackgroundColor','w','HorizontalAlignment','left','ForegroundColor',[0.3 0.3 0.3]);

% --- Tiles (slice grid) ---
[slices_per_page, tl, ax, imU, imP] = createTiles(panelLeft, default_tiles);

% Colormap + colorbar
cmap = gbhot_diverging(256);
colormap(fig, cmap);
try
    cb = colorbar(tl,'Location','southoutside');
catch
    cb = colorbar('southoutside');
end
cb.Label.String = '% change';

%% ---- Top row controls over tiles ----
btn_theme = uicontrol('Parent',panelLeft,'Style','togglebutton','Units','pixels', ...
    'String','Dark theme', 'Callback',@toggleTheme, 'TooltipString','Toggle light/dark theme');

lbl_sign = uicontrol('Parent',panelLeft,'Style','text','Units','pixels', ...
    'String','Nature of Maps:','HorizontalAlignment','right');

pop_sign = uicontrol('Parent',panelLeft,'Style','popupmenu','Units','pixels', ...
    'String',{'All Changes','Positive Changes only','Negative Changes only'}, 'Callback',@changeOverlay);

%% ---- Bottom controls (panelCtrl) widgets ----
btn_prev = uicontrol('Parent',panelCtrl,'Style','pushbutton','Units','pixels', ...
    'String','◀ Prev','FontWeight','bold','Callback',@prevPage, 'TooltipString','Previous page (←)');

btn_next = uicontrol('Parent',panelCtrl,'Style','pushbutton','Units','pixels', ...
    'String','Next ▶','FontWeight','bold','Callback',@nextPage, 'TooltipString','Next page (→)');

lbl_page = uicontrol('Parent',panelCtrl,'Style','text','Units','pixels', ...
    'String','','BackgroundColor','w','FontWeight','bold','HorizontalAlignment','center','Tag','lbl_page');

% Four sub-panels + Panel C2 (Baseline & PSC)
pnlA = uipanel('Parent',panelCtrl,'Units','pixels','BorderType','etchedin');  % Alpha / Scale (+ windowing)
pnlB = uipanel('Parent',panelCtrl,'Units','pixels','BorderType','etchedin');  % Display
pnlC = uipanel('Parent',panelCtrl,'Units','pixels','BorderType','etchedin');  % ROI (Panel D)
pnlC2 = uipanel('Parent',panelCtrl,'Units','pixels','BorderType','none');     % Baseline & PSC (Panel C)

% Title labels for Panel C (Baseline & PSC)
lblC2_title = uicontrol('Parent',pnlC2,'Style','text','Units','pixels', ...
    'String','Baseline & PSC','FontWeight','bold','BackgroundColor','w','HorizontalAlignment','center');

% Move PSC widgets into Panel C2
set([chk_psc, lbl_bmode, pop_bmode, lbl_b1, edt_b1, lbl_b2, edt_b2, lbl_bhint], 'Parent', pnlC2);

% Panel A controls
lblA_title  = uicontrol('Parent',pnlA,'Style','text','Units','pixels','String','Alpha / Scale','FontWeight','bold','BackgroundColor','w','HorizontalAlignment','center');

lblA_amode  = uicontrol('Parent',pnlA,'Style','text','Units','pixels','String','α mode:', ...
    'BackgroundColor','w','HorizontalAlignment','right','FontWeight','bold');
pop_alpha   = uicontrol('Parent',pnlA,'Style','popupmenu','Units','pixels', ...
    'String',{'Constant','|%Change|','Underlay'},'Callback',@changeOverlay);

lblA_again  = uicontrol('Parent',pnlA,'Style','text','Units','pixels','String','α gain:', ...
    'BackgroundColor','w','HorizontalAlignment','right','FontWeight','bold');
sld_again   = uicontrol('Parent',pnlA,'Style','slider','Units','pixels','Min',0,'Max',1.5,'Value',1.0, ...
    'TooltipString','Alpha gain','Callback',@changeOverlay);

lblA_scale  = uicontrol('Parent',pnlA,'Style','text','Units','pixels','String','Scale (±%):', ...
    'BackgroundColor','w','HorizontalAlignment','right','FontWeight','bold');
sld_scale   = uicontrol('Parent',pnlA,'Style','slider','Units','pixels','Min',0.1,'Max',15, ...
    'Value',min(15, max(1, pc_mag_default)),'SliderStep',[1/150 1/15],'Callback',@adjustScale, ...
    'TooltipString','Overlay scale (Ctrl + wheel)');
lbl_scale   = uicontrol('Parent',pnlA,'Style','text','Units','pixels','String',sprintf('%.1f',get(sld_scale,'Value')), ...
    'BackgroundColor','w','HorizontalAlignment','left');

% Underlay windowing controls (Panel A)
lblB_win    = uicontrol('Parent',pnlA,'Style','text','Units','pixels', ...
    'String','Underlay win:','BackgroundColor','w','HorizontalAlignment','right','FontWeight','bold');
pop_win     = uicontrol('Parent',pnlA,'Style','popupmenu','Units','pixels', ...
    'String',{'Global 2–98','Slice 2–98','Manual %'},'Callback',@changeWindowing);
sld_pLo     = uicontrol('Parent',pnlA,'Style','slider','Units','pixels','Min',0,'Max',50,'Value',2, ...
    'SliderStep',[1/50 5/50],'Callback',@changeWindowing,'TooltipString','Manual low percentile');
sld_pHi     = uicontrol('Parent',pnlA,'Style','slider','Units','pixels','Min',50,'Max',100,'Value',98, ...
    'SliderStep',[1/50 5/50],'Callback',@changeWindowing,'TooltipString','Manual high percentile');

% Panel B controls (Display)
lblB_title  = uicontrol('Parent',pnlB,'Style','text','Units','pixels','String','Display','FontWeight','bold','BackgroundColor','w','HorizontalAlignment','center');
lblB_tiles  = uicontrol('Parent',pnlB,'Style','text','Units','pixels','String','Tiles:','BackgroundColor','w','HorizontalAlignment','right','FontWeight','bold');
valTiles = 1 + (default_tiles==8);
pop_tiles   = uicontrol('Parent',pnlB,'Style','popupmenu','Units','pixels','String',{'4 (2×2)','8 (2×4)'},'Value',valTiles,'Callback',@changeTiles);

lblB_under  = uicontrol('Parent',pnlB,'Style','text','Units','pixels','String','Underlay:','BackgroundColor','w','HorizontalAlignment','right','FontWeight','bold');
pop_under   = uicontrol('Parent',pnlB,'Style','popupmenu','Units','pixels','String',{'Mean(4D)','Baseline mean','Timepoint t'},'Callback',@changeUnderlay);

% Inline timepoint controls for Underlay=Timepoint t
lbl_tp     = uicontrol('Parent',pnlB,'Style','text','Units','pixels', ...
    'String','t =','BackgroundColor','w','HorizontalAlignment','right','Visible','off');
edt_tp     = uicontrol('Parent',pnlB,'Style','edit','Units','pixels', ...
    'String','1','BackgroundColor','w','Visible','off', ...
    'Callback',@onTimepointEdit, ...
    'TooltipString','Enter timepoint (integer) and press Enter');

% Show α map and overlay (separate lines)
chk_showA   = uicontrol('Parent',pnlB,'Style','checkbox','Units','pixels','String','Show α map','BackgroundColor','w','Callback',@changeOverlay,'TooltipString','Toggle alpha map (A)');
chk_overlay = uicontrol('Parent',pnlB,'Style','checkbox','Units','pixels','String','Show overlay','Value',1,'BackgroundColor','w','Callback',@changeOverlay,'TooltipString','Toggle overlay (O)');

% Panel C controls (ROI) — includes Export ROI → NIfTI
lblC_title  = uicontrol('Parent',pnlC,'Style','text','Units','pixels','String','ROI (3D)','FontWeight','bold','BackgroundColor','w','HorizontalAlignment','center');
chk_roi     = uicontrol('Parent',pnlC,'Style','checkbox','Units','pixels','String','Use ROI (3D)','Value',1,'BackgroundColor','w','Callback',@roiChanged,'TooltipString','Toggle averaging over a 3D neighborhood');

% ROI type dropdown: placed next to checkbox
pop_roishape= uicontrol('Parent',pnlC,'Style','popupmenu','Units','pixels','String',{'3D Sphere','3D Cylinder'},'Value', 2,'Callback',@roiChanged,'TooltipString','ROI shape in 3D');

lblC_rXY    = uicontrol('Parent',pnlC,'Style','text','Units','pixels','String','XY Radius:','BackgroundColor','w','HorizontalAlignment','right');
sld_roirXY  = uicontrol('Parent',pnlC,'Style','slider','Units','pixels','Min',1,'Max',15,'Value',3,'SliderStep',[1/14 3/14],'Callback',@roiChanged,'TooltipString','ROI radius in voxels (XY)');
lbl_roirXY  = uicontrol('Parent',pnlC,'Style','text','Units','pixels','String','3','BackgroundColor','w','HorizontalAlignment','left');

lblC_rZ     = uicontrol('Parent',pnlC,'Style','text','Units','pixels','String','±Z half-depth:','BackgroundColor','w','HorizontalAlignment','right');
sld_roirZ   = uicontrol('Parent',pnlC,'Style','slider','Units','pixels','Min',0,'Max',8,'Value',2,'SliderStep',[1/8 2/8],'Callback',@roiChanged,'TooltipString','Half-depth in slices (cylinder only)');
lbl_roirZ   = uicontrol('Parent',pnlC,'Style','text','Units','pixels','String','2','BackgroundColor','w','HorizontalAlignment','left');

% --- NEW: Group assignment controls ---
lblC_group   = uicontrol('Parent',pnlC,'Style','text','Units','pixels', ...
    'String','Assign to Group:','BackgroundColor','w','HorizontalAlignment','right','FontWeight','bold');
pop_group    = uicontrol('Parent',pnlC,'Style','popupmenu','Units','pixels', ...
    'String',{'FAP-AAV + FAP','FAP-AAV + PBS','Ctrl-AAV + FAP','Ctrl-AAV + PBS'}, 'Value',1);
btn_groupAdd = uicontrol('Parent',pnlC,'Style','pushbutton','Units','pixels', ...
    'String','Save to Group','FontWeight','bold','Callback',@saveROIToGroup, ...
    'TooltipString','Append this ROI (mask + metadata + mean TS + PSC) to GroupX.mat');

% --- NEW: Group AVG PSC controls ---
lblC_gAvg     = uicontrol('Parent',pnlC,'Style','text','Units','pixels', ...
    'String','Group Avg PSC:','BackgroundColor','w','HorizontalAlignment','left','FontWeight','bold');

pop_gFile     = uicontrol('Parent',pnlC,'Style','popupmenu','Units','pixels', ...
    'String',{'<no group files found>'}, 'Value',1, 'Callback',@(~,~)[]);

btn_gRefresh  = uicontrol('Parent',pnlC,'Style','pushbutton','Units','pixels', ...
    'String','↻','FontWeight','bold','TooltipString','Rescan group folder', ...
    'Callback',@refreshGroupFileList);

btn_gPlot     = uicontrol('Parent',pnlC,'Style','pushbutton','Units','pixels', ...
    'String','Plot Avg PSC','FontWeight','bold','TooltipString','Plot mean±SEM across datasets in this group', ...
    'Callback',@plotGroupAvgPSC);


% --- ROI mask export controls (compact) ---
lblC_export   = uicontrol('Parent',pnlC,'Style','text','Units','pixels', ...
    'String','Mask (.nii.gz):','BackgroundColor','w','HorizontalAlignment','left','FontWeight','bold');
edt_roiOut    = uicontrol('Parent',pnlC,'Style','edit','Units','pixels', ...
    'String','roi_mask','BackgroundColor','w', ...
    'TooltipString','Type file name only (no extension). Will save as <name>.nii.gz');
btn_exportROI = uicontrol('Parent',pnlC,'Style','pushbutton','Units','pixels', ...
    'String','Save','FontWeight','bold','Callback',@exportROINifti, ...
    'TooltipString','Write current 3D ROI to NIfTI (uint8 mask).');

% --- ROI timeseries export (compact CSV) ---
lblC_exportTS = uicontrol('Parent',pnlC,'Style','text','Units','pixels', ...
    'String','ROI TS (.txt):','BackgroundColor','w','HorizontalAlignment','left','FontWeight','bold');
edt_tsOut     = uicontrol('Parent',pnlC,'Style','edit','Units','pixels', ...
    'String','roi_timeseries','BackgroundColor','w', ...
    'TooltipString','Type file name only (no extension). Will save as <name>.txt');
btn_exportTS  = uicontrol('Parent',pnlC,'Style','pushbutton','Units','pixels', ...
    'String','Save','FontWeight','bold','Callback',@exportROITimeSeries, ...
    'TooltipString','Export ROI mean time series (raw or PSC based on the Plot PSC toggle).');

% Status bar
lbl_status = uicontrol('Parent',panelCtrl,'Style','text','Units','pixels', ...
    'String','Ready', 'BackgroundColor','w','HorizontalAlignment','left','FontName','Consolas');

%% ---- Pack state ----
S.Y = Y4D; S.X = X; S.Ydim = Ydim; S.Z = Z; S.T = T;
S.pc = pc; S.pc_abs = pc_abs; S.mask = brainMask;
S.U_mean4D = U_mean4D; S.U_base = baseline_mean;
S.und_lo_mean=und_lo_mean; S.und_hi_mean=und_hi_mean;
S.und_lo_base=und_lo_base; S.und_hi_base=und_hi_base;
S.pc_abs_lo = pc_abs_lo; S.pc_abs_hi = pc_abs_hi;
S.alpha_base = alpha_base;
S.eps_den_global = eps_den_global;       % dynamic epsilon seed
S.global_med     = global_med;
S.PSC_YLIM = PSC_YLIM;
S.SHOW_VOXEL_PSC = SHOW_VOXEL_PSC;

S.UI_FONT = UI_FONT; S.UI_FONTSIZE = UI_FONTSIZE;

S.fig = fig;
S.panelLeft = panelLeft; S.panelRight = panelRight; S.panelCtrl = panelCtrl;
S.tl = tl; S.ax = ax; S.imU = imU; S.imP = imP; S.cb = cb; S.axTS = axTS;
S.slices_per_page = slices_per_page; S.startSlice = 1;

% Top row
S.btn_theme = btn_theme; S.isDark = false; S.lbl_sign = lbl_sign; S.pop_sign = pop_sign;

% Right-panel controls
S.chk_psc = chk_psc; S.pop_bmode = pop_bmode; S.edt_b1 = edt_b1; S.edt_b2 = edt_b2;
S.lbl_bmode = lbl_bmode; S.lbl_b1 = lbl_b1; S.lbl_b2 = lbl_b2; S.lbl_bhint = lbl_bhint;
S.baseline_idx = baseline_idx; S.signal_idx = signal_idx;

% Bottom bar widgets
S.btn_prev = btn_prev; S.btn_next = btn_next; S.lbl_page = lbl_page; S.lbl_status = lbl_status;
S.pnlA = pnlA; S.pnlB = pnlB; S.pnlC = pnlC; S.pnlC2 = pnlC2; S.lblC2_title = lblC2_title;

% Panel B (display group)
S.lblB_title = lblB_title; S.lblB_tiles = lblB_tiles; S.lblB_under = lblB_under;
S.pop_tiles  = pop_tiles;  S.pop_under  = pop_under;
S.chk_overlay = chk_overlay; S.chk_showA = chk_showA;
S.lbl_tp = lbl_tp; S.edt_tp = edt_tp;

% Panel A (alpha/scale + windowing)
S.lblA_title = lblA_title; S.lblA_amode = lblA_amode; S.pop_alpha = pop_alpha;
S.lblA_again = lblA_again; S.sld_again = sld_again;
S.lblA_scale = lblA_scale; S.sld_scale = sld_scale; S.lbl_scale = lbl_scale;
S.lblB_win = lblB_win; S.pop_win = pop_win; S.sld_pLo = sld_pLo; S.sld_pHi = sld_pHi;

% Panel D (ROI + Export + Group)
S.lblC_title = lblC_title; S.chk_roi = chk_roi; S.pop_roishape = pop_roishape;
S.lblC_rXY = lblC_rXY; S.sld_roirXY = sld_roirXY; S.lbl_roirXY = lbl_roirXY;
S.lblC_rZ = lblC_rZ; S.sld_roirZ = sld_roirZ; S.lbl_roirZ = lbl_roirZ;
S.lblC_export = lblC_export; S.edt_roiOut = edt_roiOut; S.btn_exportROI = btn_exportROI;
S.lblC_exportTS = lblC_exportTS; S.edt_tsOut = edt_tsOut; S.btn_exportTS = btn_exportTS;

% NEW: Group AVG PSC UI
S.lblC_gAvg   = lblC_gAvg;
S.pop_gFile   = pop_gFile;
S.btn_gRefresh= btn_gRefresh;
S.btn_gPlot   = btn_gPlot;

% Will hold discovered group files
S.groupFiles  = {};


% NEW: Group UI + root
S.lblC_group   = lblC_group;
S.pop_group    = pop_group;
S.btn_groupAdd = btn_groupAdd;

S.roi_center = [];                         % [rA cA z]
S.hROI  = gobjects(S.slices_per_page,1);   % per-tile ROI outlines
S.hCross = gobjects(S.slices_per_page,2);  % per-tile crosshair [horiz, vert]
S.hSliceLbl = gobjects(S.slices_per_page,1); % per-tile bottom "Slice #N" labels
S.roi3d_mask = [];                         % cached 3D mask
S.currentSlices = zeros(1,S.slices_per_page); % slice index per tile

% Baseline & TS sizes + timepoint current (inline input)
S.PSC_editW = 100;      % smaller edit width for b1/b2
S.PSC_editH = [];       % leave [] to use rowH
S.tp_current = 1;       % default timepoint when Underlay = Timepoint t

S.niiInfo   = info;     % for export geometry
S.inNiiPath = in_nii;   % full path string (original input path)
S.inDir     = fileparts(in_nii);

% Detect central Group root (…/AnalysedData), fallback to fixed path
S.groupRoot = detectGroupRoot(S.inDir);

guidata(fig,S);

S = guidata(fig);
refreshGroupFileList(fig,[]);   % populate group file list on startup
% guidata(fig,S);


% --- Global callbacks ---
set(fig,'WindowButtonMotionFcn',@onMouseMove);
set(fig,'WindowScrollWheelFcn',@onScroll);
set(fig,'WindowKeyPressFcn',@onKeyPress);
set(fig,'SizeChangedFcn',@onResize);   % responsive layout

% Initial render + layout + fonts + theme
renderPage(fig);
adjustScale(sld_scale,[]);
applyTSFont(fig);
applyTheme(S);
applyGlobalUIFont(S);
relayoutUI(S);       % first layout pass
tuneTSAxes(fig);

% Cleanup temp file (if created)
if ~isempty(tmp_nii) && exist(tmp_nii,'file')
    try, delete(tmp_nii); catch, end
end

% Ensure TS font defaults exist
S = guidata(fig);
if ~isfield(S,'tsFontName')
    S.tsFontName  = S.UI_FONT;
    S.tsFontTitle = max(16, 16);
    S.tsFontLabel = max(16, 18);
    S.tsFontXTick = max(14, 18);
    S.tsFontYTick = max(14, 18);
end
guidata(fig,S);

%% =========================
%% === CALLBACKS / FUNCS ===
%% =========================

function [spp, tl, ax, imU, imP] = createTiles(parentPanel, tiles)
    delete(findall(parentPanel,'Type','axes'));
    delete(findall(parentPanel,'Type','tiledlayout'));
    if tiles==4
        tl  = tiledlayout(parentPanel,2,2,'Padding','compact','TileSpacing','compact'); spp = 4;
    else
        tl  = tiledlayout(parentPanel,2,4,'Padding','compact','TileSpacing','compact'); spp = 8;
    end
    ax  = gobjects(spp,1); imU = gobjects(spp,1); imP = gobjects(spp,1);
    for k = 1:spp
        ax(k) = nexttile(tl,k); hold(ax(k),'on'); axis(ax(k),'image'); axis(ax(k),'off');
        set(ax(k),'XTick',[],'YTick',[],'XColor','none','YColor','none');
        set(ax(k),'ButtonDownFcn',@(h,e)onClickVoxel(h,e,NaN));
    end
end

function changeTiles(src,~)
    fig = ancestor(src,'figure'); S = guidata(fig);
    val = get(S.pop_tiles,'Value'); tiles = 4; if val==2, tiles=8; end
    [spp, tl, ax, imU, imP] = createTiles(S.panelLeft, tiles);
    S.slices_per_page = spp; S.tl = tl; S.ax = ax; S.imU = imU; S.imP = imP;
    S.hROI = gobjects(spp,1); S.hCross = gobjects(spp,2); S.hSliceLbl = gobjects(spp,1);
    delete(S.cb); try, S.cb = colorbar(S.tl,'Location','southoutside'); catch, S.cb=colorbar('southoutside'); end
    guidata(fig,S);
    renderPage(fig); applyGlobalUIFont(S); relayoutUI(S);
end

function prevPage(src,~)
    fig = ancestor(src,'figure'); S = guidata(fig);
    S.startSlice = max(S.startSlice - S.slices_per_page, 1);
    guidata(fig,S); renderPage(fig); applyGlobalUIFont(S);
end

function nextPage(src,~)
    fig = ancestor(src,'figure'); S = guidata(fig);
    S.startSlice = min(S.startSlice + S.slices_per_page, max(S.Z - S.slices_per_page + 1, 1));
    guidata(fig,S); renderPage(fig); applyGlobalUIFont(S);
end

function changeUnderlay(src,~)
    fig = ancestor(src,'figure'); S = guidata(fig);
    if get(S.pop_under,'Value')==3
        set(S.lbl_tp,'Visible','on');
        set(S.edt_tp,'Visible','on','String',num2str(S.tp_current));
    else
        set(S.lbl_tp,'Visible','off');
        set(S.edt_tp,'Visible','off');
    end
    renderPage(fig); applyGlobalUIFont(S);
end

function onTimepointEdit(h,~)
    fig = ancestor(h,'figure'); S = guidata(fig);
    tnew = str2double(get(S.edt_tp,'String'));
    if ~(isfinite(tnew) && tnew==round(tnew) && tnew>=1 && tnew<=S.T)
        tnew = max(1, min(S.T, round(S.tp_current)));
        set(S.edt_tp,'String',num2str(tnew));
        return;
    end
    S.tp_current = round(tnew);
    guidata(fig,S);
    renderPage(fig); applyGlobalUIFont(S);
end

function changeWindowing(~,~)
    fig = gcbf; 
    S = guidata(fig);
    pLo = get(S.sld_pLo,'Value'); pHi = get(S.sld_pHi,'Value');
    if pHi <= pLo, set(S.sld_pHi,'Value',min(100, pLo+1)); end
    renderPage(fig);
    S = guidata(fig); applyGlobalUIFont(S);
end

function changeOverlay(~,~)
    fig = gcbf; renderPage(fig);
    S = guidata(fig); applyGlobalUIFont(S);
end

function roiChanged(src,~)
    fig = ancestor(src,'figure'); S = guidata(fig);
    set(S.lbl_roirXY,'String',sprintf('%d',round(get(S.sld_roirXY,'Value'))));
    set(S.lbl_roirZ,'String', sprintf('%d',round(get(S.sld_roirZ,'Value'))));
    if get(S.pop_roishape,'Value')==1, set(S.sld_roirZ,'Enable','off'); else, set(S.sld_roirZ,'Enable','on'); end
    S.roi3d_mask = []; guidata(fig,S);
    if ~isempty(S.roi_center) && logical(get(S.chk_roi,'Value')), plotTSAtCenter(fig); renderPage(fig); end
    applyGlobalUIFont(S);
end

function adjustScale(src,~)
    fig = ancestor(src,'figure'); S = guidata(fig);
    showA = logical(get(S.chk_showA,'Value'));
    if showA
        for k = 1:numel(S.ax), if isgraphics(S.ax(k)), caxis(S.ax(k),[0 1]); end, end
        S.cb.Limits = [0 1]; S.cb.Ticks = [0 0.5 1]; S.cb.Label.String = 'alpha';
        set(S.lbl_scale,'String','—');
    else
        v = get(S.sld_scale,'Value'); if v<=0, v=0.1; set(S.sld_scale,'Value',v); end
        for k = 1:numel(S.ax), if isgraphics(S.ax(k)), caxis(S.ax(k),[-v v]); end, end
        S.cb.Limits = [-v v]; S.cb.Ticks = [-v 0 v];
        S.cb.TickLabels = {num2str(-v,'%.0f'),'0',num2str(v,'%.0f')}; S.cb.Label.String = '% change';
        set(S.lbl_scale,'String',sprintf('%.1f',v));
    end
    guidata(fig,S); applyGlobalUIFont(S);
end

function renderPage(fig)
    S = guidata(fig);

    % Underlay source
    uSrc = get(S.pop_under,'Value'); % 1 mean, 2 base, 3 timepoint
    switch uSrc
        case 1, U = S.U_mean4D; glo = S.und_lo_mean; ghi = S.und_hi_mean;
        case 2, U = S.U_base;    glo = S.und_lo_base; ghi = S.und_hi_base;
        otherwise
            tp = max(1, min(S.T, S.tp_current));
            U = S.Y(:,:,:,tp);
            flat = U(isfinite(U)); if isempty(flat), glo=0; ghi=1; else, [glo,ghi] = robust_window(U); end
    end

    % Controls
    sign_mode = get(S.pop_sign,'Value'); % 1 both, 2 +, 3 −
    win_mode  = get(S.pop_win,'Value');  % 1 global, 2 slice, 3 manual
    pLo = get(S.sld_pLo,'Value'); pHi = get(S.sld_pHi,'Value'); if pHi<=pLo, pHi=pLo+eps; end
    showA = logical(get(S.chk_showA,'Value'));
    showOverlay = logical(get(S.chk_overlay,'Value'));
    roi_on = logical(get(S.chk_roi,'Value'));

    first = S.startSlice; last = min(first + S.slices_per_page - 1, S.Z);
    slices = first:last; S.currentSlices(:) = 0; S.currentSlices(1:numel(slices)) = slices;
    set(S.lbl_page,'String',sprintf('Slices %d–%d of %d', first, last, S.Z));

    % Build 3D ROI mask if needed
    if roi_on && ~isempty(S.roi_center)
        if isempty(S.roi3d_mask), S.roi3d_mask = roi_mask_3d(S); guidata(fig,S); end
    else
        S.roi3d_mask = []; guidata(fig,S);
    end

    A_all = [];

    for k = 1:S.slices_per_page
        axk = S.ax(k);
        if k <= numel(slices)
            idx = slices(k);

            Uslc0  = U(:,:,idx);
            PCslc0 = S.pc(:,:,idx);
            Mslc0  = S.mask(:,:,idx);

            % Display orientation
            Uslc  = rot90(Uslc0,-1);
            PCslc = rot90(PCslc0,-1);
            Mslc  = rot90(Mslc0,-1);

            % Sign filter
            switch sign_mode
                case 2, mask_sign = PCslc > 0;
                case 3, mask_sign = PCslc < 0;
                otherwise, mask_sign = true(size(PCslc));
            end

            % Apply mask
            PCslc(~Mslc) = NaN;
            PCvis = PCslc;
            PCvis(~mask_sign) = NaN;

            % Window
            switch win_mode
                case 1, lo_use = glo; hi_use = ghi;
                case 2, [lo_use,hi_use] = robust_window(Uslc);
                otherwise
                    flat = Uslc(isfinite(Uslc));
                    if isempty(flat), lo_use=0; hi_use=1;
                    else
                        lo_use = prctile(flat, pLo);
                        hi_use = prctile(flat, pHi);
                        if hi_use<=lo_use, hi_use=lo_use+eps; end
                    end
            end

            % Alpha
            A = compute_alpha(S, PCslc, Uslc, Mslc & mask_sign);
            A_all = [A_all; A(isfinite(A) & (A>0))]; %#ok<AGROW>

            % Draw / update images
            if ~isgraphics(S.imU(k))
                axes(axk); cla(axk); hold(axk,'on'); axis(axk,'image'); axis(axk,'off');
                set(axk,'XTick',[],'YTick',[],'XColor','none','YColor','none');
                set(axk,'ButtonDownFcn',@(h,e)onClickVoxel(h,e,idx));
                S.imU(k) = image(axk, to_rgb_gray(Uslc,lo_use,hi_use), ...
                    'ButtonDownFcn',@(h,e)onClickVoxel(h,e,idx),'HitTest','on','PickableParts','all');
                S.imP(k) = imagesc(axk, PCvis, 'HitTest','off','PickableParts','none');

                % Bottom slice label
                if isgraphics(S.hSliceLbl(k)), delete(S.hSliceLbl(k)); end
                fg = tern(S.isDark, [0.95 0.95 0.95], [0 0 0]);
                S.hSliceLbl(k) = text(axk, 0.5, 0.02, sprintf('Slice #%d', idx), ...
                    'Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom', ...
                    'FontName', S.UI_FONT, 'FontSize', S.UI_FONTSIZE-1, 'FontWeight','bold', ...
                    'Color', fg, 'Interpreter','none');
            else
                set(axk,'ButtonDownFcn',@(h,e)onClickVoxel(h,e,idx));
                set(S.imU(k),'CData', to_rgb_gray(Uslc,lo_use,hi_use), ...
                    'ButtonDownFcn',@(h,e)onClickVoxel(h,e,idx),'HitTest','on','PickableParts','all');
                set(S.imP(k),'CData', PCvis, 'HitTest','off','PickableParts','none');
                set(axk,'XTick',[],'YTick',[],'XColor','none','YColor','none');

                % Update bottom label
                if isgraphics(S.hSliceLbl(k)), delete(S.hSliceLbl(k)); end
                fg = tern(S.isDark, [0.95 0.95 0.95], [0 0 0]);
                S.hSliceLbl(k) = text(axk, 0.5, 0.02, sprintf('Slice #%d', idx), ...
                    'Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom', ...
                    'FontName', S.UI_FONT, 'FontSize', S.UI_FONTSIZE-1, 'FontWeight','bold', ...
                    'Color', fg, 'Interpreter','none');
            end

            % ROI outline
            if roi_on && ~isempty(S.roi3d_mask)
                maskSlice = S.roi3d_mask(:,:,idx);
                if any(maskSlice(:))
                    maskDisp = rot90(maskSlice,-1);
                    if isgraphics(S.hROI(k)), delete(S.hROI(k)); end
                    [~, h] = contour(axk, double(maskDisp), [0.5 0.5], 'LineColor','m', 'LineWidth',1.2);
                    if ~isempty(h), S.hROI(k) = h; else, S.hROI(k) = gobjects(1); end
                else
                    if isgraphics(S.hROI(k)), delete(S.hROI(k)); S.hROI(k)=gobjects(1); end
                end
            else
                if isgraphics(S.hROI(k)), delete(S.hROI(k)); S.hROI(k)=gobjects(1); end
            end

            % Crosshair on clicked voxel
            if ~isempty(S.roi_center) && S.roi_center(3)==idx
                disp_i = S.roi_center(2);           % cA
                disp_j = S.X - S.roi_center(1) + 1; % X - rA + 1
                nx = S.X; ny = S.Ydim;
                if isgraphics(S.hCross(k,1)), delete(S.hCross(k,1)); end
                if isgraphics(S.hCross(k,2)), delete(S.hCross(k,2)); end
                S.hCross(k,1) = plot(axk, [1 nx], [disp_i disp_i], 'w-', 'LineWidth', 0.8);
                S.hCross(k,2) = plot(axk, [disp_j disp_j], [1 ny], 'w-', 'LineWidth', 0.8);
            else
                if isgraphics(S.hCross(k,1)), delete(S.hCross(k,1)); S.hCross(k,1)=gobjects(1); end
                if isgraphics(S.hCross(k,2)), delete(S.hCross(k,2)); S.hCross(k,2)=gobjects(1); end
            end

            % Overlay visibility
            if showOverlay
                set(S.imP(k),'Visible','on');
                if showA
                    set(S.imP(k),'CData', A, 'AlphaData', 1, 'CDataMapping','scaled');
                    caxis(axk,[0 1]);
                else
                    set(S.imP(k),'CData', PCvis, 'AlphaData', A, 'CDataMapping','scaled');
                end
            else
                set(S.imP(k),'Visible','off');
            end
            set(axk,'Visible','on');
        else
            cla(axk); set(axk,'Visible','off');
            if isgraphics(S.hROI(k)), delete(S.hROI(k));  S.hROI(k)=gobjects(1); end
            if isgraphics(S.hCross(k,1)), delete(S.hCross(k,1)); S.hCross(k,1)=gobjects(1); end
            if isgraphics(S.hCross(k,2)), delete(S.hCross(k,2)); S.hCross(k,2)=gobjects(1); end
            if isgraphics(S.hSliceLbl(k)), delete(S.hSliceLbl(k)); S.hSliceLbl(k)=gobjects(1); end
        end
    end

    % Alpha stats (status)
    if isempty(A_all), amin=NaN; amed=NaN; amax=NaN;
    else, amin=min(A_all); amed=median(A_all); amax=max(A_all);
    end
    set(S.lbl_status,'String',sprintf('α[min med max]=[%.2f %.2f %.2f]',amin,amed,amax));

    if showOverlay && showA
        S.cb.Label.String = 'alpha'; S.cb.Limits=[0 1]; S.cb.Ticks=[0 .5 1];
    else
        S.cb.Label.String = '% change'; adjustScale(S.sld_scale,[]);
    end
    guidata(fig,S);
end

function setStartSlice(idx)
    fig = gcbf; S = guidata(fig);
    S.startSlice = max(1, min(S.Z - S.slices_per_page + 1, idx));
    guidata(fig,S); renderPage(fig); applyGlobalUIFont(S);
end

function refreshTSLabel(axTS)
    fig = ancestor(axTS,'figure'); S = guidata(fig);
    S.tsFontName  = S.UI_FONT;
    S.tsFontTitle = max(16, 16);
    S.tsFontLabel = max(16, 18);
    S.tsFontXTick = max(14, 18);
    S.tsFontYTick = max(14, 18);
    guidata(fig,S);
    if logical(get(S.chk_psc,'Value')), ylabel(axTS,'PSC (%)'); else, ylabel(axTS,'Signal (a.u.)'); end
    applyGlobalUIFont(S);
end

function applyTSFont(fig)
    S = guidata(fig);
    if ~isfield(S,'tsFontName'),  S.tsFontName  = S.UI_FONT; end
    if ~isfield(S,'tsFontTitle'), S.tsFontTitle = max(16,16); end
    if ~isfield(S,'tsFontLabel'), S.tsFontLabel = max(16,16); end
    if ~isfield(S,'tsFontXTick'), S.tsFontXTick = max(14,16); end
    if ~isfield(S,'tsFontYTick'), S.tsFontYTick = max(14,16); end
    guidata(fig,S);

    ax = S.axTS;
    set(ax,'FontName',S.UI_FONT);
    set(ax,'FontSize',max(S.UI_FONTSIZE-2, get(ax,'FontSize')));

    set(get(ax,'Title'),  'FontName',S.UI_FONT,'FontSize',max(S.UI_FONTSIZE-2,S.tsFontTitle),'FontWeight','bold');
    set(get(ax,'XLabel'), 'FontName',S.UI_FONT,'FontSize',max(S.UI_FONTSIZE-2,S.tsFontLabel));
    set(get(ax,'YLabel'), 'FontName',S.UI_FONT,'FontSize',max(S.UI_FONTSIZE-2,S.tsFontLabel));

    try
        ax.XAxis.FontName = S.UI_FONT; ax.YAxis.FontName = S.UI_FONT;
        ax.XAxis.FontSize = max(S.UI_FONTSIZE-4, S.tsFontXTick);
        ax.YAxis.FontSize = max(S.UI_FONTSIZE-4, S.tsFontYTick);
    catch
        set(ax,'FontSize', max([S.UI_FONTSIZE-2, S.tsFontXTick, S.tsFontYTick]));
    end
end

function tuneTSAxes(fig)
    S  = guidata(fig);
    ax = S.axTS;
    ax.PositionConstraint = 'innerposition';
    ti = get(ax,'TightInset');
    padLBRT = [max(ti(1), 0.10), max(ti(2), 0.14), max(ti(3), 0.06), max(ti(4), 0.12)];
    ax.LooseInset = padLBRT;
    try
        ax.YAxis.Exponent = 0;
    catch
        ax.YRuler.Exponent = 0;
    end
    try, ax.YRuler.TickLabelFormat = '%,.0f'; end
end

function onClickVoxel(h,~, zIdx)
    fig = ancestor(h,'figure'); S = guidata(fig);
    ax = ancestor(h,'axes'); if isnan(zIdx), return; end
    pt = get(ax,'CurrentPoint'); xB = pt(1,1); yB = pt(1,2);
    xL = get(ax,'XLim'); yL = get(ax,'YLim');
    nx = S.X; ny = S.Ydim;
    jB = round( ((xB - xL(1)) / (xL(2) - xL(1))) * (nx-1) + 1 );
    iB = round( ((yB - yL(1)) / (yL(2) - yL(1))) * (ny-1) + 1 );
    if iB<1 || iB>ny || jB<1 || jB>nx, return; end
    % Convert back to original array indices (inverse of rot90(A,-1))
    rA = S.X - jB + 1;  cA = iB;
    if rA<1 || rA>S.X || cA<1 || cA>S.Ydim || zIdx<1 || zIdx>S.Z, return; end
    S.roi_center = [rA cA zIdx]; S.roi3d_mask = []; guidata(fig,S);
    plotTSAtCenter(fig); renderPage(fig); applyGlobalUIFont(S);
end

function onBaselineCtrlChanged(~,~)
    fig = gcbf; S = guidata(fig);
    % Debounce mid-typing (non-digits) when Manual
    if get(S.pop_bmode,'Value')==2
        if ~all(isstrprop(strrep(strrep(get(S.edt_b1,'String'),' ',''),'-',''),'digit')) || ...
           ~all(isstrprop(strrep(strrep(get(S.edt_b2,'String'),' ',''),'-',''),'digit'))
            refreshTSLabel(S.axTS); return;
        end
    end
    refreshTSLabel(S.axTS); plotTSAtCenter(fig); recomputeOverlayFromBaseline(fig);
end

function recomputeOverlayFromBaseline(fig)
    % Recompute PSC map with MASK-AWARE denominator and dynamic epsilon
    S = guidata(fig);
    [b1,b2] = get_baseline_window(S);
    b1 = max(1, round(b1)); b2 = min(S.T, round(b2)); if b2 <= b1, b2 = min(S.T, b1+1); end
    base_idx = b1:b2; sig_idx = S.signal_idx;
    if isempty(sig_idx), sig_idx = max(1, round(S.T*0.55)) : min(S.T, round(S.T*0.75)); end

    baseline_mean = mean(S.Y(:,:,:,base_idx),4,'omitnan');
    signal_mean   = mean(S.Y(:,:,:,sig_idx),4,'omitnan');
    den = baseline_mean;

    eps_den = S.eps_den_global;
    tiny = abs(den) < eps_den;
    den(tiny & S.mask) = eps_den;

    pc = ((signal_mean - baseline_mean) ./ den) * 100;
    pc(~S.mask) = NaN;

    pc_abs = abs(pc);
    pc_abs_flat = pc_abs(isfinite(pc_abs) & S.mask);
    if isempty(pc_abs_flat), pc_abs_flat = 0; end
    S.pc_abs_lo = prctile(pc_abs_flat,2); S.pc_abs_hi = prctile(pc_abs_flat,98);
    if S.pc_abs_hi <= S.pc_abs_lo, S.pc_abs_hi = S.pc_abs_lo + eps; end

    S.pc = pc; S.pc_abs = pc_abs; S.U_base = baseline_mean;
    [S.und_lo_base, S.und_hi_base] = robust_window(baseline_mean);

    guidata(fig,S); renderPage(fig); adjustScale(S.sld_scale,[]); applyGlobalUIFont(S);
end

function plotTSAtCenter(fig)
    S = guidata(fig); if isempty(S.roi_center), return; end
    rA = S.roi_center(1); cA = S.roi_center(2); zIdx = S.roi_center(3);
    axTS = S.axTS; cla(axTS);

    % Extract voxel trace
    ts_vox = squeeze(S.Y(rA, cA, zIdx, :));

    % ROI toggle & baseline window
    useROI  = logical(get(S.chk_roi,'Value'));
    plotPSC = logical(get(S.chk_psc,'Value'));
    [b1,b2] = get_baseline_window(S);

    % Build ROI trace if requested
    haveROI = false;
    if useROI
        if isempty(S.roi3d_mask), S.roi3d_mask = roi_mask_3d(S); guidata(fig,S); end
        lin = find(S.roi3d_mask(:) & S.mask(:));
        if ~isempty(lin)
            Y2 = reshape(S.Y, [], S.T);
            ts_roi = mean(Y2(lin,:),1)';    % ROI mean trace
            haveROI = true;
        end
    end

    if plotPSC
        % --- Baseline indices (clamped) ---
        b1c = max(1, min(S.T, round(b1)));
        b2c = max(1, min(S.T, round(b2)));
        if b2c <= b1c, b2c = min(S.T, b1c + 1); end

        if haveROI
            base_roi = mean(ts_roi(b1c:b2c), 'omitnan');
            psc_roi  = 100 * (ts_roi - base_roi) / max(S.eps_den_global, base_roi);
            plot(axTS, 1:S.T, psc_roi, 'LineWidth', 1.8);
        else
            if S.SHOW_VOXEL_PSC
                base_vox = mean(ts_vox(b1c:b2c), 'omitnan');
                psc_vox  = 100 * (ts_vox - base_vox) / max(S.eps_den_global, base_vox);
                plot(axTS, 1:S.T, psc_vox, 'LineWidth', 1.2);
            else
                % ROI-only desired; nothing to plot if no ROI yet
            end
        end

        ylabel(axTS,'PSC (%)');

        % ---- Hardcode y-limits for PSC ----
        if isfield(S,'PSC_YLIM') && numel(S.PSC_YLIM)==2 && all(isfinite(S.PSC_YLIM))
            set(axTS,'YLim', S.PSC_YLIM);
            yl = S.PSC_YLIM;
        else
            yl = ylim(axTS);
        end

        title(axTS, sprintf('PSC @ [%d,%d,%d]%s', rA, cA, zIdx, tern(haveROI,' (ROI avg)','')));

    else
        % Raw signal plotting
        if haveROI
            plot(axTS, 1:S.T, ts_roi, 'LineWidth',1.8); hold(axTS,'on');
            plot(axTS, 1:S.T, ts_vox, '--', 'LineWidth',1.0);
        else
            plot(axTS, 1:S.T, ts_vox, 'LineWidth',1.2);
        end
        ylabel(axTS,'Signal (a.u.)');
        title(axTS, sprintf('Time series @ [%d,%d,%d]%s', rA, cA, zIdx, tern(haveROI,' (ROI avg)','')));
        yl = ylim(axTS);
    end

    grid(axTS,'on'); xlabel(axTS,'Time (in sec)');
    bl = min(S.baseline_idx); bh = max(S.baseline_idx);
    sl = min(S.signal_idx);   sh = max(S.signal_idx);
    patch(axTS, [bl bh bh bl], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.90 1.00], 'FaceAlpha',0.25, 'EdgeColor','none');
    patch(axTS, [sl sh sh sl], [yl(1) yl(1) yl(2) yl(2)], [1.00 0.90 0.90], 'FaceAlpha',0.25, 'EdgeColor','none');
    uistack(findobj(axTS,'Type','line'),'top');

    applyTSFont(fig); 
    applyGlobalUIFont(S);
end

function mask3 = roi_mask_3d(S)
    if isempty(S.roi_center), mask3 = false(S.X,S.Ydim,S.Z); return; end
    rA = S.roi_center(1); cA = S.roi_center(2); zA = S.roi_center(3);
    Rxy = round(get(S.sld_roirXY,'Value'));
    shape = get(S.pop_roishape,'Value'); % 1: Sphere, 2: Cylinder
    Dz = round(get(S.sld_roirZ,'Value'));
    [RR,CC,ZZ] = ndgrid(1:S.X, 1:S.Ydim, 1:S.Z);
    switch shape
        case 1, mask3 = (RR - rA).^2 + (CC - cA).^2 + (ZZ - zA).^2 <= Rxy^2;
        otherwise, mask3 = ((RR - rA).^2 + (CC - cA).^2 <= Rxy^2) & (abs(ZZ - zA) <= Dz);
    end
    mask3 = mask3 & S.mask;
end

function exportROINifti(~,~)
    fig = gcbf; S = guidata(fig);

    if isempty(S.roi_center)
        warndlg('Click a voxel to set the ROI center first.','Export ROI'); return;
    end

    mask3 = roi_mask_3d(S);
    if ~any(mask3(:))
        warndlg('ROI mask is empty. Increase radius/depth or move center.','Export ROI'); return;
    end

    baseName = strtrim(get(S.edt_roiOut,'String'));
    if isempty(baseName), baseName = 'roi_mask'; end
    outPath = fullfile(S.inDir, [baseName '.nii.gz']);

    % Option 1: simplest & robust (no header passed)
    try
        niftiwrite(uint8(mask3), stripGZ(outPath), 'Compressed', true);
        set(S.lbl_status,'String',sprintf('ROI mask saved: %s', outPath));
        return
    catch ME1
        warning('Simple write failed (%s). Trying with sanitized header...', ME1.message);
    end

    % Option 2: use original header, but force 3D + valid time fields
    infoOut = S.niiInfo;                    
    infoOut.ImageSize       = [S.X S.Ydim S.Z];         
    if isfield(infoOut,'PixelDimensions')
        pd = infoOut.PixelDimensions;
        infoOut.PixelDimensions = pd(1:min(3,numel(pd)));
    end
    infoOut.Datatype        = 'uint8';
    infoOut.BitsPerPixel    = 8;
    infoOut.Description     = '3D ROI mask (1=inside, 0=outside)';

    validTU = {'Unknown','Seconds','Milliseconds','Microseconds'};
    if isfield(infoOut,'TimeUnits')
        if ~any(strcmp(infoOut.TimeUnits, validTU))
            infoOut.TimeUnits = 'Seconds';
        end
    end
    if isfield(infoOut,'raw') && isfield(infoOut.raw,'xyzt_units')
        infoOut.raw.xyzt_units = bitand(uint8(infoOut.raw.xyzt_units), uint8(7));
    end

    try
        niftiwrite(uint8(mask3), stripGZ(outPath), infoOut, 'Compressed', true);
        set(S.lbl_status,'String',sprintf('ROI mask saved: %s', outPath));
    catch ME2
        errordlg(sprintf(['Failed to save NIfTI after sanitizing header:\n%s\n\n' ...
                          'Tip: As a last resort, comment out the header entirely and rely on default NIfTI.' ], ME2.message), ...
                 'Export ROI');
    end

    function pathNoGZ = stripGZ(pth)
        if endsWith(lower(pth), '.nii.gz'), pathNoGZ = pth(1:end-3); else, pathNoGZ = pth; end
    end
end

function info3 = sanitizeTo3DHeader(infoIn, imgSize3)
% Create a 3-D NIfTI header from a 4-D one, preserving orientation but
% removing all temporal metadata that causes niftiwrite to error.
    info3 = infoIn;

    % Dimensions
    info3.ImageSize = imgSize3(:).';
    if isfield(info3,'PixelDimensions')
        pd = info3.PixelDimensions;
        if numel(pd) >= 3, info3.PixelDimensions = pd(1:3);
        else,              info3.PixelDimensions = [1 1 1];
        end
    end

    % Datatype for mask
    info3.Datatype     = 'uint8';
    info3.BitsPerPixel = 8;
    info3.Description  = '3D ROI mask (1=inside, 0=outside)';

    % Remove time-unit fields if present
    if isfield(info3,'TimeUnits'), info3 = rmfield(info3,'TimeUnits'); end
    if isfield(info3,'TimeScale'), info3 = rmfield(info3,'TimeScale'); end

    if isfield(info3,'raw')
        info3.raw.dim = int16([3 imgSize3(:).' 1 1 1 1]);   % ndim=3, T=1
        if isfield(info3.raw,'pixdim') && numel(info3.raw.pixdim) >= 5
            qfac = info3.raw.pixdim(1);
            dx = info3.PixelDimensions(1); dy = info3.PixelDimensions(2); dz = info3.PixelDimensions(3);
            info3.raw.pixdim = single([qfac dx dy dz 0 0 0 0]);
        end
        if isfield(info3.raw,'xyzt_units')
            info3.raw.xyzt_units = bitand(uint8(info3.raw.xyzt_units), uint8(7));
        end
        if isfield(info3.raw,'intent_code'), info3.raw.intent_code = int16(0); end
        if isfield(info3.raw,'intent_p1'),   info3.raw.intent_p1   = single(0); end
        if isfield(info3.raw,'intent_p2'),   info3.raw.intent_p2   = single(0); end
        if isfield(info3.raw,'intent_p3'),   info3.raw.intent_p3   = single(0); end
    end
end

function exportROITimeSeries(~,~)
    fig = gcbf; S = guidata(fig);
    if isempty(S.roi_center)
        warndlg('Click a voxel to set the ROI center first.', 'Export ROI TS'); return;
    end
    if isempty(S.roi3d_mask), S.roi3d_mask = roi_mask_3d(S); guidata(fig,S); end
    lin = find(S.roi3d_mask(:) & S.mask(:));
    if isempty(lin)
        warndlg('ROI mask is empty. Increase radius/depth or move center.', 'Export ROI TS'); return;
    end
    Y2 = reshape(S.Y, [], S.T);
    ts = mean(Y2(lin,:),1)';  % ROI mean

    exportPSC = logical(get(S.chk_psc,'Value'));
    if exportPSC
        [b1,b2] = get_baseline_window(S);
        b1c = max(1, min(S.T, round(b1)));
        b2c = max(1, min(S.T, round(b2))); if b2c<=b1c, b2c=min(S.T,b1c+1); end

        base_roi = mean(ts(b1c:b2c), 'omitnan');
        den_roi  = base_roi;
        eps_plot = S.eps_den_global;
        if ~isfinite(den_roi) || abs(den_roi) < eps_plot, den_roi = eps_plot; end

        ts = 100 * (ts - base_roi) / den_roi;  % PSC(t)
    end

    baseName = strtrim(get(S.edt_tsOut,'String'));
    if isempty(baseName), baseName = 'roi_timeseries'; end
    outPath = fullfile(S.inDir, [baseName '.txt']);

    try
        fid = fopen(outPath,'w');
        if fid<0, error('Cannot open file for writing'); end
        if exportPSC, fprintf(fid,'frame,psc\n'); else, fprintf(fid,'frame,signal\n'); end
        for t=1:S.T, fprintf(fid,'%d,%.6g\n', t, ts(t)); end
        fclose(fid);
        msg = sprintf('ROI time series saved: %s', outPath);
        disp(msg); set(S.lbl_status,'String',msg);
    catch ME
        errordlg(sprintf('Failed to save time series:\n%s', ME.message), 'Export ROI TS');
    end
end

% ===== NEW: Save ROI to GroupX.mat =====
function saveROIToGroup(~,~)
    fig = gcbf; S = guidata(fig);

    if isempty(S.roi_center)
        warndlg('Click a voxel to set the ROI center first.','Save to Group'); 
        return;
    end
    if isempty(S.roi3d_mask), S.roi3d_mask = roi_mask_3d(S); guidata(fig,S); end
    mask3 = S.roi3d_mask;
    if ~any(mask3(:))
        warndlg('ROI mask is empty. Increase radius/depth or move center.','Save to Group'); 
        return;
    end

    % Group index and target file
    % Group index and target file (use visible label)
    gIdx   = get(S.pop_group,'Value');
    rawStr = get(S.pop_group,'String');
    if ischar(rawStr) || isstring(rawStr)
        grpList = cellstr(rawStr);
    else
        grpList = rawStr; % already a cell array
    end
    grpName = strtrim(grpList{gIdx});                 % e.g., 'FAP-AAV + FAP'

    % Make a filesystem-safe filename from the label
    grpSlug = regexprep(lower(grpName), '[^\w-]+', '_');  % keep a-z0-9_-
    grpSlug = regexprep(grpSlug, '_+', '_');
    grpSlug = regexprep(grpSlug, '^_|_$', '');
    if isempty(grpSlug), grpSlug = sprintf('group_%d', gIdx); end

    gFile  = fullfile(S.groupRoot, [grpSlug '.mat']);     % /.../AnalysedData/fap-aav_fap.mat


    % --- Build record to append ---
    rec.timestamp      = datestr(now,'yyyy-mm-ddTHH:MM:SS');
    rec.runPath        = S.inDir;
    rec.inNiiPath      = S.inNiiPath;
    rec.imageSize      = [S.X S.Ydim S.Z];
    rec.center_rcz     = S.roi_center;                   % [r c z]
    if get(S.pop_roishape,'Value')==1
    rec.shape = 'Sphere';
else
    rec.shape = 'Cylinder';
end

    rec.rxy            = round(get(S.sld_roirXY,'Value'));
    rec.dz             = round(get(S.sld_roirZ,'Value'));
    [b1,b2]            = get_baseline_window(S);
    rec.baseline_b1b2  = [b1 b2];

    % Save mask as linear indices to be compact; also save size
    rec.mask_size      = size(mask3);
    rec.mask_idx       = find(mask3);

    % ROI mean time series (raw)
    Y2 = reshape(S.Y, [], S.T);
    lin = rec.mask_idx;
    ts_roi = mean(Y2(lin,:),1)';      % column vector
    rec.ts_raw = ts_roi;

    % ROI PSC (using current baseline window; with safe epsilon)
    b1c = max(1, min(S.T, round(b1)));
    b2c = max(1, min(S.T, round(b2))); 
    if b2c <= b1c, b2c = min(S.T, b1c + 1); end
    base_roi = mean(ts_roi(b1c:b2c), 'omitnan');
    den_roi  = base_roi;
    eps_plot = S.eps_den_global;
    if ~isfinite(den_roi) || abs(den_roi) < eps_plot, den_roi = eps_plot; end
    rec.ts_psc = 100 * (ts_roi - base_roi) / den_roi;

    % --- Append to GroupX.mat ---
    try
        if exist(gFile,'file')
            L = load(gFile); 
            if isfield(L,'ROIs')
                ROIs = L.ROIs;
            else
                ROIs = struct([]);
            end
        else
            ROIs = struct([]);
        end

        if isempty(ROIs)
            ROIs = rec;                      
        else
            ROIs(end+1) = rec;               
        end

        if ~exist(S.groupRoot,'dir')
            mkdir(S.groupRoot);
        end
        save(gFile, 'ROIs', '-v7.3');
        % keep the group-file list fresh
        refreshGroupFileList(fig,[]);


        set(S.lbl_status,'String', sprintf('Saved ROI to %s', gFile));
        disp(get(S.lbl_status,'String'));
    catch ME
        errordlg(sprintf('Failed to append ROI to %s\n\n%s', gFile, ME.message), 'Save to Group');
    end
end

% ====== Helpers ======
function out = tern(cond, a, b), if cond, out = a; else, out = b; end, end

function [b1,b2] = get_baseline_window(S)
    if get(S.pop_bmode,'Value') == 1
        b1 = min(S.baseline_idx); b2 = max(S.baseline_idx);
    else
        b1 = str2double(get(S.edt_b1,'String'));
        b2 = str2double(get(S.edt_b2,'String'));
        if isnan(b1) || isnan(b2)
            b1 = min(S.baseline_idx); b2 = max(S.baseline_idx);
        end
    end
end

function root = detectGroupRoot(runDir)
    % Try to find ".../AnalysedData" upwards from the current runDir.
    parts = split(runDir, filesep); parts = parts(:)';
    k = find(strcmp(parts,'AnalysedData'), 1, 'last');
    if ~isempty(k)
        root = strjoin(parts(1:k), filesep);
        if startsWith(runDir, filesep) && ~startsWith(root, filesep)
            root = [filesep root];
        end
    else
        % Fallback (your fixed mount path)
        root = '/Volumes/Extreme_Pro/fMRI/AnalysedData';
    end
end

%% ---- THEME / INTERACTION HELPERS ----
function A = compute_alpha(S, PCslc, Uslc, validMask)
    base  = S.alpha_base; again = get(S.sld_again,'Value'); mode  = get(S.pop_alpha,'Value');
    switch mode
        case 1, modF = ones(size(PCslc));          % Constant
        case 2
            absPC = abs(PCslc); lo = S.pc_abs_lo; hi = S.pc_abs_hi; if hi<=lo, hi=lo+eps; end
            modF = (absPC - lo) / (hi - lo); modF = max(0,min(1,modF));
        case 3
            [ulo,uhi] = robust_window(Uslc); if uhi<=ulo, uhi=ulo+eps; end
            modF = (double(Uslc) - ulo) / (uhi - ulo); modF = max(0,min(1,modF));
        otherwise, modF = ones(size(PCslc));
    end
    A = base * (again .* modF); A(~validMask) = 0; A = max(0,min(1,A));
end

function [lo,hi] = robust_window(U)
    U = U(isfinite(U));
    if isempty(U), lo=0; hi=1;
    else, lo=prctile(U,2); hi=prctile(U,98); if hi<=lo, hi=lo+eps; end
    end
end

function RGB = to_rgb_gray(slice2d,lo,hi)
    if hi<=lo, hi=lo+eps; end
    N = (double(slice2d)-lo)/(hi-lo); N = max(0,min(1,N));
    RGB = repmat(N,[1 1 3]);
end

function cmap = gbhot_diverging(n)
    if mod(n,2)==1, n=n+1; end, half = n/2;
    t = linspace(0,1,half)'; g = (1-t); c = t.*(1-t); b = t;
    neg = [0*g + 0*c, 1*g + 1*c, 0*g + 1*b]; neg = neg .* (1 - t);
    pos = hot(half); cmap = [neg; pos];
end

function toggleTheme(src,~)
    fig = ancestor(src,'figure'); S = guidata(fig);
    S.isDark = logical(get(S.btn_theme,'Value')); applyTheme(S); guidata(fig,S); applyGlobalUIFont(S);
end

function applyTheme(S)
    if S.isDark
        bg=[0.12 0.12 0.12]; fg=[0.95 0.95 0.95]; axbg=[0.10 0.10 0.10];
        ctrlBG=[0.18 0.18 0.18]; ctrlFG=fg; btnBG=[0.22 0.22 0.22]; btnFG=fg; textBG=bg; textFG=fg;
        if ismac, popupBG=[1 1 1]; popupFG=[0 0 0]; else, popupBG=ctrlBG; popupFG=ctrlFG; end
        set(S.btn_theme,'String','Light theme'); borderCol=[0.30 0.30 0.33];
    else
        bg=[1 1 1]; fg=[0 0 0]; axbg=[1 1 1];
        ctrlBG=[1 1 1]; ctrlFG=fg; btnBG=[0.94 0.94 0.94]; btnFG=fg; textBG=bg; textFG=fg;
        popupBG=ctrlBG; popupFG=ctrlFG; set(S.btn_theme,'String','Dark theme'); borderCol=[0.80 0.80 0.85];
    end

    % Root panels
    set(S.fig,'Color',bg); set(S.panelLeft,'BackgroundColor',bg);
    set(S.panelRight,'BackgroundColor',bg); set(S.panelCtrl,'BackgroundColor',bg);

    % Card panels (A/B/C/D)
    cards = {'pnlA','pnlB','pnlC2','pnlC'}; tints = { ...
        tern(S.isDark,[0.16 0.16 0.18],[0.97 0.985 1.00]), ...
        tern(S.isDark,[0.15 0.16 0.19],[0.98 0.975 1.00]), ...
        tern(S.isDark,[0.18 0.16 0.19],[0.985 0.98 1.00]), ...
        tern(S.isDark,[0.17 0.16 0.19],[0.985 0.98 0.97]) ...
    };
    for i=1:numel(cards)
        if isfield(S,cards{i}) && isgraphics(S.(cards{i}))
            set(S.(cards{i}),'BackgroundColor',tints{i}, 'BorderType','line','BorderWidth',1, ...
                              'ForegroundColor',tern(S.isDark,[0.30 0.30 0.33],[0.80 0.80 0.85]), ...
                              'HighlightColor',tern(S.isDark,[0.30 0.30 0.33],[0.80 0.80 0.85]), ...
                              'ShadowColor',tern(S.isDark,[0.30 0.30 0.33],[0.80 0.80 0.85]));
        end
    end

    % Axes & colorbar
    set(S.axTS,'Color',axbg,'XColor',fg,'YColor',fg);
    for k=1:numel(S.ax)
        if isgraphics(S.ax(k))
            set(S.ax(k),'Color',axbg,'XColor','none','YColor','none','XTick',[],'YTick',[],'Box','off');
        end
    end
    try, set(S.cb,'Color',fg); end

    % Global control styling
    allCtrls = [findall(S.panelLeft,'Type','uicontrol'); findall(S.panelRight,'Type','uicontrol')];
    for h = reshape(allCtrls,1,[])
        if ~ishandle(h), continue; end
        st = lower(get(h,'Style'));
        try
            switch st
                case 'text', set(h,'BackgroundColor',tern(S.isDark,bg,[1 1 1]),'ForegroundColor',fg);
                case 'edit', set(h,'BackgroundColor',tern(S.isDark,[0.18 0.18 0.18],[1 1 1]),'ForegroundColor',fg);
                case 'popupmenu', set(h,'ForegroundColor',tern(S.isDark,[0 0 0],fg)); try, set(h,'BackgroundColor',tern(S.isDark,[1 1 1],[1 1 1])); catch, end
                case 'slider', set(h,'BackgroundColor',tern(S.isDark,[0.18 0.18 0.18],[1 1 1]),'ForegroundColor',fg);
                case {'checkbox','radiobutton','pushbutton','togglebutton'}
                    set(h,'BackgroundColor',tern(S.isDark,[0.22 0.22 0.22],[0.94 0.94 0.94]),'ForegroundColor',fg);
                otherwise, set(h,'BackgroundColor',bg,'ForegroundColor',fg);
            end
        catch, end
    end

    % Tint children inside each card
    function tintPanelChildren(pnl, pnlBG)
        if ~isgraphics(pnl), return; end
        kids = findall(pnl,'Type','uicontrol');
        for hh = reshape(kids,1,[])
            if ~ishandle(hh), continue; end
            st = lower(get(hh,'Style')); 
            try
                switch st
                    case 'text'
                        set(hh,'BackgroundColor',pnlBG,'ForegroundColor',fg);
                    case {'checkbox','radiobutton','pushbutton','togglebutton'}
                        set(hh,'BackgroundColor',pnlBG,'ForegroundColor',fg);
                    case 'edit'
                        set(hh,'BackgroundColor',tern(S.isDark,[0.18 0.18 0.18],[1 1 1]),'ForegroundColor',fg);
                    case 'popupmenu'
                        set(hh,'ForegroundColor',tern(S.isDark,[0 0 0],fg)); if ~ismac, try, set(hh,'BackgroundColor',pnlBG); catch, end, end
                    case 'slider'
                        set(hh,'BackgroundColor',pnlBG,'ForegroundColor',fg);
                    otherwise
                        set(hh,'BackgroundColor',pnlBG,'ForegroundColor',fg);
                end
            catch, end
        end
    end
    tintPanelChildren(S.pnlA,  tints{1});
    tintPanelChildren(S.pnlB,  tints{2});
    tintPanelChildren(S.pnlC2, tints{3});
    tintPanelChildren(S.pnlC,  tints{4});

    % Keep top-row controls above tiles
    uistack(S.btn_theme,'top'); uistack(S.lbl_sign,'top'); uistack(S.pop_sign,'top');

    % Recolor slice labels
    fgTxt = tern(S.isDark, [0.95 0.95 0.95], [0 0 0]);
    if isfield(S,'hSliceLbl') && ~isempty(S.hSliceLbl)
        for k = 1:numel(S.hSliceLbl)
            if isgraphics(S.hSliceLbl(k)), try, set(S.hSliceLbl(k),'Color',fgTxt); end, end
        end
    end
end

function onMouseMove(fig,~)
    S = guidata(fig);
    if ~isstruct(S) || ~isfield(S,'ax') || isempty(S.ax) || ~all(ishandle(S.ax)), return; end
    obj = hittest(fig); ax  = ancestor(obj, 'axes');
    if isempty(ax) || ~isscalar(ax), if isfield(S,'lbl_status') && ishandle(S.lbl_status), set(S.lbl_status,'String','Ready'); end; return; end
    [tf, k] = ismember(ax, S.ax(:));
    if ~tf || k<1 || k>numel(S.currentSlices) || S.currentSlices(k)==0
        if isfield(S,'lbl_status') && ishandle(S.lbl_status), set(S.lbl_status,'String','Ready'); end, return;
    end
    idx = S.currentSlices(k);
    pt = get(ax,'CurrentPoint'); xB = pt(1,1); yB = pt(1,2);
    xL = get(ax,'XLim'); yL = get(ax,'YLim'); nx = S.X; ny = S.Ydim;
    jB = round(((xB-xL(1))/(xL(2)-xL(1)))*(nx-1)+1);
    iB = round(((yB-yL(1))/(yL(2)-yL(1)))*(ny-1)+1);
    if iB<1 || iB>ny || jB<1 || jB>nx
        set(S.lbl_status,'String',sprintf('Slice #%d (out of bounds)',idx)); return;
    end
    rA = S.X - jB + 1;  cA = iB;
    if rA<1 || rA>S.X || cA<1 || cA>S.Ydim, set(S.lbl_status,'String','Ready'); return; end
    switch get(S.pop_under,'Value')
        case 1, Uval = S.U_mean4D(rA,cA,idx);
        case 2, Uval = S.U_base(rA,cA,idx);
        otherwise, tp = max(1, min(S.T, S.tp_current)); Uval = S.Y(rA,cA,idx,tp);
    end
    PCval = S.pc(rA,cA,idx); Acur = compute_alpha(S, PCval, Uval, true); msk = S.mask(rA,cA,idx);
    set(S.lbl_status,'String', sprintf('[%d %d %d] U=%.3f %%Δ=%.3f α=%.2f mask=%d', rA, cA, idx, Uval, PCval, Acur, msk));
end

function onScroll(fig, evt)
    S = guidata(fig); step = -evt.VerticalScrollCount;  % wheel up = +1
    mods = get(fig,'CurrentModifier'); % {'control','alt','shift'}
    if any(strcmp(mods,'control'))
        v = get(S.sld_scale,'Value') + 0.5*step; v = min(max(v,0.1), get(S.sld_scale,'Max'));
        set(S.sld_scale,'Value',v); adjustScale(S.sld_scale,[]);
    else
        if step>0, nextPage(fig,[]); else, prevPage(fig,[]); end
    end
end

function onKeyPress(fig, evt)
    S = guidata(fig);
    switch lower(evt.Key)
        case 'rightarrow', nextPage(fig,[]);
        case 'leftarrow',  prevPage(fig,[]);
        case 'uparrow',    v = get(S.sld_scale,'Value') + 0.5; set(S.sld_scale,'Value',min(v,get(S.sld_scale,'Max'))); adjustScale(S.sld_scale,[]);
        case 'downarrow',  v = get(S.sld_scale,'Value') - 0.5; set(S.sld_scale,'Value',max(v,0.1)); adjustScale(S.sld_scale,[]);
        case 'o', set(S.chk_overlay,'Value',~get(S.chk_overlay,'Value')); changeOverlay(S.chk_overlay,[]);
        case 'a', set(S.chk_showA,'Value',~get(S.chk_showA,'Value')); adjustScale(S.sld_scale,[]); renderPage(fig);
        case 'p', set(S.chk_psc,'Value',~get(S.chk_psc,'Value')); refreshTSLabel(S.axTS); plotTSAtCenter(fig);
        case 'add',       v = get(S.sld_scale,'Value') + 0.5; set(S.sld_scale,'Value',min(v,get(S.sld_scale,'Max'))); adjustScale(S.sld_scale,[]);
        case 'subtract',  v = get(S.sld_scale,'Value') - 0.5; set(S.sld_scale,'Value',max(v,0.1)); adjustScale(S.sld_scale,[]);
        case 'bracketleft'
            if strcmp(get(S.sld_roirZ,'Enable'),'on'), v = max(0, get(S.sld_roirZ,'Value')-1); set(S.sld_roirZ,'Value',v); roiChanged(S.sld_roirZ,[]); end
        case 'bracketright'
            if strcmp(get(S.sld_roirZ,'Enable'),'on'), v = min(get(S.sld_roirZ,'Max'), get(S.sld_roirZ,'Value')+1); set(S.sld_roirZ,'Value',v); roiChanged(S.sld_roirZ,[]); end
    end
end

%% ===== Centralized font helper =====
function applyGlobalUIFont(S)
    targetFont = S.UI_FONT; minSize = S.UI_FONTSIZE;
    handles = [findall(S.panelLeft,'Type','uicontrol'); findall(S.panelRight,'Type','uicontrol'); findall(S.panelCtrl,'Type','uicontrol')];
    for h = reshape(handles,1,[])
        if ~ishandle(h), continue; end
        try, set(h,'FontName',targetFont); fs = get(h,'FontSize'); if isempty(fs), fs = minSize; end; set(h,'FontSize', max(minSize, fs)); catch, end
    end
    axAll = [S.ax(:); S.axTS];
    for ax = reshape(axAll,1,[])
        if ~ishandle(ax), continue; end
        try
            set(ax,'FontName',targetFont); fs = get(ax,'FontSize'); if isempty(fs), fs = minSize; end
            set(ax,'FontSize', max(minSize-2, fs));
            t = get(ax,'Title');  if isgraphics(t), set(t,'FontName',targetFont,'FontSize',max(minSize-2,get(t,'FontSize'))); end
            xl = get(ax,'XLabel');if isgraphics(xl), set(xl,'FontName',targetFont,'FontSize',max(minSize-2,get(xl,'FontSize'))); end
            yl = get(ax,'YLabel');if isgraphics(yl), set(yl,'FontName',targetFont,'FontSize',max(minSize-2,get(yl,'FontSize'))); end
            try
                ax.XAxis.FontName = targetFont; ax.YAxis.FontName = targetFont;
                ax.XAxis.FontSize = max(minSize-4, ax.XAxis.FontSize);
                ax.YAxis.FontSize = max(minSize-4, ax.YAxis.FontSize);
            catch, end
        catch, end
    end
    if isfield(S,'cb') && isgraphics(S.cb)
        try, set(S.cb,'FontName',targetFont); set(S.cb,'FontSize',max(minSize-2, get(S.cb,'FontSize'))); catch, end
    end
end

%% ===== Responsive layout manager =====
function onResize(fig,~)
    S = guidata(fig); if isempty(S) || ~isfield(S,'fig') || ~ishandle(S.fig), return; end
    relayoutUI(S);
end

function relayoutUI(S)
    % layout constants derived from font size
    F = S.UI_FONTSIZE;
    pad = round(0.6*F); gap = round(0.5*F); ctrlH = round(1.6*F); rowH  = max(ctrlH, 24);

    % Figure size
    figpos = getpixelposition(S.fig); W = figpos(3); H = figpos(4);

    % Bottom control panel height
    rowsNeeded = max([6,7,7,9]); % Increased ROI panel rows to include Group row
    statusH = rowH; navH = rowH;
    ctrlPanelH = pad + rowsNeeded*(rowH+gap) + pad + navH + pad + statusH + pad;
    ctrlPanelH = min(max(ctrlPanelH, 4.2*rowH + 4*pad), max(0.40*H, ctrlPanelH));

    % Main panels area
    topY = pad; bottomY = ctrlPanelH + pad; mainH = H - bottomY - pad;

    % Split left/right width
    extraGap = 40;  % space between images and graph

    leftW  = round(0.64*(W - 3*pad - extraGap));
    rightW = W - 3*pad - leftW - extraGap;
    leftX  = pad; 
    rightX = 2*pad + leftW + extraGap;

    % Set main panel positions
    set(S.panelLeft, 'Position', [leftX, bottomY, leftW, mainH]);
    set(S.panelRight,'Position', [rightX, bottomY, rightW, mainH]);

    % Layout top row controls on panelLeft
    layoutTopRow(S, leftX, bottomY, leftW, mainH, pad, gap, rowH);

    % Layout right panel (time series only)
    layoutRightPanel(S, rightX, bottomY, rightW, mainH, pad, gap, rowH);

    % Bottom control panel
    set(S.panelCtrl, 'Position', [pad, topY, W - 2*pad, ctrlPanelH]);
    layoutPanelCtrl(S, pad, topY, W - 2*pad, ctrlPanelH, pad, gap, rowH);
    drawnow limitrate;
end

function layoutTopRow(S, x, y, w, h, pad, gap, rowH)
    p  = getpixelposition(S.panelLeft); Lx = 1;
    Ly = p(4) - (rowH + pad) + ceil(0.30*rowH); Ly = min(max(Ly, 1), p(4) - rowH);
    Lw = p(3);
    bw = max(120, 7.5*S.UI_FONTSIZE); set(S.btn_theme,'Position',[Lx+pad, Ly, bw, rowH]);
    ddw = max(260, 12*S.UI_FONTSIZE); labw = max(180, 9*S.UI_FONTSIZE);
    set(S.pop_sign,'Position',[Lx+Lw - pad - ddw, Ly, ddw, rowH]);
    set(S.lbl_sign,'Position',[Lx+Lw - pad - ddw - gap - labw, Ly, labw, rowH]);
    uistack(S.btn_theme,'top'); uistack(S.lbl_sign,'top'); uistack(S.pop_sign,'top');
end

function layoutRightPanel(S, x, y, w, h, pad, gap, rowH)
    % Fixed gutters
    leftG   = 50; rightG  = 28; topG    = 28; bottomG = 40;
    innerX = pad + leftG; innerY = pad + bottomG;
    innerW = max(60, w - 2*pad - (leftG + rightG));
    innerH = max(60, h - 2*pad - (topG  + bottomG));
    set(S.axTS,'Units','pixels','Position', [innerX, innerY, innerW, innerH], 'PositionConstraint','innerposition');
    try, S.axTS.YAxis.Exponent = 0; catch, S.axTS.YRuler.Exponent = 0; end
end

function layoutPanelCtrl(S, x, y, w, h, pad, gap, rowH)
    statusH = rowH; navH = rowH;
    statusY = y + pad;
    navY = statusY + statusH + pad;
    areaY = navY + navH + pad;
    areaH = h - ( (areaY - y) + pad ); areaH = max(areaH, rowH);

    set(S.lbl_status,'Position',[x + pad, statusY, w - 2*pad, statusH]);

    % 4 columns for A / B / C(PSC) / D(ROI)
    colGap = gap; colW = floor( (w - 2*pad - 3*colGap) / 4 ); baseX = x + pad;
    set(S.pnlA, 'Position',[baseX + 0*(colW + colGap), areaY, colW, areaH]);
    set(S.pnlB, 'Position',[baseX + 1*(colW + colGap), areaY, colW, areaH]);
    set(S.pnlC2,'Position',[baseX + 2*(colW + colGap), areaY, colW, areaH]);
    set(S.pnlC, 'Position',[baseX + 3*(colW + colGap), areaY, colW, areaH]);

    layoutPanelA(S, S.pnlA,  pad, colGap, rowH);
    layoutPanelB(S, S.pnlB,  pad, colGap, rowH);
    layoutPanelPSC(S, S.pnlC2, pad, colGap, rowH);
    layoutPanelD_ROI(S, S.pnlC,  pad, colGap, rowH);

    % Navigation row centered
    btnW = max(140, 7.5*S.UI_FONTSIZE); pageW = max(220, 10*S.UI_FONTSIZE);
    totalW = btnW + gap + pageW + gap + btnW;
    leftX = x + (w - totalW)/2;
    set(S.btn_prev,'Position',[leftX,                     navY, btnW,  rowH]);
    set(S.lbl_page,'Position',[leftX + btnW + gap,        navY, pageW, rowH], 'HorizontalAlignment','center');
    set(S.btn_next,'Position',[leftX + btnW + gap + pageW + gap, navY, btnW,  rowH]);
end

function layoutPanelA(S, pnl, pad, gap, rowH)
    p = getpixelposition(pnl); x0 = pad; w = p(3) - 2*pad;
    y = p(4) - pad - rowH; set(S.lblA_title,'Position',[x0, y, w, rowH]);

    y = y - (rowH + gap);
    labW = max(140, 7*S.UI_FONTSIZE);
    ddW_same = max(200, 11*S.UI_FONTSIZE);
    set(S.lblA_amode,'Position',[x0, y, labW, rowH]);
    set(S.pop_alpha, 'Position',[x0 + labW + gap, y, ddW_same, rowH]);

    y = y - (rowH + gap);
    set(S.lblA_again,'Position',[x0, y, labW, rowH]);
    sldW = w - (labW + gap + 80);
    set(S.sld_again, 'Position',[x0 + labW + gap, y + round(0.25*rowH), sldW, round(0.5*rowH)]);

    y = y - (rowH + gap);
    set(S.lblA_scale,'Position',[x0, y, labW, rowH]);
    set(S.sld_scale,'Position',[x0 + labW + gap, y + round(0.25*rowH), sldW, round(0.5*rowH)]);
    set(S.lbl_scale,'Position',[x0 + labW + gap + sldW + gap, y, 80, rowH]);

    y = y - (rowH + gap);
    set(S.lblB_win,'Position',[x0, y, labW, rowH]);
    set(S.pop_win, 'Position',[x0 + labW + gap, y, ddW_same, rowH]);

    y = y - (rowH + gap);
    halfW = floor((w - gap) / 2) - gap;
    set(S.sld_pLo,'Position',[x0, y + round(0.25*rowH), halfW, round(0.5*rowH)]);
    set(S.sld_pHi,'Position',[x0 + halfW + gap, y + round(0.25*rowH), halfW, round(0.5*rowH)]);
end

function layoutPanelB(S, pnl, pad, gap, rowH)
    p = getpixelposition(pnl); x0 = pad; w = p(3) - 2*pad;
    y = p(4) - pad - rowH; set(S.lblB_title,'Position',[x0, y, w, rowH]);

    y = y - (rowH + gap);
    labW = max(120, 6*S.UI_FONTSIZE); ddW  = max(160, 9*S.UI_FONTSIZE);
    set(S.lblB_tiles,'Position',[x0, y, labW, rowH]);
    set(S.pop_tiles, 'Position',[x0 + labW + gap, y, ddW, rowH]);

    y = y - (rowH + gap);
    set(S.lblB_under,'Position',[x0, y, labW, rowH]);
    set(S.pop_under, 'Position',[x0 + labW + gap, y, ddW, rowH]);

    y = y - (rowH + gap);
    tpLabW = max(40, 2*S.UI_FONTSIZE); tpEditW = max(100, 5*S.UI_FONTSIZE);
    set(S.lbl_tp,'Position',[x0, y, tpLabW, rowH]);
    set(S.edt_tp,'Position',[x0 + tpLabW + gap, y + max(0,round(0.2*rowH)), tpEditW, max(24,round(0.8*rowH))]);

    y = y - (rowH + gap);
    set(S.chk_showA,'Position',[x0, y, w - 2*pad, rowH]);

    y = y - (rowH + gap);
    set(S.chk_overlay,'Position',[x0, y, w - 2*pad, rowH]);
end

function layoutPanelPSC(S, pnl, pad, gap, rowH)
    p = getpixelposition(pnl); x0 = pad; w = p(3) - 2*pad;
    y = p(4) - pad - rowH; set(S.lblC2_title,'Parent',pnl,'Position',[x0, y, w, rowH]);

    y = y - (rowH + gap);
    set(S.chk_psc,'Parent',pnl,'Position',[x0, y, w, rowH]);

    y = y - (rowH + gap);
    labW = max(220, 11*S.UI_FONTSIZE); ddW  = max(280, 12*S.UI_FONTSIZE);
    set(S.lbl_bmode,'Parent',pnl,'Position',[x0, y, labW, rowH], 'HorizontalAlignment','right');
    set(S.pop_bmode,'Parent',pnl,'Position',[x0 + labW + gap, y, min(ddW, w - labW - gap), rowH]);

    y = y - (rowH + gap);
    editW = max(90, 0.5*max(140, 8*S.UI_FONTSIZE));
    editH = rowH;
    yEdit1 = y + max(0, (rowH - editH)/2);
    labW2 = max(60, 3*S.UI_FONTSIZE);
    set(S.lbl_b1,'Parent',pnl,'Position',[x0, y, labW2, rowH]);
    set(S.edt_b1,'Parent',pnl,'Position',[x0 + labW2 + gap, yEdit1, editW, editH]);

    y = y - (rowH + gap);
    yEdit2 = y + max(0, (rowH - editH)/2);
    set(S.lbl_b2,'Parent',pnl,'Position',[x0, y, labW2, rowH]);
    set(S.edt_b2,'Parent',pnl,'Position',[x0 + labW2 + gap, yEdit2, editW, editH]);

    y = y - (rowH + gap);
    set(S.lbl_bhint,'Parent',pnl,'Position',[x0, y, w, rowH], 'HorizontalAlignment','left');
end

function layoutPanelD_ROI(S, pnl, pad, gap, rowH)
    p = getpixelposition(pnl); x0 = pad; w = p(3) - 2*pad;

    y = p(4) - pad - rowH; set(S.lblC_title,'Position',[x0, y, w, rowH]);

    y = y - (rowH + gap);
    chkW = max(180, 10*S.UI_FONTSIZE);
    ddW  = max(120, 9*S.UI_FONTSIZE);
    set(S.chk_roi,    'Position',[x0,               y, chkW, rowH]);
    set(S.pop_roishape,'Position',[x0 + chkW + gap, y, ddW,  rowH]);

    y = y - (rowH + gap);
    labW = max(120, 6*S.UI_FONTSIZE);
    set(S.lblC_rXY,'Position',[x0, y, labW, rowH]);
    sldW = w - (labW + gap + 50);
    set(S.sld_roirXY,'Position',[x0 + labW + gap, y + round(0.25*rowH), sldW, round(0.5*rowH)]);
    set(S.lbl_roirXY,'Position',[x0 + labW + gap + sldW + gap, y, 50, rowH]);

    y = y - (rowH + gap);
    set(S.lblC_rZ,'Position',[x0, y, labW, rowH]);
    sldW2 = w - (labW + gap + 50);
    set(S.sld_roirZ,'Position',[x0 + labW + gap, y + round(0.25*rowH), sldW2, round(0.5*rowH)]);
    set(S.lbl_roirZ,'Position',[x0 + labW + gap + sldW2 + gap, y, 50, rowH]);

    % NEW: Group row (Assign + Save)
    y = y - (rowH + gap);
    labWg = max(150, 8*S.UI_FONTSIZE);
    ddWg  = max(160, 9*S.UI_FONTSIZE);
    btnWg = max(150, 8*S.UI_FONTSIZE);
    extra = w - (labWg + gap + ddWg + gap + btnWg);
    extra = max(0, extra);
    set(S.lblC_group,  'Position',[x0, y, labWg, rowH]);
    set(S.pop_group,   'Position',[x0 + labWg + gap, y, ddWg, rowH]);
    set(S.btn_groupAdd,'Position',[x0 + labWg + gap + ddWg + gap + extra, y, btnWg, rowH]);

    % Export mask row (.nii.gz)
    y = y - (rowH + gap);
    labWexp = max(150, 8*S.UI_FONTSIZE); 
    btnW = max(110, 6*S.UI_FONTSIZE);
    editW   = max(120, w - (labWexp + btnW + 3*gap));
    set(S.lblC_export, 'String','Mask (.nii.gz):', 'Position',[x0, y, labWexp, rowH]);
    set(S.edt_roiOut,  'Position',[x0 + labWexp + gap, y, editW, rowH]);
    set(S.btn_exportROI,'String','Save', 'Position',[x0 + labWexp + gap + editW + gap, y, btnW, rowH]);

    % Group Avg PSC row
    y = y - (rowH + gap);
    labWga = max(150, 8*S.UI_FONTSIZE);
    ddWga  = max(220, 12*S.UI_FONTSIZE);
    btnWga = max(44,  3*S.UI_FONTSIZE);   % refresh is a small square
    btnWga2= max(150, 8*S.UI_FONTSIZE);   % Plot button
    extra  = w - (labWga + gap + ddWga + gap + btnWga + gap + btnWga2);
    extra  = max(0, extra);

    set(S.lblC_gAvg,   'Position',[x0, y, labWga, rowH]);
    set(S.pop_gFile,   'Position',[x0 + labWga + gap, y, ddWga, rowH]);
    set(S.btn_gRefresh,'Position',[x0 + labWga + gap + ddWga + gap, y, btnWga, rowH]);
    set(S.btn_gPlot,   'Position',[x0 + labWga + gap + ddWga + gap + btnWga + gap + extra, y, btnWga2, rowH]);

end


function refreshGroupFileList(src, ~)
    % Accepts either a UI control/figure handle (callback) or a direct figure handle.
    % Falls back to gcf if needed.

    % Resolve figure handle
    fig = [];
    if nargin >= 1 && ~isempty(src)
        if ishghandle(src)
            if strcmp(get(src,'Type'),'figure')
                fig = src;
            else
                fig = ancestor(src,'figure');
            end
        end
    end
    if isempty(fig) || ~ishandle(fig)
        fig = gcf;   % last resort
    end
    if isempty(fig) || ~ishandle(fig)
        % No valid figure available; silently bail out
        return;
    end

    S = guidata(fig);

    % Discover *.mat under S.groupRoot
    files = {};
    try
        if exist(S.groupRoot,'dir')
            d = dir(fullfile(S.groupRoot, '*.mat'));
            files = sort({d.name});
        end
    catch
        files = {};
    end

    if isempty(files)
        set(S.pop_gFile,'String', {'<no group files found>'}, 'Value',1);
        S.groupFiles = {};
    else
        curVal = get(S.pop_gFile,'Value');
        set(S.pop_gFile,'String', files, 'Value', min(max(1,curVal), numel(files)));
        S.groupFiles = files;
    end

    guidata(fig,S);
end


function plotGroupAvgPSC(~,~)
    fig = gcbf; if isempty(fig), fig = ancestor(gco,'figure'); end
    S = guidata(fig);

    if isempty(S.groupFiles)
        warndlg('No group files detected. Save ROIs to a group first or click refresh (↻).','Group Avg PSC');
        return;
    end

    gIdx = get(S.pop_gFile,'Value');
    gList = get(S.pop_gFile,'String');
    if ischar(gList), gList = cellstr(gList); end
    if gIdx<1 || gIdx>numel(gList) || startsWith(gList{gIdx},'<')
        warndlg('Select a valid group file.','Group Avg PSC'); return;
    end
    gFile = fullfile(S.groupRoot, gList{gIdx});

    % Load ROIs struct array
    [ok, ROIs, msg] = loadGroupMat(gFile);
    if ~ok, errordlg(msg,'Group Avg PSC'); return; end
    if isempty(ROIs)
        warndlg('Selected group has no ROIs saved yet.','Group Avg PSC'); return;
    end

    % Collect ts_psc column vectors; allow different lengths by truncating to min T
    Ts = cellfun(@(r) numel(r), {ROIs.ts_psc});
    Tm = min(Ts);
    if Tm < 3
        errordlg('Time series too short to compute variance/SEM.','Group Avg PSC'); return;
    end

    M = numel(ROIs);
    Mx = zeros(Tm, M);
    for i=1:M
        v = ROIs(i).ts_psc(:);
        v = v(1:Tm);
        Mx(:,i) = v;
    end

    mu  = mean(Mx, 2, 'omitnan');
    se  = std(Mx, 0, 2, 'omitnan') ./ sqrt(sum(isfinite(Mx),2));

    % Plot
    ax = S.axTS; cla(ax,'reset'); hold(ax,'on'); grid(ax,'on');
    % light individual traces
    for i=1:M
        plot(ax, 1:Tm, Mx(:,i), 'LineWidth', 0.8, 'Color', [0.6 0.6 0.6 0.5]); %#ok<*NBRAK>
    end
    % shaded mean ± SEM
    x = (1:Tm)';
    f = fill(ax, [x; flipud(x)], [mu-se; flipud(mu+se)], [0.8 0.88 1.0], ...
             'EdgeColor','none', 'FaceAlpha',0.45);
    uistack(f,'bottom');
    plot(ax, x, mu, 'b-', 'LineWidth', 2.2);

    title(ax, sprintf('Group Avg PSC (%s), n=%d', gList{gIdx}, M), 'Interpreter','none');
    xlabel(ax,'Frame'); ylabel(ax,'PSC (%)');

    % Apply hardcoded PSC limits if provided
    if isfield(S,'PSC_YLIM') && numel(S.PSC_YLIM)==2 && all(isfinite(S.PSC_YLIM))
        ylim(ax, S.PSC_YLIM);
    end

    % Keep fonts/theme consistent
    S.chk_psc.Value = 1;   % ensure PSC label
    guidata(fig,S);
    applyTSFont(fig); 
    applyGlobalUIFont(S);
    applyTheme(S);
    tuneTSAxes(fig);
end

function [ok, ROIs, msg] = loadGroupMat(pathMat)
    ok = false; ROIs = struct([]); msg = '';
    try
        L = load(pathMat, 'ROIs');
        if ~isfield(L,'ROIs') || isempty(L.ROIs)
            msg = 'File has no variable named ROIs.'; return;
        end
        ROIs = L.ROIs;
        % sanity: ensure ts_psc exists
        hasPSC = arrayfun(@(r) isfield(r,'ts_psc') && ~isempty(r.ts_psc), ROIs);
        if ~any(hasPSC)
            msg = 'No ts_psc time series found in ROIs.'; return;
        end
        ROIs = ROIs(hasPSC);
        ok = true;
    catch ME
        msg = sprintf('Failed to read %s:\n%s', pathMat, ME.message);
    end
end



end  % ===== end of main function file =====
