%% signal_change_map_matlab.m
% Axial percent-change viewer with:
% - Underlay: Mean(4D), Baseline mean, or single Timepoint t
% - Underlay windowing: Global 2–98, Slice 2–98, Manual % (pLo/pHi)
% - Overlay: %Change with sign filter, alpha modulation (Constant/|%|/Underlay) + gain
% - Toggle overlay visibility, Show α map
% - Tiles: 4 (2×2) or 8 (2×4)
% - Click a voxel to plot raw time series or PSC on the right (baseline window configurable)

clear; clc;

%% ---- SETTINGS ----
in_nii          = 'cleaned_mc_func.nii.gz';
baseline_idx    = 350:550;
signal_idx      = 1300:1500;
eps_baseline    = 1e-6;
alpha_base      = 0.6;              % base alpha (multiplied by modulation)
default_tiles   = 8;                % 4 or 8

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
U_mean4D = mean(Y,4);                                   % 3D

% Scale defaults for % change
pc_flat = pc(isfinite(pc)); if isempty(pc_flat), pc_flat = 0; end
pc_lo = prctile(pc_flat,2); pc_hi = prctile(pc_flat,98);
pc_mag_default = max(5, max(abs([pc_lo pc_hi])));

% For alpha-by-|pc| normalization
pc_abs = abs(pc); pc_abs_flat = pc_abs(isfinite(pc_abs));
if isempty(pc_abs_flat), pc_abs_flat = 0; end
pc_abs_lo = prctile(pc_abs_flat,2); pc_abs_hi = prctile(pc_abs_flat,98);
if pc_abs_hi<=pc_abs_lo, pc_abs_hi = pc_abs_lo + eps; end

% Global underlay windows
[und_lo_mean, und_hi_mean] = robust_window(U_mean4D);
[und_lo_base, und_hi_base] = robust_window(baseline_mean);

% Mask from 4D mean
try
    uflat = U_mean4D(isfinite(U_mean4D));
    u_lo = prctile(uflat,2); u_hi = prctile(uflat,98);
    thr = graythresh(mat2gray(U_mean4D));
    brainMask = U_mean4D > (u_lo + thr*(u_hi - u_lo));
catch
    brainMask = U_mean4D > 0;
end

%% ---- FIGURE LAYOUT (left: images; right: time-series) ----
fig = figure('Color','w','Name','Axial %Change Viewer', ...
             'Units','normalized','Position',[0.04 0.04 0.92 0.90]);

% Panels
panelLeft   = uipanel('Parent',fig,'Units','normalized','Position',[0.03 0.18 0.64 0.79], 'BorderType','none');
panelRight  = uipanel('Parent',fig,'Units','normalized','Position',[0.69 0.18 0.28 0.79], 'BorderType','none','BackgroundColor','w');
panelCtrl   = uipanel('Parent',fig,'Units','normalized','Position',[0.03 0.03 0.94 0.12], 'BorderType','none','BackgroundColor','w');

% Right-side time-series axes + PSC controls
axTS = axes('Parent',panelRight,'Position',[0.12 0.25 0.83 0.70]); box(axTS,'on'); grid(axTS,'on');
xlabel(axTS,'Time (frames)'); ylabel(axTS,'Signal (a.u.)'); title(axTS,'Voxel time series (click a slice)');

chk_psc = uicontrol('Parent',panelRight,'Style','checkbox','Units','normalized', ...
    'Position',[0.12 0.14 0.25 0.06], 'String','Plot PSC', 'Value',0, ...
    'BackgroundColor','w', 'Callback',@(~,~)refreshTSLabel(axTS));
uicontrol('Parent',panelRight,'Style','text','Units','normalized', ...
    'Position',[0.39 0.14 0.20 0.06], 'String','Baseline mode:', ...
    'BackgroundColor','w','HorizontalAlignment','right','FontWeight','bold');
pop_bmode = uicontrol('Parent',panelRight,'Style','popupmenu','Units','normalized', ...
    'Position',[0.60 0.145 0.26 0.06], ...
    'String',{'Global baseline_idx','Manual [b1 b2]'}, 'Callback',@(~,~)refreshTSLabel(axTS));
uicontrol('Parent',panelRight,'Style','text','Units','normalized', ...
    'Position',[0.12 0.07 0.10 0.05], 'String','b1:', ...
    'BackgroundColor','w','HorizontalAlignment','right');
edt_b1 = uicontrol('Parent',panelRight,'Style','edit','Units','normalized', ...
    'Position',[0.23 0.07 0.12 0.06], 'String',num2str(min(baseline_idx)), ...
    'Callback',@(~,~)refreshTSLabel(axTS));
uicontrol('Parent',panelRight,'Style','text','Units','normalized', ...
    'Position',[0.38 0.07 0.10 0.05], 'String','b2:', ...
    'BackgroundColor','w','HorizontalAlignment','right');
edt_b2 = uicontrol('Parent',panelRight,'Style','edit','Units','normalized', ...
    'Position',[0.49 0.07 0.12 0.06], 'String',num2str(max(baseline_idx)), ...
    'Callback',@(~,~)refreshTSLabel(axTS));
uicontrol('Parent',panelRight,'Style','text','Units','normalized', ...
    'Position',[0.64 0.07 0.22 0.05], 'String','(Used if Manual)', ...
    'BackgroundColor','w','HorizontalAlignment','left','ForegroundColor',[0.3 0.3 0.3]);

% Slices grid
[slices_per_page, tl, ax, imU, imP] = createTiles(panelLeft, default_tiles);

% Colormap + shared colorbar
cmap = gbhot_diverging(256); colormap(fig, cmap);
try, cb = colorbar(tl,'Location','southoutside'); catch, cb = colorbar('southoutside'); end
cb.Label.String = '% change';

%% ---------- Controls (bottom panel) ----------
% Paging
uicontrol('Parent',panelCtrl,'Style','pushbutton','Units','normalized','Position',[0.01 0.12 0.07 0.76], ...
    'String','◀ Prev','FontWeight','bold','Callback',@prevPage);
uicontrol('Parent',panelCtrl,'Style','pushbutton','Units','normalized','Position',[0.09 0.12 0.07 0.76], ...
    'String','Next ▶','FontWeight','bold','Callback',@nextPage);

% Page label
lbl_page = uicontrol('Parent',panelCtrl,'Style','text','Units','normalized','Position',[0.17 0.12 0.14 0.76], ...
    'String','','BackgroundColor','w','FontWeight','bold','HorizontalAlignment','center','Tag','lbl_page');

% Tiles (4 or 8)
uicontrol('Parent',panelCtrl,'Style','text','Units','normalized','Position',[0.32 0.58 0.05 0.32], ...
    'String','Tiles:','BackgroundColor','w','HorizontalAlignment','right','FontWeight','bold');
valTiles = 1 + (default_tiles==8);
pop_tiles = uicontrol('Parent',panelCtrl,'Style','popupmenu','Units','normalized','Position',[0.38 0.58 0.07 0.32], ...
    'String',{'4 (2×2)','8 (2×4)'}, 'Value', valTiles, 'Callback',@changeTiles);

% Underlay source (Mean, Baseline, Timepoint)
uicontrol('Parent',panelCtrl,'Style','text','Units','normalized','Position',[0.46 0.58 0.09 0.32], ...
    'String','Underlay:','BackgroundColor','w','HorizontalAlignment','right','FontWeight','bold');
pop_under = uicontrol('Parent',panelCtrl,'Style','popupmenu','Units','normalized','Position',[0.56 0.58 0.12 0.32], ...
    'String',{'Mean(4D)','Baseline mean','Timepoint t'}, 'Callback',@changeUnderlay);

% Timepoint slider (for underlay = timepoint)
sld_tp = uicontrol('Parent',panelCtrl,'Style','slider','Units','normalized','Position',[0.69 0.62 0.12 0.22], ...
    'Min',1, 'Max',max(1,T), 'Value',1, 'SliderStep',[1/max(1,T-1) 10/max(1,T-1)], 'Callback',@changeUnderlay, 'Enable','off');
lbl_tp = uicontrol('Parent',panelCtrl,'Style','text','Units','normalized','Position',[0.82 0.58 0.06 0.32], ...
    'String','t=1','BackgroundColor','w','HorizontalAlignment','left');

% Windowing (Global/Slice/Manual)
uicontrol('Parent',panelCtrl,'Style','text','Units','normalized','Position',[0.32 0.12 0.09 0.32], ...
    'String','Underlay win:','BackgroundColor','w','HorizontalAlignment','right','FontWeight','bold');
pop_win = uicontrol('Parent',panelCtrl,'Style','popupmenu','Units','normalized','Position',[0.42 0.12 0.11 0.32], ...
    'String',{'Global 2–98','Slice 2–98','Manual %'}, 'Callback',@changeWindowing);
sld_pLo = uicontrol('Parent',panelCtrl,'Style','slider','Units','normalized','Position',[0.54 0.16 0.10 0.22], ...
    'Min',0, 'Max',50, 'Value',2, 'SliderStep',[1/50 5/50], 'Callback',@changeWindowing);
sld_pHi = uicontrol('Parent',panelCtrl,'Style','slider','Units','normalized','Position',[0.65 0.16 0.10 0.22], ...
    'Min',50, 'Max',100, 'Value',98, 'SliderStep',[1/50 5/50], 'Callback',@changeWindowing);

% Scale slider (symmetric ±v)
uicontrol('Parent',panelCtrl,'Style','text','Units','normalized','Position',[0.76 0.12 0.07 0.32], ...
    'String','Scale (±%)','BackgroundColor','w','HorizontalAlignment','right','FontWeight','bold');
sld_scale = uicontrol('Parent',panelCtrl,'Style','slider','Units','normalized','Position',[0.84 0.16 0.10 0.22], ...
    'Min',0.1, 'Max',15, 'Value',min(15, max(1, pc_mag_default)), 'SliderStep',[1/150 1/15], 'Callback',@adjustScale);
lbl_scale = uicontrol('Parent',panelCtrl,'Style','text','Units','normalized','Position',[0.95 0.12 0.02 0.32], ...
    'String',sprintf('%.1f',get(sld_scale,'Value')),'BackgroundColor','w','HorizontalAlignment','left');

% --- Overlay controls ---
uicontrol('Parent',panelCtrl,'Style','text','Units','normalized','Position',[0.01 0.58 0.05 0.32], ...
    'String','Sign:','BackgroundColor','w','HorizontalAlignment','right','FontWeight','bold');
pop_sign = uicontrol('Parent',panelCtrl,'Style','popupmenu','Units','normalized','Position',[0.07 0.58 0.07 0.32], ...
    'String',{'±','+ only','− only'}, 'Callback',@changeOverlay);
uicontrol('Parent',panelCtrl,'Style','text','Units','normalized','Position',[0.15 0.58 0.07 0.32], ...
    'String','Alpha:','BackgroundColor','w','HorizontalAlignment','right','FontWeight','bold');
pop_alpha = uicontrol('Parent',panelCtrl,'Style','popupmenu','Units','normalized','Position',[0.23 0.58 0.07 0.32], ...
    'String',{'Constant','|%Change|','Underlay'}, 'Callback',@changeOverlay);
sld_again = uicontrol('Parent',panelCtrl,'Style','slider','Units','normalized','Position',[0.31 0.62 0.05 0.22], ...
    'Min',0, 'Max',1.5, 'Value',1.0, 'TooltipString','Alpha gain', 'Callback',@changeOverlay);
chk_overlay = uicontrol('Parent',panelCtrl,'Style','checkbox','Units','normalized','Position',[0.23 0.12 0.07 0.32], ...
    'String','Show overlay','Value',1,'BackgroundColor','w','Callback',@changeOverlay);
chk_showA = uicontrol('Parent',panelCtrl,'Style','checkbox','Units','normalized','Position',[0.14 0.12 0.08 0.32], ...
    'String','Show α map','BackgroundColor','w','Callback',@changeOverlay);
lbl_alpha = uicontrol('Parent',panelCtrl,'Style','text','Units','normalized','Position',[0.54 0.58 0.14 0.32], ...
    'String','α[min med max]=[— — —]','BackgroundColor','w','HorizontalAlignment','center');

%% ---- Pack state ----
S.Y = Y; S.X = X; S.Ydim = Ydim; S.Z = Z; S.T = T;
S.pc = pc; S.pc_abs = pc_abs; S.mask = brainMask;
S.U_mean4D = U_mean4D; S.U_base = baseline_mean;
S.und_lo_mean=und_lo_mean; S.und_hi_mean=und_hi_mean;
S.und_lo_base=und_lo_base; S.und_hi_base=und_hi_base;
S.pc_abs_lo = pc_abs_lo; S.pc_abs_hi = pc_abs_hi;
S.alpha_base = alpha_base;

S.fig = fig;
S.panelLeft = panelLeft; S.panelRight = panelRight; S.panelCtrl = panelCtrl;
S.tl = tl; S.ax = ax; S.imU = imU; S.imP = imP; S.cb = cb; S.axTS = axTS;

S.slices_per_page = slices_per_page; S.startSlice = 1;
S.lbl_page = lbl_page;
S.pop_tiles = pop_tiles;
S.pop_under = pop_under; S.sld_tp = sld_tp; S.lbl_tp = lbl_tp;
S.pop_win = pop_win; S.sld_pLo = sld_pLo; S.sld_pHi = sld_pHi;
S.sld_scale = sld_scale; S.lbl_scale = lbl_scale;
S.pop_sign = pop_sign; S.pop_alpha = pop_alpha; S.sld_again = sld_again;
S.chk_overlay = chk_overlay; S.chk_showA = chk_showA;
S.lbl_alpha = lbl_alpha;

% PSC controls
S.chk_psc = chk_psc; S.pop_bmode = pop_bmode; S.edt_b1 = edt_b1; S.edt_b2 = edt_b2;
S.baseline_idx = baseline_idx; S.signal_idx = signal_idx;

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

function [spp, tl, ax, imU, imP] = createTiles(parentPanel, tiles)
    delete(findall(parentPanel,'Type','axes'));
    delete(findall(parentPanel,'Type','tiledlayout'));
    if tiles==4
        tl = tiledlayout(parentPanel,2,2,'Padding','compact','TileSpacing','compact'); spp = 4;
    else
        tl = tiledlayout(parentPanel,2,4,'Padding','compact','TileSpacing','compact'); spp = 8;
    end
    ax  = gobjects(spp,1); imU = gobjects(spp,1); imP = gobjects(spp,1);
    for k = 1:spp
        ax(k) = nexttile(tl,k); hold(ax(k),'on'); axis(ax(k),'image'); axis(ax(k),'off');
        set(ax(k),'ButtonDownFcn',@(h,e)onClickVoxel(h,e,NaN));
    end
end

function changeTiles(src,~)
    fig = ancestor(src,'figure'); S = guidata(fig);
    val = get(S.pop_tiles,'Value'); tiles = 4; if val==2, tiles=8; end
    [spp, tl, ax, imU, imP] = createTiles(S.panelLeft, tiles);
    S.slices_per_page = spp; S.tl = tl; S.ax = ax; S.imU = imU; S.imP = imP;
    delete(S.cb); try, S.cb = colorbar(S.tl,'Location','southoutside'); catch, S.cb=colorbar('southoutside'); end
    guidata(fig,S); renderPage(fig);
end

function prevPage(src,~)
    fig = ancestor(src,'figure'); S = guidata(fig);
    S.startSlice = max(S.startSlice - S.slices_per_page, 1);
    guidata(fig,S); renderPage(fig);
end

function nextPage(src,~)
    fig = ancestor(src,'figure'); S = guidata(fig);
    S.startSlice = min(S.startSlice + S.slices_per_page, max(S.Z - S.slices_per_page + 1, 1));
    guidata(fig,S); renderPage(fig);
end

function changeUnderlay(src,~)
    fig = ancestor(src,'figure'); S = guidata(fig);
    if get(S.pop_under,'Value')==3, set(S.sld_tp,'Enable','on'); else, set(S.sld_tp,'Enable','off'); end
    set(S.lbl_tp,'String',sprintf('t=%d', round(get(S.sld_tp,'Value'))));
    renderPage(fig);
end

function changeWindowing(src,~)
    fig = ancestor(src,'figure'); renderPage(fig);
end

function changeOverlay(src,~)
    fig = ancestor(src,'figure'); renderPage(fig);
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
        S.cb.TickLabels = {num2str(-v,'%.0f'),'0',num2str(v,'%.0f')};
        S.cb.Label.String = '% change';
        set(S.lbl_scale,'String',sprintf('%.1f',v));
    end
    guidata(fig,S);
end

function renderPage(fig)
    S = guidata(fig);
    if ~isfield(S,'lbl_page') || ~ishandle(S.lbl_page)
        S.lbl_page = findobj(S.panelCtrl,'Type','uicontrol','Style','text','Tag','lbl_page');
    end

    % Underlay source 3D and global lo/hi
    uSrc = get(S.pop_under,'Value'); % 1 mean, 2 base, 3 timepoint
    switch uSrc
        case 1, U = S.U_mean4D; glo = S.und_lo_mean; ghi = S.und_hi_mean;
        case 2, U = S.U_base;   glo = S.und_lo_base; ghi = S.und_hi_base;
        otherwise
            tp = max(1, min(S.T, round(get(S.sld_tp,'Value'))));
            U = S.Y(:,:,:,tp);
            flat = U(isfinite(U)); if isempty(flat), glo=0; ghi=1; else, [glo,ghi] = robust_window(U); end
    end
    set(S.lbl_tp,'String',sprintf('t=%d', round(get(S.sld_tp,'Value'))));

    % Controls
    sign_mode   = get(S.pop_sign,'Value');   % 1: both, 2:+, 3:−
    win_mode    = get(S.pop_win,'Value');    % 1 global, 2 slice, 3 manual
    pLo         = get(S.sld_pLo,'Value');
    pHi         = get(S.sld_pHi,'Value'); if pHi<=pLo, pHi=pLo+eps; end
    showA       = logical(get(S.chk_showA,'Value'));
    showOverlay = logical(get(S.chk_overlay,'Value'));

    first = S.startSlice; last  = min(first + S.slices_per_page - 1, S.Z);
    slices = first:last;
    set(S.lbl_page,'String',sprintf('Slices %d–%d of %d', first, last, S.Z));

    A_all = [];

    for k = 1:S.slices_per_page
        axk = S.ax(k);
        if k <= numel(slices)
            idx   = slices(k);
            Uslc0 = U(:,:,idx);
            PCslc0 = S.pc(:,:,idx);
            Mslc0 = S.mask(:,:,idx);

            % Display orientation (rot90 -1)
            Uslc  = rot90(Uslc0,-1);
            PCslc = rot90(PCslc0,-1);
            Mslc  = rot90(Mslc0,-1);

            % Sign filter
            switch sign_mode
                case 2, mask_sign = PCslc > 0;
                case 3, mask_sign = PCslc < 0;
                otherwise, mask_sign = true(size(PCslc));
            end
            PCvis = PCslc; PCvis(~mask_sign) = NaN;

            % Underlay window per slice
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

            % Alpha map
            A = compute_alpha(S, PCslc, Uslc, Mslc & mask_sign);
            A_all = [A_all; A(isfinite(A) & (A>0))]; %#ok<AGROW>

            % Draw & wire clicks
            if ~isgraphics(S.imU(k))
                axes(axk); cla(axk); hold(axk,'on'); axis(axk,'image'); axis(axk,'off');
                set(axk,'ButtonDownFcn',@(h,e)onClickVoxel(h,e,idx));
                S.imU(k) = image(axk, to_rgb_gray(Uslc,lo_use,hi_use), ...
                                  'ButtonDownFcn',@(h,e)onClickVoxel(h,e,idx), ...
                                  'HitTest','on','PickableParts','all');
                S.imP(k) = imagesc(axk, PCvis, 'HitTest','off','PickableParts','none');
                title(axk, sprintf('Z = %d', idx));
            else
                set(axk,'ButtonDownFcn',@(h,e)onClickVoxel(h,e,idx));
                set(S.imU(k),'CData', to_rgb_gray(Uslc,lo_use,hi_use), ...
                             'ButtonDownFcn',@(h,e)onClickVoxel(h,e,idx), ...
                             'HitTest','on','PickableParts','all');
                set(S.imP(k),'CData', PCvis, 'HitTest','off','PickableParts','none');
                title(axk, sprintf('Z = %d', idx));
            end

            % Overlay visibility / mode
            if showOverlay
                set(S.imP(k),'Visible','on');
                if showA
                    set(S.imP(k),'CData', A, 'AlphaData', 1, 'CDataMapping','scaled'); caxis(axk,[0 1]);
                else
                    set(S.imP(k),'CData', PCvis, 'AlphaData', A, 'CDataMapping','scaled');
                end
            else
                set(S.imP(k),'Visible','off');
            end

            set(axk,'Visible','on');
        else
            cla(axk); set(axk,'Visible','off');
        end
    end

    % Alpha stats + colorbar label
    if isempty(A_all), amin=NaN; amed=NaN; amax=NaN; else, amin=min(A_all); amed=median(A_all); amax=max(A_all); end
    set(S.lbl_alpha,'String',sprintf('α[min med max]=[%.2f %.2f %.2f]',amin,amed,amax));

    if showOverlay && showA
        S.cb.Label.String = 'alpha'; S.cb.Limits=[0 1]; S.cb.Ticks=[0 .5 1];
    else
        S.cb.Label.String = '% change'; adjustScale(S.sld_scale,[]);
    end

    guidata(fig,S);
end

function refreshTSLabel(axTS)
    fig = ancestor(axTS,'figure'); S = guidata(fig);
    if logical(get(S.chk_psc,'Value')), ylabel(axTS,'PSC (%)'); else, ylabel(axTS,'Signal (a.u.)'); end
end

function onClickVoxel(h,evt, zIdx)
    fig = ancestor(h,'figure'); S = guidata(fig);
    ax = ancestor(h,'axes'); if isnan(zIdx), return; end

    % Click position in axes coords (robust to zoom/pan/YDir)
    pt = get(ax,'CurrentPoint'); xB = pt(1,1); yB = pt(1,2);
    xL = get(ax,'XLim'); yL = get(ax,'YLim');
    nx = S.X; ny = S.Ydim;

    % Map to pixel indices on displayed slice (rotated)
    jB = round( ( (xB - xL(1)) / (xL(2) - xL(1)) ) * (nx-1) + 1 );
    iB = round( ( (yB - yL(1)) / (yL(2) - yL(1)) ) * (ny-1) + 1 );
    if iB < 1 || iB > ny || jB < 1 || jB > nx, return; end

    % ✅ Inverse of rot90(A,-1):  B(i,j)=A(n - j + 1, i)
    rA = S.X - jB + 1;   % original row (X)
    cA = iB;             % original col (Y)

    if rA<1 || rA>S.X || cA<1 || cA>S.Ydim || zIdx<1 || zIdx>S.Z, return; end

    ts = squeeze(S.Y(rA, cA, zIdx, :));   % T x 1
    if isempty(ts), return; end

    plotPSC = logical(get(S.chk_psc,'Value'));
    [b1,b2] = get_baseline_window(S);

    axTS = S.axTS; cla(axTS);
    if plotPSC
        b1c = max(1, min(S.T, round(b1)));
        b2c = max(1, min(S.T, round(b2)));
        if b2c <= b1c, b2c = min(S.T, b1c+1); end
        base = mean(ts(b1c:b2c), 'omitnan');
        psc = 100 * (ts - base) / max(base, eps);
        plot(axTS, 1:S.T, psc, 'LineWidth',1.2); ylabel(axTS,'PSC (%)');
        title(axTS, sprintf('PSC @ voxel [%d,%d,%d]', rA, cA, zIdx));
    else
        plot(axTS, 1:S.T, ts, 'LineWidth',1.2); ylabel(axTS,'Signal (a.u.)');
        title(axTS, sprintf('Time series @ voxel [%d,%d,%d]', rA, cA, zIdx));
    end
    grid(axTS,'on'); xlabel(axTS,'Time (frames)');

    % Shade baseline/signal windows
    bl = min(S.baseline_idx); bh = max(S.baseline_idx);
    sl = min(S.signal_idx);   sh = max(S.signal_idx);
    yl = ylim(axTS);
    patch(axTS, [bl bh bh bl], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.90 1.00], 'FaceAlpha',0.25, 'EdgeColor','none');
    patch(axTS, [sl sh sh sl], [yl(1) yl(1) yl(2) yl(2)], [1.00 0.90 0.90], 'FaceAlpha',0.25, 'EdgeColor','none');
    uistack(findobj(axTS,'Type','line'),'top');
end

function [b1,b2] = get_baseline_window(S)
    if get(S.pop_bmode,'Value') == 1
        b1 = min(S.baseline_idx); b2 = max(S.baseline_idx);
    else
        b1 = str2double(get(S.edt_b1,'String'));
        b2 = str2double(get(S.edt_b2,'String'));
        if isnan(b1) || isnan(b2), b1 = min(S.baseline_idx); b2 = max(S.baseline_idx); end
    end
end

%% ---- HELPERS ----
function A = compute_alpha(S, PCslc, Uslc, validMask)
    base = S.alpha_base; again = get(S.sld_again,'Value'); mode = get(S.pop_alpha,'Value');
    switch mode
        case 1, modF = ones(size(PCslc));
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
    if mod(n,2)==1, n=n+1; end
    half = n/2; t = linspace(0,1,half)'; g = (1-t); c = t.*(1-t); b = t;
    neg = [0*g + 0*c, 1*g + 1*c, 0*g + 1*b]; neg = neg .* (1 - t);
    pos = hot(half); cmap = [neg; pos];
end
