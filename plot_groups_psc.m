%% 
function plot_groups_psc
% plot_groups_psc.m
% Layout: 2x2 subplots
%   (1,1) Original group mean PSC ± SEM  [x in minutes]
%   (1,2) Smoothed (movmean)             [x in minutes]
%   (2,1) Smoothed, aligned to injection (0 min), window [-5,+25] min  [x in minutes]
%   (2,2) Baseline vs Post (SMOOTHED): per-subject points + mean line + significance (no bars)
%
% Subplot (2,2) stats:
%   - Within each group: paired t-test between baseline (350:550) and post (1300:1500).
%   - Between groups (post only): Welch t-test comparing fap-aav_fap vs each other group.

    % ---- USER STYLE SETTINGS ----
    PSC_YLIM     = [-2 5];    % y-limits for PSC plots (subplots 1–3) and subplot 4
    AX_FONTSIZE  = 16;        % axis labels, tick labels, titles
    LEG_FONTSIZE = 12;        % legend text (used in subplots 1–3)

    % ---- Group comparison windows (indices; 1 Hz) ----
    BASE_WIN = 350:550;        % baseline window (in samples)
    POST_WIN = 1300:1500;      % post-injection window (in samples)

    % ---- Choose group files ----
    [fn, fp] = uigetfile('*.mat','Select group .mat files','MultiSelect','on');
    if isequal(fn,0), disp('No files selected.'); return; end
    if ischar(fn), fn = {fn}; end

    % ---- Ask smoothing window ----
    defWin = '21';
    answ = inputdlg({'Moving mean window (samples):'}, 'Smoothing', 1, {defWin});
    if isempty(answ), return; end
    win = max(1, round(str2double(answ{1})));

    % ---- Ask injection start ----
    defInj = '600';
    answ2 = inputdlg({'Injection start timepoint (samples; 1 sample = 1 sec):'}, ...
                     'Injection Start', 1, {defInj});
    if isempty(answ2), return; end
    inj_idx_user = round(str2double(answ2{1}));
    if ~isfinite(inj_idx_user) || inj_idx_user < 1
        errordlg('Injection index must be a positive integer.','Input Error'); return;
    end

    % ---- Colors ----
    cmap = lines(numel(fn));                 % for subplots 1–3 (per-group curves)
    baselineColor = [0.20 0.55 0.90];        % fixed Baseline color (subplot 4)
    postColor     = [0.90 0.40 0.25];        % fixed Post-injection color (subplot 4)

    legends = cell(1,numel(fn));

    % ---- Figure ----
    figure('Color','w');
    tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    % For storage across panels
    group_M_list  = cell(1,numel(fn));   % original matrices (time x subjects)
    group_Ms_list = cell(1,numel(fn));   % smoothed matrices
    group_mu_list = cell(1,numel(fn));   % original group mean curves
    Lmins         = zeros(1,numel(fn));  % lengths after trimming
    group_names   = cell(1,numel(fn));
    plotted_mask  = false(1,numel(fn));

    % ---------- (1,1) ORIGINAL ----------
    ax1 = nexttile(tl,1); hold(ax1,'on'); grid(ax1,'on');
    xlabel(ax1,'Time (min)','FontSize',AX_FONTSIZE);
    ylabel(ax1,'PSC (%)','FontSize',AX_FONTSIZE);
    title(ax1,'Original mean PSC \pm SEM','FontSize',AX_FONTSIZE+1);
    set(ax1,'LineWidth',1.2,'FontName','Calibri','FontSize',AX_FONTSIZE);

    h1_lines = gobjects(1,numel(fn));

    for gi = 1:numel(fn)
        S = load(fullfile(fp, fn{gi}));

        % --- Extract PSC per subject ---
        if isfield(S,'ROIs') && isstruct(S.ROIs)
            ROIs = S.ROIs; Ts = {};
            for k=1:numel(ROIs)
                if isfield(ROIs(k),'ts_psc') && isnumeric(ROIs(k).ts_psc)
                    v = ROIs(k).ts_psc(:);
                    if ~isempty(v) && all(isfinite(v)), Ts{end+1} = v; end %#ok<AGROW>
                end
            end
            if isempty(Ts)
                warning('No ts_psc in %s', fn{gi}); continue;
            end
            Lmin = min(cellfun(@numel, Ts));
            M = zeros(Lmin, numel(Ts));
            for c=1:numel(Ts), M(:,c) = Ts{c}(1:Lmin); end
        else
            % Fallback: largest 2D numeric matrix assumed time x subjects
            names = fieldnames(S); bestName=''; bestSize=0;
            for q=1:numel(names)
                v=S.(names{q});
                if isnumeric(v)&&ismatrix(v)
                    if numel(v)>bestSize, bestName=names{q}; bestSize=numel(v); end
                end
            end
            if isempty(bestName), warning('No suitable matrix in %s', fn{gi}); continue; end
            M = S.(bestName); if size(M,1)<size(M,2), M = M.'; end; Lmin=size(M,1);
        end

        mu  = mean(M,2,'omitnan');
        sem = std(M,0,2,'omitnan') ./ sqrt(size(M,2));
        tmin = (1:numel(mu))'/60;

        c = cmap(gi,:);
        hfill = fill(ax1,[tmin;flipud(tmin)],[mu-sem;flipud(mu+sem)],c,'FaceAlpha',0.15,'EdgeColor','none');
        % Hide fill from legend
        set(get(get(hfill,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

        [~,base,~] = fileparts(fn{gi});
        h1_lines(gi) = plot(ax1,tmin,mu,'Color',c,'LineWidth',2,'DisplayName',base);

        legends{gi}      = base;
        group_names{gi}  = base;
        group_M_list{gi} = M;
        group_mu_list{gi}= mu;
        Lmins(gi)        = Lmin;
        plotted_mask(gi) = true;
    end
    legend(ax1,h1_lines(plotted_mask),legends(plotted_mask), ...
           'Interpreter','none','FontSize',LEG_FONTSIZE,'Location','best');
    ylim(ax1,PSC_YLIM); box(ax1,'on');

    % ---------- (1,2) SMOOTHED ----------
    ax2 = nexttile(tl,2); hold(ax2,'on'); grid(ax2,'on');
    xlabel(ax2,'Time (min)','FontSize',AX_FONTSIZE);
    ylabel(ax2,'PSC (%)','FontSize',AX_FONTSIZE);
    title(ax2,sprintf('Smoothed (movmean=%d)',win),'FontSize',AX_FONTSIZE+1);
    set(ax2,'LineWidth',1.2,'FontName','Calibri','FontSize',AX_FONTSIZE);

    h2_lines = gobjects(1,numel(fn));

    for gi=1:numel(fn)
        M = group_M_list{gi}; if isempty(M), continue; end
        Ms = movmean(M,win,1,'Endpoints','shrink');
        mu_s  = mean(Ms,2,'omitnan');
        sem_s = std(Ms,0,2,'omitnan')./sqrt(size(Ms,2));
        tmin = (1:numel(mu_s))'/60;
        c = cmap(gi,:);
        hfill = fill(ax2,[tmin;flipud(tmin)],[mu_s-sem_s;flipud(mu_s+sem_s)],c,'FaceAlpha',0.15,'EdgeColor','none');
        set(get(get(hfill,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

        h2_lines(gi) = plot(ax2,tmin,mu_s,'Color',c,'LineWidth',2,'DisplayName',group_names{gi});
        group_Ms_list{gi} = Ms;
    end
    legend(ax2,h2_lines(plotted_mask),legends(plotted_mask), ...
           'Interpreter','none','FontSize',LEG_FONTSIZE,'Location','best');
    ylim(ax2,PSC_YLIM); box(ax2,'on');

    % ---------- (2,1) INJECTION-ALIGNED (smoothed) ----------
    ax3 = nexttile(tl,3); hold(ax3,'on'); grid(ax3,'on');
    xlabel(ax3,'Time (min)','FontSize',AX_FONTSIZE);
    ylabel(ax3,'PSC (%)','FontSize',AX_FONTSIZE);
    title(ax3,'Smoothed, aligned to injection (0 min)','FontSize',AX_FONTSIZE+1);
    set(ax3,'LineWidth',1.2,'FontName','Calibri','FontSize',AX_FONTSIZE);

    h3_lines = gobjects(1,numel(fn));

    rel_sec = (-300:1500)';  % -5..+25 minutes (in seconds)
    x_min   = rel_sec/60;
    for gi=1:numel(fn)
        Ms   = group_Ms_list{gi};
        Lmin = Lmins(gi);
        if isempty(Ms) || Lmin<=0, continue; end
        abs_idx = inj_idx_user + rel_sec;              % vector of indices
        valid   = (abs_idx >= 1) & (abs_idx <= Lmin);  % element-wise AND
        Ms_win  = NaN(numel(rel_sec), size(Ms,2));
        Ms_win(valid,:) = Ms(abs_idx(valid),:);
        mu_a  = mean(Ms_win,2,'omitnan');
        sem_a = std(Ms_win,0,2,'omitnan') ./ sqrt(sum(~isnan(Ms_win),2));
        c = cmap(gi,:);
        hfill = fill(ax3,[x_min;flipud(x_min)],[mu_a-sem_a;flipud(mu_a+sem_a)],c,'FaceAlpha',0.15,'EdgeColor','none');
        set(get(get(hfill,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

        h3_lines(gi) = plot(ax3,x_min,mu_a,'Color',c,'LineWidth',2,'DisplayName',group_names{gi});
    end
    yL = PSC_YLIM; ylim(ax3,yL);
    plot(ax3,[0 0],yL,'--','Color',[0 0 0],'LineWidth',1);
    xlim(ax3,[-5 25]); xticks(ax3,-5:5:25);
    legend(ax3,h3_lines(plotted_mask),legends(plotted_mask), ...
           'Interpreter','none','FontSize',LEG_FONTSIZE,'Location','best');
    box(ax3,'on');

    % ---------- (2,2) POINTS + MEAN LINE + SIGNIFICANCE (use SMOOTHED PSC) ----------
    ax4 = nexttile(tl,4); hold(ax4,'on'); grid(ax4,'on');
    xlabel(ax4,'Groups','FontSize',AX_FONTSIZE);
    ylabel(ax4,'PSC (%)','FontSize',AX_FONTSIZE);
    title(ax4,'Baseline vs Post (SMOOTHED): subject max + mean line + significance','FontSize',AX_FONTSIZE);
    set(ax4,'LineWidth',1.2,'FontName','Calibri','FontSize',AX_FONTSIZE);

    nG = numel(fn);
    xg = 1:nG;             % group positions
    halfWidth = 0.35/2;    % cluster half-spacing
    line_half = 0.12;      % half-length for the mean line segments
    dy = diff(PSC_YLIM);

    % Store post values for between-group tests (SMOOTHED)
    post_vals_all = cell(1,nG);
    base_vals_all = cell(1,nG);

    % Find pivot group (fap-aav_fap)
    pivot_idx = NaN;
    for gi = 1:nG
        [~,b,~] = fileparts(fn{gi});
        if contains(lower(b), 'fap-aav_fap')
            pivot_idx = gi; break;
        end
    end

    % Plot points + mean lines; do within-group paired t-test (SMOOTHED)
    for gi = 1:nG
        if ~plotted_mask(gi), continue; end

        % Use SMOOTHED PSC
        Ms = group_Ms_list{gi};
        % Fallback (defensive)
        if isempty(Ms)
            Mraw = group_M_list{gi};
            if isempty(Mraw), continue; end
            Ms = movmean(Mraw,win,1,'Endpoints','shrink');
        end

        Lmin = size(Ms,1);
        if Lmin <= 0, continue; end

        % Windows inside data range (same BASE/POST windows, trimmed)
        b_idx = max(1,BASE_WIN(1)) : min(Lmin, BASE_WIN(end));
        p_idx = max(1,POST_WIN(1)) : min(Lmin, POST_WIN(end));
        if isempty(b_idx) || isempty(p_idx)
            warning('Group %s: window(s) exceed data length (%d). Skipping.', group_names{gi}, Lmin);
            continue;
        end

        % --- Per-subject MAX over time in each window (SMOOTHED) ---
        base_vals = max(Ms(b_idx,:), [], 1, 'omitnan');   % 1×Nsubjects
        post_vals = max(Ms(p_idx,:), [], 1, 'omitnan');

        base_vals_all{gi} = base_vals;
        post_vals_all{gi} = post_vals;

        % --- Scatter points (jittered) ---
        xb = (xg(gi) - halfWidth) + (rand(size(base_vals)) - 0.5) * (halfWidth*1.2);
        xp = (xg(gi) + halfWidth) + (rand(size(post_vals)) - 0.5) * (halfWidth*1.2);
        scatter(ax4, xb, base_vals, 25, 'MarkerFaceColor', baselineColor, 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.8);
        scatter(ax4, xp, post_vals, 25, 'MarkerFaceColor', postColor,     'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.8);

        % --- Horizontal mean line across subjects for each condition ---
        base_mean = mean(base_vals, 'omitnan');
        post_mean = mean(post_vals, 'omitnan');
        plot(ax4, [xg(gi)-halfWidth - line_half, xg(gi)-halfWidth + line_half], ...
                  [base_mean,                  base_mean], '-', 'Color', baselineColor, 'LineWidth', 2);
        plot(ax4, [xg(gi)+halfWidth - line_half, xg(gi)+halfWidth + line_half], ...
                  [post_mean,                  post_mean], '-', 'Color', postColor,     'LineWidth', 2);

        % --- Within-group paired t-test (baseline vs post, same subjects; SMOOTHED) ---
        try
            [~, p_within] = ttest(base_vals, post_vals);   % paired by default
        catch
            p_within = NaN;
        end
        stars = p_to_stars(p_within);
        if ~isempty(stars)
            % bracket from baseline center to post center, at a bit above local max
            y_local = max([base_vals(:); post_vals(:)], [], 'omitnan');
            y_level = y_local + 0.06*dy;
            draw_sig(ax4, xg(gi)-halfWidth, xg(gi)+halfWidth, y_level, stars);
        end
    end

    % --- Between-group: fap-aav_fap SMOOTHED post vs others (Welch's t-test) ---
    if ~isnan(pivot_idx) && ~isempty(post_vals_all{pivot_idx})
        pivot_post = post_vals_all{pivot_idx}(:);
        % global top for stacking the between-group brackets
        all_post = cellfun(@(v) v(:), post_vals_all(plotted_mask), 'uni', 0);
        all_post = vertcat(all_post{:});
        y_top = max(all_post, [], 'omitnan');
        if isempty(y_top) || isnan(y_top), y_top = PSC_YLIM(2); end
        baseY = y_top + 0.10*dy;
        stepY = 0.08*dy;  % vertical spacing between brackets

        layer = 0;
        for gi = 1:nG
            if gi == pivot_idx || ~plotted_mask(gi) || isempty(post_vals_all{gi}), continue; end
            other_post = post_vals_all{gi}(:);

            try
                [~, p_between] = ttest2(pivot_post, other_post, 'Vartype','unequal');
            catch
                p_between = NaN;
            end
            stars = p_to_stars(p_between);
            if ~isempty(stars)
                y_level = baseY + layer*stepY;
                draw_sig(ax4, xg(pivot_idx)+halfWidth, xg(gi)+halfWidth, y_level, stars);
                layer = layer + 1;
            end
        end
    end

    set(ax4,'XTick',xg(plotted_mask),'XTickLabel',group_names(plotted_mask));
    ylim(ax4, PSC_YLIM);
    % -- No legend on subplot 4 --
    box(ax4,'on');

    % Overall title
    sgtitle(tl,'Group PSC Plots','FontName','Calibri','FontSize',AX_FONTSIZE+4);

    % ===== Local helper functions =====
    function stars = p_to_stars(p)
        if ~isfinite(p) || isempty(p), stars = ''; return; end
        if p < 1e-3, stars = '***';
        elseif p < 1e-2, stars = '**';
        elseif p < 5e-2, stars = '*';
        else, stars = '';
        end
    end

    function draw_sig(ax, x1, x2, y, stars)
        % draws a bracket from x1 to x2 at height y with star label
        cap = 0.015*diff(PSC_YLIM);  % small vertical cap
        plot(ax, [x1 x1 x2 x2], [y-cap y y y-cap], 'k-', 'LineWidth', 1);
        text(ax, mean([x1 x2]), y + 0.01*diff(PSC_YLIM), stars, ...
            'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold');
    end
end
