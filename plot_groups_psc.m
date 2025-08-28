function plot_groups_psc
% plot_groups_psc.m
% Layout: 2x2 subplots, all x-axes in minutes
%   (1,1) Original group mean PSC ± SEM
%   (1,2) Smoothed (movmean)
%   (2,1) Smoothed, aligned to injection start (0 min), window [-5,+25] min
%   (2,2) Reserved

    % ---- USER STYLE SETTINGS ----
    PSC_YLIM   = [-5 20];   % y-limits for all PSC subplots
    AX_FONTSIZE = 16;       % axis labels, tick labels, titles
    LEG_FONTSIZE = 12;      % legend text

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

    % ---- Colors & legends ----
    cmap = lines(numel(fn));
    legends = cell(1,numel(fn));

    % ---- Figure ----
    figure('Color','w');
    tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    % For storage
    group_M_list = cell(1,numel(fn));
    group_t_list = cell(1,numel(fn));
    group_Ms_list= cell(1,numel(fn));
    Lmins        = zeros(1,numel(fn));

    % ---------- (1,1) ORIGINAL ----------
    ax1 = nexttile(tl,1); hold(ax1,'on'); grid(ax1,'on');
    xlabel(ax1,'Time (min)','FontSize',AX_FONTSIZE);
    ylabel(ax1,'PSC (%)','FontSize',AX_FONTSIZE);
    title(ax1,'Original mean PSC \pm SEM','FontSize',AX_FONTSIZE+1);
    set(ax1,'LineWidth',1.2,'FontName','Calibri','FontSize',AX_FONTSIZE);

    for gi = 1:numel(fn)
        S = load(fullfile(fp, fn{gi}));
        % Extract PSC
        if isfield(S,'ROIs') && isstruct(S.ROIs)
            Ts={}; ROIs=S.ROIs;
            for k=1:numel(ROIs)
                if isfield(ROIs(k),'ts_psc'), v=ROIs(k).ts_psc(:); Ts{end+1}=v; end
            end
            Lmin=min(cellfun(@numel,Ts));
            M=cell2mat(cellfun(@(v)v(1:Lmin),Ts,'uni',0));
        else
            names=fieldnames(S); best=''; bestSize=0;
            for q=1:numel(names)
                v=S.(names{q});
                if isnumeric(v)&&ismatrix(v)
                    if numel(v)>bestSize, best=v; bestSize=numel(v); end
                end
            end
            M=S.(best); if size(M,1)<size(M,2), M=M.'; end; Lmin=size(M,1);
        end
        mu=mean(M,2,'omitnan'); sem=std(M,0,2,'omitnan')/sqrt(size(M,2));
        t=(1:numel(mu))'/60;
        c=cmap(gi,:);
        fill(ax1,[t;flipud(t)],[mu-sem;flipud(mu+sem)],c,'FaceAlpha',0.15,'EdgeColor','none');
        plot(ax1,t,mu,'Color',c,'LineWidth',2);
        [~,base,~]=fileparts(fn{gi}); legends{gi}=base;
        group_M_list{gi}=M; group_t_list{gi}=(1:Lmin)'; Lmins(gi)=Lmin;
    end
    legend(ax1,legends,'Interpreter','none','FontSize',LEG_FONTSIZE,'Location','best');
    ylim(ax1,PSC_YLIM);

    % ---------- (1,2) SMOOTHED ----------
    ax2 = nexttile(tl,2); hold(ax2,'on'); grid(ax2,'on');
    xlabel(ax2,'Time (min)','FontSize',AX_FONTSIZE);
    ylabel(ax2,'PSC (%)','FontSize',AX_FONTSIZE);
    title(ax2,sprintf('Smoothed (movmean=%d)',win),'FontSize',AX_FONTSIZE+1);
    set(ax2,'LineWidth',1.2,'FontName','Calibri','FontSize',AX_FONTSIZE);

    for gi=1:numel(fn)
        M=group_M_list{gi}; if isempty(M), continue; end
        Ms=movmean(M,win,1,'Endpoints','shrink');
        mu=mean(Ms,2,'omitnan'); sem=std(Ms,0,2,'omitnan')./sqrt(size(Ms,2));
        t=(1:numel(mu))'/60;
        c=cmap(gi,:);
        fill(ax2,[t;flipud(t)],[mu-sem;flipud(mu+sem)],c,'FaceAlpha',0.15,'EdgeColor','none');
        plot(ax2,t,mu,'Color',c,'LineWidth',2);
        group_Ms_list{gi}=Ms;
    end
    legend(ax2,legends,'Interpreter','none','FontSize',LEG_FONTSIZE,'Location','best');
    ylim(ax2,PSC_YLIM);

    % ---------- (2,1) INJECTION-ALIGNED ----------
    ax3 = nexttile(tl,3); hold(ax3,'on'); grid(ax3,'on');
    xlabel(ax3,'Time (min)','FontSize',AX_FONTSIZE);
    ylabel(ax3,'PSC (%)','FontSize',AX_FONTSIZE);
    title(ax3,'Smoothed, aligned to injection (0 min)','FontSize',AX_FONTSIZE+1);
    set(ax3,'LineWidth',1.2,'FontName','Calibri','FontSize',AX_FONTSIZE);

    rel_sec=(-300:1500)'; x_min=rel_sec/60;
    for gi=1:numel(fn)
        Ms=group_Ms_list{gi}; Lmin=Lmins(gi);
        if isempty(Ms), continue; end
        abs_idx=inj_idx_user+rel_sec;
        valid=abs_idx>=1 & abs_idx<=Lmin;
        Ms_win=NaN(numel(rel_sec),size(Ms,2));
        Ms_win(valid,:)=Ms(abs_idx(valid),:);
        mu=mean(Ms_win,2,'omitnan');
        sem=std(Ms_win,0,2,'omitnan')./sqrt(sum(~isnan(Ms_win),2));
        c=cmap(gi,:);
        fill(ax3,[x_min;flipud(x_min)],[mu-sem;flipud(mu+sem)],c,'FaceAlpha',0.15,'EdgeColor','none');
        plot(ax3,x_min,mu,'Color',c,'LineWidth',2);
    end
    yL=PSC_YLIM; ylim(ax3,yL);
    plot(ax3,[0 0],yL,'--','Color',[0 0 0],'LineWidth',1);
    xlim(ax3,[-5 25]); xticks(ax3,-5:5:25);
    legend(ax3,legends,'Interpreter','none','FontSize',LEG_FONTSIZE,'Location','best');

    % ---------- (2,2) Reserved ----------
    ax4 = nexttile(tl,4); axis(ax4,'off');
    text(0.5,0.5,'Reserved','Parent',ax4,'FontSize',AX_FONTSIZE,'FontName','Calibri','HorizontalAlignment','center');

    % Overall title
    sgtitle(tl,'Group PSC Plots','FontName','Calibri','FontSize',AX_FONTSIZE+4);
end
