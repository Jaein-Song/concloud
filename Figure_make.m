%%VAR PRESET
% for i=1:2
% for j=1:4
% out{i,j}.var=cfad{i,j}.var./cfad_num{i,j}.var;
% end
% end
clf    
close all
    
    varset=out;
    QCf=1; %1for unQCed, 2 for QCed
    abcd={'a)','b)','c)','d)'};
%%FIGURE PRESET
%colormap setting

%font setting
    font_name='times new roman'; %font tame
        fs_ax=20; %font size of axis
        fs_cbt=18;%font size of color bar title
        fs_cb=15; %font size of colorbar
        fs_tl=20; %font size of title
        fs_ab=30; %font size of abcd
%line setting
        lw_ax=2; %line width of the axis
        lw_cb=1; %line width of the colorbar
%figure size setting
    SCS              = get(0,'screensize');
    height           = 1080*0.9;
    width            = 1920*0.9;
    ncol             = 2;
    nrow             = 2;
    figcolor='w';
    MenuBar          = 'figure';
    ToolBar          = 'figure';
    marginl=0.02;
    marginr=0.02;
    margind=0.02;
    marginu=0.02;
    margini=0.01;
    figoh=(1-margind-marginu)/ncol;
    figow=(1-marginr-marginl)/nrow;
%%PRESET END
    if ~exist('ymax')
        ymax=15;
        yintv=2.5;
    end
%%FOR-LOOP PLOTTING
% cont_fig = figure;
for QCf=1:2
        cont_fig = figure('color',figcolor,'PaperSize',[height*4 width*4],'units','pixels','position',[0 0 width height],'menubar',MenuBar,'toolbar',ToolBar);

    for figi=1:4
        

        clear var
    var=varset{QCf,figi}.var;
    x=varset{QCf,figi}.xax;
    xlb=varset{QCf,figi}.xlb;
    xlims=varset{QCf,figi}.xlim;
        if ~isfield(varset{QCf,figi},'cmin')
            cmin=-2;cmax=log10(1);cint=200;cb_int=1;cticks={'0.01','0.1','1','10'};
        else
            cmin=varset{QCf,figi}.cmin;
            cmax=varset{QCf,figi}.cmax;
            cint=varset{QCf,figi}.cint;
            cb_int=varset{QCf,figi}.cbint;
            cticks=varset{QCf,figi}.cticks;
        end
        if ~isfield(varset{QCf,figi},'yax')
            y=h/1000;
            ylb='Height (km)';
            ylims=[0 15];
        else
            y=varset{QCf,figi}.yax;
            ylb=varset{QCf,figi}.ylb;
            ylims=varset{QCf,figi}.ylim;
        end
                colormap(jet(cint));

        ax = axes('Parent',cont_fig);
        hold(ax,'on');
        fposdim=[marginl+figow*(mod(figi+1,2)) margind+(ncol-floor((figi-1)/ncol)-1)*(figoh+margini) figow figoh];
        aposdim=[marginl/2+figow*(mod(figi+1,2)) margind+(ncol-floor((figi-1)/ncol)-1)*(figoh+margini) figow figoh-2*margini];
        set(ax,'OuterPosition',fposdim)
        if exist('contourflag')
            if contourflag==1
                contourf(x,y,var,'linestyle','none','levelstep',(cmax-cmin)/cint);
            else
                                pcolor(x,y,var);shading interp
            end
        else
            if ~exist('CFADflag')
                pcolor(x,y,var);shading interp
            elseif CFADflag==1
                var=var*100;
                var(var==0)=1e-6;

                contourf(x,y,log10(var),'linestyle','none','levelstep',(cmax-cmin)/cint);
            else
                pcolor(x,y,var);shading interp
            end
        end
        if exist('contourline')&&contourline>0
            for ci=1:contourline
                if ci==1
                    clinecolor=[0 0 0];
%                 elseif ci==contourline
%                     clinecolor=[0 0 0];
                else
                    clinecolor=[ci/contourline ci/contourline ci/contourline];
                end
                contour(x,y,var,varset{QCf,figi}.ctl(ci,:),'linecolor',clinecolor,'linewidth',2);
            end
        end

        xlabel(xlb);
        ylabel(ylb)
        box(ax,'on');
                colormap(jet(cint));
        annotation('textbox',aposdim,'String',abcd{figi},'FontName',font_name,'FontSize',fs_ab,'FontWeight','bold','linestyle','none','HorizontalAlignment','left','VerticalAlignment','top')
        % axis(ax,'tight');
        % 나머지 axes 속성 설정
        set(ax,'BoxStyle','full','Layer','top');
        % colorbar 생성
        cb=colorbar('peer',ax,'Ticks',[cmin:cb_int:cmax],...
            'TickLabels',cticks,...
            'LineWidth',lw_cb,...
            'FontName','minion pro','FontWeight','bold',...
            'FontSize',fs_cb);
        %     'TickLabels',{'0.0001','0.001','0.01','0.1','1'},...
        if isfield(varset{QCf,figi},'name')
            title(varset{QCf,figi}.name,'FontWeight','bold','FontName',font_name,'FontSize',fs_tl);
        end
        if ~isfield(varset{QCf,figi},'cb_name')
            title(cb,'Frequency (%)','FontSize',fs_cbt)
        else
            title(cb,varset{QCf,figi}.cb_name,'FontSize',fs_cbt)
        end
        set(ax,'BoxStyle','full','CLim',[cmin cmax],'FontSize',fs_ax,'FontWeight','bold',...
            'FontName',font_name,'Layer','top','LineWidth',lw_ax,'YMinorTick','on','YTick',...
            [0:yintv:ymax],'xlim',xlims,'ylim',ylims);
%             [0 2.5 5 7.5 10 12.5 15],'xlim',xlims,'ylim',ylims);
        if ~exist('CFADflag')
            set(ax,'XTick',[0:4:24]);
        end
    end
    if exist('prefix_filename') 
        printfilename='fig';
        printfilename=strcat(prefix_filename{QCf},printfilename);
    else    
        if QCf==1
            printfilename='unQCed';
        else
            printfilename='QCed';
        end
    end
    print(printfilename,'-dtiffn','-r100')
    print(printfilename,'-depsc','-r100')
    clf
    close all
 end