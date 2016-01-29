function fig = alternating_psych_analysis()

            
    inDir = './alternation_data/';
    

    MyFiles = dir([inDir,'*.csv']);


    NumFiles = size(MyFiles); NumFiles = NumFiles(1);
    CustomDisp = [-0.25, -0.125, -0.0625, 0.0625, ...
                0.125, 0.25];
    MergeFiles = {
        [4,5]
        };

    Exclude = [6,7];
    Merged = [];
    Count = NumFiles - length(Exclude) - length(MergeFiles);
    IncludeCount = 0;
    initials = {};
    %Colours = {'Red';'Green';'Blue'; [1,0,1]; 'Black'; 'Cyan'; [1,1,0]};
    Colours = {[1,1,0],'Black','Red','Cyan','Green', 'Blue',[1,0,1]};
    first = 1;
    for file = 1:NumFiles;
        if ~any(Merged == file) && ~any(Exclude == file);
            merge = 0;
            for merger = 1:length(MergeFiles);
                if any(file == MergeFiles{merger});
                    merge = 1;
                    mergewith = setdiff(MergeFiles{merger},file);
                end
            end

            if merge
                CurrentData = [];
                FilesToMerge = [file,mergewith];
                for mergefile = FilesToMerge;
                    CurrentFile = fopen([inDir,MyFiles(file).name],'r');
                    ReadData = cell2mat(textscan(CurrentFile, '%f%f%f%f', 'Delimiter',',', 'HeaderLines',2));
                    CurrentData = [CurrentData; ReadData];
                    Merged = [Merged,mergefile];
                end

            else

                CurrentFile = fopen([inDir,MyFiles(file).name],'r');
                CurrentData = cell2mat(textscan(CurrentFile, '%f%f%f%f', 'Delimiter',',', 'HeaderLines',2));
                
            end


            Disparities = unique(CurrentData(:,1));
            AlternationRates = unique(CurrentData(:,4));
            
            %IsFar = zeros(1,length(Disparities));
            if first
                IsNear = zeros(length(AlternationRates),length(Disparities),Count);
                first = 0;
            end

            for disp = 1:length(Disparities);
                for AR = 1:length(AlternationRates);
                    IsCurrent = (CurrentData(:,1) == Disparities(disp)) .* ...
                        (CurrentData(:,4) == AlternationRates(AR));
                    CurrentResponses = CurrentData(logical(IsCurrent),2);
                    IsNear(AR,disp,IncludeCount+1) = sum(CurrentResponses == 1)/length(CurrentResponses);
                end
            end
            initials{IncludeCount+1} = MyFiles(file).name(1:2);
            IncludeCount = IncludeCount +1;
            
            
            fclose(CurrentFile);
        end
    end
    AlternationRates = AlternationRates/2;
    isCorrect = (IsNear(:,1:3,:)+(1-IsNear(:,4:6,:)))/2;

    fig=figure();
    % So we want two figures;
    % 1) Performance as a function of alt rate for diff subs
    % 2) Performance as a function of alt rate for diff disps
    perfSubs = squeeze(mean(isCorrect,2));
    perfDisps = squeeze(mean(isCorrect,3));


    subColor = {'red','blue','black','green'};
    dispColor = {rand([1,3]),rand([1,3]),rand([1,3])};
    subplot(1,2,2); hold on
    set(gca)
    for j = 1:size(perfSubs,2);
        N = length(CurrentResponses)*size(perfDisps,2);
        x = AlternationRates(3:7); y = perfSubs(3:7,j);
        [L,U] = BinoConf_Score(y*N,N);
        L = y-L; U = U-y;
        
        E = errorbar(log(x),y,L,U);
        set(E,'markersize',10,'markerfacecolor',subColor{j},'markeredgecolor','k',...
            'linewidth',2,'linestyle','-','marker','o','color',subColor{j});
    end
    set(gca,'XTick',log(x),'XTicklabel',x,'YTick',0:0.25:1);
    xlim([min(log(x))-0.25,max(log(x))+0.25])
    ylim([0,1])
    line([1,6],[0.5,0.5],'linewidth',2,'color','k','linestyle','--');
    leg = legend(initials,'location','southeast'); set(leg,...
        'box','off');
    xlabel('Alternation Rate (Hz)')
    ylabel('Proportion correct')
    
    
    
    subplot(1,2,1); hold on
    set(gca);
    disps = abs(Disparities)*0.045;
    disps_leg = {};
    
    
    for j = 1:size(perfDisps,2);
        N = length(CurrentResponses)*size(perfSubs,2);
        x = AlternationRates(3:7); y = perfDisps(3:7,j);
        [L,U] = BinoConf_Score(y*N,N);
        L = y-L; U = U-y;
        
        E = errorbar(log(x),y,L,U);
        set(E,'markersize',10,'markerfacecolor',dispColor{j},'markeredgecolor','k',...
            'linewidth',2,'linestyle','-','marker','o','color',dispColor{j});
        
        disps_leg{j} = ['\pm ',num2str(disps(j)), '^\circ'];
    end
    
    leg2 = legend(disps_leg,'location','southeast'); set(leg2,...
        'box','off');
    set(gca,'XTick',log(x),'XTickLabel',x,'YTick',0:0.25:1);
    xlim([min(log(x))-0.25,max(log(x))+0.25])
    ylim([0,1])
    line([1,6],[0.5,0.5],'linewidth',2,'color','k','linestyle','--');
    xlabel('Alternation Rate (Hz');
    ylabel('Proportion correct')
    
    set_plot_params(gcf,'reverse_labels',1);
   
end
