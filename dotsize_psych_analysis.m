function BigAnswerMatrix = dotsize_psych_analysis();
    inDir = 'dotsize_data/';
    
    MyNames = dir([inDir,'*.csv']);
    for fname = 1:length(MyNames);
        MyFiles{fname} = MyNames(fname).name;
    end
    NumFiles = length(MyNames);% NumFiles = NumFiles(1);
    CustomDisp = [-0.1,0.1];
    MergeFiles = {
        [5,6,7]
        };

    Exclude = [4];
    Merged = [];



    NumMerge = 0;
    for j = 1:length(MergeFiles);
        NumMerge = NumMerge + length(MergeFiles{j});
    end

    Count = length(MyFiles) - length(Exclude) - (NumMerge-1)*(NumMerge>0);
    IncludeCount = 0;


    %AnswerMatrix = zeros(length(Disparities),length(RefreshRates),length(DotMatchLevels),Count);
    %AnswerMatrix = zeros(4,4,3,1);
    %AnswerMatrix = zeros(2,1,3,2); 

    Colours = {'Red';'Green';'Blue';'Black';[1,0,1]; 'Cyan'; [1,1,0]};

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
                    CurrentFile = fopen([inDir,MyFiles{mergefile}],'r');
                    ReadData = cell2mat(textscan(CurrentFile, '%f%f%f%f%f%f%f', 'Delimiter',',', 'HeaderLines',2));
                    CurrentData = [CurrentData; ReadData];
                    Merged = [Merged,mergefile];
                end

            else

                CurrentFile = fopen([inDir,MyFiles{file}],'r');
                CurrentData = cell2mat(textscan(CurrentFile, '%f%f%f%f%f%f%f', 'Delimiter',',', 'HeaderLines',2));
            end


            Disparities = unique(CurrentData(:,1));
            RefreshRates = unique(CurrentData(:,4));
            DotMatchLevels = unique(CurrentData(:,5));
            Densities = unique(CurrentData(:,6));
            DotSizes = unique(CurrentData(:,7));

            if ~IncludeCount
                AnswerMatrix = zeros(length(RefreshRates), length(DotMatchLevels), ...
                    length(Densities),length(DotSizes),Count);
            end

            CorrectAnswer = (CurrentData(:,1) > 0) *2 +1;
            IsCorrect = CurrentData(:,2) == CorrectAnswer;
            for RefRate = 1:length(RefreshRates);
                for DM = 1:length(DotMatchLevels);
                    for density = 1:length(Densities);
                        for dotsize = 1:length(DotSizes);
                            CurrentIndices = ...
                                (CurrentData(:,4) == RefreshRates(RefRate)) .* ...
                                (CurrentData(:,5) == DotMatchLevels(DM)) .* ...
                                (CurrentData(:,6) == Densities(density)) .* ...
                                (CurrentData(:,7) == DotSizes(dotsize));

                            CurrentIndices = logical(CurrentIndices);

                            PropCorrect = sum(IsCorrect(CurrentIndices))/sum(CurrentIndices);

                            AnswerMatrix(RefRate,DM,density,dotsize,IncludeCount+1) = PropCorrect;
                        end
                    end
                end
            end
            IncludeCount = IncludeCount +1;

            %figure(); hold on;
            Initials = MyFiles{file}; AllInitials{IncludeCount} = Initials(1:2);

            Colours = {[1,0,1], [0,0,1], [1,0.6,0],[1,0,0]};
            Legend = {};
            %set(gcf, 'Color','White');


            %xlim([-0.1, 1.1]);
            %ylim([0,1]);
        end
    end
    BigAnswerMatrix = squeeze(AnswerMatrix);
end