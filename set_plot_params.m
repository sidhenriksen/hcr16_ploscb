function set_plot_params(myFig,varargin);
    % Sets default properties on a figure, labels the subplots, etc.
    % Usage: set_plot_params(myFig,...)
    % myFig: a figure handle.
    % Optional arguments:
    % 'reverse_labels': if set to 1 labels b,a instead of a,b
    % 'labels': toggle on, off
    
    reverse_labels = 0;
    toggle_labels = 'on';
    for j = 1:2:length(varargin);
        
        switch varargin{j};
            case 'reverse_labels';
                reverse_labels = varargin{j+1};
                
            case 'labels'
                toggle_labels = varargin{j+1};
                
                
        end

    end

    children = get(myFig,'children');
    
    % This first bit just gets all the axes (i.e. excludes legends, etc.)
    axis_index = [];
    legend_index = [];
    for j = 1:length(children);
        
        if strcmp(children(j).Type,'axes');
            axis_index = [axis_index,j];
        elseif strcmp(children(j).Type,'legend');
            legend_index = [legend_index,j];
       
        end
        
    end
    if ~reverse_labels;
        axis_index = axis_index(end:-1:1);
        legend_index = legend_index(end:-1:1);
    end
    
    labels = {'a','b','c','d','e','f','g','h','i','j','k','l','m'};
    
    nAxes = length(axis_index);
    nRows = floor(nAxes/2); 
    nRows(nRows == 0) = 1;
    nCols = nAxes/nRows;
    
    labelsize = 20;
    axis_fontsize = 18;
    position = [200,200,nCols*600,400*nRows];
    set(myFig,'Position',position,'color','white');
    
    % Loop over axes, set the appropriate axis parameters
    for ax = 1:nAxes;
        current_axis = children(axis_index(ax));
        
        set(current_axis,'fontsize',axis_fontsize);
        xlab = get(current_axis,'xlabel');
        ylab = get(current_axis,'ylabel');
        set(xlab,'fontsize',labelsize);
        set(ylab,'fontsize',labelsize);
        
        xlim = get(current_axis,'xlim'); ylim = get(current_axis,'ylim');
            
        xrange = range(xlim);
        yrange = range(ylim);
        
        subplot(current_axis);
        
        if nAxes > 1 && strcmpi(toggle_labels,'on');
            text(xlim(1)+xrange*0.05,ylim(2)+yrange*0.05,labels{ax},'fontsize',labelsize);
        end
        
    end
    
    legend_fontsize = 16;
    
    nLegends = length(legend_index);
    for leg = 1:nLegends;
        current_legend = children(legend_index(leg));
        
        set(current_legend,'fontsize',legend_fontsize);
    end
    
    


end