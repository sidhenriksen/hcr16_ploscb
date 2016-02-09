function fig2()


    figure();
    subplot(1,2,1);
    hmrds = imread('hmrds.png');
    image(flipud(hmrds));
    set(gca,'xtick',[],'ytick',[],'ydir','norm','box','off','visible','off',...
        'ydir','norm');
    
    subplot(1,2,2); hold on;
    
    x = linspace(-1,1,501);
    a = 1; b = 5; theta = 0;
    matching_computation = sigmoid(x,a,b,theta)*0.5+0.5;
    correlation_computation = sigmoid(x,a,b,theta);
    plot(x,matching_computation,'-','linewidth',4,'color',[0.2,0.8,0.2]);
    plot(x,correlation_computation,'-','linewidth',4,'color',[0.2,0.2,0.8]);
    legend('Matching','Correlation','location','southeast');
    set(gca,'xtick',-1:0.5:1,'ytick',0:0.25:1);
    xlabel('Binocular correlation');
    ylabel('Proportion correct');
    set_plot_params(gcf)
    plot([-1,1],[0.5,0.5],'k --','linewidth',2);
    plot([0,0],[0,1],'k --','linewidth',2);
    
    

end

function y = sigmoid(x,a,b,theta);

    y = 1./(1+a*exp(-(x-theta)*b));

end