function fig1();

fig_dims = [800,500];
x = linspace(-0.5,0.5,501);
x2 = linspace(-0.5,0.5,21);

lw=3;
ms=9;

gray = [0.5,0.5,0.5];

[c1,a1] = tcs(x);
[c2,a2] = tcs(x2);

figure(); hold on;
plot(x,c1,'r -','linewidth',lw);
plot(x2,c2,'k o','markersize',ms,'linewidth',2);

plot(x,a1,'k --','linewidth',lw)
plot(x2,a2,'k o','markersize',ms,'markerfacecolor','k','linewidth',2);

plot(x,x*0+0.5,':','color',gray,'linewidth',lw);

set(gca,'ytick',0:0.25:1,'xtick',-0.5:0.25:0.5,'fontsize',16);
xlabel('Disparity','fontsize',20);
ylabel('Firing rate','fontsize',20);
set(gcf,'color','white','position',[100,100,100+fig_dims(1),100+fig_dims(2)])

text(0.15,0.8,sprintf('Correlated\namplitude'),'color','r','fontsize',18);
text(0.1,0.2,sprintf('Anticorrelated\n   amplitude'),'color','k','fontsize',18);
text(-0.45,0.25,sprintf('Uncorrelated\n   baseline'),'color',gray,'fontsize',18)
end

function [c,a] = tcs(x);
    sigma=0.15;
    f = 0.3/sigma;
    g = @(x)(exp(-x.^2./(2*sigma^2)) .* cos(2*pi*f*x));
    
    c0 = g(x);
    a0 = -g(x);

    c = (c0-min(a0))./(max(c0)-min(a0));
    a = (a0-min(a0))./(max(c0)-min(a0));
end
