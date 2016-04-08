function fig1();


x = linspace(-0.5,0.5,501);
x2 = linspace(-0.5,0.5,21);

lw=3;
ms=9;

gray = [0.5,0.5,0.5];

[c1,a1] = tcs(x);
[c2,a2] = tcs(x2);

figure(); 
subplot(1,3,1);
hold on;
plot(x,c1,'r -','linewidth',lw);
plot(x2,c2,'k o','markersize',ms,'linewidth',2);

plot(x,a1,'k --','linewidth',lw)
plot(x2,a2,'k o','markersize',ms,'markerfacecolor','k','linewidth',2);
plot(x,x*0+0.5,':','color',gray,'linewidth',lw);

quiver(0,0.735,0,0.225,'color','red','linewidth',3);
quiver(0,0.735,0,-0.225,'color','red','linewidth',3);

quiver(0,0.265,0,0.225,'color','k','linewidth',3);
quiver(0,0.265,0,-0.225,'color','k','linewidth',3);

set(gca,'ytick',0:0.25:1,'xtick',-0.5:0.25:0.5)
xlabel('Disparity')
ylabel('Firing rate')


text(0.15,0.8,sprintf('Correlated\namplitude'),'color','r','fontsize',15);
text(0.1,0.2,sprintf('Anticorrelated\n   amplitude'),'color','k','fontsize',15);
text(-0.45,0.25,sprintf('Uncorrelated\n   baseline'),'color',gray,'fontsize',15)



subplot(1,3,2);

hold on;
cs = linspace(-1,1,51);
N = 1e6;

C = cs;
M = C; 
M(M<0)=0;

C2 = (((cs+1) * 0.5 ).^2 -0.25)* 1/0.75;

plot(cs,C,'b -','linewidth',4);
plot(cs,M,'--','color',[0.2,0.8,0.3],'linewidth',4)
xlabel('Binocular correlation');
ylabel('Response')
%plot(cs,C2,'k -','linewidth',2); 

leg1=legend('Pure correlation','Matching');
set(leg1,'location','southeast');

ps = linspace(0,1,501);
SSE = zeros(1,length(ps));
for j = 1:length(SSE);
    cm=C*ps(j) + M*(1-ps(j));
    ss=sum((cm-C2).^2);
    SSE(j) = ss;
end
[~,idx] = min(SSE);


subplot(1,3,3); hold on;
plot(cs,(ps(idx)*C + (1-ps(idx))*M),'k -','linewidth',4);
plot(cs,C2,'- m','linewidth',4)
ylim([-1,1]);
xlim([-1,1]);
xlabel('Binocular correlation');
leg2 = legend('Correlation+matching','Simple squaring');
set(leg2,'location','southeast');

set_plot_params(gcf,'fs_scalar',1.1)


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


