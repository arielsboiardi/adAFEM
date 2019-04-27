%% Esporta figure successione

gif('BasiAnimate.gif')

for k=2:14
    figure(k)
%     nome=strcat('imgs\base',num2str(k),'.png');
%     print('-dpng','-r600',nome)
    gif
end