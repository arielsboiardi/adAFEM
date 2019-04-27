function plotBase(space)
% Rappresenta la base lagrangiana dello spazio space

T=space.T;
hold on
for kdx=1:space.dim
    B=zeros(size(T));
    B(kdx)=1;
    plot(T,B,'-b','LineWidth',2)
end

box on
pbaspect([2.5,1,1])
set(gca,'FontName','Calibri Light','FontSize',12);
end