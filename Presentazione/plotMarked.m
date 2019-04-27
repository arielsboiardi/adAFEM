function plotMarked(marked, space)
% evidenzia graficamente gli elementi marcati

idMarked=marked.id;
T=space.T;

for j=idMarked
    plot([T(j),T(j+1)],[0,0],'rs-','LineWidth',3);
end
end