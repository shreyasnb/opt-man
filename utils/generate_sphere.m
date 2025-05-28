function [x,y,z] = generate_sphere(p)
    figure
    [x,y,z]=sphere(p-1);
    h = surf(x,y,z,'FaceAlpha',0.3,'EdgeColor','none', 'FaceColor','yellow','FaceLighting','gouraud');
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    axis equal
    camlight
end