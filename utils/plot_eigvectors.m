function p = plot_eigvectors(A)
    [V,D] = eig(A);
    [~,ind] = sort(diag(D));
    Vs = V(:,ind);
    c = {'r','g','b'};
    for i=1:3
        vi = Vs(:,i);
        [az,el,~] = cart2sph(vi(1), vi(2), vi(3));
        [xs,ys,zs] = sph2cart(az, el, linspace(0,1,100));
        plot3(xs,ys,zs, 'LineWidth',1.5 ,'Color', c{i});
        s = scatter3(xs(end),ys(end),zs(end), 'k', 'filled');
        set(get(get(s,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    legend({'$v_1$', '$v_2$', '$v_3 = v^*$'}, 'Interpreter','latex');
end