[X,Y,Z] = cylinder([0 .3], 50);
axis([-1 1,-1 1,-.5 .5])
M=makehgtform('translate',[0,0,0],'xrotate',pi/2,'yrotate',pi/2);
Mp = makehgtform('translate', [0,0,0], 'xrotate', pi/2, 'yrotate',-pi/2);
h=surf(X,Y,Z,'Parent',hgtransform('Matrix',M),'LineStyle','none','FaceAlpha',0.4);
hold on
hp = surf(X,Y,Z, 'Parent',hgtransform('Matrix', Mp), 'LineStyle','none', 'FaceAlpha',0.4);
view([30,35])
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on
light

% 
% alpha_store = [alpha_store, alpha];
% alpha = alpha_0;
% while true
%     if((f(x_old,Y, Y_hat, b, lambda)-f((x_old-alpha*v), Y, Y_hat, b, lambda))< -r*alpha*norm(v)^2)
%         break
%     else
%         alpha = tau*alpha;
%     end
% end
