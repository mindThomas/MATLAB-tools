function plotTransform(tf, varargin)
   
    origin = tf(1:3,4);
    x = origin + tf(1:3,1);
    y = origin + tf(1:3,2);
    z = origin + tf(1:3,3);
    
    hold on;
    plot3([origin(1),x(1)],[origin(2),x(2)],[origin(3),x(3)],'r','linewidth',2); 
    plot3([origin(1),y(1)],[origin(2),y(2)],[origin(3),y(3)],'g','linewidth',2); 
    plot3([origin(1),z(1)],[origin(2),z(2)],[origin(3),z(3)],'b','linewidth',2); 
    hold off;

    if (size(varargin) == 1)
        text(origin(1),origin(2),origin(3),['   ' varargin{1}],'HorizontalAlignment','left','FontSize',8,'color','w');
    end

end