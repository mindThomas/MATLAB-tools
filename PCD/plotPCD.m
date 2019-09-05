function plotPCD(pcd, varargin)    
    if (length(varargin) == 0)
        idx = 1:length(pcd.x);
        color = 'g';
    elseif (length(varargin) == 1)
        if (isnumeric(varargin{1}))
            idx = varargin{1};
            color = 'g';
        else
            idx = 1:length(pcd.x);
            color = varargin{1};
        end
    elseif (length(varargin) == 2)
        color = 'g';
        idx = 1:length(pcd.x);
        if (isnumeric(varargin{1}))
            idx = varargin{1};            
        else            
            color = varargin{1};
        end
        if (isnumeric(varargin{2}))
            idx = varargin{2};            
        else            
            color = varargin{2};
        end        
    else
        error('Incorrect arguments');
    end     
    
    [AZ,EL] = view(gca);
    target = camtarget;
    pos = campos;
    newPlot = (AZ==0 && EL==90);
    
    scatter3(pcd.x(idx), pcd.y(idx), pcd.z(idx), 1, color);       
    
    % Configure plot
    axis equal;
    xlabel('x')
    ylabel('y')
    zlabel('z')
    xlim([-100, 100]);
    ylim([-100, 100]);
    zlim([0, 20]);
    set(gcf,'renderer','opengl');
    set(gca, 'Clipping', 'off');
    set(gcf,'Color','k')
    color = get(gcf,'Color');
    set(gca,'Color',color)  
    set(gca,'XColor',color,'YColor',color,'ZColor',color,'TickDir','out')
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    ax.ZAxis.Visible = 'off';             
    
    if (newPlot)
        view([-90, 22]);
    else
        view([AZ,EL]);
    end                    
    
    zoom(5.0);  
end