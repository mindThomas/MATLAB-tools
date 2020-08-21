function plotValues(pcd, valueLabel, rangeMin, rangeMax, varargin)
    if (nargin < 3)
        rangeMin = 0;
        rangeMax = 0;
    end
    
    if (size(varargin) == 0)
        idx = 1:length(pcd.x);
    elseif (size(varargin) == 1)
        idx = varargin{1};
    else
        error('Incorrect arguments');
    end        
        
    [AZ,EL] = view(gca);
    newPlot = (AZ==0 && EL==90);
    eval(['values = pcd.' valueLabel ';']);
    
    if (rangeMin == rangeMax)
        rangeMin = min(values(idx));
        rangeMax = max(values(idx));
    end
    
    if (rangeMin ~= rangeMax)
        [r,g,b] = jetColor(values(idx), rangeMin, rangeMax);
    else
        r = 0 * ones(size(values(idx)));
        g = 0 * ones(size(values(idx)));
        b = 1 * ones(size(values(idx)));
    end
    
    scatter3(pcd.x(idx), pcd.y(idx), pcd.z(idx), 1, [r,g,b]);
    
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