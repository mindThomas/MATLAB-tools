function saveImage(folder, filename) 
    fig = gcf;
    %fig = figure('Position', get(0, 'Screensize'));
    set(fig, 'Position', get(0, 'Screensize'))
    frame = getframe(fig);
    imwrite(frame.cdata, [folder, filename, '.png'], 'png')
end