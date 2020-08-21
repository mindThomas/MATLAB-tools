function enablePCDtooltip(pcd)    
    dcm_obj = datacursormode(gcf);
    set(dcm_obj, 'UpdateFcn', @(a,b) myupdatefcn(a,b,pcd))

    function txt = myupdatefcn(~,event_obj,pcd)
        % ~            Currently not used (empty)
        % event_obj    Object containing event data structure
        % output_txt   Data cursor text        
        pos = get(event_obj,'Position');
        idx = get(event_obj,'DataIndex');
        txt = {   ['Index: ',num2str(idx)],...
                  ['X: ',num2str(pos(1))],...
                  ['Y: ',num2str(pos(2))],...
                  ['Z: ',num2str(pos(3))],...        
                  ['intensity: ',num2str(pcd.intensity(idx))]...
                  };                
    end
end
