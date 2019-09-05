function qj_init = initJointsFromNameList(jnt_name_list, ndof, jnt_names_body, qj_init_body)
    if (size(jnt_name_list,1) ~= ndof)
        error('initJointsFromNameList: %s', WBM.wbmErrorMsg.DIM_MISMATCH);
    end

    nJnts_body = size(jnt_names_body,1);
    if (nJnts_body ~= size(qj_init_body,1))
        % the data-arrays have not the same length ...
        error('initJointsFromNameList: %s', WBM.wbmErrorMsg.DIM_MISMATCH);
    end

    % Set the initial joint position values from the name and value lists of the
    % initial body configuration. The order of the position values will be set
    % in the same order as given in the joint-name list of the URDF-file.
    qj_init = zeros(ndof,1);
    for i = 1:nJnts_body
        idx = find(strcmp(jnt_name_list, jnt_names_body{i,1}), 1);
        if ~isempty(idx)
            qj_init(idx,1) = qj_init_body(i,1);
        end
    end
end
