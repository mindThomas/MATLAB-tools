function [jnt_lnk_names, nLnks, jnt_pair_idx, nJntPairs, shp_size_sf, foot_jnt_idx] = getSimConfigFromURDF(urdf_file_name, lnk_sf_list, sf)
    if isempty(urdf_file_name)
        error('getSimConfigFromURDF: %s', WBM.wbmErrorMsg.EMPTY_STRING);
    end

    if ~exist('sf', 'var')
        % use the default scale factor for the shape sizes ...
        sf = 0.03;
    end

    % Try to convert the URDF-file of the robot model into a Matlab structure for easy access to the data:
    % Note: This function uses an open 3rd-party method of Wouter Falkena from the
    %       Delft University of Technology.
    %       Source: <https://www.mathworks.com/matlabcentral/fileexchange/28518-xml2struct>
    try
        xml_tree = xml2struct(urdf_file_name);
    catch
        strError = sprintf('Failed to read the URDF file \"%s\"!', urdf_file_name);
        error('getSimConfigFromURDF: %s', strError);
    end

    if ~strcmp(xml_tree.robot.link{1,1}.Attributes.name, 'base_link')
       error('getSimConfigFromURDF: The order of the links is different! Verify the URDF-file!');
    end

    nLnks    = size(xml_tree.robot.link,2);
    nJnts    = size(xml_tree.robot.joint,2);
    nJntLnks = nLnks - 1;

    % create the name list for the joint-links:
    jnt_lnk_names = cell(nJntLnks,1);
    for i = 2:nLnks
        jnt_lnk_names{i-1,1} = xml_tree.robot.link{1,i}.Attributes.name;
    end
    if isempty( find(strcmp(jnt_lnk_names, 'com'), 1) )
        jnt_lnk_names{nLnks,1} = 'com'; % add the CoM to the list ...
    end

    % create an adjacency matrix of the connected joints from the name list (without 'com'):
    jnt_adj = zeros(nJntLnks,nJntLnks);
    for i = 1:nJntLnks % for each row
        for j = 1:nJntLnks % for each column
            child_lnk = xml_tree.robot.joint{1,j}.child.Attributes.link;

            if strcmp(child_lnk, jnt_lnk_names{i,1})
                parent_lnk = xml_tree.robot.joint{1,j}.parent.Attributes.link;
                idx = find(strcmp(jnt_lnk_names, parent_lnk), 1);

                jnt_type = xml_tree.robot.joint{1,j}.Attributes.type;
                switch jnt_type
                    case 'revolute'
                        jnt_adj(i,idx) = 1;
                    case 'prismatic'
                        jnt_adj(i,idx) = 2;
                    case 'fixed'
                        jnt_adj(i,idx) = 3;
                    otherwise
                        error('getSimConfigFromURDF: %s', WBM.wbmErrorMsg.UNKNOWN_JNT_TYPE);
                end
            end
        end
    end

    % create the joint-pair index table for the xyz-positions of the joints in the 3D-space:
    % note: this table is needed to build up the skeleton of the robot.
    jnt_pair_idx = uint8(zeros(nJnts,6));
    for i = 1:nJnts
        [res, idx] = find((jnt_adj(i,1:nJntLnks) > 0), 1);
        if res
            % store the indices for the joint-positions ...
            % joints 1 & 2:               x1  x2  y1  y2  z1  z2
            jnt_pair_idx(i,1:6) = horzcat(i, idx, i, idx, i, idx);
            %jnt_pair_idx(i,1:6) = horzcat(idx, i, idx, i, idx, i);
        end
    end

    % create the list with the scale factors for the link-shapes of the robot's body:
    shp_size_sf = ones(nJnts,2) * sf;
    for i = 1:size(lnk_sf_list,1)
        % replace the current values with them from the list ...
        idx = find(strcmp(jnt_lnk_names(1:nLnks-1,1), lnk_sf_list{i,1}), 1);
        shp_size_sf(idx,1:2) = lnk_sf_list{i,2};
    end

    % search and delete those rows from the index-matrix which are completely zero,
    % because these joints are not connected with the tree (skeleton):
    i = 1;
    nJntPairs = nJnts;
    while (i ~= nJntPairs)
        if (sum(jnt_pair_idx(i,1:6) > 0) == 0)
            % the i-th row has only zeros: --> delete this row from the matrix ...
            jnt_pair_idx = jnt_pair_idx([1:i-1,i+1:nJntPairs],1:6);
            % delete the same row also in the shape-size matrix ...
            shp_size_sf  = shp_size_sf([1:i-1,i+1:nJntPairs],1:2);
            % update the number of index-pairs ...
            nJntPairs = size(jnt_pair_idx,1);
        end
        i = i + 1;
    end
    % check the number joint-pairs (a tree has n-1 edges + without 'com'):
    if (nJntPairs ~= (nLnks - 2))
        error('getSimConfigFromURDF: The model has more edges than vertices! Verify the URDF-file!');
    end

    % get the indices of those joints, where the feet are connected to them:
    foot_jnt_idx = uint8(zeros(1,2));
    foot_jnt_idx(1,1) = find(strcmp(jnt_lnk_names, 'l_sole'), 1);
    foot_jnt_idx(1,2) = find(strcmp(jnt_lnk_names, 'r_sole'), 1);
end
