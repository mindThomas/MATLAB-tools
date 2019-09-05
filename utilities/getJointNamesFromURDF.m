function [jnt_names, ndof] = getJointNamesFromURDF(urdf_file_name, jnt_type)
    if isempty(urdf_file_name)
        error('getJointNamesFromURDF: %s', WBM.wbmErrorMsg.EMPTY_STRING);
    end

    if ~exist('jnt_type', 'var')
        % set 'revolute' as default type ...
        jnt_type = 'revolute';
    elseif ~find(strcmp(jnt_type, {'revolute', 'prismatic', 'fixed'}), 1)
        error('getJointNamesFromURDF: %s', WBM.wbmErrorMsg.UNKNOWN_JNT_TYPE);
    end

    % Try to convert the URDF-file of the robot model into a Matlab structure for easy access to the data:
    % Note: This function uses an open 3rd-party method of Wouter Falkena from the
    %       Delft University of Technology.
    %       Source: <https://www.mathworks.com/matlabcentral/fileexchange/28518-xml2struct>
    try
        xml_tree = xml2struct(urdf_file_name);
    catch
        strError = sprintf('Failed to read the URDF file \"%s\"!', urdf_file_name);
        error('getJointNamesFromURDF: %s', strError);
    end

    nJnts = size(xml_tree.robot.joint,2);
    jnt_names = cell(nJnts,1);
    ndof = 0;

    % get the joint names of the specified type in the same order as listed in the URDF-file:
    idx = 1;
    for i = 1:nJnts
        jnt = xml_tree.robot.joint{1,i};

        if strcmp(jnt.Attributes.type, jnt_type)
            jnt_names{idx,1} = jnt.Attributes.name;
            idx = idx + 1;
        end
    end
    len = idx - 1;
    if (len > 0)
        % remove the remaining empty elements from the array ...
        jnt_names = jnt_names(1:len,1);
        if ~strcmp(jnt_type, 'fixed')
            ndof = len; % update the DOF ...
        end
        return
    end
    % else ...
    jnt_names = cell(0); % return an empty cell-array ...
end
