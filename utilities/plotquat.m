function plotquat(vqT, twait)
    if ~exist('twait', 'var')
        % use the default time duration of waiting (in seconds) ...
        twait = 0.0005;
    end

    if ~ismatrix(vqT)
        error('plotquat: %s', WBM.wbmErrorMsg.WRONG_DATA_TYPE);
    end
    % the matrix-dimension must be of m-by-7:
    [nRows, nCols] = size(vqT);
    if (nCols ~= 7)
        error('plotquat: %s', WBM.wbmErrorMsg.WRONG_MAT_DIM);
    end
    lin_width = 1;
    mkr_size  = 8.5;

    figure('Name', 'Quaternions:', 'NumberTitle', 'off');
    cidx = nRows - 1; % color-index ...
    cmap = colormap( copper(cidx) );
    %cmap = colormap( summer(cidx) ); % optional

    %% Translations:
    % get the position values of each element:
    pos_x = vqT(1:nRows,1);
    pos_y = vqT(1:nRows,2);
    pos_z = vqT(1:nRows,3);

    %% Orientations:
    % convert the quaternions into axis-angle rotations:
    ang_t = acos(vqT(1:nRows,4))*2; % compute the axis-angle theta ...
    % compute the rotation axis elements ...
    rax_x = vqT(1:nRows,5)./sqrt(1 - vqT(1:nRows,4).*vqT(1:nRows,4));
    rax_y = vqT(1:nRows,6)./sqrt(1 - vqT(1:nRows,4).*vqT(1:nRows,4));
    rax_z = vqT(1:nRows,7)./sqrt(1 - vqT(1:nRows,4).*vqT(1:nRows,4));
    % convert the axis-angle rotation of the 1st quaternion element
    % into a rotation matrix representation:
    R = WBM.utilities.axang2rotm( vertcat(rax_x(1,1), rax_y(1,1), rax_z(1,1), ang_t(1,1)) );

    %% Plot the orientation axes (xyz-axes) of all quaternions:

    axis( repmat([-2 2], 1, 3) );

    hold on;

    xlabel('pos_x');
    ylabel('pos_y');
    zlabel('pos_z');
    legend;

    grid on;

    strTitle = 'Quaternion (#%d)';
    title(sprintf(strTitle, 1));

    % 3D-position of the first quaternion ...
    hPos3D = plot3(pos_x(1,1), pos_y(1,1), pos_z(1,1), 'Marker', 'o', 'MarkerEdgeColor', 'black');
    if verLessThan('matlab', '8.4.0')
        % Matlab <= R2014a:
        set(hPos3D, 'EraseMode', 'normal');
    end

    if ~verLessThan('matlab', '8.4.0')
        % for Matlab R2014b and later:
        set(gca, 'SortMethod', 'childorder');

        % draw the orientation axes (lines) of the 1st quaternion ...
        hLin_orX = animatedline(horzcat(pos_x(1,1), pos_x(1,1)+R(1,1)), horzcat(pos_y(1,1), pos_y(1,1)+R(2,1)), ...
                                horzcat(pos_z(1,1), pos_z(1,1)+R(3,1)), 'LineWidth', lin_width, 'Color', 'red');
        hLin_orY = animatedline(horzcat(pos_x(1,1), pos_x(1,1)+R(1,2)), horzcat(pos_y(1,1), pos_y(1,1)+R(2,2)), ...
                                horzcat(pos_z(1,1), pos_z(1,1)+R(3,2)), 'LineWidth', lin_width, 'Color', 'green');
        hLin_orZ = animatedline(horzcat(pos_x(1,1), pos_x(1,1)+R(1,3)), horzcat(pos_y(1,1), pos_y(1,1)+R(2,3)), ...
                                horzcat(pos_z(1,1), pos_z(1,1)+R(3,3)), 'LineWidth', lin_width, 'Color', 'blue');

        pause; % stop the execution and wait for the user to press any key ...

        % draw the orientation axes (xyz-axes) of all remaining quaternions in the array list ...
        for i = 2:nRows
            title(sprintf(strTitle, i));

            % draw the position of the previous quaternion ...
            prev = i-1;
            plot3(pos_x(prev,1), pos_y(prev,1), pos_z(prev,1), 'Marker', '.', 'MarkerSize', mkr_size, 'MarkerEdgeColor', cmap(cidx,1:3));
            cidx = cidx - 1;
            % set the 3D-position of the new quaternion ...
            set(hPos3D, 'XData', pos_x(i,1), 'YData', pos_y(i,1), 'ZData', pos_z(i,1));
            % calculate the rotation matrix (orientation) from the axis-angle representation of the new quaternion ...
            R = WBM.utilities.axang2rotm( vertcat(rax_x(i,1), rax_y(i,1), rax_z(i,1), ang_t(i,1)) );

            % remove the old orientation axes from the previous quaternion ...
            clearpoints(hLin_orX);
            clearpoints(hLin_orY);
            clearpoints(hLin_orZ);
            % draw the orientation axes (lines) of the new quaterion element ...
            addpoints(hLin_orX, horzcat(pos_x(i,1), pos_x(i,1)+R(1,1)), horzcat(pos_y(i,1), pos_y(i,1)+R(2,1)), horzcat(pos_z(i,1), pos_z(i,1)+R(3,1)));
            addpoints(hLin_orY, horzcat(pos_x(i,1), pos_x(i,1)+R(1,2)), horzcat(pos_y(i,1), pos_y(i,1)+R(2,2)), horzcat(pos_z(i,1), pos_z(i,1)+R(3,2)));
            addpoints(hLin_orZ, horzcat(pos_x(i,1), pos_x(i,1)+R(1,3)), horzcat(pos_y(i,1), pos_y(i,1)+R(2,3)), horzcat(pos_z(i,1), pos_z(i,1)+R(3,3)));

            pause(twait);
        end
    else
        % for older Matlab versions (<= R2014a):
        set(gca, 'DrawMode', 'fast');

        hLin_orX = line('XData', horzcat(pos_x(1,1), pos_x(1,1)+R(1,1)), 'YData', horzcat(pos_y(1,1), pos_y(1,1)+R(2,1)), ...
                        'ZData', horzcat(pos_z(1,1), pos_z(1,1)+R(3,1)), 'EraseMode', 'normal', ...
                        'LineWidth', lin_width, 'Color', 'red');
        hLin_orY = line('XData', horzcat(pos_x(1,1), pos_x(1,1)+R(1,2)), 'YData', horzcat(pos_y(1,1), pos_y(1,1)+R(2,2)), ...
                        'ZData', horzcat(pos_z(1,1), pos_z(1,1)+R(3,2)), 'EraseMode', 'normal', ...
                        'LineWidth', lin_width, 'Color', 'green');
        hLin_orZ = line('XData', horzcat(pos_x(1,1), pos_x(1,1)+R(1,3)), 'YData', horzcat(pos_y(1,1), pos_y(1,1)+R(2,3)), ...
                        'ZData', horzcat(pos_z(1,1), pos_z(1)+R(3,3)), 'EraseMode', 'normal', ...
                        'LineWidth', lin_width, 'Color', 'blue');
        pause;

        for i = 2:nRows
            title(sprintf(strTitle, i));

            prev = i-1;
            plot3(pos_x(prev,1), pos_y(prev,1), pos_z(prev,1), 'Marker', '.', 'MarkerSize', mkr_size, 'MarkerEdgeColor', cmap(cidx,1:3));
            cidx = cidx - 1;

            set(hPos3D, 'XData', pos_x(i,1), 'YData', pos_y(i,1), 'ZData', pos_z(i,1));

            R = WBM.utilities.axang2rotm( vertcat(rax_x(i,1), rax_y(i,1), rax_z(i,1), ang_t(i,1)) );

            set(hLin_orX, 'XData', horzcat(pos_x(i,1), pos_x(i,1)+R(1,1)), 'YData', horzcat(pos_y(i,1), pos_y(i,1)+R(2,1)), ...
                          'ZData', horzcat(pos_z(i,1), pos_z(i,1)+R(3,1)), 'EraseMode', 'normal', ...
                          'LineWidth', lin_width, 'Color', 'blue');
            set(hLin_orY, 'XData', horzcat(pos_x(i,1), pos_x(i,1)+R(1,2)), 'YData', horzcat(pos_y(i,1), pos_y(i,1)+R(2,2)), ...
                          'ZData', horzcat(pos_z(i,1), pos_z(i,1)+R(3,2)), 'EraseMode', 'normal', ...
                          'LineWidth', lin_width, 'Color', 'green');
            set(hLin_orZ, 'XData', horzcat(pos_x(i,1), pos_x(i,1)+R(1,3)), 'YData', horzcat(pos_y(i,1), pos_y(i,1)+R(2,3)), ...
                          'ZData', horzcat(pos_z(i,1), pos_z(i,1)+R(3,3)), 'EraseMode', 'normal', ...
                          'LineWidth', lin_width, 'Color', 'red');
            pause(twait);
        end
    end
end
