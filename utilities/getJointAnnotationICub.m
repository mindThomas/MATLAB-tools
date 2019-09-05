function jnt_name = getJointAnnotationICub(lnk_name_group, idx)
    switch lnk_name_group
        case 'torso'
            switch idx
                case 1
                    jnt_name = 'torso pitch';
                case 2
                    jnt_name = 'torso roll';
                case 3
                    jnt_name = 'torso yaw';
                otherwise
                    error('getJointAnnotationICub:%s: %s', lnk_name_group, WBM.wbmErrorMsg.IDX_OUT_OF_BOUNDS);
            end
        case 'right_arm'
            switch idx
                case 1
                    jnt_name = 'r. shoulder pitch';
                case 2
                    jnt_name = 'r. shoulder roll';
                case 3
                    jnt_name = 'r. shoulder yaw';
                case 4
                    jnt_name = 'r. elbow';
                case 5
                    jnt_name = 'r. wrist prosup';
                otherwise
                    error('getJointAnnotationICub:%s: %s', lnk_name_group, WBM.wbmErrorMsg.IDX_OUT_OF_BOUNDS);
            end
        case 'left_arm'
            switch idx
                case 1
                    jnt_name = 'l. shoulder pitch';
                case 2
                    jnt_name = 'l. shoulder roll';
                case 3
                    jnt_name = 'l. shoulder yaw';
                case 4
                    jnt_name = 'l. elbow';
                case 5
                    jnt_name = 'l. wrist prosup';
                otherwise
                    error('getJointAnnotationICub:%s: %s', lnk_name_group, WBM.wbmErrorMsg.IDX_OUT_OF_BOUNDS);
            end
        case 'right_leg'
            switch idx
                case 1
                    jnt_name = 'r. hip pitch';
                case 2
                    jnt_name = 'r. hip roll';
                case 3
                    jnt_name = 'r. hip yaw';
                case 4
                    jnt_name = 'r. knee';
                case 5
                    jnt_name = 'r. ankle pitch';
                case 6
                    jnt_name = 'r. ankle roll';
                otherwise
                    error('getJointAnnotationICub:%s: %s', lnk_name_group, WBM.wbmErrorMsg.IDX_OUT_OF_BOUNDS);
            end
        case 'left_leg'
            switch idx
                case 1
                    jnt_name = 'l. hip pitch';
                case 2
                    jnt_name = 'l. hip roll';
                case 3
                    jnt_name = 'l. hip yaw';
                case 4
                    jnt_name = 'l. knee';
                case 5
                    jnt_name = 'l. ankle pitch';
                case 6
                    jnt_name = 'l. ankle roll';
                otherwise
                    error('getJointAnnotationICub:%s: %s', lnk_name_group, WBM.wbmErrorMsg.IDX_OUT_OF_BOUNDS);
            end
        otherwise
            error('getJointAnnotationICub: %s', WBM.wbmErrorMsg.UNKNOWN_LNK_NAME);
    end
end
