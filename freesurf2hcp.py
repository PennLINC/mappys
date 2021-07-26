import os
import sys

lh_name = '{0}.32k_fs_LR.lh.shape.gii'.format(sys.argv[3])
rh_name = '{0}.32k_fs_LR.rh.shape.gii'.format(sys.argv[3])

lh_cmd = 'wb_command -metric-resample {0} resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii \
    resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii ADAP_BARY_AREA {1} \
    -area-metrics resample_fsaverage/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii \
    resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii'.format(sys.argv[1],lh_name)

os.system(lh_cmd)

rh_cmd = 'wb_command -metric-resample {0} resample_fsaverage/fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii \
    resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii ADAP_BARY_AREA {1} \
    -area-metrics resample_fsaverage/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii \
    resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii'.format(sys.argv[2],rh_name)

os.system(rh_cmd)