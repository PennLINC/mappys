import nibabel as nib
import os
import numpy as np
import pandas as pd

# n_nodes = sys.argv[1]
n_nodes = 400

parcel_path = 'https://github.com/ThomasYeoLab/CBIG/raw/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/HCP/fslr32k/cifti/Schaefer2018_{0}Parcels_17Networks_order.dlabel.nii'.format(n_nodes)
parcel_file = 'Schaefer2018_{0}Parcels_17Networks_order.dlabel.nii'.format(n_nodes)
os.system('rm {0}'.format(parcel_file))
os.system('wget {0}'.format(parcel_path))
parcels= nib.load(parcel_file)
p = parcels.get_fdata()
df = pd.DataFrame(columns=['Schaefer{0}_17'.format(n_nodes)])
df['Schaefer{0}_17'.format(n_nodes)] = np.arange(n_nodes).astype(int)

#princple gradients
sa_file = 'hcp.gradients.dscalar.nii'
if os.path.exists(sa_file) == False:
    sa_path = 'https://github.com/PennLINC/Brain_Organization/raw/master/PrincipleGradient/hcp.gradients.dscalar.nii'
    os.system('wget {0}'.format(sa_path))
out = "sa_{0}.pscalar.nii".format(n_nodes)
cmd = 'wb_command -cifti-parcellate {0} {1} 2 {2}'.format(sa_file,parcel_file,out)
os.system(cmd)
df['SA_1'] = nib.load(out).get_fdata()[0]

    
# mylein
mfile = 'ave_MyelinMap_BC_MSMAll.32k_fs_LR.dscalar.nii'
if os.path.exists(mfile) == False:
    os.system('wget https://github.com/PennLINC/Brain_Organization/raw/master/Myelin/ave_MyelinMap_BC_MSMAll.32k_fs_LR.dscalar.nii')
out = "myelin_{0}.pscalar.nii".format(n_nodes)
cmd = 'wb_command -cifti-parcellate {0} {1} 2 {2}'.format(mfile,parcel_file,out)
os.system(cmd)
df['myelin'] = nib.load(out).get_fdata()[0]


#thickness
tfile = 'Q1-Q6_RelatedValidation210.corrThickness_MSMAll_2_d41_WRN_DeDrift_grad_s0.32k_fs_LR.dscalar.nii'
if os.path.exists(tfile) == False:
    os.system('wget https://github.com/PennLINC/Brain_Organization/raw/master/CorrThickness/Q1-Q6_RelatedValidation210.corrThickness_MSMAll_2_d41_WRN_DeDrift_grad_s0.32k_fs_LR.dscalar.nii')
out = "corrthick_{0}.pscalar.nii".format(n_nodes)
cmd = 'wb_command -cifti-parcellate {0} {1} 2 {2}'.format(tfile,parcel_file,out)
os.system(cmd)
df['CorrThickness'] = nib.load(out).get_fdata()[0]

# evo
lfile = 'lh.macaque-human.RelativeAreaExpansion.human.32k_fs_LR.nii.gz'
rfile = 'rh.macaque-human.RelativeAreaExpansion.human.32k_fs_LR.nii.gz'
if os.path.exists(lfile) == False:
    os.system('wget https://github.com/TingsterX/alignment_macaque-human/raw/main/area_expansion/lh.macaque-human.RelativeAreaExpansion.human.32k_fs_LR.nii.gz')
if os.path.exists(rfile) == False:
    os.system('wget https://github.com/TingsterX/alignment_macaque-human/raw/main/area_expansion/rh.macaque-human.RelativeAreaExpansion.human.32k_fs_LR.nii.gz')
df['evo_expansion'] = np.zeros((400))
evo = np.zeros((p.shape[1]))
evo[:32492] = nib.load(lfile).get_fdata().flatten()
evo[32492:] = nib.load(rfile).get_fdata().flatten()

for i in range(n_nodes):
    df['evo_expansion'][i] = evo[p[0]==i+1].mean()

out = "evo_expansion_{0}.pscalar.nii".format(n_nodes)
nib.save(nib.Cifti2Image(df['evo_expansion'].values.reshape(1,n_nodes),header=nib.load("corrthick_{0}.pscalar.nii".format(n_nodes)).header),out)


# allo
os.system('wget https://github.com/PennLINC/Brain_Organization/raw/master/AllometricScaling/rh.AllometricScaling_fsaverage5.func.gii')
os.system('wget https://github.com/PennLINC/Brain_Organization/raw/master/AllometricScaling/lh.AllometricScaling_fsaverage5.func.gii')
cmd = 'python freesurf2hcp.py lh.AllometricScaling_fsaverage5.func.gii rh.AllometricScaling_fsaverage5.func.gii alloscale'
os.system(cmd)


lfile = 'alloscale.32k_fs_LR.lh.shape.gii'
rfile = 'alloscale.32k_fs_LR.rh.shape.gii'
df['alloscale'] = np.zeros((400))
evo = np.zeros((p.shape[1]))
evo[:32492] = nib.load(lfile).agg_data().flatten()
evo[32492:] = nib.load(rfile).agg_data().flatten()

for i in range(n_nodes):
    df['alloscale'][i] = evo[p[0]==i+1].mean()

out = "alloscale_{0}.pscalar.nii".format(n_nodes)
nib.save(nib.Cifti2Image(df['alloscale'].values.reshape(1,n_nodes),header=nib.load("corrthick_{0}.pscalar.nii".format(n_nodes)).header),out)

# cbf
os.system('wget https://github.com/PennLINC/Brain_Organization/raw/master/MeanCBF/lh.MeanCBF.fsaverage5.func.gii')
os.system('wget https://github.com/PennLINC/Brain_Organization/raw/master/MeanCBF/rh.MeanCBF.fsaverage5.func.gii')

lfile = 'lh.MeanCBF.fsaverage5.func.gii'
rfile = 'rh.MeanCBF.fsaverage5.func.gii'

cmd = 'python freesurf2hcp.py lh.MeanCBF.fsaverage5.func.gii rh.MeanCBF.fsaverage5.func.gii cbf'
os.system(cmd)

lfile = 'cbf.32k_fs_LR.lh.shape.gii'
rfile = 'cbf.32k_fs_LR.rh.shape.gii'

df['cbf'] = np.zeros((400))
evo = np.zeros((p.shape[1]))
evo[:32492] = nib.load(lfile).agg_data().flatten()
evo[32492:] = nib.load(rfile).agg_data().flatten()

for i in range(n_nodes):
    df['cbf'][i] = evo[p[0]==i+1].mean()

out = "cbf_{0}.pscalar.nii".format(n_nodes)
nib.save(nib.Cifti2Image(df['cbf'].values.reshape(1,n_nodes),header=nib.load("corrthick_{0}.pscalar.nii".format(n_nodes)).header),out)



df.to_csv('mappings_{0}.csv'.format(n_nodes))
