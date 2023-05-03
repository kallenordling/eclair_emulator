# Python code for filtering ECHAM netCDF outputs.
#
# Required Puhti modules:
#	python-data/3.10-22.09
#
import netCDF4 as netcdf
import numpy as np

# Output variables - to be saved to a text file
var_list=["h2o_inv","tpot_inv","emu_lwp","tpot_pbl","pbl_h","emu_cdnc","cos_mu","emu_w","emu_precip","autoconv_int","accretion1_int","accretion2_int","W"]

def filter_extract_data(emu_fname,act_fname,mask_fname):
    # Open the target NetCDF files
    # a) 3D variables from the emulator outputs
    ncid_emu = netcdf.Dataset(emu_fname,"r")
    vars_emu = ncid_emu.variables
    if "emu_mask" not in vars_emu:
        print("Mask not found from "+emu_fname+"!")
        ncid_emu.close()
        return np.zeros((1,))
    # b) 4D variables from active outputs (optional)
    if len(act_fname)>1:
        ncid_act = netcdf.Dataset(act_fname,"r")
        vars_act = ncid_act.variables
    else:
        vars_act=[""]
    # c) 4D  emulator mask (emu_mask_3d) either from a separate netCDF file or from the other files
    mask3d_flag = True
    if len(mask_fname)>1:
        ncid_mask = netcdf.Dataset(mask_fname,"r")
        if "emu_mask_3d" not in ncid_mask.variables:
            print("3D mask not found!")
            ncid_emu.close()
            if len(act_fname)>1: ncid_act.close()
            ncid_mask.close()
            return np.zeros((1,))
    elif "emu_mask_3d" in vars_emu:
            ncid_mask = ncid_emu
    elif "emu_mask_3d" in vars_act:
            ncid_mask = ncid_act
    else:
        mask3d_flag = False
    #
    # 3D mask (time,x,y)
    # Vector indexes for selecting the data
    i,j,k, = np.where( np.array(ncid_emu.variables['emu_mask'])>0.5 )
    #
    # 4D mask
    if mask3d_flag:
        # The mask wave (from float to bool)
        mask_4d = np.array( np.array(ncid_mask.variables['emu_mask_3d']) > 0.5 )
        # Count the number of valid entries (dimension 1 = z) and then
        # reduce to vector based on indexies from the 3D mask
        counts = np.sum(mask_4d, axis=1, dtype=float)[i,j,k]
    #
    out = np.zeros( (i.shape[0],len(var_list)) )
    col = 0
    for var in var_list:
        if var in vars_emu:
            # Dimensions
            n=len(ncid_emu.variables[var].shape)
            if n==3:
                # 3D variables: directly to vector
                out[:,col]=np.array(ncid_emu.variables[var])[i,j,k]
            elif n==4 and mask3d_flag:
                # 4D variables: average over z
                # First sum over z and then to vector; this is divided by the number of selected points
                out[:,col]=np.sum( np.array(ncid_emu.variables[var]), axis=1, where=mask_4d )[i,j,k]/counts
            else:
                print("Variable "+var+" dimensions not valid!")
        elif var in vars_act:
            # Dimensions
            n=len(ncid_act.variables[var].shape)
            if n==3:
                # 3D variables: directly to vector
                out[:,col]=np.array(ncid_act.variables[var])[i,j,k]
            elif n==4 and mask3d_flag:
                # 4D variables: average over z
                # First sum over z and then to vector; this is divided by the number of selected points
                out[:,col]=np.sum( np.array(ncid_act.variables[var]), axis=1, where=mask_4d )[i,j,k]/counts
            else:
                print("Variable "+var+" dimensions not valid!")
        else:
            print("Variable "+var+" not found!")
        col+=1
    #
    # Close files
    ncid_emu.close()
    if len(act_fname)>1: ncid_act.close()
    if len(mask_fname)>1: ncid_mask.close()
    #
    return out

# Collect data from the original ECHAM outputs: year 2000, months 1-12
emu_fmt = "/fmi/scratch/project_2001927/nordlin1/emu-kontrolli-r6284/emu-kontrolli-r6284_2000%02u.01_emulato.nc"
act_fmt = "/fmi/scratch/project_2001927/nordlin1/emu-kontrolli-r6284/W_emu-kontrolli-r6284_2000%02u.nc"

# New file: header
f=open("echam_filtered.dat", "w")
for var in var_list: f.write(var+" ")
f.write("\n")
f.close()
# Data - append mode
f=open("echam_filtered.dat", "a")
for mo in range(1,13):
    print(("Prosessing month %u")%(mo))
    emu_fname = (emu_fmt) % (mo)
    act_fname = (act_fmt) % (mo)
    out=filter_extract_data(emu_fname,act_fname,"")
    if len(out.shape)==2:
        np.savetxt(f,out,fmt='%.6e')
    else:
        print("Skipped!")
f.close()
