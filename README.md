# LES-based emulators for ECHAM

Code and data used in manuscript:

Nordling, K., Keskinen, J.-P., Romakkaniemi, S., Kokkola, H., Räisänen, P., Lipponen, A., Partanen, A.-I., Ahola, J., Tonttila, J., Alper, M. E., Korhonen, H., and Raatikainen, T.: Technical note: Emulation of a large-eddy simulator for stratocumulus clouds in a general circulation model, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2023-912, 2023.

## Directory gpf_master_git
Fortran implementation of Gaussian Process Emulator (GPE) is based on the GPF Fortran library (https://github.com/ots22/gpf). The original code was modified by implementing new covariance function called LINSQEXP (squared exponential plus a linear term), which is defined in the source directory as cov_linsqexp.f90. Linear covariance function cov_lin.f90 was implemented, but it not used here. The prior for the log-likelihood was also modified for out purpose (see gp.f90).

These codes are used first in offline emulator training where the emulators (or their parameters) are saved for further use. Then these codes are used to compile ECHAM with the GPE implementation. This allow reading in the saved emulators and computing the outputs. The third use case of the code is emulator validation based on the leave-one-out tests.

## Directory train
HUOM: Mikäs tämän hakemiston tarkoitus on? Pitäisikö noiden ./lib ja ./obj hakemistojen olla gpf_master_git:issä? Toki käännetyt koodit, lib:it, *.mod, yms. voi jättää pois, kun ne on konekohtaisia.

## Directory traininng_and_testing
This directory contains Fortran source codes and training inputs that are needed to produce the Gaussian Process Emulators for updraft velocity and precipitation tendencies for both day and night.

Directories data_day.precip, data_day.w, data_night.precip and data_night.w contain different versions of the input data files (emulator inputs based on Binary Space Partitioning (BSP) and corresponding LES outputs) for emulator training.

HUOM: onko data_night.auto4 käytössä?

Files gp_train_w_night.f90, gp_train_w_day.f90, gp_train_precip_day.f90, and gp_train_precip_night.f90 are the Fortran source codes for training the corresponding emulators.

HUOM: gp_train_precip_night.f90 on tyhjä? Onkohan gp_train_precip_night_plotting.f90 oikea? Tosin input “./data_tmp_2” on erilainen. Mitä muuten ovat nuo *.w, *.f90_backup, ja data_tmp_2 filut?

File test_emu.f90 is the Fortran source code for running the leave-one-out emulator validation tests. The results from these tests are in text files lou_day_precip.dat, lou_day_w.dat, lou_night_precip.dat, and lou_night_w.dat.

HUOM: Noita res*.dat filujahan ei käytetä missään (siis voinee poistaa)?

HUOM: Olisikohan emulaattoriden hyvä olla täällä? Siis ehkä, mutta ei ole pakko.

## Code mo_emulator.f90

This is the Fortran source code to be used with general circulation model ECHAM. The code is used to find out if an ECHAM column is suitable for emulation, i.e. contains suitable low cloud. If the column is suitable, the code uses the pre-trained emulators for calculating updraft velocity and precipitation tendencies. The code also takes care of additional ECHAM outputs.


## Code pp_echam.py
Python code used to extract data from the ECHAM outputs for calculating the distributions of rain water production rates and updraft velocities from the default ECHAM and the emulators implemented in ECHAM.

HUOM! Tätä dataa ei ole nyt missään tallessa, mutta on ehkä tarpeettoman suuri githubiin talletettavaksi?

## Citation
Please refer to https://doi.org/10.5281/zenodo.8405068 and the manuscript Nordling at al. (2023).
