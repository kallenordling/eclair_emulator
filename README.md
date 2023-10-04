# eclair_emulator
This code is for ECLAIR updraft and precipitation emulator

mo_emulator.f90 contains fortran code for ECHAM GCM for calculating emulated updraft valocity and precipitation tendencies for stratocumulus cloud regions
traininng_and_testing folder contains fortran code to produce gaussian process emulator based on LES model
gp_train_w_night.f90 gp_train_w_day.f90 gp_train_precip_day.f90 gp_train_precip_night.f90 files are fortran codes for produceing corresponding emulators
