@echo OFF
FOR /L %%i IN (0,1,7) DO (
    (START octave -qf driveCN.m %%i)
)