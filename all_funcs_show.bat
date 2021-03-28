setlocal enabledelayedexpansion
set func1=c_p
set func2=der_phi_Gamma
set func3=der_phi_gamma
set func4=der_phi_plume
set rows=50
set cols=100
set steps=300
set speed=0.1
set path_tmp=output\tmp\r%rows%_c%cols%_s%steps%
set /A count=1
for %%a in (%func1%, %func2%, %func3%, %func4%) do (
set fname_tmp=%%a!count!
set /A count=!count!+1
set full_tmp=%path_tmp%\!fname_tmp!.txt
python ScalarFieldAnimation\ScalarFieldAnimation.py -i !full_tmp! -s %speed% -t %%a --dont_save
)
pause