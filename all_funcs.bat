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
set path_animation=output\animation\r%rows%_c%cols%_s%steps%_speed%speed%
md %path_tmp%
md %path_animation%
set /A count=1
for %%a in (%func1%, %func2%, %func3%, %func4%) do (
set fname_tmp=%%a!count!
set fname_out=!fname_tmp!
set /A count=!count!+1
set full_tmp=%path_tmp%\!fname_tmp!.txt
set full_out=%path_animation%\!fname_out!.gif
x64\Release\Flow.exe -o !full_tmp! -r %rows% -c %cols% -s %steps% -f %%a
python ScalarFieldAnimation\ScalarFieldAnimation.py -i !full_tmp! -o !full_out! -s %speed% --dont_show -t %%a
)
pause