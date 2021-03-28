setlocal enabledelayedexpansion
set func1=c_p
set func2=der_phi_Gamma
set func3=der_phi_gamma
set func4=der_phi_plume
set rows=50
set cols=100
set steps=100
set speed=0.1
set /A count=1
for %%a in (%func1% %func2% %func3% %func4%) do (
set fname_tmp=%%a_%rows%_%cols%_%steps%
set fname_out=!fname_tmp!_speed%speed%
echo !count!
set /A count=!count!+1
echo !fname_tmp!
echo !fname_out!
)
pause