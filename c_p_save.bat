set func=c_p
set "title=$c_p$"
set rows=50
set cols=100
set steps=5
set speed=0.1
set path_tmp=output\tmp\r%rows%_c%cols%_s%steps%
set path_animation=output\animation\r%rows%_c%cols%_s%steps%_speed%speed%
set fname_tmp=%func%1
set fname_out=%fname_tmp%
set full_tmp=%path_tmp%\%fname_tmp%.txt
set full_out=%path_animation%\%fname_out%.gif
python ScalarFieldAnimation\ScalarFieldAnimation.py -i %full_tmp% -o %full_out% -s %speed% -t %title% --dont_show
pause