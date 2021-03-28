set func=c_p
set func_num=1
set rows=50
set cols=100
set steps=300
set save_steps=50 150 250
set speed=0.1
set path_tmp=output\tmp\r%rows%_c%cols%_s%steps%
set path_animation=output\animation\r%rows%_c%cols%_s%steps%_speed%speed%
set fname_tmp=%func%%func_num%
set fname_out=%fname_tmp%
set full_tmp=%path_tmp%\%fname_tmp%.txt
set full_out=%path_animation%\%fname_out%.gif
python ScalarFieldAnimation\ScalarFieldAnimation.py -i %full_tmp% -o %full_out% -s %speed% -t %func% --dont_save --dont_show --save_steps %save_steps%
pause