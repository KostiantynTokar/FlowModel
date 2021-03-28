set func=c_p
set "title=$c_p$"
set rows=50
set cols=100
set steps=5
set speed=0.1
set path_tmp=output\tmp\r%rows%_c%cols%_s%steps%
set fname_tmp=%func%1
set full_tmp=%path_tmp%\%fname_tmp%.txt
python ScalarFieldAnimation\ScalarFieldAnimation.py -i %full_tmp% -s %speed% -t %title% --dont_save
pause