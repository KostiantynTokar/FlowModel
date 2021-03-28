set func=null
set title=""
set steps=1001
set speed=1
set fps=10
rem set alpha=4.71238898038
set alpha=1.57079632679
set param_tmp=s%steps%_alpha%alpha%
set "obstacle_region={(-0.5,-0.5),(0.5,0.5)}"
set "extent=-2 -2 2 2"
set param_tmp=s%steps%
set param_animation=%param_tmp%_speed%speed%_fps%fps%
set path_tmp=output\tmp\+\%param_tmp%
set path_animation=output\animation\+\%param_animation%
md %path_tmp%
md %path_animation%
set fname_tmp=%func%
set fname_out=%fname_tmp%
set full_tmp=%path_tmp%\%fname_tmp%.txt
set full_out=%path_animation%\%fname_out%.gif
x64\Release\Flow.exe -o %full_tmp% -s %steps% -f %func% --obstacle_region %obstacle_region% --exec_policy par_unseq --alpha %alpha%
python ScalarFieldAnimation\ScalarFieldAnimation.py -i %full_tmp% -o %full_out% -s %speed% -t %title% --fps %fps% --dont_show --no_colorbar --extent %extent%
pause