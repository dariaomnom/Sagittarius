
# Файл plot3d.plot

# Установка параметров графика

#set term wxt size 1200,900 # Установка размера окна графика
set xlabel "X" # Название оси X
set ylabel "Y" # Название оси Y
set zlabel "Z" # Название оси Z
set title "Orbits" # Заголовок графика

# Отображение 3D-графика


splot   'S2_debug.txt' using 2:3:4 with lines title "S2_debug", \
		'S38_debug.txt' using 2:3:4 with lines title "S38_debug", \
		'S55_debug.txt' using 2:3:4 with lines title "S55_debug", \
		'-' using 1:2:3 with points pointtype 7 pointsize 1.5 title 'Points'
		0 0 0
		e
pause -1