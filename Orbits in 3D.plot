
# Файл plot3d.plot

# Установка параметров графика

set term wxt size 1200,900 # Установка размера окна графика
set xlabel "X" # Название оси X
set ylabel "Y" # Название оси Y
set zlabel "Z" # Название оси Z
set title "Orbits" # Заголовок графика

# Отображение 3D-графика
splot 'S2.txt' using 2:3:4 with lines title "S2", \
			'S38.txt' using 2:3:4 with lines title "S38", \
			'S55.txt' using 2:3:4 with lines title "S55"


# Пауза, чтобы пользователь мог просмотреть график
pause -1