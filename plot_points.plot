# plot_points.plot
set terminal png size 800,600
set output 'points_plot.png'

# Настройки графика
set title "2D Representation of Points"
set xlabel "X-axis"
set ylabel "Y-axis"

# Диапазоны (можно изменить по необходимости)
set xrange [-7e+11:7e+13]
set yrange [-2.5e+13:2e+14]

# Рисуем точки
plot '-' using 1:2 with points pointtype 7 pointsize 1.5 title 'Points'
4.3426582623434055e+13 1.5889472469555943e+13
6.8048019575759688e+13 -2.1779086455138918e+13
-5.7972749134143286e+11 2.9004098939704912e+12
e