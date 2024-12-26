set title "Asrometric"
set xlabel "R.A."
set ylabel "Decl"
set key outside
set grid 
#set xrange [-0.35: 0.15]
#set yrange [-0.2: 0.2]
#plot 	"S2_Equatorial.txt" using 2:3 with line title "S2 calc", \
#			"S38_Equatorial.txt" using 2:3 with line title "S38 calc", \
#			"S55_Equatorial.txt" using 2:3 with line title "S55 calc", \
#			"S2_Article.txt" using 2:3 with points title "S2 Article", \
#			"S38_Article.txt" using 2:3 with points title "S2 Article", \
#			"S55_Article.txt" using 2:3 with points title "S55 Article"

plot 	"S2_debug.txt" using 2:3 with line title "S2_debug", \
		"S38_debug.txt" using 2:3 with line title "S38_debug", \
		"S55_debug.txt" using 2:3 with line title "S55_debug", \
			
pause -1