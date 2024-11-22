set title "Asrometric"
set xlabel "R.A."
set ylabel "Decl"
plot 	"S2_Equatorial.txt" using 2:3 with line title "S2", \
			"S38_Equatorial.txt" using 2:3 with line title "S38", \
			"S55_Equatorial.txt" using 2:3 with line title "S55"


			
pause -1