set xlabel "Processes" 
set ylabel "Time, ms"

plot 	"small1proc" using 1:2 title "small" with lines , "mid1proc" using 1:2 title "middle" with lines, "big1proc" using 1:2 title "big" with lines
