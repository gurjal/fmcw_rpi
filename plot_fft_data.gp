set grid
plot \
  'data.dat' using 1:2 with lines,\
  'data.dat' using 1:3 with lines
pause 0.5
reread
