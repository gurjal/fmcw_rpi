set grid
plot [1:2048][0:6]\
  'data.dat' using 1:2 with lines,\
  'data.dat' using 1:3 with lines
pause 0.5
reread
