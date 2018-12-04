#Fitting to a line

g(x)=a*x+b
fit g(x) 'long_ave_num.dat' every ::30000 using 1:2 via a,b

and then to plot

plot "long_ave_num.dat" every using 1:2 w lp, f(x),g(x)
