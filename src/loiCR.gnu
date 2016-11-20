set terminal x11
set encoding utf8
set title "Gaussienne"
set xlabel "Points xi "
set ylabel "Images des Fi ou fi"
plot "loinormaleCR.dat" using 1:2 title "Fonction de repartition" with lines
replot "loinormaleCR.dat" using 1:3 title "Densite de probabilite loi normale centree reduite" with lines
replot "loinormaleMuSigma.dat" using 1:3 title "Densite de probabilite mu = -2 sigma**2 = 0.5" with lines
pause -1 
quit

