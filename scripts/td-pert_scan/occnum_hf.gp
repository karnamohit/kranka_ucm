set terminal postscript enhanced color "Helvetica-Bold" 12;
set output 'time_occup.dens_nofield_h_gs_20fs_y_rhf.eps'; 
set title "Real time TDHF/STO-3G: H_{3}C-CH=CH-N(CH_{3})_{2}\n{/=16 time vs MO occupations}" font "Helvetica Bold, 20"; 
set xlabel "t (fs)" font "Helvetica Bold, 14"; 
# set xrange [0.0:2.4];
set ylabel 'occupation number' font "Helvetica Bold, 14"; 
set yrange [0.0:2.0]; 
set key right top Right;
set xtics font "Helvetica, 18"; 
# set bmargin 5.0; 
set ytics font "Helvetica, 18"; 
# set mytics 0.0005; 
plot './time_occup.dens_nofield_h_gs_20fs_y_rhf.txt' u 1:($16*1) w l lw 4 title '(HOMO-4)',\
     './time_occup.dens_nofield_h_gs_20fs_y_rhf.txt' u 1:($17*1) w l lw 4 title '(HOMO-3)',\
     './time_occup.dens_nofield_h_gs_20fs_y_rhf.txt' u 1:($18*1) w l lw 4 title '(HOMO-2)',\
     './time_occup.dens_nofield_h_gs_20fs_y_rhf.txt' u 1:($19*1) w l lw 4 title '(HOMO-1)',\
     './time_occup.dens_nofield_h_gs_20fs_y_rhf.txt' u 1:($20*1) w l linecolor rgb "black" lw 4 title 'HOMO',\
     './time_occup.dens_nofield_h_gs_20fs_y_rhf.txt' u 1:($21*1) w l linecolor rgb "red" lw 4 title 'LUMO',\
     './time_occup.dens_nofield_h_gs_20fs_y_rhf.txt' u 1:($22*1) w l lw 4 title '(LUMO+1)',\
     './time_occup.dens_nofield_h_gs_20fs_y_rhf.txt' u 1:($23*1) w l lw 4 title '(LUMO+2)',\
     './time_occup.dens_nofield_h_gs_20fs_y_rhf.txt' u 1:($24*1) w l lw 4 title '(LUMO+3)',\
     './time_occup.dens_nofield_h_gs_20fs_y_rhf.txt' u 1:($25*1) w l lw 4 title '(LUMO+4)'
