
rm -rf particle_info
mkdir particle_info

rm -rf particle_vel_info
mkdir particle_vel_info

rm -rf images
mkdir images

split -d -a 4 -l 7290 velocity.txt ./particle_vel_info/part_vel_data_

split -d -a 4 -l 7290 data.txt ./particle_info/part_data_

echo set term png font "'Times New Roman, 20'" size 2000,2000 fontscale 3.0 > gnu_plotting.gp
echo unset xtics >> gnu_plotting.gp
echo unset ytics >> gnu_plotting.gp
echo unset xlabel >> gnu_plotting.gp
echo unset ylabel >> gnu_plotting.gp
echo set xrange [0.15:0.65] >> gnu_plotting.gp
echo set yrange [0.15:0.65] >> gnu_plotting.gp


i=0
for file in `ls ./particle_info/`
do
echo
echo set out "'"./images/$file.png"'";
time=$(echo "scale=4; 0.05*$i" | bc)
echo set title "'"t = $time s"'"

echo plot "'"./particle_info/$file"'" u 2:3:4 every ::0::3000 w circle lc rgb "'"#800000"'" fill solid noborder title "''","'"./particle_info/$file"'" u 2:3:4 every ::3001::6119 w circle lc rgb "'"#F1C40F"'" fill solid noborder title "''","'"./particle_info/$file"'" u 2:3:4 every ::6120 w circle lc rgb "'"#34495E"'" fill solid noborder title "''"

#echo plot "'"./particle_info/$file"'" u 2:3:4 every ::0::30000 w circle lc rgb "'"#800000"'" fill solid noborder title "''","'"./particle_info/$file"'" u 2:3:4 every ::30001::60000 w circle lc rgb "'"#008080"'" fill solid noborder title "''","'"./particle_info/$file"'" u 2:3:4 every ::60001::80503 w circle lc rgb "'"#F1C40F"'" fill solid noborder title "''", "'"./particle_info/$file"'" u 2:3:4 every ::80504 w circle lc rgb "'"#34495E"'" fill solid noborder title "''"

echo unset out
echo
i=$(($i+1))
done >> gnu_plotting.gp


gnuplot gnu_plotting.gp


avconv -r 20 -i ./images/part_data_%4d.png video.avi

#ffmpeg -i video_particles.avi -vf "movie=lol.avi [in1]; [in]pad=2000*2:2000[in0]; [in0][in1] overlay=2000:0 [out]" out.avi

