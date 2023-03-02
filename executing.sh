
clear

cv=3.3
rm cv.txt
echo "$cv" > cv.txt

g=0.5
rm gamma_input.txt
echo "$g" > gamma_input.txt

matlab -nosplash -nodisplay<nucleation_dispersed_forward.m

#matlab -nosplash -nodisplay<nucleation_pocket_forward.m

#matlab -nosplash -nodisplay<nucleation_stripe_forward.m

rm init_positions.txt
cut -f 2 part_data_2000>x.txt
cut -f 3 part_data_2000>y.txt
paste line.txt x.txt y.txt k.txt>init_positions.txt
rm x.txt y.txt

rm init_velocities.txt
cut -f 2 part_vel_data_2000>vx.txt
cut -f 3 part_vel_data_2000>vy.txt
paste line.txt vx.txt vy.txt>init_velocities.txt
rm vx.txt vy.txt


gcc universal_code_gamma.c -lm -O3
time ./a.out
sh plot_regular_rand.sh

matlab -nosplash -nodisplay<vavg.m

rm -rf particle_info
rm -rf particle_vel_info
rm -rf images
rm data.txt velocity.txt gnu_plotting.gp
rm scatter_forward.eps


exit 0

