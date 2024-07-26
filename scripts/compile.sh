
#openmp
# qcc -source -fopenmp -grid=octree -autolink blur3D.c -lm

module load gcc/9.4.0

gcc -Wall -std=c99 -D_XOPEN_SOURCE=700 -O2 _blur3D.c -o blur3Do -L$BASILISK/gl -lglutils -lfb_tiny -lm -fopenmp


