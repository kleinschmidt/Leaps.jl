cd(joinpath(Pkg.dir("Leaps"), "deps"))
pic = @windows ? "" : "-fPIC"

## normal: 
## run(`gfortran -ffixed-form $pic -O3 -shared leaps.f leapshdr.f -o libleaps.so`)

## extra helpful debugging stuff: 
run(`gfortran -g -fbounds-check -Wall -fbacktrace -finit-real=nan -ffixed-form $pic -O3 -shared leaps.f leapshdr.f -o libleaps.so`)
