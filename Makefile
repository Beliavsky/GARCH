exec   = gjr_gfort.exe
obj    = kind.o util.o constants.o random.o basic_stats.o obj_fun_garch.o nelder_mead.o gjr_garch.o xgjr_garch.o
FC     = gfortran
FFLAGS = -O0 -Wall -Werror=unused-parameter -Werror=unused-variable -Werror=unused-function -Wno-maybe-uninitialized -Wno-surprising -fbounds-check -static -g -fmodule-private

all: $(exec)

# Compile .f90 to .o
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

$(exec): $(obj)
	$(FC) -o $(exec) $(obj) $(FFLAGS)

run: $(exec)
	./$(exec)

clean:
	rm -f $(exec) $(obj)

