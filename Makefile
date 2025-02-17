executables = xew_arch_data_gfort.exe xew_arch_sim_gfort.exe xgarch_sim_gfort.exe xgarch_data_gfort.exe xgeneral_garch_data_gfort.exe xgjr_garch_loop_gfort.exe xgjr_garch_data_gfort.exe xgjr_garch_gfort.exe
FC     = gfortran
FFLAGS = -O0 -Wall -Werror=unused-parameter -Werror=unused-variable -Werror=unused-function -Wno-maybe-uninitialized -Wno-surprising -fbounds-check -static -g -fmodule-private
obj    = kind.o util.o watch.o constants.o density.o random.o basic_stats.o table.o table_stats.o dataframe.o dataframe_stats.o obj_fun_garch.o nelder_mead.o garch.o xew_arch_data.o xew_arch_sim.o xgarch_sim.o xgarch_data.o xgeneral_garch_data.o gjr_garch.o xgjr_garch_loop.o xgjr_garch_data.o xgjr_garch.o

all: $(executables)

# Compile .f90 to .o
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

xew_arch_data_gfort.exe: kind.o util.o watch.o constants.o density.o random.o basic_stats.o table.o table_stats.o dataframe.o dataframe_stats.o obj_fun_garch.o nelder_mead.o garch.o xew_arch_data.o
	$(FC) -o xew_arch_data_gfort.exe kind.o util.o watch.o constants.o density.o random.o basic_stats.o table.o table_stats.o dataframe.o dataframe_stats.o obj_fun_garch.o nelder_mead.o garch.o xew_arch_data.o $(FFLAGS)

xew_arch_sim_gfort.exe: kind.o util.o constants.o density.o random.o basic_stats.o obj_fun_garch.o nelder_mead.o garch.o xew_arch_sim.o
	$(FC) -o xew_arch_sim_gfort.exe kind.o util.o constants.o density.o random.o basic_stats.o obj_fun_garch.o nelder_mead.o garch.o xew_arch_sim.o $(FFLAGS)

xgarch_sim_gfort.exe: kind.o util.o constants.o density.o random.o basic_stats.o obj_fun_garch.o nelder_mead.o garch.o xgarch_sim.o
	$(FC) -o xgarch_sim_gfort.exe kind.o util.o constants.o density.o random.o basic_stats.o obj_fun_garch.o nelder_mead.o garch.o xgarch_sim.o $(FFLAGS)

xgarch_data_gfort.exe: kind.o util.o watch.o constants.o density.o random.o basic_stats.o table.o table_stats.o dataframe.o dataframe_stats.o obj_fun_garch.o nelder_mead.o garch.o xgarch_data.o
	$(FC) -o xgarch_data_gfort.exe kind.o util.o watch.o constants.o density.o random.o basic_stats.o table.o table_stats.o dataframe.o dataframe_stats.o obj_fun_garch.o nelder_mead.o garch.o xgarch_data.o $(FFLAGS)

xgeneral_garch_data_gfort.exe: kind.o util.o watch.o constants.o density.o random.o basic_stats.o table.o table_stats.o dataframe.o dataframe_stats.o obj_fun_garch.o nelder_mead.o garch.o xgeneral_garch_data.o
	$(FC) -o xgeneral_garch_data_gfort.exe kind.o util.o watch.o constants.o density.o random.o basic_stats.o table.o table_stats.o dataframe.o dataframe_stats.o obj_fun_garch.o nelder_mead.o garch.o xgeneral_garch_data.o $(FFLAGS)

xgjr_garch_loop_gfort.exe: kind.o util.o constants.o density.o random.o basic_stats.o obj_fun_garch.o nelder_mead.o gjr_garch.o xgjr_garch_loop.o
	$(FC) -o xgjr_garch_loop_gfort.exe kind.o util.o constants.o density.o random.o basic_stats.o obj_fun_garch.o nelder_mead.o gjr_garch.o xgjr_garch_loop.o $(FFLAGS)

xgjr_garch_data_gfort.exe: kind.o util.o watch.o constants.o density.o random.o basic_stats.o table.o table_stats.o dataframe.o dataframe_stats.o obj_fun_garch.o nelder_mead.o gjr_garch.o xgjr_garch_data.o
	$(FC) -o xgjr_garch_data_gfort.exe kind.o util.o watch.o constants.o density.o random.o basic_stats.o table.o table_stats.o dataframe.o dataframe_stats.o obj_fun_garch.o nelder_mead.o gjr_garch.o xgjr_garch_data.o $(FFLAGS)

xgjr_garch_gfort.exe: kind.o util.o constants.o density.o random.o basic_stats.o obj_fun_garch.o nelder_mead.o gjr_garch.o xgjr_garch.o
	$(FC) -o xgjr_garch_gfort.exe kind.o util.o constants.o density.o random.o basic_stats.o obj_fun_garch.o nelder_mead.o gjr_garch.o xgjr_garch.o $(FFLAGS)

run: $(executables)
	./xew_arch_data_gfort.exe
	./xew_arch_sim_gfort.exe
	./xgarch_sim_gfort.exe
	./xgarch_data_gfort.exe
	./xgeneral_garch_data_gfort.exe
	./xgjr_garch_loop_gfort.exe
	./xgjr_garch_data_gfort.exe
	./xgjr_garch_gfort.exe

clean:
	rm -f $(executables) $(obj)

