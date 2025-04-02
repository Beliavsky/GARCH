executables = xasymm_garch_sim_gfort.exe xgarch_sim_gfort.exe xgarch_data_gfort.exe xgeneral_garch_data_gfort.exe xgjr_garch_loop_gfort.exe xgjr_garch_data_gfort.exe xgjr_garch_gfort.exe
FC     = gfortran
FFLAGS = -O0 -Wall -Werror=unused-parameter -Werror=unused-variable -Werror=unused-function -Wno-maybe-uninitialized -Wno-surprising -fbounds-check -static -g -fmodule-private
obj    = kind.o constants.o util.o random.o basic_stats.o asymm_garch_sim.o xasymm_garch_sim.o density.o obj_fun_garch.o nelder_mead.o garch.o xgarch_sim.o watch.o table.o table_stats.o dataframe.o dataframe_stats.o xgarch_data.o xgeneral_garch_data.o optim_methods.o uobyqa.o gjr_garch.o xgjr_garch_loop.o qsort.o gradient.o xgjr_garch_data.o xgjr_garch.o

all: $(executables)

# Compile .f90 to .o
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

xasymm_garch_sim_gfort.exe: kind.o constants.o util.o random.o basic_stats.o asymm_garch_sim.o xasymm_garch_sim.o
	$(FC) -o xasymm_garch_sim_gfort.exe kind.o constants.o util.o random.o basic_stats.o asymm_garch_sim.o xasymm_garch_sim.o $(FFLAGS)

xgarch_sim_gfort.exe: kind.o constants.o util.o random.o basic_stats.o density.o obj_fun_garch.o nelder_mead.o garch.o xgarch_sim.o
	$(FC) -o xgarch_sim_gfort.exe kind.o constants.o util.o random.o basic_stats.o density.o obj_fun_garch.o nelder_mead.o garch.o xgarch_sim.o $(FFLAGS)

xgarch_data_gfort.exe: kind.o constants.o util.o random.o basic_stats.o density.o obj_fun_garch.o nelder_mead.o garch.o watch.o table.o table_stats.o dataframe.o dataframe_stats.o xgarch_data.o
	$(FC) -o xgarch_data_gfort.exe kind.o constants.o util.o random.o basic_stats.o density.o obj_fun_garch.o nelder_mead.o garch.o watch.o table.o table_stats.o dataframe.o dataframe_stats.o xgarch_data.o $(FFLAGS)

xgeneral_garch_data_gfort.exe: kind.o constants.o util.o random.o basic_stats.o density.o obj_fun_garch.o nelder_mead.o garch.o watch.o table.o table_stats.o dataframe.o dataframe_stats.o xgeneral_garch_data.o
	$(FC) -o xgeneral_garch_data_gfort.exe kind.o constants.o util.o random.o basic_stats.o density.o obj_fun_garch.o nelder_mead.o garch.o watch.o table.o table_stats.o dataframe.o dataframe_stats.o xgeneral_garch_data.o $(FFLAGS)

xgjr_garch_loop_gfort.exe: kind.o constants.o util.o random.o basic_stats.o density.o obj_fun_garch.o nelder_mead.o optim_methods.o uobyqa.o gjr_garch.o xgjr_garch_loop.o
	$(FC) -o xgjr_garch_loop_gfort.exe kind.o constants.o util.o random.o basic_stats.o density.o obj_fun_garch.o nelder_mead.o optim_methods.o uobyqa.o gjr_garch.o xgjr_garch_loop.o $(FFLAGS)

xgjr_garch_data_gfort.exe: kind.o constants.o util.o random.o basic_stats.o density.o obj_fun_garch.o nelder_mead.o watch.o table.o table_stats.o dataframe.o dataframe_stats.o optim_methods.o uobyqa.o gjr_garch.o qsort.o gradient.o xgjr_garch_data.o
	$(FC) -o xgjr_garch_data_gfort.exe kind.o constants.o util.o random.o basic_stats.o density.o obj_fun_garch.o nelder_mead.o watch.o table.o table_stats.o dataframe.o dataframe_stats.o optim_methods.o uobyqa.o gjr_garch.o qsort.o gradient.o xgjr_garch_data.o $(FFLAGS)

xgjr_garch_gfort.exe: kind.o constants.o util.o random.o basic_stats.o density.o obj_fun_garch.o nelder_mead.o optim_methods.o uobyqa.o gjr_garch.o xgjr_garch.o
	$(FC) -o xgjr_garch_gfort.exe kind.o constants.o util.o random.o basic_stats.o density.o obj_fun_garch.o nelder_mead.o optim_methods.o uobyqa.o gjr_garch.o xgjr_garch.o $(FFLAGS)

run: $(executables)
	./xasymm_garch_sim_gfort.exe
	./xgarch_sim_gfort.exe
	./xgarch_data_gfort.exe
	./xgeneral_garch_data_gfort.exe
	./xgjr_garch_loop_gfort.exe
	./xgjr_garch_data_gfort.exe
	./xgjr_garch_gfort.exe

clean:
	rm -f $(executables) $(obj)

