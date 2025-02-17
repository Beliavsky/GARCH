program xgeneral_garch_data
! read asset prices, compute returns and fit various GARCH models to them
use            kind_mod, only: dp
use           garch_mod, only: fit_garch, nparam_garch, nparam_gjr_garch, &
                               print_garch_results, fit_gjr_garch, garch_str, &
                               gjr_garch_str, print_garch_lik, print_garch_lik
use         obj_fun_mod, only: garch_model
use            util_mod, only: write_merge, read_words_line, str
use     basic_stats_mod, only: mean, sd, kurtosis, basic_stats, &
                               correl, acf
use       dataframe_mod, only: DataFrame, print_summary, operator(*), &
                               nrow, ncol
use dataframe_stats_mod, only: simple_ret
use           table_mod, only: Table, display
use     table_stats_mod, only: basic_stats_table, corr_table
use           watch_mod, only: watch, print_elapsed_times
implicit none
character(len=*), parameter :: fmt_ci = "(a25,':',*(1x,i20))", &
   fmt_cr = "(a25,':',*(1x,f20.6))", fmt_cc = "(a25,':',*(1x,a))"
integer :: info, niter_nm, isym, nsym, nret, iu, max_col, nacfsq, &
   nacfabs, max_iter_nm, idist, ndist, nfits
real(kind=dp), allocatable :: sigma_est(:)
character (len=*), parameter :: infile = "general_garch_data.dat"
type(DataFrame) :: df, df_ret
type(Table) :: ret_stats
real(kind=dp) :: alpha0, gamma0, beta0 ! initial guesses for GARCH parameters
real(kind=dp) :: lik_tol ! convergence criterion for GARCH likelihood
real(kind=dp) :: scale_ret
real(kind=dp) :: logL ! log-likelihood
real(kind=dp) :: garch_par(nparam_garch), gjr_garch_par(nparam_gjr_garch)  ! estimated parameters
character (len=20), allocatable :: dist(:), garch_models(:), symbols(:)
real(kind=dp), allocatable :: ret(:), logl_garch(:,:), logl_gjr_garch(:,:)
character (len=100) :: prices_file
open (newunit=iu, file=infile, action="read", status="old")
print fmt_cc, "infile", infile
read (iu,*) prices_file
read (iu,*) max_col ! max # of columns of prices to read
read (iu,*) scale_ret
read (iu,*) nacfsq  ! # of autocorrelations of squared returns to print
read (iu,*) nacfabs ! # of autocorrelations of absolute returns to print
call read_words_line(iu, dist) ! conditional distributions
read (iu,*) max_iter_nm
read (iu,*) lik_tol ! convergence criterion for GARCH likelihood
read (iu,*) alpha0  ! guess for weight on past squared return
read (iu,*) gamma0   ! guess for weight on past squared return when return is negative
read (iu,*) beta0   ! guess for weight on previous variance
call read_words_line(iu, garch_models)  !list of GARCH models to fit
print fmt_cc, "prices file", trim(prices_file)
print fmt_cr, "return scaling", scale_ret
print fmt_ci, "max_iter_nm", max_iter_nm
print fmt_cr, "log10(lik_tol)", log10(lik_tol)
print fmt_cr, "alpha0", alpha0
print fmt_cr, "beta0", beta0
call watch("init")
call df%read_csv(prices_file, max_col=max_col)
nsym = ncol(df)
print fmt_ci, "#assets", nsym
call print_summary(df, fmt_trailer="()")
df_ret = scale_ret*simple_ret(df)
nret = nrow(df_ret)
allocate (sigma_est(nret))
ret_stats = basic_stats_table(df_ret%values, df_ret%columns, &
   index_name="returns")
symbols = df_ret%columns
call display(ret_stats, fmt_trailer="()")
if (nsym > 1) call display(corr_table(df_ret%values,df_ret%columns), &
   title="return correlations")
call watch("read prices and computed return stats")
ndist = size(dist)
nfits = 0
allocate (logl_garch(nsym, ndist), logl_gjr_garch(nsym, ndist))
do isym=1,nsym
do idist=1,ndist
   write (*, "(/,a25,':',*(1x,a))") "symbol", trim(df_ret%columns(isym))
   ret = df_ret%values(:,isym)
   if (any(garch_models == garch_str)) then
      call fit_garch(ret, dist(idist), garch_par, logL, info, niter=niter_nm, &
         sigma=sigma_est, max_iter=max_iter_nm, tol=lik_tol, alpha0=alpha0, &
         beta0=beta0)
      nfits = nfits + 1
      logl_garch(isym, idist) = logl
      call print_garch_results(garch_model, dist(idist), garch_par, &
         logl, niter=niter_nm, info=info, sigma_est=sigma_est, ret=ret, &
         nacfsq=nacfsq, nacfabs=nacfabs, fmt_header="()")
   end if
   if (any(garch_models == gjr_garch_str)) then
      call fit_gjr_garch(ret, dist(idist), gjr_garch_par, logL, info, niter=niter_nm, &
         sigma=sigma_est, max_iter=max_iter_nm, tol=lik_tol, alpha0=alpha0, &
         gamma0=gamma0, beta0=beta0)
      nfits = nfits + 1
      logl_gjr_garch(isym, idist) = logl
      call print_garch_results(garch_model, dist(idist), gjr_garch_par, &
         logl, niter=niter_nm, info=info, sigma_est=sigma_est, ret=ret, &
         nacfsq=nacfsq, nacfabs=nacfabs, fmt_header="()")
   end if
end do
end do
call watch("fit " // trim(str(nfits)) // " GARCH models")
if (any(garch_models == garch_str)) &
   call print_garch_lik(garch_str, symbols, dist, logl_garch)
if (any(garch_models == gjr_garch_str)) &
   call print_garch_lik(gjr_garch_str, symbols, dist, logl_gjr_garch)
call print_elapsed_times()
print "(/,a)", "(8) finished xgeneral_garch_data"
end program xgeneral_garch_data
