program xgarch_data
! read asset prices, compute returns and fit GARCH models to them
use            kind_mod, only: dp
use           garch_mod, only: fit_garch, nparam_garch
use         obj_fun_mod, only: garch_model
use          random_mod, only: random_student_t
use            util_mod, only: write_merge, read_words_line, str
use     basic_stats_mod, only: mean, sd, kurtosis, basic_stats, &
                               basic_stats_names, correl, acf
use       dataframe_mod, only: DataFrame, print_summary, operator(*), &
                               nrow, ncol
use dataframe_stats_mod, only: simple_ret
use           table_mod, only: Table, display
use     table_stats_mod, only: basic_stats_table, corr_table
use           watch_mod, only: watch, print_elapsed_times
implicit none
character(len=*), parameter :: fmt_ci = "(a25,':',*(1x,i20))", &
   fmt_cr = "(a25,':',*(1x,f20.6))", fmt_cc = "(a25,':',*(1x,a))",&
   fmt_par = "(a25,':',*(f10.6))"
integer :: info, niter_nm, isym, nsym, nret, iu, max_col, nacf, &
   max_iter_nm, idist, ndist, nfits
real(kind=dp), allocatable :: sigma_est(:)
character (len=*), parameter :: infile = "garch_data.dat"
type(DataFrame) :: df, df_ret
type(Table) :: ret_stats
real(kind=dp) :: alpha0, beta0 ! initial guesses for GARCH parameters
real(kind=dp) :: lik_tol ! convergence criterion for GARCH likelihood
real(kind=dp) :: scale_ret
real(kind=dp) :: logL ! log-likelihood
real(kind=dp) :: par_out(nparam_garch) ! estimated parameters
character (len=20), allocatable :: dist(:)
! Conditional return distribution
real(kind=dp), allocatable :: ret(:)
character (len=100) :: prices_file
open (newunit=iu, file=infile, action="read", status="old")
print fmt_cc, "infile", infile
read (iu,*) prices_file
read (iu,*) max_col ! max # of columns of prices to read
read (iu,*) scale_ret
read (iu,*) nacf ! # of autocorrelations of squared returns to print
call read_words_line(iu, dist)
! read (iu,*) dist ! conditional distribution
read (iu,*) max_iter_nm
read (iu,*) lik_tol ! convergence criterion for GARCH likelihood
read (iu,*) alpha0  ! guess for weight on past squared return
read (iu,*) beta0   ! guess for weight on previous variance
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
call display(ret_stats, fmt_trailer="()")
if (nsym > 1) call display(corr_table(df_ret%values,df_ret%columns), &
   title="return correlations")
call watch("read prices and computed return stats")
! garch_model = "garch"
ndist = size(dist)
nfits = 0
do isym=1,nsym
do idist=1,ndist
   ! Fit the GARCH model to the returns.
   write (*, "(/,a25,':',*(1x,a))") "symbol", trim(df_ret%columns(isym))
   ret = df_ret%values(:,isym)
   call fit_garch(ret, dist(idist), par_out, logL, info, niter=niter_nm, &
      sigma=sigma_est, max_iter=max_iter_nm, tol=lik_tol, alpha0=alpha0, &
      beta0=beta0)
   nfits = nfits + 1
   write (*, fmt_cc) "GARCH model", garch_model
   write (*, fmt_cc) "Conditional distribution", trim(dist(idist))
   write (*, "(/,a25,':',*(a10))") "parameters", "mu", "omega", "alpha", &
                                   "beta"
   write (*, fmt_par) "estimated", par_out
   write (*, fmt_cr) "Log-likelihood", logL
   write (*, fmt_ci) "#Nelder-Mead iterations", niter_nm
   call write_merge(info==0, "Convergence achieved.", &
        "Maximum iterations reached without full convergence.")
   write (*,"(/,a25,':',*(a10))") "stat", basic_stats_names
   write (*,fmt_par) "sigma_est", basic_stats(sigma_est)
   write (*, "(/, a, *(1x,f10.4))") "kurtosis of ret, ret/sigma_est", &
      kurtosis(ret), kurtosis(ret/sigma_est)
   if (nacf > 0) then
      write (*,fmt_par) "ACF(ret^2)", acf(ret**2, nacf)
      write (*,fmt_par) "ACF((ret/sigma_est)^2)", acf((ret/sigma_est)**2, nacf)
   end if
end do
end do
call watch("fit " // trim(str(nfits)) // " GARCH models")
call print_elapsed_times()
print "(/,a)", "(2) finished xgarch_data"
end program xgarch_data
