module optim_methods_mod
implicit none
private
public :: nelder_mead_str, uobyqa_str
character (len=*), parameter :: nelder_mead_str = "nelder-mead", uobyqa_str = "uobyqa"
end module optim_methods_mod