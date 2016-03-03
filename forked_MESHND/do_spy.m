function do_spy (A)
%DO_SPY use cspy(A) to plot a matrix, or spy(A) if cspy not installed.
try
    % This function is in CSparse.  It generates better looking plots than spy.
    spy (A) ;
catch
    spy (A) ;
end
