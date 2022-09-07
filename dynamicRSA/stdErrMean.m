function sem = stdErrMean(mat,dim)
% computes SEM of 2D mat across columns (dim=1) or rows (dim=2)
if dim==2
    mat = mat';
end

%sem = std(mat)/sqrt(length(mat));
sem = std(mat)./sqrt(size(mat,1));

end
