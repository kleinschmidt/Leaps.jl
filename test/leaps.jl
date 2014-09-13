using DataFrames

## generate data (and export so R can read)

srand(1)
X = rand(100, 6)
y = rand(100)

d = convert(DataFrame, X)
d[:y] = y


writetable(joinpath(Pkg.dir("Leaps"), "test", "test.csv"), d)

using Leaps

## test normal case
rs = RegSubsets(X, y, intercept=false)
xhaust!(rs)

## test case with extra linear dependency
rs2 = RegSubsets(hcat(X[:,1], X), y, intercept=false)
xhaust!(rs2)
