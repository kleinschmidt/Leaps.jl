module Leaps

import Base: show

export RegSubsets, QRLeaps, makeqr, tolset, ssleaps, sing, initr, xhaust!


const libleaps = joinpath(Pkg.dir("Leaps"), "deps", "libleaps.so")
## const libleaps = "/Library/Frameworks/R.framework/Resources/library/leaps/libs/leaps.so"

##       SUBROUTINE MAKEQR(NP,NN,WEIGHTS,TXMAT,YVEC,D,RBAR,THETAB,
##      $     SSERR,IER)
## C     Calls INCLUD to construct Banachiewicz factorisation
## C
## C     
##       INTEGER NP, NN, IER
##       DOUBLE PRECISION WEIGHTS(NN), TXMAT(*), YVEC(NN), D(NP), RBAR(*), 
##      +     THETAB(NP), SSERR
## C     local variables
##       INTEGER I, NRBAR

## (pretty sure that INTEGER --> Int32 and DOUBLE PRECISION --> Float64)
## (makes sense given that all Rs numbers are double)


X = rand(100, 6)
y = rand(100)
wts = ones(size(y))


type RegSubsets
    np::Int64                           # Number of predictors
    nn::Int64                           # Number of cases?
    nrbar::Int64                        # Length of upper triangular of R from
                                        # QR decomp of X
    rbar::Vector{Float64}               # upper triangular part of R
    d::Vector{Float64}                  # sqrt of row magnitudes from QR
    thetab::Vector{Float64}             # ??
    sserr::Float64                      # ??
    tol::Vector{Float64}
    lindep::Vector{Bool}
    rss::Vector{Float64}
    bound::Vector{Float64}
    ress::Matrix{Float64}
    lopt::Matrix{Int64}
    nvmax::Int64
    nbest::Int64
    il::Int64
    ir::Int64
    vorder::Vector{Int64}
    intercept::Bool
    nullrss::Float64
end


function RegSubsets(X::Matrix{Float64}, y::Vector{Float64},
                    wts::Vector{Float64} = ones(size(y));
                    intercept::Bool = true,
                    nvmax::Int64 = 8, nbest::Int64 = 1)
    nn, np = size(X)
    nn == size(y, 1) ||
        error(Base.LinAlg.DimensionMismatch("length of y must match rows in X"))
    nn == length(wts) ||
        error(Base.LinAlg.DimensionMismatch("length of wts must match rows in X"))

    nvmax = min(nvmax, np)
    
    if intercept
        np += 1
        X = hcat(ones(nn,1), X)
        nvmax += 1
    end

    vorder = Int64[1:np]
    il = int(nvmax * (nvmax+1) / 2)
    nrbar = int(np * (np-1) / 2)

    qrleap = makeqr(X, y, wts, nn, np)
    tol = tolset(qrleap)
    rss = ssleaps(qrleap)

    ## check for linear dependency between variables
    lindep = sing(qrleap, tol, rss)
    if any(lindep)
        ## figure out how many non-linearly-dependent variables we have
        n_dep = sum(lindep)
        warn("Linear dependencies found: $n_dep")
        
        ## re-order the variables so that all the dependent are at the end
        ## and then try the whole thing again.
        first_dep = findfirst(lindep)
        if first_dep < (length(lindep) - n_dep + 1)
            ## dependent variables are NOT all at the end.  reorder to put them
            ## there
            warn("  Re-ordering and re-trying...")
            dep_at_end_order = sortperm(lindep)
            return RegSubsets(X[:, dep_at_end_order], y, wts,
                              intercept=intercept, nvmax=nvmax, nbest=nbest)
        end

        ## no re-ordering needed, reduce the number of variables to consider
        nvmax = first_dep - 1
        warn("Max number of variables reduced to $(nvmax - intercept)")
        ## recompute residual sum of squares (works because sing modifies the
        ## QRLeaps object in place because it's fields are passed by reference
        rss = ssleaps(qrleap)
    end

    bound, ress, lopt = initr(np, nvmax, nbest, il, vorder, rss)
    
    nullrss = intercept ? rss[1] : sum(y.^2)

    ir = nvmax
    
    RegSubsets(np, nn, nrbar, qrleap.rbar, qrleap.d, qrleap.thetab,
               qrleap.sserr[1], tol, lindep, rss, bound, ress, lopt,
               nvmax, nbest, il, ir, vorder, intercept, nullrss)
end


function xhaust!(rs::RegSubsets)
    first = int(1)
    last = int(rs.np)

    dimwk = 3*last
    dimiwk = rs.nvmax
    wk = Array(Float64, dimwk)
    iwk = Array(Int32, dimiwk)

    ier = Int32[0]

    vorder = convert(Array{Int32}, rs.vorder)
    lopt = convert(Array{Int32}, rs.lopt)
    
    ccall((:xhaust_, libleaps), Void,
          (Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
           Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64},
           Ptr{Float64}, Ptr{Int32}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Ptr{Int32}),
          &rs.np, &rs.nrbar, rs.d, rs.rbar, rs.thetab,
          &first, &last, vorder, rs.tol, rs.rss,
          rs.bound, &rs.nvmax, rs.ress, &rs.ir, &rs.nbest,
          lopt, &rs.il, wk, &dimwk, iwk,
          &dimiwk, ier)

    rs.vorder = vorder
    rs.lopt = lopt

    
    
    if ier[1] == 0
        return rs
    else
        error("leaps: error in xhaust ($(ier[1]))")
    end
end


## Some kind of QR decomposition.
## c.f. ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
## http://lib.stat.cmu.edu/apstat/274
type QRLeaps
    d::Array{Float64,1}
    rbar::Array{Float64,1}
    thetab::Array{Float64,1}
    sserr::Array{Float64,1}             # Represented as an array of length
                                        # one to allow modification in place
    np::Int64
    nn::Int64
end

function show(io::IO, qr::QRLeaps)
    println("QRLeaps ($(qr.nn) x $(qr.np)): ")
    println("  d:      $(qr.d)")
    println("  thetab: $(qr.thetab)")
    println("  rbar:   $(qr.rbar)")
    println("  sserr:  $(qr.sserr)")
end

function makeqr(X::Matrix{Float64}, y::Vector{Float64},
                wts::Vector{Float64}, nn::Int64, np::Int64)

    d = Array(Float64, np)
    n_rbar = div(np * (np-1), 2)
    rbar = Array(Float64, n_rbar)
    thetab = Array(Float64, np)
    sserr = Float64[0]
    ier = Int32[0]
    
    ccall((:makeqr_, libleaps), Void,
          (Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, # input types
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}), # output types
          &np, &nn, wts, X', y,         # input values
          d, rbar, thetab, sserr, ier) # output values

    if ier[1]==0
        ## no error, return d, rbar, tehtab, sserr
        return QRLeaps(d, rbar, thetab, sserr, np, nn)
    else
        error("leaps: error in qrmake ($(ier[1]))")
    end
end



function tolset(qr::QRLeaps)
    workspace = Array(Float64, qr.np)
    tol = Array(Float64, qr.np)
    ier = Int32[0]

    ccall((:tolset_, libleaps), Void,
          (Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64},
           Ptr{Float64}, Ptr{Float64}, Ptr{Int32}),
          &qr.np, &length(qr.rbar), qr.d, qr.rbar,
          tol, workspace, ier)

    if ier[1]==0
        ## no error, return tol
        return tol
    else
        error("leaps: error in tolset ($(ier[1]))")
    end
    
end

function ssleaps(qr::QRLeaps)
    ier = Int32[0]
    rss = Array(Float64, qr.np)

    ccall((:ssleaps_, libleaps), Void,
          (Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
           Ptr{Float64}, Ptr{Int32}),
          &qr.np, qr.d, qr.thetab, qr.sserr,
          rss, ier)

    if ier[1]==0
        ## no error, return rss
        return rss
    else
        error("leaps: error in ssleaps ($(ier[1]))")
    end
end

function sing(qr::QRLeaps, tol::Vector{Float64}, rss::Vector{Float64})
    ier = Int32[0]
    lindep = Array(Int32, qr.np)
    workspace = Array(Float64, qr.np)

    ccall((:sing_, libleaps), Void,
          (Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64},
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
           Ptr{Int32}, Ptr{Float64}, Ptr{Int32}),
          &qr.np, &length(qr.rbar), qr.d, qr.rbar,
          qr.thetab, qr.sserr, tol,
          lindep, workspace, ier)

    ## following leaps.R's lead, only errors >0 are a problem
    if ier[1]<=0
        ## no error, return lindep
        return convert(Vector{Bool}, lindep)
    else
        error("leaps: error in sing ($(ier[1])). (lindep returned as $lindep )")
    end
end

function initr(np::Int64, nvmax::Int64, nbest::Int64, il::Int64,
               vorder::Vector{Int64}, rss::Vector{Float64})
    ir = nvmax
    bound = Array(Float64, np)
    ress = Array(Float64, ir, nbest)   # make me a matrix?
    lopt = Array(Int32, il, nbest)
    ier = Int32[0]
    ## vorder_32 = convert(Vector{Int32}, vorder)

    ## ccall((:initr_, libleaps), Void,
    ##       (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64},
    ##        Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
    ##        Ptr{Int32}, Ptr{Float64}, Ptr{Int32}),
    ##       &np, &nvmax, &nbest, bound,
    ##       ress, &nvmax, lopt, il,
    ##       vorder_32, rss, ier)

    ## vorder[:] = convert(Vector{Int32}, vorder_32)[:]
    
    for best in 1:nbest
        pos = 1
        for nvar in 1:nvmax
            if best == 1
                ress[nvar, best] = rss[nvar]
            else
                ress[nvar, best] = typemax(eltype(ress))
            end
            if best == nbest
                bound[nvar] = ress[nvar, nbest]
            end
            for i in 1:nvar
                if best == 1
                    lopt[pos, best] = vorder[i]
                else
                    lopt[pos, best] = 0
                end
                pos += 1
            end
        end
    end

    
    if ier[1]==0
        return (bound, ress, lopt)
    else
        error("leaps: error in initr ($(ier[1])).")
    end
    
end

end


