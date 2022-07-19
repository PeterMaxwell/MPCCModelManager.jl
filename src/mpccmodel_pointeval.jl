

function mpccmodel_pointeval_basic(model::AbstractMPCCModel, ptevalreq::MPCCPointEvalReq, x0::Vector{S}, pr0::Vector{T}, ps0::Vector{Int64}) where {S <: Real, T <: Real}

    # std entries
    if ( ptevalreq.f )
        f = model.f(x0, pr0, ps0)
    else
        f = missing
    end

    if ( ptevalreq.ce )
        ce = model.ce(x0, pr0, ps0)
    else
        ce = missing
    end

    if ( ptevalreq.ci )
        ci = model.ci(x0, pr0, ps0)
    else
        ci = missing
    end

    if ( ptevalreq.F )
        F = model.F(x0, pr0, ps0)
    else
        F = missing
    end

    # grad entries
    if ( ptevalreq.gradf )
        gradf = model.gradf(x0, pr0, ps0)
    else
        gradf = missing
    end

    if ( ptevalreq.gradce )
        gradce = model.gradce(x0, pr0, ps0)
    else
        gradce = missing
    end

    if ( ptevalreq.gradci )
        gradci = model.gradci(x0, pr0, ps0)
    else
        gradci = missing
    end

    if ( ptevalreq.gradF )
        gradF = model.gradF(x0, pr0, ps0)
    else
        gradF = missing
    end

    # hess entries
    if ( ptevalreq.hessf )
        hessf = model.hessf(x0, pr0, ps0)
    else
        hessf = missing
    end

    if ( ptevalreq.hessce )
        hessce = model.hessce(x0, pr0, ps0)
    else
        hessce = missing
    end

    if ( ptevalreq.hessci )
        hessci = model.hessci(x0, pr0, ps0)
    else
        hessci = missing
    end

    if ( ptevalreq.hessF )
        hessF = model.hessF(x0, pr0, ps0)
    else
        hessF = missing
    end

    # dp entries
    if ( ptevalreq.fdp )
        fdp = model.fdp(x0, pr0, ps0)
    else
        fdp = missing
    end

    if ( ptevalreq.cedp )
        cedp = model.cedp(x0, pr0, ps0)
    else
        cedp = missing
    end

    if ( ptevalreq.cidp )
        cidp = model.cidp(x0, pr0, ps0)
    else
        cidp = missing
    end

    if ( ptevalreq.Fdp )
        Fdp = model.Fdp(x0, pr0, ps0)
    else
        Fdp = missing
    end

    # graddp entries
    if ( ptevalreq.gradfdp )
        gradfdp = model.gradfdp(x0, pr0, ps0)
    else
        gradfdp = missing
    end

    if ( ptevalreq.gradcedp )
        gradcedp = model.gradcedp(x0, pr0, ps0)
    else
        gradcedp = missing
    end

    if ( ptevalreq.gradcidp )
        gradcidp = model.gradcidp(x0, pr0, ps0)
    else
        gradcidp = missing
    end

    if ( ptevalreq.gradFdp )
        gradFdp = model.gradFdp(x0, pr0, ps0)
    else
        gradFdp = missing
    end

    return MPCCPointEval{promote_type(S, T)}(
        f, ce, ci, F,
        gradf, gradce, gradci, gradF,
        hessf, hessce, hessci, hessF,
        fdp, cedp, cidp, Fdp,
        gradfdp, gradcedp, gradcidp, gradFdp
    )

end