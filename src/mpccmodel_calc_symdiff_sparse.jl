

function mpccmodel_setup_symdiff_sparse(config::MPCCModelConfig)

    @unpack dimspec, defn, fns = config
    @unpack n, q, l, me, mi, r, s = config.dimspec

    s_x = config.x
    s_pr = config.pr
    s_ps = config.ps

    # We need to setup the num arrays for the sparse jacobian and hessians
    # Note that although SparseMatrixCSC is used consistently below, most of these
    # are vectors; Symbolics spits out the mtrx type and can't be bothered converting

    # grad entries
    gradf_sparse_num = transpose( Symbolics.sparsejacobian([ defn.f ], s_x) )
    gradce_sparse_num = transpose( Symbolics.sparsejacobian( defn.ce, s_x ) )
    gradci_sparse_num = transpose( Symbolics.sparsejacobian( defn.ci, s_x ) )
    gradF_sparse_num = Matrix{SparseMatrixCSC}(undef, l, q)
    for lp_q=1:q
        for lp_l=1:l
            gradF_sparse_num[lp_l, lp_q] = transpose( Symbolics.sparsejacobian( [ defn.F[lp_l, lp_q] ], s_x ) )
        end
    end

    println("gradf: ", gradf_sparse_num)
    println("gradce: ", gradce_sparse_num)
    println("gradci: ", gradci_sparse_num)
    println("gradF: ", gradF_sparse_num)


    # hess entries
    hessf_sparse_num = Symbolics.sparsehessian(defn.f, s_x)
    hessce_sparse_num = Vector{SparseMatrixCSC}(undef, me)
    for lp_ce=1:me
        hessce_sparse_num[lp_ce] = Symbolics.sparsehessian(defn.ce[lp_ce], s_x)
    end
    hessci_sparse_num = Vector{SparseMatrixCSC}(undef, mi)
    for lp_ci=1:mi
        hessci_sparse_num[lp_ci] = Symbolics.sparsehessian(defn.ci[lp_ci], s_x)
    end
    hessF_sparse_num = Matrix{SparseMatrixCSC}(undef, l, q)
    for lp_q=1:q
        for lp_l=1:l
            hessF_sparse_num[lp_l, lp_q] = Symbolics.sparsehessian(defn.F[lp_l, lp_q], s_x)
        end
    end

    println("hessf: ", hessf_sparse_num)
    println("hessce: ", hessce_sparse_num)
    println("hessci: ", hessci_sparse_num)
    println("hessF: ", hessF_sparse_num)


    # dp entries
    fdp_sparse_num = Vector{Num}(undef, r)
    for lp_r=1:r
        fdp_sparse_num[lp_r] = Symbolics.derivative(defn.f, s_pr[lp_r]; simplify=true)
    end
    cedp_sparse_num = Vector{Vector{Num}}(undef, r)
    for lp_r=1:r        
        cedp_sparse_num[lp_r] = Symbolics.derivative(defn.ce, s_pr[lp_r]; simplify=true)
    end
    cidp_sparse_num = Vector{Vector{Num}}(undef, r)
    for lp_r=1:r        
        cidp_sparse_num[lp_r] = Symbolics.derivative(defn.ci, s_pr[lp_r]; simplify=true)
    end
    Fdp_sparse_num = Vector{Matrix{Num}}(undef, r)
    for lp_r=1:r        
        Fdp_sparse_num[lp_r] = Symbolics.derivative(defn.F, s_pr[lp_r]; simplify=true)
    end

    println("fdp: ", fdp_sparse_num)
    println("cedp: ", cedp_sparse_num)
    println("cidp: ", cidp_sparse_num)
    println("Fdp: ", Fdp_sparse_num)


    # graddp entries
    gradfdp_sparse_num = Vector{SparseMatrixCSC}(undef, r)
    for lp_r=1:r
        gradfdp_sparse_num[lp_r] = transpose( Symbolics.sparsejacobian([fdp_sparse_num[lp_r]], s_x) )        
    end
    gradcedp_sparse_num = Vector{SparseMatrixCSC}(undef, r)
    for lp_r=1:r
        gradcedp_sparse_num[lp_r] = transpose( Symbolics.sparsejacobian( cedp_sparse_num[lp_r], s_x ) )
    end
    gradcidp_sparse_num = Vector{SparseMatrixCSC}(undef, r)
    for lp_r=1:r
        gradcidp_sparse_num[lp_r] = transpose( Symbolics.sparsejacobian( cidp_sparse_num[lp_r], s_x ) )
    end
    gradFdp_sparse_num = Vector{Matrix{SparseMatrixCSC}}(undef, r)
    for lp_r=1:r
        gradFdp_sparse_num[lp_r] = Matrix{SparseMatrixCSC}(undef, l, q)
        for lp_q=1:q
            for lp_l=1:l
                gradFdp_sparse_num[lp_r][lp_l, lp_q] = transpose( Symbolics.sparsejacobian( [ Fdp_sparse_num[lp_r][lp_l, lp_q] ], s_x ) )
            end
        end
    end

    println("gradfdp: ", gradfdp_sparse_num)
    println("gradcedp: ", gradcedp_sparse_num)
    println("gradcidp: ", gradcidp_sparse_num)
    println("gradFdp: ", gradFdp_sparse_num)

    # f, ce, ci, F all return usual dense vectors/matrices because they are, well, dense.

    # Full f functions (really just converting a 1 element vector to a scalar)
    function local_f(x0::Vector{S}, pr0::Vector{T}, ps0::Vector{Int64}) where {S <: Real, T <: Real}
        return fns.f(x, pr, ps)[1]
    end



    # Full ce functions
    function local_ce(x0::Vector{S}, pr0::Vector{T}, ps0::Vector{Int64}) where {S <: Real, T <: Real}   # ::Vector{S}
        return fns.ce(x, pr, ps)
    end



    # Indexed ce functions
    function local_ce(x0::Vector{S}, pr0::Vector{T}, ps0::Vector{Int64}, idxs_ce::Vector{Int64}) where {S <: Real, T <: Real} # ::Vector{S}
        return mm_fd_dn_ce_i(dimspec, fns.ce_i, x, pr, ps, idxs_ce)
    end
    


    # Full ci functions
    function local_ci(x0::Vector{S}, pr0::Vector{T}, ps0::Vector{Int64}) where {S <: Real, T <: Real}  # ::Vector{S}
        return fns.ci(x, pr, ps)
    end


    # Indexed ci functions
    function local_ci(x0::Vector{S}, pr0::Vector{T}, ps0::Vector{Int64}, idxs_ci::Vector{Int64}) where {S <: Real, T <: Real}  # ::Vector{S}
        return mm_fd_dn_ci_i(dimspec, fns.ci_i, x, pr, ps, idxs_ci)
    end
    



    # Full F functions
    function local_F(x0::Vector{S}, pr0::Vector{T}, ps0::Vector{Int64}) where {S <: Real, T <: Real}  # ::Matrix{S}
        return fns.F(x, pr, ps)
    end




    # Indexed F functions; NOTE, these yeild a vector in the order that the (l,q) index tuples were in!
    function local_F(x0::Vector{S}, pr0::Vector{T}, ps0::Vector{Int64}, idxs_F::Vector{Tuple{Int64,Int64}}) where {S <: Real, T <: Real}  # ::Vector{S}
        return mm_fd_dn_F_i(dimspec, fns.F_i, x, pr, ps, idxs_F)
    end
    



    # gradf functions
    gradf_sparse_fn = build_function(gradf_sparse_num, x, pr, ps; expression=Val{false})

    function local_gradf(x0::Vector{S}, pr0::Vector{T}, ps0::Vector{Int64}) where {S <: Real, T <: Real}
        return gradf_sparse_fn[1](x0, pr0, ps0)
    end



    # Development placeholder functions
    # local_f(x0::Vector) = sum(x)
    # function local_f!(out_f::AbstractVector, x0::Vector) 
    #     out_f[1] = sum(x)
    #     return nothing
    # end

    # local_ce(x0::Vector) = sum(x)
    # function local_ce!(out_ce::AbstractVector, x0::Vector) 
    #     out_ce[1] = sum(x)
    #     return nothing
    # end

    # local_ci(x0::Vector) = sum(x)
    # function local_ci!(out_ci::AbstractVector, x0::Vector) 
    #     out_ci[1] = sum(x)
    #     return nothing
    # end

    # local_F(x0::Vector) = sum(x)
    # function local_F!(out_F::AbstractVector, x0::Vector) 
    #     out_F[1] = sum(x)
    #     return nothing
    # end

    # local_gradf(x0::Vector) = sum(x)
    # function local_gradf!(out_gradf::AbstractVector, x0::Vector) 
    #     out_gradf[1] = sum(x)
    #     return nothing
    # end
    
    local_gradce(x0::Vector) = sum(x)
    function local_gradce!(out_gradce::AbstractVector, x0::Vector) 
        out_gradce[1] = sum(x)
        return nothing
    end
    
    local_gradci(x0::Vector) = sum(x)
    function local_gradci!(out_gradci::AbstractVector, x0::Vector) 
        out_gradci[1] = sum(x)
        return nothing
    end
    
    local_gradF(x0::Vector) = sum(x)
    function local_gradF!(out_gradF::AbstractVector, x0::Vector) 
        out_gradF[1] = sum(x)
        return nothing
    end
    
    local_hessf(x0::Vector) = sum(x)
    function local_hessf!(out_hessf::AbstractVector, x0::Vector) 
        out_hessf[1] = sum(x)
        return nothing
    end
    
    local_hessce(x0::Vector) = sum(x)
    function local_hessce!(out_hessce::AbstractVector, x0::Vector) 
        out_hessce[1] = sum(x)
        return nothing
    end
    
    local_hessci(x0::Vector) = sum(x)
    function local_hessci!(out_hessci::AbstractVector, x0::Vector) 
        out_hessci[1] = sum(x)
        return nothing
    end
    
    local_hessF(x0::Vector) = sum(x)
    function local_hessF!(out_hessF::AbstractVector, x0::Vector) 
        out_hessF[1] = sum(x)
        return nothing
    end
    
    local_fdp(x0::Vector) = sum(x)
    function local_fdp!(out_fdp::AbstractVector, x0::Vector) 
        out_fdp[1] = sum(x)
        return nothing
    end
    
    local_cedp(x0::Vector) = sum(x)
    function local_cedp!(out_cedp::AbstractVector, x0::Vector) 
        out_cedp[1] = sum(x)
        return nothing
    end
    
    local_cidp(x0::Vector) = sum(x)
    function local_cidp!(out_cidp::AbstractVector, x0::Vector) 
        out_cidp[1] = sum(x)
        return nothing
    end
    
    local_Fdp(x0::Vector) = sum(x)
    function local_Fdp!(out_Fdp::AbstractVector, x0::Vector) 
        out_Fdp[1] = sum(x)
        return nothing
    end

    local_gradfdp(x0::Vector) = sum(x)
    function local_gradfdp!(out_gradfdp::AbstractVector, x0::Vector) 
        out_gradfdp[1] = sum(x)
        return nothing
    end
    
    local_gradcedp(x0::Vector) = sum(x)
    function local_gradcedp!(out_gradcedp::AbstractVector, x0::Vector) 
        out_gradcedp[1] = sum(x)
        return nothing
    end
    
    local_gradcidp(x0::Vector) = sum(x)
    function local_gradcidp!(out_gradcidp::AbstractVector, x0::Vector) 
        out_gradcidp[1] = sum(x)
        return nothing
    end
    
    local_gradFdp(x0::Vector) = sum(x)
    function local_gradFdp!(out_gradFdp::AbstractVector, x0::Vector) 
        out_gradFdp[1] = sum(x)
        return nothing
    end

    # TODO need to use the [1] format
    return MPCCModel(
            config,
            local_f, missing,
            local_ce, missing,
            local_ci, missing,
            local_F, missing,
            local_gradf, missing,
            local_gradce, missing,
            local_gradci, missing,
            local_gradF, missing,
            local_hessf, missing,
            local_hessce, missing,
            local_hessci, missing,
            local_hessF, missing,
            local_fdp, missing,
            local_cedp, missing,
            local_cidp, missing,
            local_Fdp, missing,
            local_gradfdp, missing,
            local_gradcedp, missing,
            local_gradcidp, missing,
            local_gradFdp, missing
        )
end



