

function mpccmodel_setup_newton_penalty(
        config::MPCCModelConfig        
    )

    @unpack dimspec = config
    @unpack n, q, l, me, mi, r, s = dimspec

    # Std variables
    @variables x[1:n] pr[1:r] ps[1:s]

    # The masking variables (these will either be 0 or 1 when called)
    @variables m_ci_n[1:mi] m_ci_p[1:mi] m_F_n[1:l, 1:q] m_Fq[1:q]

    # The penalty coefficients
    @variables α_ci α_F β_ci γ

    println(" ** Doing max formulation **")

    # Construct the penalised objective function
    ϕ = config.defn.f
    for lp_ci=1:mi
        ϕ += α_ci * m_ci_n[lp_ci] * max(0, -config.defn.ci[lp_ci])
    end
    for lp_ci=1:mi
        ϕ += β_ci * m_ci_p[lp_ci] * config.defn.ci[lp_ci]
    end
    for lp_q=1:q
        for lp_l=1:l
            ϕ += α_F * m_F_n[lp_l, lp_q] * max(0, -config.defn.F[lp_l, lp_q])
        end
    end
    temp_Fq = Vector{Num}(undef, q)
    for lp_q=1:q
        temp_Fq[lp_q] = Num(1)
        for lp_l=1:l
            temp_Fq[lp_q] *= config.defn.F[lp_l, lp_q]
        end
        ϕ += γ * m_Fq[lp_q] * temp_Fq[lp_q]
    end

    ϕ_fn = build_function([ϕ], x, pr, ps, α_ci, α_F, β_ci, γ, m_ci_n, m_ci_p, m_F_n, m_Fq; expression=Val{false})

    function local_gradϕ(
            x::Vector{S}, pr::Vector{T}, ps::Vector{Int64},
            α_ci::S, α_F::S, β_ci::S, γ::S,
            m_ci_n::Vector{Float64}, m_ci_p::Vector{Float64}, m_F_n::Matrix{Float64}, m_Fq::Vector{Float64}           
        ) where {S <: Real, T <: Real}  # ::Vector{S}

        return mm_fd_ntpen_gradϕ(dimspec, ϕ_fn[1], x, pr, ps, α_ci, α_F, β_ci, γ, m_ci_n, m_ci_p, m_F_n, m_Fq)
    end

    function local_hessϕ(
            x::Vector{S}, pr::Vector{T}, ps::Vector{Int64},
            α_ci::S, α_F::S, β_ci::S, γ::S,
            m_ci_n::Vector{Float64}, m_ci_p::Vector{Float64}, m_F_n::Matrix{Float64}, m_Fq::Vector{Float64}           
        ) where {S <: Real, T <: Real}  # ::Vector{S}

        return mm_fd_ntpen_hessϕ(dimspec, ϕ_fn[1], x, pr, ps, α_ci, α_F, β_ci, γ, m_ci_n, m_ci_p, m_F_n, m_Fq)
    end

    return MPCCNewtonPenaltyFunctions(
            ϕ_fn[1],
            local_gradϕ,
            local_hessϕ
        )
end





# Define gradient wrt x
function mm_fd_ntpen_gradϕ(
        dimspec::MPCCDimSpec, ϕ_fn::Function,
        x::Vector{S}, pr::Vector{T}, ps::Vector{Int64},
        α_ci::S, α_F::S, β_ci::S, γ::S,
        m_ci_n::Vector{Float64}, m_ci_p::Vector{Float64}, m_F_n::Matrix{Float64}, m_Fq::Vector{Float64}
    ) where {S <: Real, T <: Real}

    local_ϕ(z) = ϕ_fn(z, pr, ps, α_ci, α_F, β_ci, γ, m_ci_n, m_ci_p, m_F_n, m_Fq)[1]
    return ForwardDiff.gradient(local_ϕ, x)
end


# Define gradient wrt x
function mm_fd_ntpen_hessϕ(
        dimspec::MPCCDimSpec, ϕ_fn::Function,
        x::Vector{S}, pr::Vector{T}, ps::Vector{Int64},
        α_ci::S, α_F::S, β_ci::S, γ::S,
        m_ci_n::Vector{Float64}, m_ci_p::Vector{Float64}, m_F_n::Matrix{Float64}, m_Fq::Vector{Float64}
    ) where {S <: Real, T <: Real}

    local_ϕ(z) = ϕ_fn(z, pr, ps, α_ci, α_F, β_ci, γ, m_ci_n, m_ci_p, m_F_n, m_Fq)[1]
    return ForwardDiff.hessian(local_ϕ, x)
end

