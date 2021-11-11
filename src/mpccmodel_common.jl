

# TODO Change Vector{Num} to either Symbolics.Arr{Num, 1} or something more general. May require some finessing


export  MPCCDimSpec,
        MPCCFunctions,
        MPCCPointEvalReq,
        MPCCPointEval,
        MPCCModelTestVector,
        MPCCModelNZMask,        
        MPCCModelConfig,
        MPCCModel,
        MPCCJumpFixedFunctions,
        MPCCNewtonPenaltyFunctions

# Trick from https://discourse.julialang.org/t/aliases-for-union-t-nothing-and-union-t-missing/15402/15
const Opt{T} = Union{Missing,T}



# Information on the dimensions of a problem
struct MPCCDimSpec
    n::Int64		# Dimension of spatial variables
    q::Int64		# Columns in F
    l::Int64        # Number of F functions per complementarity variable (usually 2)
    me::Int64		# Number of equality constraints
    mi::Int64		# Number of inequality constraints
	r::Int64		# Number of continuous parameters for PMPCC
	s::Int64		# Number of discrete parameters for PMPCC
end
function MPCCDimSpec(dict::Dict)
    return MPCCDimSpec(
        dict["n"],
        dict["q"],
        dict["l"],
        dict["me"],
        dict["mi"],
        dict["r"],
        dict["s"]
    )
end






struct MPCCDefinition
    f::Num
    ce::Vector{Num}
    ci::Vector{Num}
    F::Matrix{Num}
end



# Core functions defining a model, returns (possibly indexed subselection of vectors) of functions 
struct MPCCFunctions
    f::Function
    f!::Function
    ce::Function                   # Returns vector of length me
    ce!::Function
    ce_i::AbstractVector           # Vector so we can write scalar values in to mutating function argument
    ce_i!::AbstractVector
    ci::Function                   # Returns vector of length mi
    ci!::Function
    ci_i::AbstractVector
    ci_i!::AbstractVector
    F::Function                    # Returns vatrix of dim l x q
    F!::Function
    F_i::AbstractMatrix
    F_i!::AbstractMatrix
    Fq::Function                    # Returns vector of length q
    Fq!::Function
    Fq_i::AbstractVector
    Fq_i!::AbstractVector
end



# Struct to state what should be evaluated
struct MPCCPointEvalReq
    f::Bool
    ce::Bool
    ci::Bool
    F::Bool
    gradf::Bool
    gradce::Bool
    gradci::Bool
    gradF::Bool
    hessf::Bool
    hessce::Bool
    hessci::Bool
    hessF::Bool
    fdp::Bool
    cedp::Bool
    cidp::Bool
    Fdp::Bool
    gradfdp::Bool
    gradcedp::Bool
    gradcidp::Bool
    gradFdp::Bool                           
end
function MPCCPointEvalReq()
    return MPCCPointEvalReq(    true, true, true, true,
                                true, true, true, true,
                                true, true, true, true,
                                true, true, true, true,
                                true, true, true, true
                            )
end




# Evaluation results of a problem at a point
struct MPCCPointEval{T <: AbstractFloat}
    f::Opt{T}
    ce::Opt{Vector{T}}                          # Vector of length me
    ci::Opt{Vector{T}}                          # Vector of length mi
    F::Opt{Matrix{T}}                           # Matrix of dim l x q
    gradf::Opt{Vector{T}}                       # Vector of length n
    gradce::Opt{Matrix{T}}                      # Changed 20211029: now matrix of me by n
    gradci::Opt{Matrix{T}}                      # Changed 20211029: now matrix of mi by n
    gradF::Opt{Matrix{Vector{T}}}               # Matrix of dim l x q of n length vectors (the gradients wrt x)
    hessf::Opt{Matrix{T}}                       # Matrix dim n x n
    hessce::Opt{Vector{Matrix{T}}}              # Vector of length me of 2d Hessian matrices
    hessci::Opt{Vector{Matrix{T}}}              # Vector of length mi of 2d Hessian matrices
    hessF::Opt{Matrix{Matrix{T}}}               # Matrix of dim l x q of 2d Hessian matrices
    fdp::Opt{Vector{T}}                         # Vector of length r
    cedp::Opt{Vector{Vector{T}}}                # Vector of length r of vectors of length me
    cidp::Opt{Vector{Vector{T}}}                # Vector of length r of vectors of length mi
    Fdp::Opt{Vector{Matrix{T}}}                 # Vector of length r of matrix of dim l x q
    gradfdp::Opt{Vector{Vector{T}}}             # Vector of length r of gradients of f/dp (vectors of length n)
    gradcedp::Opt{Vector{Matrix{T}}}            # Vector of length r of gradce, each /dp
    gradcidp::Opt{Vector{Matrix{T}}}            # Vector of length r of gradci, each /dp
    gradFdp::Opt{Vector{Matrix{Vector{T}}}}     # Vector of length r of gradF, each /dp
end
function MPCCPointEval()  # Quick way to create empty/zero struct
	return MPCCPointEval(
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing
    )
end
Base.:(==)(x::MPCCPointEval, y::MPCCPointEval) = (
	x.f == y.f &&
    x.ce == y.ce &&
    x.ci == y.ci &&
	x.F == y.F &&
	x.gradf == y.gradf &&
    x.gradce == y.gradce &&
    x.gradci == y.gradci &&
    x.gradF == y.gradF &&
    x.hessf == y.hessf &&
    x.hessce == y.hessce &&
    x.hessci == y.hessci &&
    x.hessF == y.hessF &&
    x.fdp == y.fdp &&
    x.cedp == y.cedp &&
    x.cidp == y.cidp &&
    x.Fdp == y.Fdp &&
    x.gradfdp == y.gradfdp &&
    x.gradcedp == y.gradcedp &&
    x.gradcidp == y.gradcidp &&
    x.gradFdp == y.gradFdp )
Base.:(≈)(x::MPCCPointEval, y::MPCCPointEval) = (
    x.f ≈ y.f &&
    x.ce ≈ y.ce &&
    x.ci ≈ y.ci &&
    x.F ≈ y.F &&
    x.gradf ≈ y.gradf &&
    x.gradce ≈ y.gradce &&
    x.gradci ≈ y.gradci &&
    x.gradF ≈ y.gradF &&
    x.hessf ≈ y.hessf &&
    x.hessce ≈ y.hessce &&
    x.hessci ≈ y.hessci &&
    x.hessF ≈ y.hessF &&
    x.fdp ≈ y.fdp &&
    x.cedp ≈ y.cedp &&
    x.cidp ≈ y.cidp &&
    x.Fdp ≈ y.Fdp &&
    x.gradfdp ≈ y.gradfdp &&
    x.gradcedp ≈ y.gradcedp &&
    x.gradcidp ≈ y.gradcidp &&
    x.gradFdp ≈ y.gradFdp )





# An aid for plotting only, determines which variables are used in which constraints
struct MPCCModelNZMask
    ce::Vector{Set{Int64}}
    ci::Vector{Set{Int64}}
    F::Vector{Vector{Set{Int64}}}
end


struct MPCCModelTestVector{R <: Real, S <: Real, T <: Real}
    x_val::Vector{S}
    pr_val::Vector{T}
    ps_val::Vector{Int64}
    eval_test_pt::MPCCPointEval{R}
end



struct MPCCParameterisationDefn{R <: Real}
    pr::Vector{Num}
    # prdt::Num
    tspan::Tuple{R, R}
    descr::String
end


struct MPCCParameterisationFunctions
    pr::Function            # t -> Vector pr
    pr!::Function
    prdt::Function          # t -> Vector of d(pr) / dt
    prdt!::Function
end



struct MPCCParameterisations
    t::Num  # Symbolics variable
    defns::Vector{MPCCParameterisationDefn}
    fns::Vector{MPCCParameterisationFunctions}
end




struct MPCCModelConfig
    x::Vector{Num}      # Symbolics variables
    pr::Vector{Num}     # Symbolics variables
    ps::Vector{Num}     # Symbolics variables
    dimspec::MPCCDimSpec
    defn::MPCCDefinition
    fns::MPCCFunctions
    nzmask::Opt{MPCCModelNZMask}
    testvectors::Vector{MPCCModelTestVector}
    knownsols::Vector{Function}
    parameterisations::MPCCParameterisations
end



# Notes: add additional methods to handle calls to MPCCModel functions with only t parameter.

struct MPCCModel
    config::MPCCModelConfig
    f::Function
    f!::Opt{Function}
    ce::Function
    ce!::Opt{Function}
    ci::Function
    ci!::Opt{Function}
    F::Function
    F!::Opt{Function}
    Fq::Function
    Fq!::Opt{Function}
    gradf::Function
    gradf!::Opt{Function}
    gradce::Function
    gradce!::Opt{Function}
    gradci::Function
    gradci!::Opt{Function}
    gradF::Function
    gradF!::Opt{Function}
    gradFq::Function
    gradFq!::Opt{Function}
    hessf::Function
    hessf!::Opt{Function}
    hessce::Function
    hessce!::Opt{Function}
    hessci::Function
    hessci!::Opt{Function}
    hessF::Function
    hessF!::Opt{Function}
    hessFq::Function
    hessFq!::Opt{Function}
    fdp::Function
    fdp!::Opt{Function}
    cedp::Function
    cedp!::Opt{Function}
    cidp::Function
    cidp!::Opt{Function}
    Fdp::Function
    Fdp!::Opt{Function}
    gradfdp::Function
    gradfdp!::Opt{Function}
    gradcedp::Function
    gradcedp!::Opt{Function}
    gradcidp::Function
    gradcidp!::Opt{Function}
    gradFdp::Function
    gradFdp!::Opt{Function}
end



struct MPCCJumpFixedFunctions
    jump_f::Function
    jump_f_pen::Function
    jump_ce::Vector{Function}
    jump_ci::Vector{Function}
    jump_F::Matrix{Function}
end


# Only std return functions for the moment
struct MPCCNewtonPenaltyFunctions
    ϕ::Function
    gradϕ::Function
    hessϕ::Function
end
