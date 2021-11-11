using MPCCModelManager, Test, UnPack

# @test Int64 <: Integer
# @test Float64 <: AbstractFloat


model_cfg = mpccmodel_load_defn_from_file("dev1");

@unpack n, q, l, me, mi, r, s = model_cfg.dimspec


x = [ 1.0, 2, 3, 4 ]
pr0 = [ -0.5, 1.0 ]
ps0 = [ 1 ]

#--------------------------------------------------------------------------------
# Check mutating functions match standard return values
#--------------------------------------------------------------------------------

# f
res_f = model_cfg.fns.f(x, pr0, ps0)
res_f_mut = similar(res_f) * 0.0
res_f_mut *= 0.0
model_cfg.fns.f!(res_f_mut, x, pr0, ps0)
@test res_f ≈ res_f_mut

# ce
res_ce = model_cfg.fns.ce(x, pr0, ps0)
res_ce_mut = similar(res_ce) * 0.0
res_ce_mut *= 0.0
model_cfg.fns.ce!(res_ce_mut, x, pr0, ps0)
@test res_ce ≈ res_ce_mut

# ci
res_ci = model_cfg.fns.ci(x, pr0, ps0)
res_ci_mut = similar(res_ci) * 0.0
res_ci_mut *= 0.0
model_cfg.fns.ci!(res_ci_mut, x, pr0, ps0)
@test res_ci ≈ res_ci_mut

# F
res_F = model_cfg.fns.F(x, pr0, ps0)
res_F_mut = similar(res_F) * 0.0
model_cfg.fns.F!(res_F_mut, x, pr0, ps0)
@test res_F ≈ res_F_mut

# ce indexed
for lp_ce=1:me
    res_ce_i = model_cfg.fns.ce_i[lp_ce](x, pr0, ps0)
    res_ce_i_mut = similar(res_ce_i) * 0.0
    model_cfg.fns.ce_i![lp_ce](res_ce_i_mut, x, pr0, ps0)
    @test res_ce_i ≈ res_ce_i_mut
end

# ci indexed
for lp_ci=1:mi
    res_ci_i = model_cfg.fns.ci_i[lp_ci](x, pr0, ps0)
    res_ci_i_mut = similar(res_ci_i) * 0.0
    model_cfg.fns.ci_i![lp_ci](res_ci_i_mut, x, pr0, ps0)
    @test res_ci_i ≈ res_ci_i_mut
end

# F indexed
for lp_l=1:l
    for lp_q=1:q
        res_F_i = model_cfg.fns.F_i[lp_l, lp_q](x, pr0, ps0)
        res_F_i_mut = similar(res_F_i) * 0.0
        model_cfg.fns.F_i![lp_l, lp_q](res_F_i_mut, x, pr0, ps0)
        @test res_F_i ≈ res_F_i_mut
    end
end