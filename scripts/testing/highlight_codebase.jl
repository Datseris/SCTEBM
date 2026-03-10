using DynamicalSystems, ConceptualClimateModels
CTMLM = ConceptualClimateModels.CloudToppedMixedLayerModel

eqs = [
    # clouds
    CTMLM.cf_dynamic(:sigmoid),
    CTMLM.decoupling_variable(),
    # mixed layer
    CTMLM.sst_dynamic(),
    CTMLM.mlm_dynamic(),
    CTMLM.mlm_q₊(:relative),
    ⋮
]

ds = processes_to_coupledodes(eqs, CTMLM)
u0 = [:C => 0.9, :SST => 290, :q_b => 9, :s_b => 291, :z_b => 900]
X, t = trajectory(ds, 100, u0) # time in days

# additional components
@parameters N_d = 250 [description = "cloud droplet number"]
@parameters g_C = 0.9 [description = "backscattering factor"]
@variables tau(t) [description = "cloud optical depth"]'
@variables LWP(t) [description = "liquid water path in g/m2"]
# only have to modify `eqs` vector
eqs = [
    ⋮
    LWP ~ CTMLM.liquid_water_path_linear(),
    tau ~ (0.19*(N_d*10^-6)^1/3*(LWP)^5/6), # aerosol-dependent
    CTMLM.α_C ~ tau/(2/(sqrt(3)*(1 - g_C)) + tau) * CTMLM.C,
    ⋮
]

