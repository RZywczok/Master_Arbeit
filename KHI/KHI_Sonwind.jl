using OrdinaryDiffEq, Trixi, Trixi2Vtk, HOHQMesh, Plots

tspan = (0.0,5.)

coordinates_min = (-1.0, -1.0)
coordinates_max = (1.0, 1.0)

name_of_file = "KHI_Merkur/"

workin_dir = "Master/Code/"

path_to_save = workin_dir * name_of_file


ν  = 0 #1.26*10^(-6) / Rescale_Const[8] # WICHTIG DAS IST NICHT MU SONDER kinematische Viskosität ν 
eta() = 0#0.0001
prandtl_number() = 0.72

equations_hyperbolic = IdealGlmMhdEquations2D(1.8)

μ_0 = 1.256637 * big(10)^(-6)                        # Magnetische Feldkonstante
R_s = 8.314462618                                    # Universelle Gaskonstante
k_B = 1.380649 * big(10)^(-23)                       # Boltzmann konstante
ato_u = 1.6605402 * big(10)^(-27)                   # Atommasse u in kg 
m_SW = ato_u * 0.9 + ato_u * 4.002602 * 0.1

R_M_D = 4879.9/2 * 10^3

struct Sunwind
    B::Vector{Float64}
    T::Float64
    n::Int64
    ρ::Float64
    v::Vector{Float64}
    β::Float64
    p::Float64
    c_A::Float64
    Ma::Float64
end

R_MP(R_M, g01, rho1, u1) = Float64(((2*(g01 * 10^(-9))^2) / ( μ_0 * rho1 * u1^2))^(1/6), RoundDown)
R_BS(R_MP) = R_MP / 0.8

function initial_condition_defintion(B,T,v,n)
    l_R  = 1000/2 * 10^3 
    v_R  = sqrt(sum(v .^2)) 
    ρ = (ato_u * 0.9 + ato_u * 4.002602 * 0.1) * n        # Dichte als kg/(m³) Berechnung: (n Teilchen)/cm³ * 10^6 [m³/cm³] * (u[als Atommasse] * 90% + Helium * 10%)  <- Massenzusammensetung des Sonnenwindes
    ρ_R  = ρ 
    B_R  = sqrt(v_R^2 * ρ_R * μ_0)

    t_R  = l_R/v_R
    p_R  = ρ_R* (v_R ^2)
    μ0_R = B_R^2 /(v_R ^2 * ρ_R)
    μR_R = (B_R^2 * l_R) /(v_R * ρ_R)
    μNS_R = v_R * l_R
    V_R  = B_R * v_R * ρ_R
    E_R  = v_R^2
    
    β = (n * k_B * T)/(sum(B.^2)/(2*μ_0))                 # Plasma Betha
    p =  n * k_B * T                                      # Druck im Sonnenwind nach idealer Gasgleichung

    c_A = sqrt(sum(B.^2))/sqrt(μ_0 * ρ)                   # Alvén-Geschindigkeit
    Ma = sqrt(sum(v.^2)) / c_A                            # Alvénische Machzahl   
    return Sunwind(Float64.(B./B_R, RoundDown), Float64(T, RoundDown), Float64(n, RoundDown), Float64(ρ / ρ_R, RoundDown), Float64.(v./v_R, RoundDown), Float64(β, RoundDown), Float64(p/p_R, RoundDown), Float64(c_A, RoundDown), Float64(Ma, RoundDown)), [Float64(B_R, RoundDown),Float64(v_R, RoundDown), Float64(ρ_R, RoundDown),Float64(p_R, RoundDown), Float64(t_R, RoundDown), Float64(μ0_R, RoundDown), Float64(μR_R, RoundDown), Float64(μNS_R, RoundDown)]
end

function Def_B(x)
    slope       = 50
    amplitude  = 0.03
    B          = tanh.(slope .* x[2] .+ 15) .- tanh.(slope .* x[2] .- 15) 
    rho        = initial_sunwind.ρ + initial_sunwind.ρ * 1.5 * B
    V1,V2,V3   = initial_sunwind.v
    v1         = V1 * 0.5 * (B - 1)
    v2         = V1 * amplitude * sin(2 * pi * x[1])
    v3         = V3
    p          = (rho/m_SW) * k_B * T 
    B1, B2, B3 = initial_sunwind.B
    psi        = 0
    return prim2cons(SVector(rho, v1, v2, v3, p, B1, B2, B3, psi), equations_hyperbolic)
end

function Def_Boundery(x)
    slope       = 50
    amplitude  = 0.03
    B          = tanh.(slope .* x[2] .+ 15) .- tanh.(slope .* x[2] .- 15) 
    rho        = initial_sunwind.ρ + initial_sunwind.ρ * 1.5 * B
    V1,V2,V3   = initial_sunwind.v
    v1         = V1 * 0.5 * (B - 1)
    v2         = V2
    v3         = V3
    p          = (rho/m_SW) * k_B * T 
    B1, B2, B3 = initial_sunwind.B
    psi        = 0
    return prim2cons(SVector(rho, v1, v2, v3, p, B1, B2, B3, psi), equations_hyperbolic)
end

B1 = [26.61, 14.86, 10.54] .* 10^(-9)                  # Tesla
T = 15 * big(10)^(4)                                   # Kelvin
v = [423,0,0]  * 10 ^ 3                                # Sonnenwind Geschwindigkeit
n = 50* big(10)^(6)                                    # Teilchenzahl [1/m³]
initial_sunwind, Rescale_Const = initial_condition_defintion(B1,T,v,n)

k = 1 # the wave number of the perturbation
rescaling = Float64.([Rescale_Const[3],Rescale_Const[2],Rescale_Const[2],Rescale_Const[2], Rescale_Const[4],Rescale_Const[1],Rescale_Const[1],Rescale_Const[1],1], RoundDown)


println("Massendichte: ", round(initial_sunwind.ρ; digits=4)," | Geschwindigkeit: ",  round.(initial_sunwind.v; digits=4)," | Druck: ", round(initial_sunwind.p; digits=7)," | Magnetfeld: ", round.(initial_sunwind.B; digits=3)," | B_R: ",  Rescale_Const[1]," | v_R: ",  Rescale_Const[2]," | ρ_R: ",  Rescale_Const[3]," | p_R: ",  Rescale_Const[4], " | Zeit: ",  round(Rescale_Const[5]; digits=3)," | μ_0: ", round(μ_0/Rescale_Const[6])," | μ_R: ", round(Rescale_Const[7])," | μ_NS: ", Rescale_Const[8])

k = 1 # the wave number of the perturbation
rescaling = [Rescale_Const[3],Rescale_Const[2],Rescale_Const[2],Rescale_Const[2], Rescale_Const[4],Rescale_Const[1],Rescale_Const[1],Rescale_Const[1],1]

B_inside = cons2prim(Def_B([0,0]), equations_hyperbolic) .* rescaling
B_outside = cons2prim(Def_B([0,1]), equations_hyperbolic).* rescaling

N_fluid = ((B_inside[1]/ m_SW) + (B_outside[1]/ m_SW)) / ((B_inside[1]/ m_SW) * (B_outside[1]/ m_SW)*μ_0*m_SW);

if k * (B_outside[2] - B_inside[2] )^2 > ((k * B_inside[6])^2 + (k * B_outside[6])^2)* N_fluid
    println("Es gibt eine KHI da: ", k * (B_inside[2] - B_outside[2])^2  , " > ", ((k * B_inside[6])^2 + (k * B_outside[6])^2)* N_fluid)
else 
    println("Es gibt KEINE KHI da: ", k * (B_inside[2] - B_outside[2])^2 , " < ", ((k * B_inside[6])^2 + (k * B_outside[6])^2)* N_fluid)
end

c_A = sqrt(sum(initial_sunwind.B.^2))/sqrt(1 * initial_sunwind.ρ)                   # Alvén-Geschindigkeit
Ma = sqrt(sum(initial_sunwind.v.^2)) / c_A    

@inline function boundary_condition_outflow(u_inner, normal_direction, x, t,
                                              surface_flux_function,
                                              equations::IdealGlmMhdEquations2D)
        flux = Trixi.flux(u_inner, normal_direction, equations)
    return flux
end
function initial_condition_kelvin_helmholtz_instability(x, t, equations::IdealGlmMhdEquations2D)
    return Def_B(x)
end
initial_condition = initial_condition_kelvin_helmholtz_instability

function boundary_condition_Wall(x, t, equations::IdealGlmMhdEquations2D)
    return Def_Boundery(x)
end

boundary_conditions_hyperbolic =   Dict(#:x_neg => boundary_condition_do_nothing,
                                        :y_neg => BoundaryConditionDirichlet(boundary_condition_Wall),
                                        :y_pos => BoundaryConditionDirichlet(boundary_condition_Wall))#,
                                        #:x_pos => boundary_condition_do_nothing)

boundary_conditions_parabolic =   Dict(#:x_neg => boundary_condition_do_nothing,
                                       :y_neg => boundary_condition_do_nothing,
                                       :y_pos => boundary_condition_do_nothing)#,
                                       #:x_pos => boundary_condition_do_nothing)


trees_per_dimension = (19, 21)
mesh = P4estMesh(trees_per_dimension, polydeg = 3,
                 coordinates_min = coordinates_min, coordinates_max = coordinates_max,
                 periodicity = (true,true))

surface_flux = (flux_lax_friedrichs, flux_nonconservative_powell)
volume_flux  = (flux_central, flux_nonconservative_powell)

polydeg = 4
basis = LobattoLegendreBasis(polydeg)
indicator_sc = IndicatorHennemannGassner(equations_hyperbolic, basis,
                                         alpha_max=0.002,
                                         alpha_min=0.0001,
                                         alpha_smooth=true,
                                         variable=Trixi.density)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg=volume_flux,
                                                 volume_flux_fv=surface_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

initial_condition = initial_condition_kelvin_helmholtz_instability

semi = SemidiscretizationHyperbolicParabolic(mesh,
                                             (equations_hyperbolic, equations_parabolic),
                                             initial_condition, solver;
                                            # boundary_conditions=(boundary_conditions_hyperbolic,boundary_conditions_parabolic)
)


###############################################################################
# ODE solvers, callbacks etc.

ode = semidiscretize(semi, tspan);

summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

alive_callback = AliveCallback(analysis_interval=analysis_interval)

save_solution = SaveSolutionCallback(interval=100,
                                     save_initial_solution=true,
                                     output_directory= path_to_save,
                                     save_final_solution=true,
                                     solution_variables = cons2prim)

amr_indicator = IndicatorLöhner(semi,variable=density_pressure)


amr_controller = ControllerThreeLevel(semi, amr_indicator;
                                      base_level=1,
                                      med_level=5, med_threshold=0.01,
                                      max_level=6, max_threshold=0.04)

amr_callback = AMRCallback(semi, amr_controller,
                            interval=1,
                            adapt_initial_condition=true,
                            adapt_initial_condition_only_refine=true)

cfl = 0.4

stepsize_callback = StepsizeCallback(cfl=cfl)

glm_speed_callback = GlmSpeedCallback(glm_scale=0.3, cfl=cfl)

callbacks = CallbackSet(summary_callback,
                        #analysis_callback,
                        alive_callback,
                        save_solution,
                        amr_callback,
                        stepsize_callback,
                        glm_speed_callback
      );

stage_limiter! = PositivityPreservingLimiterZhangShu(thresholds=(5.0e-5, 1.0e-5),
                                                      variables=(pressure, Trixi.density))

try
    rm(workin_dir * name_of_file, recursive=true)
catch
    @warn "Could not read file."
end


sol = solve(ode, CarpenterKennedy2N54(stage_limiter!, williamson_condition=false),
            dt=0.001, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep=false, callback=callbacks)
summary_callback()

plot(sol)

try
    for i in readdir(workin_dir * name_of_file * "Restart_MHD_Files/")
        cp(workin_dir * name_of_file * "Restart_MHD_Files/" * i, path_to_save * i)#;force=true)
    end
    cp(workin_dir * name_of_file * "Restart_MHD_Files/p4est_data" , path_to_save * "p4est_data")#;force=true)
catch
    println("No File")
end

try
    rm(workin_dir * "vtk_" * name_of_file, recursive=true)
catch
    @warn "Could not read file."
end

saved_solutions = filter(x->contains(x,".h5"), readdir(path_to_save))

saved_solutions_path = joinpath.(workin_dir,name_of_file .* saved_solutions)

for i in saved_solutions_path
    trixi2vtk(i, output_directory= workin_dir * "vtk_" * name_of_file)
end


