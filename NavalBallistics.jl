include("ReferenceFrames.jl")


using Plots
using LinearAlgebra
using DifferentialEquations
using Interpolations
using AstroTime
using .ReferenceFrames


function Saturation_Vapour_Pressure(T, Antoine)
    if Antoine == false
        if T > 647.096
            p_sat = 0
        elseif T > 273.15
            p_sat =  611.21 * exp((18.678 - ((T - 273.15) / 234.5)) * ((T - 273.15) / (-16.01 + T)))
        else
            p_sat = 611.15 * exp((23.036 - ((T - 273.15) / 333.7)) * ((T - 273.15) / (6.67 + T)))
        end
    elseif Antoine == true
        if T > 647.096
            p_sat = 0
        elseif T > 501.486
            p_sat = 133.322387415 * 10^(8.14019 - 1810.94 / (-28.665 + T))
        elseif T > 379
            p_sat = 10^(8.55959 - (643.748 / (-198.043 + T)))
        elseif T > 363
            p_sat = 10^(10.08354 - (1663.125 / (-45.622 + T)))
        elseif T > 334
            p_sat = 10^(10.0768 - (1659.793 / (-45.854 + T)))
        elseif T > 273
            p_sat = 10^(10.40221 - (1838.675 / (-31.737 + T)))
        else
            p_sat = 10^(9.6543 - (1435.264 / (-64.848 + T)))
        end
    else
        p_sat = 0
    end
end


function Atmosphere_Initialization(μ, h_ref, h_array, T_array, Γ_array, molm_ref, p_ref, ρ_ref, atm_size)
    p_array = zeros(atm_size)
    ρ_array = zeros(atm_size)
    p_array[1] = p_ref
    ρ_array[1] = ρ_ref
    for i::UInt8 = 1:(atm_size - 1)
        g_0 = μ * (h_ref)^-2
        if Γ_array[i] == 0.0
            p_array[i + 1] = p_array[i] * exp((-g_0 * molm_ref * (h_array[i + 1] - h_array[i])) / (8.31446261815324 * T_array[i]))
        else
            p_array[i + 1] = p_array[i] * (T_array[i + 1] / T_array[i])^((g_0 * molm_ref) / (8.31446261815324 * Γ_array[i]))
        end
        ρ_array[i + 1] = ρ_array[i] * (p_array[i + 1] / p_array[i]) * (T_array[i] / T_array[i + 1])
        i += 1
    end
    return p_array, ρ_array
end


function Met_Con(r, μ, h_ref, h_array, T_array, Γ_array, p_array, ρ_array, molm_ref, γ_ref)
    h = norm(r) - h_ref
    h_interpolant = extrapolate(interpolate((h_array, ), h_array, Gridded(Constant{Previous}())), Flat())
    h_floor = h_interpolant(h)
    g_0 = μ * (h_ref)^-2
    T_interpolant = extrapolate(interpolate((h_array, ), T_array, Gridded(Linear())), Linear())
    T = T_interpolant(h)
    T_floor = T_interpolant(h_floor)
    Γ_interpolant = extrapolate(interpolate((h_array, ), Γ_array, Gridded(Constant{Previous}())), Flat())
    Γ = Γ_interpolant(h)
    p_interpolant = extrapolate(interpolate((h_array, ), p_array, Gridded(Constant{Previous}())), Flat())
    p_floor = p_interpolant(h) 
    ρ_interpolant = extrapolate(interpolate((h_array, ), ρ_array, Gridded(Constant{Previous}())), Flat())
    ρ_floor = ρ_interpolant(h)
    if Γ == 0.0
        p = p_floor * exp((-g_0 * molm_ref * (h - h_floor)) / (8.31446261815324 * T_floor))
    else
        p = p_floor * (T / T_floor)^((g_0 * molm_ref) / (8.31446261815324 * Γ))
    end
    ρ = ρ_floor * (p / p_floor) * (T_floor / T)
    v_sound = sqrt((γ_ref * 8.31446261815324 * T) / molm_ref)
    return T, p, ρ, v_sound
end


function Fleeman_Drag(t, Mach, q, S_ref, S_tip, S_nozzle, l_proj, l_nose, d_proj, t_burn)
    if Mach < 1
        C_dbasecoast = 0.12 + (0.13 * (Mach^2))
    else
        C_dbasecoast = 0.25 / Mach
    end
    if t <= t_burn
        C_dbase = C_dbasecoast * (1 - (S_nozzle / S_ref))
    else
        C_dbase = C_dbasecoast
    end
    C_dfriction = 0.053 * (l_proj / d_proj) * ((Mach / (q * l_proj))^0.2)
    C_dnose = ((atan(0.5 / (l_nose / d_proj)))^1.69) * (1.586 + (1.834 / (Mach^2)))
    C_dtip = ((π / 4)^1.69) * (1.586 + (1.834 / (Mach^2)))
    C_dwave = (C_dnose * ((S_ref - S_tip) / S_ref)) + (C_dtip * (S_tip / S_ref))
    C_d = C_dbase + C_dfriction + C_dwave
    return C_d
end


function Ballistics_Initial_State(v_muzzle, firing_azimuth, firing_elevation, ϕ_shooter, λ_shooter, a_body, b_body, e_first)
    v_horizontal = v_muzzle * cos(firing_elevation)
    v_ENU = [v_horizontal * sin(firing_azimuth); v_horizontal * cos(firing_azimuth); v_muzzle * sin(firing_elevation)]
    if ϕ_shooter == π / 2
        z_ini = b_body
    elseif ϕ_shooter == -π / 2
        z_ini = -b_body
    elseif ϕ_shooter == 0.0
        z_ini = 0.0
    else
        z_ini = (a_body * (1 - e_first^2) * sin(ϕ_shooter)) / (sqrt(1 - (e_first^2 * ((sin(ϕ_shooter))^2))))
    end
    if ϕ_shooter == π / 2
        p = 0.0
    elseif ϕ_shooter == -π / 2
        p = 0.0
    elseif ϕ_shooter == 0.0
        p = a_body
    else
        p = (a_body * cos(ϕ_shooter)) / (sqrt(1 - (e_first^2 * ((sin(ϕ_shooter))^2))))
    end
    x_ini = p * cos(λ_shooter)
    y_ini = p * sin(λ_shooter)
    v_transform = [-sin(λ_shooter) -sin(ϕ_shooter)*cos(λ_shooter) cos(ϕ_shooter)*cos(λ_shooter); cos(λ_shooter) -sin(ϕ_shooter)*sin(λ_shooter) cos(ϕ_shooter)*sin(λ_shooter); 0.0 cos(ϕ_shooter) sin(ϕ_shooter)]
    v_ECEF = v_transform * v_ENU
    S_ini = [x_ini; y_ini; z_ini; v_ECEF]
    return S_ini
end


function Naval_Ballistics!(∂S∂t, S, Ballistics_Parameters, t)
    μ, t_iJ2000, t_iYear, T_syn, r_syn, ω_0, ψ_sec, Θ_sec, ψ_ip, ψ_oop, Θ_ip, Θ_oop, ω_nut, χ_nut, h_ref, h_array, T_array, Γ_array, p_array, ρ_array, molm_ref, γ_ref, m_proj, S_ref, S_tip, S_nozzle, l_proj, l_nose, d_proj, t_burn = Ballistics_Parameters
    r = S[1:3]
    v = S[4:6]
    a_g = -μ  * ((norm(r))^-3) * r
    ψ, Θ, (), ψ̇, Θ̇, ω̇_temp, ψ̇̇, Θ̇̇ = ReferenceFrames.Body_Rotation(t, t_iJ2000, t_iYear, T_syn, r_syn, ω_0, ψ_sec, Θ_sec, ψ_ip, ψ_oop, Θ_ip, Θ_oop, ω_nut, χ_nut)
    ω = ReferenceFrames.Body_Angular_Velocity(ψ, Θ, ψ̇, Θ̇, ω̇_temp)
    α = ReferenceFrames.Body_Angular_Acceleration(ψ, Θ, ψ̇, Θ̇, ω̇_temp, ψ̇̇, Θ̇̇ )
    a_Coriolis = -2 * cross(ω, v)
    a_centrifugal = -cross(ω, cross(ω, r))
    a_Euler = -cross(α, r)
    (), (), ρ, v_sound = Met_Con(r, μ, h_ref, h_array, T_array, Γ_array, p_array, ρ_array, molm_ref, γ_ref)
    u_mag = norm(v)
    q = 0.5 * ρ * (u_mag^2)
    Mach = u_mag / v_sound
    C_d = Fleeman_Drag(t, Mach, q, S_ref, S_tip, S_nozzle, l_proj, l_nose, d_proj, t_burn)
    a_drag = -((q * C_d * S_ref) / m_proj) * (v / u_mag)
    a = a_g + a_Coriolis + a_centrifugal + a_Euler + a_drag
    ∂S∂t[1:3] = v
    ∂S∂t[4:6] = a
end


function Halting_Condition(r, a_body, b_body, e_first)
    x = r[1]
    y = r[2]
    p = sqrt(x^2 + y^2)
    z = r[3]
    if p == 0.0
        if z == 0.0
            ϕ = 0.0
            λ = 0.0
            geoalt = -b_body
        else
            ϕ = π * sign(z)
            λ = 0.0
            geoalt = abs(z) - b_body
        end
    else
        c = a_body * (e_first)^2
        e′ = sqrt(1.0 - (e_first)^2)
        z′ = e′ * abs(z)
        u = 2.0 * (z′ - c)
        v = 2.0 * (z′ + c)
        t_M = (c - z′) / p
        f_M = (p * (t_M)^4) + (u * (t_M)^3) + (v * t_M) - p
        t_1 = (p - c + z′) / (p - c + 2.0 * z′)
        t_0 = p / (z′ + c)
        if t_M <= 0.0
            t_temp = t_1
        elseif t_M >= 1.0 
            t_temp = t_0
        elseif f_M >= 0.0
            t_temp = t_0
        else
            t_temp = t_1
        end
        while true
            Δ_t_temp = (p - ((p * t_temp^4) + (u * t_temp^3) + (v * t_temp))) / ((4.0 * p * t_temp^3) + (3.0 * u * t_temp^2) + v)
            last_t_temp = t_temp
            t_temp += Δ_t_temp
            TE = abs(t_temp - last_t_temp)
            if TE <= 10^-12
                break
            end
        end
        ϕ = atan((1.0 - t_temp^2), (2.0 * e′ * t_temp)) * sign(z)
        λ = atan(y, x)
        geoalt = ((2.0 * p * e′ * t_temp) + (abs(z) * (1.0 - t_temp^2)) - (a_body * e′ * (1.0 + t_temp^2))) / (sqrt(((1.0 + t_temp^2)^2) - (4.0 * e_first^2 * t_temp^2)))
    end
    N = a_body / (sqrt(1.0 - (e_first^2 * (sin(ϕ))^2)))
    Δ_z = N * e_first^2 * sin(ϕ)
    return ϕ, λ, Δ_z, geoalt
end


function Vincenty(ϕ2, λ2, ϕ1, λ1, a_body, b_body, f_body)
    ϕ2_reduced = atan((1.0 - f_body) * tan(ϕ2))
    ϕ1_reduced = atan((1.0 - f_body) * tan(ϕ1))
    Δ_λ = λ2 - λ1
    Δ_λ_reduced = Δ_λ
    α_cos2 = 0.0
    σ_sin = 0.0
    σ_cos = 0.0
    σ = 0.0
    σ_midpoint = 0.0
    for i = 1:100
        σ_sin = sqrt(((cos(ϕ2_reduced) * sin(Δ_λ_reduced))^2) + (((cos(ϕ1_reduced) * sin(ϕ2_reduced)) - (sin(ϕ1_reduced) * cos(ϕ2_reduced) * cos(Δ_λ_reduced)))^2))
        σ_cos = (sin(ϕ1_reduced) * sin(ϕ2_reduced)) + (cos(ϕ1_reduced) * cos(ϕ2_reduced) * cos(Δ_λ_reduced))
        σ = atan(σ_sin, σ_cos)
        α_sin = (cos(ϕ1_reduced) * cos(ϕ2_reduced) * sin(Δ_λ_reduced)) / (sin(σ))
        α_cos2 = 1.0 - ((α_sin)^2)
        σ_midpoint = σ_cos - ((2.0 * sin(ϕ1_reduced) * sin(ϕ2_reduced)) / (α_cos2))
        Clairaut = (f_body / 16.0) * α_cos2 * (4.0 + (f_body * (4.0 - (3.0 * α_cos2))))
        lastguess = Δ_λ_reduced
        Δ_λ_reduced = Δ_λ + (((1.0 - Clairaut) * f_body * α_sin) * (σ + ((Clairaut * σ_sin) * (σ_midpoint + (Clairaut * σ_cos * (-1.0 + (2.0 * (σ_midpoint^2))))))))
        TE = abs(Δ_λ_reduced - lastguess)
        if TE <= 10^-12
            break
        end
        i += 1
    end
    u2 = α_cos2 * (((a_body^2) - (b_body^2)) / (b_body^2))
    A = 1.0 + ((u2 / 16384.0) * (4096.0 + (u2 * (-768.0 + (u2 * (320.0 - (175.0 * u2)))))))
    B = (u2 / 1024.0) * (256.0 + (u2 * (-128.0 + (u2 * (74.0 - (47.0 * u2))))))
    Δ_σ = (B * σ_sin) * (σ_midpoint + (((1.0 / 4.0) * B) * ((σ_cos * (-1.0 + (2.0 * (σ_midpoint^2)))) - ((1.0 / 6.0) * B * σ_midpoint * (-3.0 + (4.0 * (σ_sin^2))) * (-3.0 + (4.0 * (σ_midpoint^2)))))))
    s_ellipsoid = b_body * A * (σ - Δ_σ)
    return s_ellipsoid
end


function Final_Time_Interpolant(u, Final_Time_Parameters)
    a_body, b_body, e_first, Naval_Ballistics_Integrator = Final_Time_Parameters
    S_out = Naval_Ballistics_Integrator.sol(u)
    (), (), (), geoalt = Halting_Condition(S_out[1:3], a_body, b_body, e_first)
    h_proj = geoalt
    return h_proj
end



function Single_Shot(firing_solution, shooter_coordinates, Ballistics_Parameters, a_body, b_body, e_first, f_body)
    firing_azimuth = firing_solution[1]
    firing_elevation = firing_solution[2]
    ϕ_shooter = shooter_coordinates[1]
    λ_shooter = shooter_coordinates[2]
    Naval_Ballistics_ODE = ODEFunction(Naval_Ballistics!)
    S_ini = Ballistics_Initial_State(v_muzzle, firing_azimuth, firing_elevation, ϕ_shooter, λ_shooter, a_body, b_body, e_first)
    ODE_timespan = (0.0, 10.0)
    Naval_Ballistics_IVP = ODEProblem(Naval_Ballistics_ODE, S_ini, ODE_timespan, Ballistics_Parameters)
    Naval_Ballistics_Integrator = init(Naval_Ballistics_IVP, AutoVern9(Rodas5P()), abstol = 10^-12, reltol = 10^-12)
    ϕ_new, λ_new = ϕ_shooter, λ_shooter
    impact_range = 0.0
    while true
        ϕ_old, λ_old = ϕ_new, λ_new
        step!(Naval_Ballistics_Integrator)
        S_out = Naval_Ballistics_Integrator.sol(Naval_Ballistics_Integrator.t)
        r_out = S_out[1:3]
        v_out = S_out[4:6]
        ϕ_new, λ_new, Δ_z, geoalt = Halting_Condition(r_out, a_body, b_body, e_first)
        impact_range += Vincenty(ϕ_new, λ_new, ϕ_old, λ_old, a_body, b_body, f_body)
        surfnorm = r_out + [0; 0; Δ_z]
        if (geoalt <= 0.0) & (dot(surfnorm, v_out) <= 0.0)
            break
        end
    end
    Final_Time_Parameters = a_body, b_body, e_first, Naval_Ballistics_Integrator
    Final_Time_Span = (Naval_Ballistics_Integrator.tprev, Naval_Ballistics_Integrator.t)
    Final_Time_Problem = IntervalNonlinearProblem(IntervalNonlinearFunction(Final_Time_Interpolant), Final_Time_Span, Final_Time_Parameters)
    Final_Time_Solution = solve(Final_Time_Problem, ITP(), abstol = 10^-12, reltol = 10^-12)
    t_f = Final_Time_Solution.u
    S_f = Naval_Ballistics_Integrator.sol(t_f)
    return impact_range, t_f, S_f
    end


μ::Float64 = 6.67430 * 5.9722 * 10^13
a_body::Float64 = 6378136.3
f_inverse::Float64 = 298.257
f_body::Float64 = f_inverse^-1
b_body::Float64 = a_body * (1 - f_body)
e_first::Float64 = sqrt(1 - ((b_body^2) / (a_body^2)))


ω_0::Float64 = 7.292115 * 10^-5
r_syn::Float64 = 1.00273781191135448
T_syn::Float64 = 86400.002
# Secular terms from Laskar's NGT
ψ_sec = [0.0; ((502909.66 * (1 / 3600)) * (pi / 180)); ((11119.71 * (1 / 3600)) * (pi / 180)); ((77.32 * (1 / 3600)) * (pi / 180)); (((-2353.16) * (1 / 3600)) * (pi / 180)); (((-180.55) * (1 / 3600)) * (pi / 180)); ((174.51 * (1 / 3600)) * (pi / 180)); ((130.95 * (1 / 3600)) * (pi / 180)); ((24.24 * (1 / 3600)) * (pi / 180)); (((-47.59) * (1 / 3600)) * (pi / 180)); ((-8.66 * (1 / 3600)) * (pi / 180))]
Θ_sec = [0.0; (((-4680.93) * (1 / 3600)) * (pi / 180)); (((-1.55) * (1 / 3600)) * (pi / 180)); ((1999.25 * (1 / 3600)) * (pi / 180)); (((-51.38) * (1 / 3600)) * (pi / 180)); (((-249.67) * (1 / 3600)) * (pi / 180)); (((-39.05) * (1 / 3600)) * (pi / 180)); ((7.12 * (1 / 3600)) * (pi / 180)); ((27.87 * (1 / 3600)) * (pi / 180)); ((5.79 * (1 / 3600)) * (pi / 180)); ((2.45 * (1 / 3600)) * (pi / 180))]
# Spectral terms from Dehant & Matthews
ψ_ip = [((-17208) * (10^-3) * (1 / 3600) * (pi / 180)); (207 * (10^-3) * (1 / 3600) * (pi / 180)); (148 * (10^-3) * (1 / 3600) * (pi / 180)); ((-1317) * (10^-3) * (1 / 3600) * (pi / 180)); ((-228) * (10^-3) * (1 / 3600) * (pi / 180))]
ψ_oop = [(3.3 * (10^-3) * (1 / 3600) * (pi / 180)); (0.1 * (10^-3) * (1 / 3600) * (pi / 180)); (1.2 * (10^-3) * (1 / 3600) * (pi / 180)); ((-1.4) * (10^-3) * (1 / 3600) * (pi / 180)); (0.3 * (10^-3) * (1 / 3600) * (pi / 180))]
Θ_ip = [(9205 * (10^-3) * (1 / 3600) * (pi / 180)); ((-90) * (10^-3) * (1 / 3600) * (pi / 180)); (7 * (10^-3) * (1 / 3600) * (pi / 180)); (573 * (10^-3) * (1 / 3600) * (pi / 180)); (98 * (10^-3) * (1 / 3600) * (pi / 180))]
Θ_oop = [(1.5 * (10^-3) * (1 / 3600) * (pi / 180)); (0.0 * (10^-3) * (1 / 3600) * (pi / 180)); ((-0.2) * (10^-3) * (1 / 3600) * (pi / 180)); ((-0.5) * (10^-3) * (1 / 3600) * (pi / 180)); (0.1 * (10^-3) * (1 / 3600) * (pi / 180))]
ω_nut = [((-18.6) * 365.25 * 86400)^(-1); ((-9.3) * 365.25 * 86400)^(-1); (1 * 365.25 * 86400)^(-1); (0.5 * 365.25 * 86400)^(-1); (13.66 * 86400)^(-10)]
χ_nut = [0.0; 0.0; 0.0; 0.0; 0.0]


Date = [2008 2 6 20 45]
CurrentEpoch = Epoch{UniversalTime}(Date[1], Date[2], Date[3], Date[4], Date[5], 0.0)
ReferenceEpoch = UT1Epoch(0.0seconds, origin = :j2000)
YearEpoch = Epoch{UniversalTime}(Date[1], 1, 1, 0, 0, 0.0)
t_iJ2000 = value(seconds(CurrentEpoch - ReferenceEpoch))
t_iYear = value(seconds(CurrentEpoch - YearEpoch))


p_dry::Float64 = 101325.0
molm_dry::Float64 = 0.028946
γ_dry::Float64 = 1.4
hum_rel::Float64 = 0.0
Antoine::Bool = false
T_ref::Float64 = 288.15
p_water::Float64 = hum_rel * Saturation_Vapour_Pressure(T_ref, Antoine)
p_ref::Float64 = p_dry + p_water;
molm_ref::Float64 = (p_dry / p_ref) * molm_dry + (p_water / p_ref) * 0.01801528
γ_ref::Float64 = (p_dry / p_ref) * γ_dry + (p_water / p_ref) * 8/6
ρ_ref::Float64 = (molm_ref * p_ref) / (8.31446261815324 * T_ref)
h_ref::Float64 = sqrt(μ / 9.80665)
h_array = [0.0; 11000.0; 20000.0; 32000.0; 47000.0; 51000.0; 71000.0; 84852.0; 84853.0]
atm_size::UInt8 = size(h_array, 1)
T_array = [T_ref; 216.65; 216.65; 228.65; 270.65; 270.65; 214.65; 186.946; 186.946]
Γ_array = [0.0065; 0.0; -0.001; -0.0028; 0.0; 0.0028; 0.002; 0.0; 0.0]
p_array, ρ_array = Atmosphere_Initialization(μ, h_ref, h_array, T_array, Γ_array, molm_ref, p_ref, ρ_ref, atm_size)


m_pounds::Float64 = 2700.0
m_proj::Float64 = m_pounds * 0.45359237
d_proj::Float64 = 0.406
S_ref::Float64 = π * ((0.5 * d_proj)^2)
l_proj::Float64 = 1.829
l_nose::Float64 = 0.978
d_tip::Float64 = 0.0697
S_tip::Float64 = π * ((0.5 * d_tip)^2)
d_nozzle::Float64 = 0.0
S_nozzle::Float64 = π * ((0.5 * d_nozzle)^2)
t_burn::Float64 = 0.0
v_fps::Float64 = 2500
v_muzzle::Float64 = v_fps * 0.3048


Ballistics_Parameters = μ, t_iJ2000, t_iYear, T_syn, r_syn, ω_0, ψ_sec, Θ_sec, ψ_ip, ψ_oop, Θ_ip, Θ_oop, ω_nut, χ_nut, h_ref, h_array, T_array, Γ_array, p_array, ρ_array, molm_ref, γ_ref, m_proj, S_ref, S_tip, S_nozzle, l_proj, l_nose, d_proj, t_burn
firing_azimuth = 90.0 * π / 180
firing_elevation = -2.0 * π / 180
ϕ_shooter::Float64 = 30.0 * π / 180
λ_shooter::Float64 = 30.0 * π / 180
firing_solution = [firing_azimuth; firing_elevation]
shooter_coordinates = [ϕ_shooter; λ_shooter]


test_size = 53
sampling_rate = 1.0 * π / 180
elevation_matrix = zeros(test_size, 1)
range_matrix = zeros(test_size, 1)
time_matrix = zeros(test_size, 1)
velocity_matrix = zeros(test_size, 1)
deck_obliquity_matrix = zeros(test_size, 1)
for i = 1:test_size
    impact_range, time_matrix[i], S_f = Single_Shot(firing_solution, shooter_coordinates, Ballistics_Parameters, a_body, b_body, e_first, f_body)
    elevation_matrix[i] = firing_solution[2]
    range_matrix[i] = impact_range / 1852.0
    r_f = S_f[1:3]
    v_f = S_f[4:6]
    velocity_matrix[i] = norm(v_f) / 0.3048
    (), (), Δ_z, () = Halting_Condition(r_f, a_body, b_body, e_first)
    decknorm = r_f + [zeros(2, 1); Δ_z]
    deck_obliquity_matrix[i] = acos(dot((-v_f ./ (sqrt.(sum(abs2.(v_f), dims = 1)))), (decknorm ./ (sqrt.(sum(abs2.(decknorm), dims = 1)))))) * (180.0 / π)
    firing_solution += [0.0; sampling_rate]
    i += 1
end
elevation_matrix *= 180 / π
belt_obliquity_matrix = 90.0 .- deck_obliquity_matrix
range_plot = plot(range_matrix, elevation_matrix, title = "Firing Angle vs Range", label = "Simulated Results", xlabel = "Range (nautical miles)", ylabel = "Elevation (degrees)")
time_plot = plot(range_matrix, time_matrix, title = "Time of Flight vs Range", label = "Simulated Results", xlabel = "Range (nautical miles)", ylabel = "Time of Flight (seconds)")
velocity_plot = plot(range_matrix, velocity_matrix, title = "Striking Velocity vs Range", label = "Simulated Results", ylims = (0.0, v_fps), xlabel = "Range (nautical miles)", ylabel = "Striking Velocity (feet per second)")
obliquity_plot = plot(range_matrix, [belt_obliquity_matrix deck_obliquity_matrix], title = "Obliquity vs Range", label = ["Simulation Against Belt" "Simulation Against Deck"], ylims = (0.0, 90.0), xlabel = "Range (nautical miles)", ylabel = "Obliquity (degrees)")
elevation_real = [2.36; 5.05; 8.16; 11.77; 16.03; 21.11; 27.37; 36.08; 45.00]
range_real = [4572.0/1852.0; 9140.0/1852.0; 13716.0/1852.0; 18290.0/1852.0; 22860.0/1852.0; 27430.0/1852.0; 32000.0/1852.0; 36580.0/1852.0; 38720.0/1852.0]
time_real = [6.29; 13.24; 20.98; 29.59; 39.30; 50.32; 63.22; 79.96; 95.32]
velocity_real = [2279.0; 2074.0; 1892.0; 1740.0; 1632.0; 1567.0; 1556.0; 1607.0; 1686.0]
belt_obliquity_real = [2.50; 5.01; 9.78; 14.92; 21.12; 28.25; 36.27; 47.73; 51.23]
deck_obliquity_real = 90.0 .- belt_obliquity_real
plot!(range_plot, range_real, elevation_real, seriestype = :scatter, markershape = :star4, label = "US AP Mk 8 Shell (16 inch/50 calibre Mk 7)")
plot!(time_plot, range_real, time_real, seriestype = :scatter, markershape = :star4, label = "US AP Mk 8 Shell (16 inch/50 calibre Mk 7)")
plot!(velocity_plot, range_real, velocity_real, seriestype = :scatter, markershape = :star4, label = "US AP Mk 8 Shell (16 inch/50 calibre Mk 7)")
plot!(obliquity_plot, range_real, [belt_obliquity_real deck_obliquity_real], seriestype = :scatter, markershape = :star4, label = ["US AP Mk 8 Against Belt" "US AP Mk 8 Against Deck"])
plot!(obliquity_plot, legend = :outerbottom, legendcolumns = 2)
savefig(range_plot, "Range_Plot.png")
savefig(time_plot, "Time_Plot.png")
savefig(velocity_plot, "Velocity_Plot.png")
savefig(obliquity_plot, "Obliquity_Plot.png")

