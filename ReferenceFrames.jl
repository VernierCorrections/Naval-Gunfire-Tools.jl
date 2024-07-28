module ReferenceFrames


export Body_Rotation
export Body_Rotation_Transform
export Body_Angular_Velocity
export Body_Angular_Acceleration


using AstroTime


function Sidereal_Phase(t, t_iYear, T_syn, r_syn)
    Sidereal_Phase = (t + t_iYear) - (floor((t + t_iYear) / ((T_syn) / (r_syn))) * ((T_syn) / (r_syn)))
    return Sidereal_Phase
end


function Body_Rotation(t, t_iJ2000, t_iYear, T_syn, r_syn, ω_0, ψ_sec, Θ_sec, ψ_ip, ψ_oop, Θ_ip, Θ_oop, ω_nut, χ_nut)
    t_J2000 = t + t_iJ2000
    t_Phase = Sidereal_Phase(t, t_iYear, T_syn, r_syn)
    t_reduced = t_J2000 / ((10^4) * 365.25 * 86400)
    secpol = [1; t_reduced; (t_reduced)^2; (t_reduced)^3; (t_reduced)^4; (t_reduced)^5; (t_reduced)^6; (t_reduced)^7; (t_reduced)^8; (t_reduced)^9; (t_reduced)^10]
    ψ_a = sum((ψ_sec .* secpol), dims = 1)
    Θ_a = sum((Θ_sec .* secpol), dims = 1)
    Ξ_vec = (t_J2000 * ω_nut) + χ_nut
    Δ_ψ = sum(((ψ_ip .* (sin.(Ξ_vec))) .+ (ψ_oop .* (cos.(Ξ_vec)))), dims = 1)
    Δ_Θ = sum(((Θ_ip .* (cos.(Ξ_vec))) .+ (Θ_oop .* (sin.(Ξ_vec)))), dims = 1)
    ψ = ψ_a + Δ_ψ
    Θ = Θ_a + Δ_Θ
    ω_temp = t_Phase * ω_0
    secpol′ = [0; 1; 2*(t_reduced); 3*(t_reduced)^2; 4*(t_reduced)^3; 5*(t_reduced)^4; 6*(t_reduced)^5; 7*(t_reduced)^6; 8*(t_reduced)^7; 9*(t_reduced)^8; 10*(t_reduced)^9]
    ψ̇_a = (sum((ψ_sec .* secpol′), dims = 1)) / ((10^4) * 365.25 * 86400)
    Θ̇_a = (sum((Θ_sec .* secpol′), dims = 1)) / ((10^4) * 365.25 * 86400)
    Δ̇_ψ = sum((((ω_nut .* ψ_ip) .* (cos.(Ξ_vec))) .+ (-(ω_nut .* ψ_oop) .* (sin.(Ξ_vec)))), dims = 1)
    Δ̇_Θ = sum(((-(ω_nut .* Θ_ip) .* (sin.(Ξ_vec))) .+ ((ω_nut .* Θ_oop) .* (cos.(Ξ_vec)))), dims = 1)
    ψ̇ = ψ̇_a + Δ̇_ψ
    Θ̇ = Θ̇_a + Δ̇_Θ
    ω̇_temp = ω_0
    secpol′′ = [0; 0; 2; 6*(t_reduced); 12*(t_reduced)^2; 20*(t_reduced)^3; 30*(t_reduced)^4; 42*(t_reduced)^5; 56*(t_reduced)^6; 72*(t_reduced)^7; 90*(t_reduced)^8]
    ψ̇̇_a = sum((ψ_sec .* secpol′′), dims = 1) / (((10^4) * 365.25 * 86400)^2)
    Θ̇̇_a = sum((Θ_sec .* secpol′′), dims = 1) / (((10^4) * 365.25 * 86400)^2)
    Δ̇̇_ψ = sum(((-((ω_nut .^2) .* ψ_ip) .* sin.(Ξ_vec)) .+ (-((ω_nut .^2) .* ψ_oop) .* cos.(Ξ_vec))), dims = 1)
    Δ̇̇_Θ = sum(((-((ω_nut .^2) .* Θ_ip) .* cos.(Ξ_vec)) .+ (-((ω_nut .^2) .* Θ_oop) .* sin.(Ξ_vec))), dims = 1)
    ψ̇̇ = ψ̇̇_a + Δ̇̇_ψ
    Θ̇̇ = Θ̇̇_a + Δ̇̇_Θ
    return ψ, Θ, ω_temp, ψ̇, Θ̇, ω̇_temp, ψ̇̇, Θ̇̇
end


function Body_Rotation_Transform(ψ, Θ, ω_temp)
precession = [1.0 0.0 0.0; 0.0 cos.(ψ) -sin.(ψ); 0.0 sin.(ψ) cos.(ψ)]
nutation = [cos.(Θ) 0.0 sin.(Θ); 0.0 1.0 0.0; -sin.(Θ) 0.0 cos.(Θ)]
rotation = [cos.(ω_temp) -sin.(ω_temp) 0.0; sin.(ω_temp) cos.(ω_temp) 0.0; 0.0 0.0 1.0]
Rotation_Transform = (rotation) * ((nutation) * (precession))
return Rotation_Transform
end


function Body_Angular_Velocity(ψ, Θ, ψ̇, Θ̇, ω̇_temp)
    ω_x = (ω̇_temp .* sin.(Θ) .* sin.(ψ)) + (Θ̇ .* cos.(ψ))
    ω_y = (-ω̇_temp .* sin.(Θ) .* cos.(ψ)) + (Θ̇ .* sin.(ψ))
    ω_z = (ω̇_temp .* cos.(Θ)) + ψ̇
    ω = [ω_x; ω_y; ω_z];
    return ω
end


function Body_Angular_Acceleration(ψ, Θ, ψ̇, Θ̇, ω̇_temp, ψ̇̇, Θ̇̇ )
    α_x = (Θ̇ .* ω̇_temp .* cos.(Θ) .* sin.(ψ)) + (Θ̇̇ .* cos.(ψ)) + (ψ̇ .* ((ω̇_temp .* sin.(Θ) .* cos.(ψ)) - (Θ̇ .* sin.(ψ))))
    α_y = (Θ̇ .* (-ω̇_temp) .* cos.(Θ) .* cos.(ψ)) + (Θ̇̇ .* sin.(ψ)) + (ψ̇ .* ((ω̇_temp .* sin.(Θ) .* sin.(ψ)) + (Θ̇ .* cos.(ψ))))
    α_z = (Θ̇ .* (-ω̇_temp) .* sin.(Θ)) + (ψ̇̇ )
    α = [α_x; α_y; α_z]
    return α
end


end