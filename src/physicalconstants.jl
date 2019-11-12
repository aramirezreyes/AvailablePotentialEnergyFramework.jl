const R             = 287 # Ideal gas constants
const heat_capacity = 1004; #Heat capacity of air
const L             = 2.5*1e6 #Enthalpy of phase change
const epsilon       = 29/18-1
const g             = 10 #acceleration of gravity

@with_kw struct substance
    cp = nothing
    cv = nothing
    R = nothing
    Lv = nothing
    Lf = nothing
end


dryair = substance(
    cp = 1006 ,#J/kg/k at 1013 hPa
    cv = 718,
    R  = 287.05 # J/kg/k
)

liquidwater  = substance(
     Lv = 2.5e6, #J/kg
     Lf = 3.33e5,
     cp = 4200 #j/kg/k
)


