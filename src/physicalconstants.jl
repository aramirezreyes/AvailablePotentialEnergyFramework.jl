const R             = 287 # Ideal gas constants
const epsilon       = 29/18-1
const g             = 10 #acceleration of gravity

@with_kw struct Substance
    cp = nothing
    cv = nothing
    R = nothing
    Lv = nothing
    Lf = nothing
end


Dryair = Substance(
    cp = 1006 ,#J/kg/k at 1013 hPa
    cv = 718,
    R  = 287.05 # J/kg/k
)

Liquidwater = Substance(
     Lv = 2.5e6, #J/kg
     Lf = 3.33e5,
     cp = 4200 #j/kg/k
)

Watervapor = Substance(
    R = 461.52 #j/kg/K
)
