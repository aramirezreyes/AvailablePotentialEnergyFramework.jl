const epsilon       = 18.016/28.966
const g             = 10u"m/s/s" #acceleration of gravity

@with_kw struct Substance
    cp = nothing
    cv = nothing
    R = nothing
    Lv = nothing
    Lf = nothing
end


Dryair = Substance(
    cp = 1006u"J/kg/K" ,#J/kg/k at 1013 hPa
    cv = 718u"J/kg/K",
    R  = 287.05u"J/kg/K" # J/kg/k
)

Liquidwater = Substance(
     Lv = 2.5e6u"J/kg", #J/kg
     Lf = 3.33e5u"J/kg",
     cp = 4200u"J/kg/K" #j/kg/k
)

Watervapor = Substance(
    R = 461.52u"J/kg/K" #j/kg/K
)
