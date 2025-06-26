using Pinching4Model, CairoMakie 


#Consider this OpenSees example:
#https://opensees.berkeley.edu/wiki/index.php/Pinching4MaterialExample

# set peakpts [list 0.0001 0.001 0.002 0.003 0.005 0.006 0.007 0.009 0.01 0.011 0.012 0.013 ]


half_cycles = [[0.0, 0.0001], 
           [0.0001, -0.0001], 
           [-0.0001, 0.001],
           [0.001, -0.001], 
            [-0.001, 0.002],
            [0.002, -0.003],
            [-0.003, 0.003],
            [0.003, -0.005],
            [-0.005, 0.006],
             [0.006, -0.006],
              [-0.006, 0.007],
              [0.007, -0.007],
              [-0.007, 0.009],
              [0.009, -0.009],
              [-0.009, 0.010],
              [0.010, -0.010],
               [-0.010, 0.011],
               [0.011, -0.011],
               [-0.011, 0.012],
               [0.012, -0.012],
               [-0.012, 0.013],
               [0.013, 0.0]]

incr = 300

d = []

for i in eachindex(half_cycles)
        
    if i != size(half_cycles)[1]
        segment_range = collect(range(half_cycles[i][1], half_cycles[i][2], incr))[1:end-1]
    else
        segment_range = collect(range(half_cycles[i][1], half_cycles[i][2], incr))
    end

    d = [d; segment_range]

end

scatterlines(1:length(d), d)



### no damage plot
dmgtype = "energy"

strain1p = 0.0001
strain2p = 0.0055
strain3p = 0.0188
strain4p = 0.0189
strain1n = -0.0001
strain2n = -0.0055
strain3n = -0.0188
strain4n = -0.0189

stress1p = 2.0
stress2p = 6.0
stress3p = 7.0
stress4p = 0.2
stress1n = -2.0
stress2n = -6.0
stress3n = -7.0
stress4n = -0.2

rDispP = 0.5
rForceP = 0.25
uForceP = 0.05
rDispN = 0.5
rForceN = 0.25
uForceN = 0.05


gammaK1 = 1.0
gammaK2 = 0.2
gammaK3 = 0.3
gammaK4 = 0.2
gammaKLimit = 0.9

gammaK1 = 0.0
gammaK2 = 0.0
gammaK3 = 0.0
gammaK4 = 0.0
gammaKLimit = 0.0


gammaD1 = 0.5
gammaD2 = 0.5
gammaD3 = 2.0
gammaD4 = 2.0
gammaDLimit = 0.5

gammaD1 = 0.0
gammaD2 = 0.0
gammaD3 = 0.0
gammaD4 = 0.0
gammaDLimit = 0.0

gammaF1 = 1.0
gammaF2 = 0.0
gammaF3 = 1.0
gammaF4 = 1.0
gammaFLimit = 0.9

gammaF1 = 0.0
gammaF2 = 0.0
gammaF3 = 0.0
gammaF4 = 0.0
gammaFLimit = 0.0


gE = 10



inputs = Pinching4Model.Inputs(
    d,
    dmgtype,


    strain1p,
    strain2p,
    strain3p,
    strain4p,
    strain1n,
    strain2n,
    strain3n,
    strain4n,

    stress1p,
    stress2p,
    stress3p,
    stress4p,
    stress1n,
    stress2n,
    stress3n,
    stress4n,


    rDispP,
    rForceP,
    uForceP,
    rDispN,
    rForceN,
    uForceN,


    gammaK1,
    gammaK2,
    gammaK3,
    gammaK4,
    gammaKLimit,

    gammaD1,
    gammaD2,
    gammaD3,
    gammaD4,
    gammaDLimit,

    gammaF1,
    gammaF2,
    gammaF3,
    gammaF4,
    gammaFLimit,

    gE
)



force = Pinching4Model.response(inputs)


f = Figure()
ax = Axis(f[1, 1])


lines!(ax, d, force, color=:red)
ylims!(-8, 8)
xlims!(-0.015, 0.015)
f


###




### damage plot
dmgtype = "energy"

strain1p = 0.0001
strain2p = 0.0055
strain3p = 0.0188
strain4p = 0.0189
strain1n = -0.0001
strain2n = -0.0055
strain3n = -0.0188
strain4n = -0.0189

stress1p = 2.0
stress2p = 6.0
stress3p = 7.0
stress4p = 0.2
stress1n = -2.0
stress2n = -6.0
stress3n = -7.0
stress4n = -0.2

rDispP = 0.5
rForceP = 0.25
uForceP = 0.05
rDispN = 0.5
rForceN = 0.25
uForceN = 0.05


gammaK1 = 1.0
gammaK2 = 0.2
gammaK3 = 0.3
gammaK4 = 0.2
gammaKLimit = 0.9

# gammaK1 = 0.0
# gammaK2 = 0.0
# gammaK3 = 0.0
# gammaK4 = 0.0
# gammaKLimit = 0.0


gammaD1 = 0.5
gammaD2 = 0.5
gammaD3 = 2.0
gammaD4 = 2.0
gammaDLimit = 0.5

# gammaD1 = 0.0
# gammaD2 = 0.0
# gammaD3 = 0.0
# gammaD4 = 0.0
# gammaDLimit = 0.0

gammaF1 = 1.0
gammaF2 = 0.0
gammaF3 = 1.0
gammaF4 = 1.0
gammaFLimit = 0.9

# gammaF1 = 0.0
# gammaF2 = 0.0
# gammaF3 = 0.0
# gammaF4 = 0.0
# gammaFLimit = 0.0


gE = 10



inputs = Pinching4Model.Inputs(
    d,
    dmgtype,


    strain1p,
    strain2p,
    strain3p,
    strain4p,
    strain1n,
    strain2n,
    strain3n,
    strain4n,

    stress1p,
    stress2p,
    stress3p,
    stress4p,
    stress1n,
    stress2n,
    stress3n,
    stress4n,


    rDispP,
    rForceP,
    uForceP,
    rDispN,
    rForceN,
    uForceN,


    gammaK1,
    gammaK2,
    gammaK3,
    gammaK4,
    gammaKLimit,

    gammaD1,
    gammaD2,
    gammaD3,
    gammaD4,
    gammaDLimit,

    gammaF1,
    gammaF2,
    gammaF3,
    gammaF4,
    gammaFLimit,

    gE
)



force_damage = Pinching4Model.response(inputs)





lines!(ax, d, force_damage, color=:blue)
ylims!(-8, 8)
xlims!(-0.015, 0.015)
f
