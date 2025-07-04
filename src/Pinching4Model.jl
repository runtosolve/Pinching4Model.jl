module Pinching4Model


mutable struct Inputs
    d::Vector{Real}         # Displacement/strain history
    dmgtype::String            # "energy" or "cycle"

    # Backbone envelope points (strain/stress)
    strain1p::Real
    strain2p::Real
    strain3p::Real
    strain4p::Real
    strain1n::Real
    strain2n::Real
    strain3n::Real
    strain4n::Real

    stress1p::Real
    stress2p::Real
    stress3p::Real
    stress4p::Real
    stress1n::Real
    stress2n::Real
    stress3n::Real
    stress4n::Real

    # Degradation parameters (pinching & damage)
    rDispP::Real
    rForceP::Real
    uForceP::Real
    rDispN::Real
    rForceN::Real
    uForceN::Real

    # Damage model coefficients
    gammaK1::Real; gammaK2::Real; gammaK3::Real; gammaK4::Real; gammaKLimit::Real
    gammaD1::Real; gammaD2::Real; gammaD3::Real; gammaD4::Real; gammaDLimit::Real
    gammaF1::Real; gammaF2::Real; gammaF3::Real; gammaF4::Real; gammaFLimit::Real

    # gK::Vector{Real}        # Stiffness degradation (gK1 to gK4)
    # gKLim::Real
    # gD::Vector{Real}        # Deformation degradation (gD1 to gD4)
    # gDLim::Real
    # gF::Vector{Real}        # Strength degradation (gF1 to gF4)
    # gFLim::Real
    gE::Real                # Energy-based degradation multiplier

    
    # gammaK1::Real; gammaK2::Real; gammaK3::Real; gammaK4::Real; gammaKLimit::Real
    # gammaD1::Real; gammaD2::Real; gammaD3::Real; gammaD4::Real; gammaDLimit::Real
    # gammaF1::Real; gammaF2::Real; gammaF3::Real; gammaF4::Real; gammaFLimit::Real


end







function response(MDL::Inputs)
    # Initialize arrays and variables
    envlpPosStrain = zeros(Real, 6)
    envlpPosStress = zeros(Real, 6)
    envlpNegStrain = zeros(Real, 6)
    envlpNegStress = zeros(Real, 6)

    state3Stress = zeros(Real, 6)
    state3Strain = zeros(Real, 6)
    state4Stress = zeros(Real, 6)
    state4Strain = zeros(Real, 6)
    envlpPosDamgdStress = zeros(Real, 6)
    envlpNegDamgdStress = zeros(Real, 6)

    # Damage type flag
    DmgCyc = MDL.dmgtype == "energy" ? 0 : 1

    TnCycle = 0.0
    CnCycle = 0.0

    Tstress = 0.0
    Tstrain = 0.0
    Ttangent = 0.0

    Cstate = 0
    Cstrain = 0.0
    Cstress = 0.0
    CstrainRate = 0.0
    lowCstateStrain = 0.0
    lowCstateStress = 0.0
    hghCstateStrain = 0.0
    hghCstateStress = 0.0
    CminStrainDmnd = 0.0
    CmaxStrainDmnd = 0.0
    Cenergy = 0.0
    CgammaK = 0.0
    CgammaD = 0.0
    CgammaF = 0.0
    gammaKUsed = 0.0
    gammaFUsed = 0.0

    Tstate = 0
    dstrain = 0.0
    TstrainRate = 0.0
    lowTstateStrain = 0.0
    lowTstateStress = 0.0
    hghTstateStrain = 0.0
    hghTstateStress = 0.0
    TminStrainDmnd = 0.0
    TmaxStrainDmnd = 0.0
    Tenergy = 0.0
    TgammaK = 0.0
    TgammaD = 0.0
    TgammaF = 0.0

    kElasticPos = 0.0
    kElasticNeg = 0.0
    kElasticPosDamgd = 0.0
    kElasticNegDamgd = 0.0
    uMaxDamgd = 0.0
    uMinDamgd = 0.0

    energyCapacity = 0.0
    kunload = 0.0
    elasticStrainEnergy = 0.0

    force = zeros(Real, length(MDL.d))

    # Envelope calculation
    envlpPosStrain, envlpPosStress, envlpNegStrain, envlpNegStress, kElasticPos, kElasticNeg, energyCapacity =
        SetEnvelop(MDL)

    for i in 1:length(MDL.d)
        if i == 1
            Cstate, Cstrain, Cstress, CstrainRate, lowCstateStrain, lowCstateStress, hghCstateStrain, hghCstateStress,
            CminStrainDmnd, CmaxStrainDmnd, Cenergy, CgammaK, CgammaD, CgammaF, CnCycle,
            Ttangent, dstrain, gammaKUsed, gammaFUsed, kElasticPosDamgd, kElasticNegDamgd,
            uMaxDamgd, uMinDamgd =
                revertToStart(envlpPosStrain, envlpPosStress, envlpNegStrain, envlpNegStress, kElasticPos, kElasticNeg)
        else
            Tstate, TstrainRate, lowTstateStrain, lowTstateStress, hghTstateStrain, hghTstateStress,
            TminStrainDmnd, TmaxStrainDmnd, Tenergy, Tstrain, Tstress, TgammaD, TgammaK, TgammaF, TnCycle =
                revertToLastCommit(Cstate, CstrainRate, lowCstateStrain, lowCstateStress, hghCstateStrain,
                                   hghCstateStress, CminStrainDmnd, CmaxStrainDmnd, Cenergy, Cstrain,
                                   Cstress, CgammaD, CgammaK, CgammaF, CnCycle)
        end

        Tstate, Tenergy, Tstrain, lowTstateStrain, hghTstateStrain, lowTstateStress, hghTstateStress,
        TgammaF, TgammaK, TgammaD, dstrain, Ttangent, Tstress, state3Strain, state3Stress,
        state4Strain, state4Stress, kunload, elasticStrainEnergy, TminStrainDmnd, TmaxStrainDmnd, gammaKUsed, gammaFUsed,
        kElasticPosDamgd, kElasticNegDamgd, envlpPosDamgdStress, envlpNegDamgdStress, TnCycle =
            setTrialStrain(MDL.d[i], CstrainRate, Cstate, Cenergy, lowCstateStrain, hghCstateStrain, lowCstateStress,
                           hghCstateStress, CminStrainDmnd, CmaxStrainDmnd, CgammaF, CgammaK, CgammaD,
                           envlpPosStress, envlpPosStrain, kElasticPosDamgd, kElasticNegDamgd,
                           state3Strain, state3Stress, kunload, state4Strain, state4Stress,
                           Cstrain, uMaxDamgd, uMinDamgd, envlpNegStrain, envlpNegStress, kElasticPos, kElasticNeg,
                           Cstress, DmgCyc, CnCycle, energyCapacity, MDL, Tstate, Tenergy, Tstrain,
                           lowTstateStrain, hghTstateStrain, lowTstateStress, hghTstateStress, TgammaF, TgammaK,
                           TgammaD, dstrain, Ttangent, Tstress, elasticStrainEnergy,
                           TminStrainDmnd, TmaxStrainDmnd, gammaKUsed, gammaFUsed,
                           envlpPosDamgdStress, envlpNegDamgdStress, TnCycle)

        Cstate, CstrainRate, lowCstateStrain, lowCstateStress, hghCstateStrain, hghCstateStress,
        CminStrainDmnd, CmaxStrainDmnd, Cenergy, Cstress, Cstrain, CgammaK, CgammaD, CgammaF,
        kElasticPosDamgd, kElasticNegDamgd, uMaxDamgd, uMinDamgd,
        envlpPosDamgdStress, envlpNegDamgdStress, CnCycle =
            commitState(Tstate, dstrain, TstrainRate, lowTstateStrain, lowTstateStress, hghTstateStrain,
                        hghTstateStress, TminStrainDmnd, TmaxStrainDmnd, Tenergy, Tstress, Tstrain,
                        TgammaK, TgammaD, TgammaF, kElasticPos, kElasticNeg, gammaKUsed, gammaFUsed,
                        envlpPosStress, envlpNegStress, TnCycle)

        force[i] = Cstress
    end

    return force
end



function SetEnvelop(MDL::Inputs)
    # Initialize outputs
    envlpPosStrain = zeros(Real, 6)
    envlpPosStress = zeros(Real, 6)
    envlpNegStrain = zeros(Real, 6)
    envlpNegStress = zeros(Real, 6)

    # Initial tangent stiffnesses
    kPos = MDL.stress1p / MDL.strain1p
    kNeg = MDL.stress1n / MDL.strain1n
    k = max(kPos, kNeg)

    if MDL.strain1p > -1.0 * MDL.strain1n
        u = 1.0e-4 * MDL.strain1p
    else
        u = -1.0e-4 * MDL.strain1n
    end

    envlpPosStrain[1] = u
    envlpPosStress[1] = u * k
    envlpNegStrain[1] = -u
    envlpNegStress[1] = -u * k

    envlpPosStrain[2:5] .= [MDL.strain1p, MDL.strain2p, MDL.strain3p, MDL.strain4p]
    envlpNegStrain[2:5] .= [MDL.strain1n, MDL.strain2n, MDL.strain3n, MDL.strain4n]

    envlpPosStress[2:5] .= [MDL.stress1p, MDL.stress2p, MDL.stress3p, MDL.stress4p]
    envlpNegStress[2:5] .= [MDL.stress1n, MDL.stress2n, MDL.stress3n, MDL.stress4n]

    # Post-peak branch extension
    k1 = (MDL.stress4p - MDL.stress3p) / (MDL.strain4p - MDL.strain3p)
    k2 = (MDL.stress4n - MDL.stress3n) / (MDL.strain4n - MDL.strain3n)

    envlpPosStrain[6] = 1e6 * MDL.strain4p
    envlpPosStress[6] = k1 > 0.0 ? MDL.stress4p + k1 * (envlpPosStrain[6] - MDL.strain4p) : 1.1 * MDL.stress4p

    envlpNegStrain[6] = 1e6 * MDL.strain4n
    envlpNegStress[6] = k2 > 0.0 ? MDL.stress4n + k2 * (envlpNegStrain[6] - MDL.strain4n) : 1.1 * MDL.stress4n

    # Elastic stiffnesses
    kElasticPos = envlpPosStress[2] / envlpPosStrain[2]
    kElasticNeg = envlpNegStress[2] / envlpNegStrain[2]

    # Energy capacity (trapezoidal rule)
    energypos = 0.5 * envlpPosStrain[1] * envlpPosStress[1]
    for j in 1:4
        energypos += 0.5 * (envlpPosStress[j] + envlpPosStress[j + 1]) * (envlpPosStrain[j + 1] - envlpPosStrain[j])
    end

    energyneg = 0.5 * envlpNegStrain[1] * envlpNegStress[1]
    for j in 1:4
        energyneg += 0.5 * (envlpNegStress[j] + envlpNegStress[j + 1]) * (envlpNegStrain[j + 1] - envlpNegStrain[j])
    end

    energyCapacity = MDL.gE * max(energypos, energyneg)

    return envlpPosStrain, envlpPosStress, envlpNegStrain, envlpNegStress, kElasticPos, kElasticNeg, energyCapacity
end


function revertToStart(
    envlpPosStrain::Vector{Real},
    envlpPosStress::Vector{Real},
    envlpNegStrain::Vector{Real},
    envlpNegStress::Vector{Real},
    kElasticPos::Real,
    kElasticNeg::Real
)
    # Initial converged state variables
    Cstate = 0
    Cstrain = 0.0
    Cstress = 0.0
    CstrainRate = 0.0

    lowCstateStrain = envlpNegStrain[1]
    lowCstateStress = envlpNegStress[1]
    hghCstateStrain = envlpPosStrain[1]
    hghCstateStress = envlpPosStress[1]

    CminStrainDmnd = envlpNegStrain[2]
    CmaxStrainDmnd = envlpPosStrain[2]

    Cenergy = 0.0
    CgammaK = 0.0
    CgammaD = 0.0
    CgammaF = 0.0
    CnCycle = 0.0

    # Trial variables
    Ttangent = envlpPosStress[1] / envlpPosStrain[1]
    dstrain = 0.0
    gammaKUsed = 0.0
    gammaFUsed = 0.0

    # Damage-influenced elastic stiffness and limits
    kElasticPosDamgd = kElasticPos
    kElasticNegDamgd = kElasticNeg
    uMaxDamgd = CmaxStrainDmnd
    uMinDamgd = CminStrainDmnd

    return Cstate, Cstrain, Cstress, CstrainRate, lowCstateStrain, lowCstateStress, hghCstateStrain, hghCstateStress,
           CminStrainDmnd, CmaxStrainDmnd, Cenergy, CgammaK, CgammaD, CgammaF, CnCycle,
           Ttangent, dstrain, gammaKUsed, gammaFUsed, kElasticPosDamgd, kElasticNegDamgd,
           uMaxDamgd, uMinDamgd
end


function revertToLastCommit(
    Cstate::Int,
    CstrainRate::Real,
    lowCstateStrain::Real,
    lowCstateStress::Real,
    hghCstateStrain::Real,
    hghCstateStress::Real,
    CminStrainDmnd::Real,
    CmaxStrainDmnd::Real,
    Cenergy::Real,
    Cstrain::Real,
    Cstress::Real,
    CgammaD::Real,
    CgammaK::Real,
    CgammaF::Real,
    CnCycle::Real
)
    Tstate = Cstate
    TstrainRate = CstrainRate

    lowTstateStrain = lowCstateStrain
    lowTstateStress = lowCstateStress
    hghTstateStrain = hghCstateStrain
    hghTstateStress = hghCstateStress
    TminStrainDmnd = CminStrainDmnd
    TmaxStrainDmnd = CmaxStrainDmnd
    Tenergy = Cenergy

    Tstrain = Cstrain
    Tstress = Cstress

    TgammaD = CgammaD
    TgammaK = CgammaK
    TgammaF = CgammaF

    TnCycle = CnCycle

    return Tstate, TstrainRate, lowTstateStrain, lowTstateStress, hghTstateStrain, hghTstateStress,
           TminStrainDmnd, TmaxStrainDmnd, Tenergy, Tstrain, Tstress,
           TgammaD, TgammaK, TgammaF, TnCycle
end


function setTrialStrain(
    strain::Real, CstrainRate::Real, Cstate::Int, Cenergy::Real,
    lowCstateStrain::Real, hghCstateStrain::Real,
    lowCstateStress::Real, hghCstateStress::Real,
    CminStrainDmnd::Real, CmaxStrainDmnd::Real,
    CgammaF::Real, CgammaK::Real, CgammaD::Real,
    envlpPosStress::Vector{Real}, envlpPosStrain::Vector{Real},
    kElasticPosDamgd::Real, kElasticNegDamgd::Real,
    state3Strain::Vector{Real}, state3Stress::Vector{Real},
    kunload::Real, state4Strain::Vector{Real}, state4Stress::Vector{Real},
    Cstrain::Real, uMaxDamgd::Real, uMinDamgd::Real,
    envlpNegStrain::Vector{Real}, envlpNegStress::Vector{Real},
    kElasticPos::Real, kElasticNeg::Real, Cstress::Real,
    DmgCyc::Int, CnCycle::Real, energyCapacity::Real,
    MDL::Inputs,
    Tstate::Int, Tenergy::Real, Tstrain::Real,
    lowTstateStrain::Real, hghTstateStrain::Real,
    lowTstateStress::Real, hghTstateStress::Real,
    TgammaF::Real, TgammaK::Real, TgammaD::Real,
    dstrain::Real, Ttangent::Real, Tstress::Real,
    elasticStrainEnergy::Real,
    TminStrainDmnd::Real, TmaxStrainDmnd::Real,
    gammaKUsed::Real, gammaFUsed::Real,
    envlpPosDamgdStress::Vector{Real}, envlpNegDamgdStress::Vector{Real},
    TnCycle::Real
)
    # Trial state initialization
    Tstate = Cstate
    Tenergy = Cenergy
    Tstrain = strain
    lowTstateStrain = lowCstateStrain
    hghTstateStrain = hghCstateStrain
    lowTstateStress = lowCstateStress
    hghTstateStress = hghCstateStress
    TminStrainDmnd = CminStrainDmnd
    TmaxStrainDmnd = CmaxStrainDmnd
    TgammaF = CgammaF
    TgammaK = CgammaK
    TgammaD = CgammaD

    # Strain increment
    dstrain = Tstrain - Cstrain
    if abs(dstrain) < 1e-12
        dstrain = 0.0
    end

    # State categorization
    lowTstateStrain, lowTstateStress, hghTstateStrain, hghTstateStress,
    TmaxStrainDmnd, gammaFUsed, envlpNegDamgdStress, gammaKUsed,
    kElasticPosDamgd, TminStrainDmnd, envlpPosDamgdStress, kElasticNegDamgd,
    Tstate =
        getstate(Tstrain, dstrain, CstrainRate, Tstate,
                 lowTstateStrain, lowTstateStress, hghTstateStrain, hghTstateStress,
                 envlpPosStrain, envlpPosStress, Cstrain, TmaxStrainDmnd,
                 uMaxDamgd, uMinDamgd, CgammaF, envlpNegStrain, envlpNegStress,
                 Cstress, CgammaK, kElasticPos, TminStrainDmnd, kElasticNeg,
                 gammaFUsed, envlpNegDamgdStress, gammaKUsed, kElasticPosDamgd,
                 envlpPosDamgdStress, kElasticNegDamgd)

    # Evaluate tangent and stress
    if Tstate == 0
        Ttangent = envlpPosStress[1] / envlpPosStrain[1]
        Tstress = Ttangent * Tstrain
    elseif Tstate == 1
        Tstress = posEnvlpStress(Tstrain, envlpPosDamgdStress, envlpPosStrain)
        Ttangent = posEnvlpTangent(Tstrain, envlpPosDamgdStress, envlpPosStrain)
    elseif Tstate == 2
        Ttangent = negEnvlpTangent(Tstrain, envlpNegDamgdStress, envlpNegStrain)
        Tstress = negEnvlpStress(Tstrain, envlpNegDamgdStress, envlpNegStrain)
    elseif Tstate == 3
        kunload = hghTstateStrain < 0.0 ? kElasticNegDamgd : kElasticPosDamgd
        state3Strain[1] = lowTstateStrain
        state3Strain[4] = hghTstateStrain
        state3Stress[1] = lowTstateStress
        state3Stress[4] = hghTstateStress

        state3Strain, state3Stress =
            getstate3(state3Strain, state3Stress, kunload, kElasticNegDamgd,
                      lowTstateStrain, lowTstateStress, TminStrainDmnd,
                      envlpNegStrain, envlpNegDamgdStress, hghTstateStrain,
                      hghTstateStress, MDL)

        Ttangent = Envlp3Tangent(state3Strain, state3Stress, Tstrain)
        Tstress = Envlp3Stress(state3Strain, state3Stress, Tstrain)
    elseif Tstate == 4
        kunload = lowTstateStrain < 0.0 ? kElasticNegDamgd : kElasticPosDamgd
        state4Strain[1] = lowTstateStrain
        state4Strain[4] = hghTstateStrain
        state4Stress[1] = lowTstateStress
        state4Stress[4] = hghTstateStress

        state4Strain, state4Stress =
            getstate4(state4Strain, state4Stress, kunload, kElasticPosDamgd,
                      lowTstateStrain, lowTstateStress, TmaxStrainDmnd,
                      envlpPosStrain, envlpPosDamgdStress, hghTstateStrain,
                      hghTstateStress, MDL)

        Ttangent = Envlp4Tangent(state4Strain, state4Stress, Tstrain)
        Tstress = Envlp4Stress(state4Strain, state4Stress, Tstrain)
    end

    # Energy calculation
    denergy = 0.5 * (Tstress + Cstress) * dstrain
    elasticStrainEnergy = Tstrain > 0.0 ?
        0.5 * Tstress / kElasticPosDamgd * Tstress :
        0.5 * Tstress / kElasticNegDamgd * Tstress

    Tenergy = Cenergy + denergy

    # Update degradation variables
    TnCycle, TgammaK, TgammaD, TgammaF =
        updateDmg(Tstrain, dstrain, TmaxStrainDmnd, TminStrainDmnd,
                  envlpPosStrain, envlpNegStrain, CnCycle, Tenergy,
                  energyCapacity, DmgCyc, elasticStrainEnergy,
                  envlpPosDamgdStress, envlpNegDamgdStress,
                  kElasticPos, kElasticNeg, TnCycle, TgammaK, TgammaD, TgammaF,
                  MDL)

    return Tstate, Tenergy, Tstrain, lowTstateStrain, hghTstateStrain,
           lowTstateStress, hghTstateStress, TgammaF, TgammaK, TgammaD,
           dstrain, Ttangent, Tstress, state3Strain, state3Stress,
           state4Strain, state4Stress, kunload, elasticStrainEnergy,
           TminStrainDmnd, TmaxStrainDmnd, gammaKUsed, gammaFUsed,
           kElasticPosDamgd, kElasticNegDamgd, envlpPosDamgdStress,
           envlpNegDamgdStress, TnCycle
end


function commitState(
    Tstate::Int, dstrain::Real, TstrainRate::Real,
    lowTstateStrain::Real, lowTstateStress::Real,
    hghTstateStrain::Real, hghTstateStress::Real,
    TminStrainDmnd::Real, TmaxStrainDmnd::Real,
    Tenergy::Real, Tstress::Real, Tstrain::Real,
    TgammaK::Real, TgammaD::Real, TgammaF::Real,
    kElasticPos::Real, kElasticNeg::Real,
    gammaKUsed::Real, gammaFUsed::Real,
    envlpPosStress::Vector{Real}, envlpNegStress::Vector{Real},
    TnCycle::Real
)
    # Commit state
    Cstate = Tstate

    CstrainRate = abs(dstrain) > 1e-12 ? dstrain : TstrainRate

    lowCstateStrain = lowTstateStrain
    lowCstateStress = lowTstateStress
    hghCstateStrain = hghTstateStrain
    hghCstateStress = hghTstateStress
    CminStrainDmnd = TminStrainDmnd
    CmaxStrainDmnd = TmaxStrainDmnd
    Cenergy = Tenergy

    Cstress = Tstress
    Cstrain = Tstrain

    CgammaK = TgammaK
    CgammaD = TgammaD
    CgammaF = TgammaF

    # Damage-modified elastic stiffness
    kElasticPosDamgd = kElasticPos * (1 - gammaKUsed)
    kElasticNegDamgd = kElasticNeg * (1 - gammaKUsed)

    # Damage-modified deformation limits
    uMaxDamgd = TmaxStrainDmnd * (1 + CgammaD)
    uMinDamgd = TminStrainDmnd * (1 + CgammaD)

    # Damage-modified strength envelopes
    envlpPosDamgdStress = (1 - gammaFUsed) .* envlpPosStress
    envlpNegDamgdStress = (1 - gammaFUsed) .* envlpNegStress

    CnCycle = TnCycle

    return Cstate, CstrainRate, lowCstateStrain, lowCstateStress,
           hghCstateStrain, hghCstateStress, CminStrainDmnd, CmaxStrainDmnd,
           Cenergy, Cstress, Cstrain, CgammaK, CgammaD, CgammaF,
           kElasticPosDamgd, kElasticNegDamgd, uMaxDamgd, uMinDamgd,
           envlpPosDamgdStress, envlpNegDamgdStress, CnCycle
end


function getstate(
    u::Real, du::Real, CstrainRate::Real, Tstate::Int,
    lowTstateStrain::Real, lowTstateStress::Real,
    hghTstateStrain::Real, hghTstateStress::Real,
    envlpPosStrain::Vector{Real}, envlpPosStress::Vector{Real},
    Cstrain::Real, TmaxStrainDmnd::Real,
    uMaxDamgd::Real, uMinDamgd::Real, CgammaF::Real,
    envlpNegStrain::Vector{Real}, envlpNegStress::Vector{Real},
    Cstress::Real, CgammaK::Real, kElasticPos::Real,
    TminStrainDmnd::Real, kElasticNeg::Real,
    gammaFUsed::Real, envlpNegDamgdStress::Vector{Real},
    gammaKUsed::Real, kElasticPosDamgd::Real,
    envlpPosDamgdStress::Vector{Real}, kElasticNegDamgd::Real
)

    cid = du * CstrainRate <= 0.0
    cis = false
    newState = 0

    if u < lowTstateStrain || u > hghTstateStrain || cid
        if Tstate == 0
            if u > hghTstateStrain
                cis = true
                newState = 1
                lowTstateStrain = envlpPosStrain[1]
                lowTstateStress = envlpPosStress[1]
                hghTstateStrain = envlpPosStrain[6]
                hghTstateStress = envlpPosStress[6]
            elseif u < lowTstateStrain
                cis = true
                newState = 2
                lowTstateStrain = envlpNegStrain[6]
                lowTstateStress = envlpNegStress[6]
                hghTstateStrain = envlpNegStrain[1]
                hghTstateStress = envlpNegStress[1]
            end
        elseif Tstate == 1 && du < 0.0
            cis = true
            TmaxStrainDmnd = max(u - du, TmaxStrainDmnd, uMaxDamgd)
            if u < uMinDamgd
                newState = 2
                gammaFUsed = CgammaF
                envlpNegDamgdStress .= envlpNegStress .* (1.0 - gammaFUsed)
                lowTstateStrain = envlpNegStrain[6]
                lowTstateStress = envlpNegStress[6]
                hghTstateStrain = envlpNegStrain[1]
                hghTstateStress = envlpNegStress[1]
            else
                newState = 3
                lowTstateStrain = uMinDamgd
                gammaFUsed = CgammaF
                envlpNegDamgdStress .= envlpNegStress .* (1.0 - gammaFUsed)
                lowTstateStress = negEnvlpStress(uMinDamgd, envlpNegDamgdStress, envlpNegStrain)
                hghTstateStrain = Cstrain
                hghTstateStress = Cstress
            end
            gammaKUsed = CgammaK
            kElasticPosDamgd = kElasticPos * (1.0 - gammaKUsed)
        elseif Tstate == 2 && du > 0.0
            cis = true
            TminStrainDmnd = min(Cstrain, TminStrainDmnd, uMinDamgd)
            if u > uMaxDamgd
                newState = 1
                gammaFUsed = CgammaF
                envlpPosDamgdStress .= envlpPosStress .* (1.0 - gammaFUsed)
                lowTstateStrain = envlpPosStrain[1]
                lowTstateStress = envlpPosStress[1]
                hghTstateStrain = envlpPosStrain[6]
                hghTstateStress = envlpPosStress[6]
            else
                newState = 4
                lowTstateStrain = Cstrain
                lowTstateStress = Cstress
                hghTstateStrain = uMaxDamgd
                gammaFUsed = CgammaF
                envlpPosDamgdStress .= envlpPosStress .* (1.0 - gammaFUsed)
                hghTstateStress = posEnvlpStress(uMaxDamgd, envlpPosDamgdStress, envlpPosStrain)
            end
            gammaKUsed = CgammaK
            kElasticNegDamgd = kElasticNeg * (1.0 - gammaKUsed)
        elseif Tstate == 3
            if u < lowTstateStrain
                cis = true
                newState = 2
                lowTstateStrain = envlpNegStrain[6]
                hghTstateStrain = envlpNegStrain[1]
                lowTstateStress = envlpNegDamgdStress[6]
                hghTstateStress = envlpNegDamgdStress[1]
            elseif u > uMaxDamgd && du > 0.0
                cis = true
                newState = 1
                lowTstateStrain = envlpPosStrain[1]
                lowTstateStress = envlpPosStress[1]
                hghTstateStrain = envlpPosStrain[6]
                hghTstateStress = envlpPosStress[6]
            elseif du > 0.0
                cis = true
                newState = 4
                lowTstateStrain = Cstrain
                lowTstateStress = Cstress
                hghTstateStrain = uMaxDamgd
                gammaFUsed = CgammaF
                envlpPosDamgdStress .= envlpPosStress .* (1.0 - gammaFUsed)
                hghTstateStress = posEnvlpStress(uMaxDamgd, envlpPosDamgdStress, envlpPosStrain)
                gammaKUsed = CgammaK
                kElasticNegDamgd = kElasticNeg * (1.0 - gammaKUsed)
            end
        elseif Tstate == 4
            if u > hghTstateStrain
                cis = true
                newState = 1
                lowTstateStrain = envlpPosStrain[1]
                lowTstateStress = envlpPosDamgdStress[1]
                hghTstateStrain = envlpPosStrain[6]
                hghTstateStress = envlpPosDamgdStress[6]
            elseif u < uMinDamgd && du < 0.0
                cis = true
                newState = 2
                lowTstateStrain = envlpNegStrain[6]
                lowTstateStress = envlpNegDamgdStress[6]
                hghTstateStrain = envlpNegStrain[1]
                hghTstateStress = envlpNegDamgdStress[1]
            elseif du < 0.0
                cis = true
                newState = 3
                lowTstateStrain = uMinDamgd
                gammaFUsed = CgammaF
                envlpNegDamgdStress .= envlpNegStress .* (1.0 - gammaFUsed)
                lowTstateStress = negEnvlpStress(uMinDamgd, envlpNegDamgdStress, envlpNegStrain)
                hghTstateStrain = Cstrain
                hghTstateStress = Cstress
                gammaKUsed = CgammaK
                kElasticPosDamgd = kElasticPos * (1.0 - gammaKUsed)
            end
        end
    end

    if cis
        Tstate = newState
    end

    return lowTstateStrain, lowTstateStress, hghTstateStrain, hghTstateStress,
           TmaxStrainDmnd, gammaFUsed, envlpNegDamgdStress, gammaKUsed,
           kElasticPosDamgd, TminStrainDmnd, envlpPosDamgdStress,
           kElasticNegDamgd, Tstate
end


function posEnvlpStress(u::Real, envlpPosDamgdStress::Vector{Real}, envlpPosStrain::Vector{Real})
    k = 0.0
    f = 0.0
    i = 1

    # Loop to find correct envelope segment
    while k == 0.0 && i <= 5
        if u <= envlpPosStrain[i + 1]
            k = (envlpPosDamgdStress[i + 1] - envlpPosDamgdStress[i]) /
                (envlpPosStrain[i + 1] - envlpPosStrain[i])
            f = envlpPosDamgdStress[i] + (u - envlpPosStrain[i]) * k
        end
        i += 1
    end

    # If not found in any of the first five segments, extrapolate at end
    if k == 0.0
        k = (envlpPosDamgdStress[6] - envlpPosDamgdStress[5]) /
            (envlpPosStrain[6] - envlpPosStrain[5])
        f = envlpPosDamgdStress[6] + k * (u - envlpPosStrain[6])
    end

    return f
end


function negEnvlpStress(u::Real, envlpNegDamgdStress::Vector{Real}, envlpNegStrain::Vector{Real})
    k = 0.0
    f = 0.0
    i = 1

    # Loop to find appropriate segment
    while k == 0.0 && i <= 5
        if u >= envlpNegStrain[i + 1]
            k = (envlpNegDamgdStress[i] - envlpNegDamgdStress[i + 1]) /
                (envlpNegStrain[i] - envlpNegStrain[i + 1])
            f = envlpNegDamgdStress[i + 1] + (u - envlpNegStrain[i + 1]) * k
        end
        i += 1
    end

    # Extrapolate at the far end if no segment matched
    if k == 0.0
        k = (envlpNegDamgdStress[5] - envlpNegDamgdStress[6]) /
            (envlpNegStrain[5] - envlpNegStrain[6])
        f = envlpNegDamgdStress[6] + k * (u - envlpNegStrain[6])
    end

    return f
end



function updateDmg(
    strain::Real, dstrain::Real,
    TmaxStrainDmnd::Real, TminStrainDmnd::Real,
    envlpPosStrain::Vector{Real}, envlpNegStrain::Vector{Real},
    CnCycle::Real, Tenergy::Real, energyCapacity::Real,
    DmgCyc::Int, elasticStrainEnergy::Real,
    envlpPosDamgdStress::Vector{Real}, envlpNegDamgdStress::Vector{Real},
    kElasticPos::Real, kElasticNeg::Real,
    TnCycle::Real, TgammaK::Real, TgammaD::Real, TgammaF::Real,
    MDL::Inputs
)
    tes = 0.0
    umaxAbs = max(TmaxStrainDmnd, -TminStrainDmnd)
    uultAbs = max(envlpPosStrain[5], -envlpNegStrain[5])

    TnCycle = CnCycle + abs(dstrain) / (4 * umaxAbs)

    if (abs(strain) < uultAbs) && (Tenergy < energyCapacity)
        ratio = umaxAbs / uultAbs
        TgammaK = MDL.gammaK1 * ratio^MDL.gammaK3
        TgammaD = MDL.gammaD1 * ratio^MDL.gammaD3
        TgammaF = MDL.gammaF1 * ratio^MDL.gammaF3

        if Tenergy > elasticStrainEnergy && DmgCyc == 0
            tes = (Tenergy - elasticStrainEnergy) / energyCapacity
            TgammaK += MDL.gammaK2 * tes^MDL.gammaK4
            TgammaD += MDL.gammaD2 * tes^MDL.gammaD4
            TgammaF += MDL.gammaF2 * tes^MDL.gammaF4
        elseif DmgCyc == 1
            TgammaK += MDL.gammaK2 * TnCycle^MDL.gammaK4
            TgammaD += MDL.gammaD2 * TnCycle^MDL.gammaD4
            TgammaF += MDL.gammaF2 * TnCycle^MDL.gammaF4
        end

        kminP = posEnvlpStress(TmaxStrainDmnd, envlpPosDamgdStress, envlpPosStrain) / TmaxStrainDmnd
        kminN = negEnvlpStress(TminStrainDmnd, envlpNegDamgdStress, envlpNegStrain) / TminStrainDmnd
        kmin = max(kminP / kElasticPos, kminN / kElasticNeg)
        gammaKLimEnv = max(0.0, 1.0 - kmin)

        k1 = min(TgammaK, MDL.gammaKLimit)
        TgammaK = min(k1, gammaKLimEnv)
        TgammaD = min(TgammaD, MDL.gammaDLimit)
        TgammaF = min(TgammaF, MDL.gammaFLimit)

    elseif abs(strain) < uultAbs
        kminP = posEnvlpStress(TmaxStrainDmnd, envlpPosDamgdStress, envlpPosStrain) / TmaxStrainDmnd
        kminN = negEnvlpStress(TminStrainDmnd, envlpNegDamgdStress, envlpNegStrain) / TminStrainDmnd
        kmin = max(kminP / kElasticPos, kminN / kElasticNeg)
        gammaKLimEnv = max(0.0, 1.0 - kmin)

        TgammaK = min(MDL.gammaKLimit, gammaKLimEnv)
        TgammaD = min(TgammaD, MDL.gammaDLimit)
        TgammaF = min(TgammaF, MDL.gammaFLimit)
    end

    return TnCycle, TgammaK, TgammaD, TgammaF
end


function Envlp3Tangent(s3Strain::Vector{Real}, s3Stress::Vector{Real}, u::Real)
    k = 0.0
    i = 1

    while (k == 0.0 || i <= 3) && i <= 3
        if u >= s3Strain[i]
            k = (s3Stress[i + 1] - s3Stress[i]) / (s3Strain[i + 1] - s3Strain[i])
        end
        i += 1
    end

    if k == 0.0
        i = u < s3Strain[1] ? 1 : 3
        k = (s3Stress[i + 1] - s3Stress[i]) / (s3Strain[i + 1] - s3Strain[i])
    end

    return k
end


function Envlp3Stress(s3Strain::Vector{Real}, s3Stress::Vector{Real}, u::Real)
    k = 0.0
    f = 0.0
    i = 1

    while (k == 0.0 || i <= 3) && i <= 3
        if u >= s3Strain[i]
            k = (s3Stress[i + 1] - s3Stress[i]) / (s3Strain[i + 1] - s3Strain[i])
            f = s3Stress[i] + (u - s3Strain[i]) * k
        end
        i += 1
    end

    if k == 0.0
        i = u < s3Strain[1] ? 1 : 3
        k = (s3Stress[i + 1] - s3Stress[i]) / (s3Strain[i + 1] - s3Strain[i])
        f = s3Stress[i] + (u - s3Strain[i]) * k
    end

    return f
end


function Envlp4Tangent(s4Strain::Vector{Real}, s4Stress::Vector{Real}, u::Real)
    k = 0.0
    i = 1

    while (k == 0.0 || i <= 3) && i <= 3
        if u >= s4Strain[i]
            k = (s4Stress[i + 1] - s4Stress[i]) / (s4Strain[i + 1] - s4Strain[i])
        end
        i += 1
    end

    if k == 0.0
        i = u < s4Strain[1] ? 1 : 3
        k = (s4Stress[i + 1] - s4Stress[i]) / (s4Strain[i + 1] - s4Strain[i])
    end

    return k
end


function Envlp4Stress(s4Strain::Vector{Real}, s4Stress::Vector{Real}, u::Real)
    k = 0.0
    f = 0.0
    i = 1

    while (k == 0.0 || i <= 3) && i <= 3
        if u >= s4Strain[i]
            k = (s4Stress[i + 1] - s4Stress[i]) / (s4Strain[i + 1] - s4Strain[i])
            f = s4Stress[i] + (u - s4Strain[i]) * k
        end
        i += 1
    end

    if k == 0.0
        i = u < s4Strain[1] ? 1 : 3
        k = (s4Stress[i + 1] - s4Stress[i]) / (s4Strain[i + 1] - s4Strain[i])
        f = s4Stress[i] + (u - s4Strain[i]) * k
    end

    return f
end


function getstate3!(
    state3Strain::Vector{Real},
    state3Stress::Vector{Real},
    kunload::Real,
    kElasticNegDamgd::Real,
    lowTstateStrain::Real,
    lowTstateStress::Real,
    TminStrainDmnd::Real,
    envlpNegStrain::Vector{Real},
    envlpNegDamgdStress::Vector{Real},
    hghTstateStrain::Real,
    hghTstateStress::Real,
    MDL
)
    kmax = max(kunload, kElasticNegDamgd)

    if state3Strain[1] * state3Strain[4] < 0.0
        state3Strain[2] = lowTstateStrain * MDL.rDispN
        if MDL.rForceN - MDL.uForceN > 1e-8
            state3Stress[2] = lowTstateStress * MDL.rForceN
        else
            if TminStrainDmnd < envlpNegStrain[4]
                st1 = lowTstateStress * MDL.uForceN * (1 + 1e-6)
                st2 = envlpNegDamgdStress[5] * (1 + 1e-6)
                state3Stress[2] = min(st1, st2)
            else
                st1 = envlpNegDamgdStress[4] * MDL.uForceN * (1 + 1e-6)
                st2 = envlpNegDamgdStress[5] * (1 + 1e-6)
                state3Stress[2] = min(st1, st2)
            end
        end

        if (state3Stress[2] - state3Stress[1]) / (state3Strain[2] - state3Strain[1]) > kElasticNegDamgd
            state3Strain[2] = state3Strain[1] + (state3Stress[2] - state3Stress[1]) / kElasticNegDamgd
        end

        if state3Strain[2] > state3Strain[4]
            du = state3Strain[4] - state3Strain[1]
            df = state3Stress[4] - state3Stress[1]
            state3Strain[2] = state3Strain[1] + 0.33 * du
            state3Strain[3] = state3Strain[1] + 0.67 * du
            state3Stress[2] = state3Stress[1] + 0.33 * df
            state3Stress[3] = state3Stress[1] + 0.67 * df
        else
            state3Stress[3] = if TminStrainDmnd < envlpNegStrain[4]
                MDL.uForceN * envlpNegDamgdStress[5]
            else
                MDL.uForceN * envlpNegDamgdStress[4]
            end

            state3Strain[3] = hghTstateStrain - (hghTstateStress - state3Stress[3]) / kunload

            if state3Strain[3] > state3Strain[4]
                du = state3Strain[4] - state3Strain[2]
                df = state3Stress[4] - state3Stress[2]
                state3Strain[3] = state3Strain[2] + 0.5 * du
                state3Stress[3] = state3Stress[2] + 0.5 * df
            elseif (state3Stress[3] - state3Stress[2]) / (state3Strain[3] - state3Strain[2]) > kmax
                du = state3Strain[4] - state3Strain[1]
                df = state3Stress[4] - state3Stress[1]
                state3Strain[2] = state3Strain[1] + 0.33 * du
                state3Strain[3] = state3Strain[1] + 0.67 * du
                state3Stress[2] = state3Stress[1] + 0.33 * df
                state3Stress[3] = state3Stress[1] + 0.67 * df
            elseif state3Strain[3] < state3Strain[2] || (state3Stress[3] - state3Stress[2]) / (state3Strain[3] - state3Strain[2]) < 0
                if state3Strain[3] < 0.0
                    du = state3Strain[4] - state3Strain[2]
                    df = state3Stress[4] - state3Stress[2]
                    state3Strain[3] = state3Strain[2] + 0.5 * du
                    state3Stress[3] = state3Stress[2] + 0.5 * df
                elseif state3Strain[2] > 0.0
                    du = state3Strain[3] - state3Strain[1]
                    df = state3Stress[3] - state3Stress[1]
                    state3Strain[2] = state3Strain[1] + 0.5 * du
                    state3Stress[2] = state3Stress[1] + 0.5 * df
                else
                    avgforce = 0.5 * (state3Stress[3] + state3Stress[2])
                    dfr = avgforce < 0.0 ? -avgforce / 100 : avgforce / 100
                    slope12 = (state3Stress[2] - state3Stress[1]) / (state3Strain[2] - state3Strain[1])
                    slope34 = (state3Stress[4] - state3Stress[3]) / (state3Strain[4] - state3Strain[3])
                    state3Stress[2] = avgforce - dfr
                    state3Stress[3] = avgforce + dfr
                    state3Strain[2] = state3Strain[1] + (state3Stress[2] - state3Stress[1]) / slope12
                    state3Strain[3] = state3Strain[4] - (state3Stress[4] - state3Stress[3]) / slope34
                end
            end
        end
    else
        du = state3Strain[4] - state3Strain[1]
        df = state3Stress[4] - state3Stress[1]
        state3Strain[2] = state3Strain[1] + 0.33 * du
        state3Strain[3] = state3Strain[1] + 0.67 * du
        state3Stress[2] = state3Stress[1] + 0.33 * df
        state3Stress[3] = state3Stress[1] + 0.67 * df
    end

    checkSlope = state3Stress[1] / state3Strain[1]
    slope = 0.0

    for i in 1:3
        du = state3Strain[i+1] - state3Strain[i]
        df = state3Stress[i+1] - state3Stress[i]
        if du < 0.0 || df < 0.0
            du = state3Strain[4] - state3Strain[1]
            df = state3Stress[4] - state3Stress[1]
            state3Strain[2] = state3Strain[1] + 0.33 * du
            state3Strain[3] = state3Strain[1] + 0.67 * du
            state3Stress[2] = state3Stress[1] + 0.33 * df
            state3Stress[3] = state3Stress[1] + 0.67 * df
            slope = df / du
            break
        end
        if slope > 1e-8 && slope < checkSlope
            state3Strain[2] = 0.0
            state3Stress[2] = 0.0
            state3Strain[3] = state3Strain[4] / 2
            state3Stress[3] = state3Stress[4] / 2
        end
    end

    return state3Strain, state3Stress
end


function getstate4(state4Strain, state4Stress, kunload, kElasticPosDamgd,
                    lowTstateStrain, lowTstateStress, TmaxStrainDmnd,
                    envlpPosStrain, envlpPosDamgdStress,
                    hghTstateStrain, hghTstateStress, MDL)

    kmax = max(kunload, kElasticPosDamgd)

    if state4Strain[1] * state4Strain[4] < 0.0
        state4Strain[3] = hghTstateStrain * MDL.rDispP

        if MDL.uForceP == 0.0 || (MDL.rForceP - MDL.uForceP > 1e-8)
            state4Stress[3] = hghTstateStress * MDL.rForceP
        else
            if TmaxStrainDmnd > envlpPosStrain[4]
                st1 = hghTstateStress * MDL.uForceP * (1.0 + 1e-6)
                st2 = envlpPosDamgdStress[5] * (1.0 + 1e-6)
                state4Stress[3] = max(st1, st2)
            else
                st1 = envlpPosDamgdStress[4] * MDL.uForceP * (1.0 + 1e-6)
                st2 = envlpPosDamgdStress[5] * (1.0 + 1e-6)
                state4Stress[3] = max(st1, st2)
            end
        end

        if (state4Stress[4] - state4Stress[3]) / (state4Strain[4] - state4Strain[3]) > kElasticPosDamgd
            state4Strain[3] = hghTstateStrain - (state4Stress[4] - state4Stress[3]) / kElasticPosDamgd
        end

        if state4Strain[3] < state4Strain[1]
            du = state4Strain[4] - state4Strain[1]
            df = state4Stress[4] - state4Stress[1]
            state4Strain[2] = state4Strain[1] + 0.33 * du
            state4Strain[3] = state4Strain[1] + 0.67 * du
            state4Stress[2] = state4Stress[1] + 0.33 * df
            state4Stress[3] = state4Stress[1] + 0.67 * df
        else
            if TmaxStrainDmnd > envlpPosStrain[4]
                state4Stress[2] = MDL.uForceP * envlpPosDamgdStress[5]
            else
                state4Stress[2] = MDL.uForceP * envlpPosDamgdStress[4]
            end

            state4Strain[2] = lowTstateStrain + (-lowTstateStress + state4Stress[2]) / kunload

            if state4Strain[2] < state4Strain[1] ||
               (state4Stress[3] - state4Stress[2]) / (state4Strain[3] - state4Strain[2]) > kmax
                du = state4Strain[4] - state4Strain[1]
                df = state4Stress[4] - state4Stress[1]
                state4Strain[2] = state4Strain[1] + 0.33 * du
                state4Strain[3] = state4Strain[1] + 0.67 * du
                state4Stress[2] = state4Stress[1] + 0.33 * df
                state4Stress[3] = state4Stress[1] + 0.67 * df
            elseif (state4Strain[3] < state4Strain[2]) ||
                   ((state4Stress[3] - state4Stress[2]) / (state4Strain[3] - state4Strain[2]) < 0)
                if state4Strain[2] > 0.0
                    du = state4Strain[3] - state4Strain[1]
                    df = state4Stress[3] - state4Stress[1]
                    state4Strain[2] = state4Strain[1] + 0.5 * du
                    state4Stress[2] = state4Stress[1] + 0.5 * df
                elseif state4Strain[3] < 0.0
                    du = state4Strain[4] - state4Strain[2]
                    df = state4Stress[4] - state4Stress[2]
                    state4Strain[3] = state4Strain[2] + 0.5 * du
                    state4Stress[3] = state4Stress[2] + 0.5 * df
                else
                    avgforce = 0.5 * (state4Stress[3] + state4Stress[2])
                    dfr = abs(avgforce / 100)
                    slope12 = (state4Stress[2] - state4Stress[1]) / (state4Strain[2] - state4Strain[1])
                    slope34 = (state4Stress[4] - state4Stress[3]) / (state4Strain[4] - state4Strain[3])
                    state4Stress[2] = avgforce - dfr
                    state4Stress[3] = avgforce + dfr
                    state4Strain[2] = state4Strain[1] + (state4Stress[2] - state4Stress[1]) / slope12
                    state4Strain[3] = state4Strain[4] - (state4Stress[4] - state4Stress[3]) / slope34
                end
            end
        end
    else
        du = state4Strain[4] - state4Strain[1]
        df = state4Stress[4] - state4Stress[1]
        state4Strain[2] = state4Strain[1] + 0.33 * du
        state4Strain[3] = state4Strain[1] + 0.67 * du
        state4Stress[2] = state4Stress[1] + 0.33 * df
        state4Stress[3] = state4Stress[1] + 0.67 * df
    end

    checkSlope = state4Stress[1] / state4Strain[1]
    slope = 0.0

    for i in 1:3
        du = state4Strain[i+1] - state4Strain[i]
        df = state4Stress[i+1] - state4Stress[i]
        if du < 0.0 || df < 0.0
            du = state4Strain[4] - state4Strain[1]
            df = state4Stress[4] - state4Stress[1]
            state4Strain[2] = state4Strain[1] + 0.33 * du
            state4Strain[3] = state4Strain[1] + 0.67 * du
            state4Stress[2] = state4Stress[1] + 0.33 * df
            state4Stress[3] = state4Stress[1] + 0.67 * df
            slope = df / du
            break
        end
    end

    if slope > 1e-8 && slope < checkSlope
        state4Strain[2] = 0.0
        state4Stress[2] = 0.0
        state4Strain[3] = state4Strain[4] / 2
        state4Stress[3] = state4Stress[4] / 2
    end

    return state4Strain, state4Stress

end


function posEnvlpTangent(u::Real, envlpPosDamgdStress::Vector{Real}, envlpPosStrain::Vector{Real})::Real
    k = 0.0
    i = 1

    while k == 0.0 && i <= 5
        if u <= envlpPosStrain[i+1]
            k = (envlpPosDamgdStress[i+1] - envlpPosDamgdStress[i]) / (envlpPosStrain[i+1] - envlpPosStrain[i])
        end
        i += 1
    end

    if k == 0.0
        k = (envlpPosDamgdStress[6] - envlpPosDamgdStress[5]) / (envlpPosStrain[6] - envlpPosStrain[5])
    end

    return k
end


function getstate3(state3Strain, state3Stress, kunload, kElasticNegDamgd,
                    lowTstateStrain, lowTstateStress, TminStrainDmnd,
                    envlpNegStrain, envlpNegDamgdStress,
                    hghTstateStrain, hghTstateStress, MDL)

    kmax = max(kunload, kElasticNegDamgd)

    if state3Strain[1] * state3Strain[4] < 0.0
        state3Strain[2] = lowTstateStrain * MDL.rDispN

        if MDL.rForceN - MDL.uForceN > 1e-8
            state3Stress[2] = lowTstateStress * MDL.rForceN
        else
            if TminStrainDmnd < envlpNegStrain[4]
                st1 = lowTstateStress * MDL.uForceN * (1.0 + 1e-6)
                st2 = envlpNegDamgdStress[5] * (1.0 + 1e-6)
                state3Stress[2] = min(st1, st2)
            else
                st1 = envlpNegDamgdStress[4] * MDL.uForceN * (1.0 + 1e-6)
                st2 = envlpNegDamgdStress[5] * (1.0 + 1e-6)
                state3Stress[2] = min(st1, st2)
            end
        end

        if (state3Stress[2] - state3Stress[1]) / (state3Strain[2] - state3Strain[1]) > kElasticNegDamgd
            state3Strain[2] = lowTstateStrain + (state3Stress[2] - state3Stress[1]) / kElasticNegDamgd
        end

        if state3Strain[2] > state3Strain[4]
            du = state3Strain[4] - state3Strain[1]
            df = state3Stress[4] - state3Stress[1]
            state3Strain[2] = state3Strain[1] + 0.33 * du
            state3Strain[3] = state3Strain[1] + 0.67 * du
            state3Stress[2] = state3Stress[1] + 0.33 * df
            state3Stress[3] = state3Stress[1] + 0.67 * df
        else
            state3Stress[3] = TminStrainDmnd < envlpNegStrain[4] ?
                              MDL.uForceN * envlpNegDamgdStress[5] :
                              MDL.uForceN * envlpNegDamgdStress[4]

            state3Strain[3] = hghTstateStrain - (hghTstateStress - state3Stress[3]) / kunload

            if state3Strain[3] > state3Strain[4]
                du = state3Strain[4] - state3Strain[2]
                df = state3Stress[4] - state3Stress[2]
                state3Strain[3] = state3Strain[2] + 0.5 * du
                state3Stress[3] = state3Stress[2] + 0.5 * df
            elseif (state3Stress[3] - state3Stress[2]) / (state3Strain[3] - state3Strain[2]) > kmax
                du = state3Strain[4] - state3Strain[1]
                df = state3Stress[4] - state3Stress[1]
                state3Strain[2] = state3Strain[1] + 0.33 * du
                state3Strain[3] = state3Strain[1] + 0.67 * du
                state3Stress[2] = state3Stress[1] + 0.33 * df
                state3Stress[3] = state3Stress[1] + 0.67 * df
            elseif state3Strain[3] < state3Strain[2] ||
                   (state3Stress[3] - state3Stress[2]) / (state3Strain[3] - state3Strain[2]) < 0
                if state3Strain[3] < 0.0
                    du = state3Strain[4] - state3Strain[2]
                    df = state3Stress[4] - state3Stress[2]
                    state3Strain[3] = state3Strain[2] + 0.5 * du
                    state3Stress[3] = state3Stress[2] + 0.5 * df
                elseif state3Strain[2] > 0.0
                    du = state3Strain[3] - state3Strain[1]
                    df = state3Stress[3] - state3Stress[1]
                    state3Strain[2] = state3Strain[1] + 0.5 * du
                    state3Stress[2] = state3Stress[1] + 0.5 * df
                else
                    avgforce = 0.5 * (state3Stress[3] + state3Stress[2])
                    dfr = abs(avgforce / 100)
                    slope12 = (state3Stress[2] - state3Stress[1]) / (state3Strain[2] - state3Strain[1])
                    slope34 = (state3Stress[4] - state3Stress[3]) / (state3Strain[4] - state3Strain[3])
                    state3Stress[2] = avgforce - dfr
                    state3Stress[3] = avgforce + dfr
                    state3Strain[2] = state3Strain[1] + (state3Stress[2] - state3Stress[1]) / slope12
                    state3Strain[3] = state3Strain[4] - (state3Stress[4] - state3Stress[3]) / slope34
                end
            end
        end
    else
        du = state3Strain[4] - state3Strain[1]
        df = state3Stress[4] - state3Stress[1]
        state3Strain[2] = state3Strain[1] + 0.33 * du
        state3Strain[3] = state3Strain[1] + 0.67 * du
        state3Stress[2] = state3Stress[1] + 0.33 * df
        state3Stress[3] = state3Stress[1] + 0.67 * df
    end

    checkSlope = state3Stress[1] / state3Strain[1]
    slope = 0.0

    for i in 1:3
        du = state3Strain[i+1] - state3Strain[i]
        df = state3Stress[i+1] - state3Stress[i]
        if du < 0.0 || df < 0.0
            du = state3Strain[4] - state3Strain[1]
            df = state3Stress[4] - state3Stress[1]
            state3Strain[2] = state3Strain[1] + 0.33 * du
            state3Strain[3] = state3Strain[1] + 0.67 * du
            state3Stress[2] = state3Stress[1] + 0.33 * df
            state3Stress[3] = state3Stress[1] + 0.67 * df
            slope = df / du
            break
        end
    end

    if slope > 1e-8 && slope < checkSlope
        state3Strain[2] = 0.0
        state3Stress[2] = 0.0
        state3Strain[3] = state3Strain[4] / 2
        state3Stress[3] = state3Stress[4] / 2
    end

    return state3Strain, state3Stress

end


function negEnvlpTangent(u, envlpNegDamgdStress, envlpNegStrain)
    k = 0.0
    i = 1

    while k == 0.0 && i <= 5
        if u >= envlpNegStrain[i+1]
            k = (envlpNegDamgdStress[i] - envlpNegDamgdStress[i+1]) /
                (envlpNegStrain[i] - envlpNegStrain[i+1])
        end
        i += 1
    end

    if k == 0.0
        k = (envlpNegDamgdStress[5] - envlpNegDamgdStress[6]) /
            (envlpNegStrain[5] - envlpNegStrain[6])
    end

    return k
end



end # module Pinching4Model
