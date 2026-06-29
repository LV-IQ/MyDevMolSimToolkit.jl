module Reweighting

using ProgressMeter
using LinearAlgebra: diag
using CellListMap: ParticleSystem, map_pairwise, map_pairwise!
using ..MolSimToolkit: Simulation, positions, unitcell
using ..MolSimToolkit.MolecularMinimumDistances
import OrderedCollections
import PDBTools

export reweight, lennard_jones_perturbation, Perturbation, SystemPerturbations, SystemPerturbationsOneGroup, gauss

const testdir = "$(@__DIR__)/test"

include("./Reweight.jl")
include("./perturbation_examples.jl")

end #Module Reweighting

@testitem "Reweight with small trajectory using minimum distances and atoms from the same molecule" begin
    import PDBTools
    import OrderedCollections
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory.xtc")

    g1 = PDBTools.selindex(simulation.atoms, "residue 4981 and name F12")

    g2 = PDBTools.selindex(simulation.atoms, at -> at.residue == 4981 && at.name in ["H", "H21", "H22"])

    c1 = at -> at = true

    c2 = at -> at.residue == 4981 && at.name in ["H21", "H22"]

    dist(r, δ) = r*δ

    Dict = OrderedCollections.OrderedDict(1 => Perturbation(simulation.atoms, c1, c2, dist, [1]))

    input = SystemPerturbations(g1, 1, g2, 3, Dict)

    res = reweight(simulation,
                    input;
                    all_distances = false,
                    k = 1.0,
                    T = 1.0,
                    cutoff = 12.0,
            )

    @test res[1].energies[1] ≈ [
        2.66700,
        0.0,
        0.0,
        2.64641,
        2.42596,
        2.61084,
        2.40289,
        2.48584,
        2.86653,
        2.86769      
    ] atol = 1.e-4

    @test res[1].distances ≈ [
        1,
        0,
        0,
        1,
        1,
        1,
        1,
        1,
        1,
        1,  
    ] atol = 1.e-4
end

@testitem "Reweight with small trajectory using all distances and atoms from the same molecule" begin
    import PDBTools
    import OrderedCollections
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory.xtc")

    g1 = PDBTools.selindex(simulation.atoms, "residue 4981 and name F12")

    g2 = PDBTools.selindex(simulation.atoms, at -> at.residue == 4981 && at.name in ["H", "H21", "H22"])

    c1 = at -> at = true

    c2 = at -> at.residue == 4981 && at.name in ["H21", "H22"]

    dist(r, δ) = r*δ

    Dict = OrderedCollections.OrderedDict(1 => Perturbation(simulation.atoms, c1, c2, dist, [-2, -1, 0, 1, 2]))

    input = SystemPerturbations(g1, 1, g2, 3, Dict)

    res = reweight(simulation,
                    input;
                    all_distances = true,
                    k = 1.0,
                    T = 1.0,
                    cutoff = 12.0,
            )
    @test res[1].energies ≈ [
        -2*[5.450897,6.018059,6.05581,5.926243,5.067425,5.816867,5.229802,5.382748,6.127764, 6.212621],
        -[5.450897,6.018059,6.05581,5.926243,5.067425,5.816867,5.229802,5.382748,6.127764, 6.212621],
        zeros(10),
        [5.450897,6.018059,6.05581,5.926243,5.067425,5.816867,5.229802,5.382748,6.127764, 6.212621],
        2*[5.450897,6.018059,6.05581,5.926243,5.067425,5.816867,5.229802,5.382748,6.127764, 6.212621],
    ] atol = 1.e-4
end

@testitem "Reweight with small trajectory using minimum distances and contributions from different residues" begin
    import PDBTools
    import OrderedCollections
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory.xtc")

    g1 = PDBTools.selindex(simulation.atoms, at -> at.residue == 15 && at.name in ["H", "N", "C", "OC1", "OC2"])

    g2 = PDBTools.selindex(simulation.atoms,at -> at.residue == 11 && at.name in ["CB", "HB1", "HB2", "HB3"])

    c1 = at -> at.name in ["H"]

    c2 = at -> at.name in ["HB3"]

    c11 = at -> at.name in ["H"]

    c12 = at -> at = true

    dist(r, δ) = r*δ

    Dict = OrderedCollections.OrderedDict(
        "a" => Perturbation(simulation.atoms, c1, c2, dist, [-1, 0, 1]), 
        "b" => Perturbation(simulation.atoms, c11, c12, dist, [1.0])
    )

    input = SystemPerturbations(g1, 5, g2, 4, Dict)

    res = reweight(simulation,
                    input;
                    all_distances = false,
                    k = 1.0,
                    T = 1.0,
                    cutoff = 15.0,
            )

    @test res["a"].energies ≈ [
        -[6.568418, 0, 6.064767, 0, 0, 0, 9.888358, 0, 7.498538, 0],
        zeros(10),
        [6.568418, 0, 6.064767, 0, 0, 0, 9.888358, 0, 7.498538, 0],
    ] atol = 1.e-3

    @test res["a"].probabilities ≈ [
        [0.031440, 4.41428e-5, 0.0190000, 4.41428e-5, 4.41428e-5, 4.41428e-5, 0.869599, 4.41428e-5, 0.079695, 4.41428e-5],
        ones(10)/10,
        [0.000234, 0.166546, 0.000387, 0.166546, 0.166546, 0.166546, 8.4542e-6, 0.166546, 9.224899e-5, 0.166546],
    ] atol = 1.e-5

    @test res["a"].distances ≈ [
        1, 0, 1, 0, 0, 0, 1, 0, 1, 0,
    ] atol = 1.e-1

    @test res["b"].energies[1] ≈ [
        6.568418,
        5.697282,
        6.064767,
        0,
        6.03079,
        6.464387,
        9.888358,
        9.800954,
        7.498538,
        6.192204      
    ] atol = 1.e-4
end

@testitem "Reweight with small trajectory using all distances and contributions from different residues" begin
    import PDBTools
    import OrderedCollections
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory.xtc")

    i1 = PDBTools.selindex(atoms(simulation), "resname TFE and name O")

    i2 = PDBTools.selindex(atoms(simulation), "residue 11")

    sum_of_dist = reweight(simulation, (i,j,r) -> r, [i1[239]], i2; cutoff = 25.0)
    @test sum_of_dist.energy ≈ [7.4295543149]
end

@testitem "Reweight with small trajectory using minimum distances and contributions from one group" begin #ADICIONAR TESTE COM CONTRIBUIÇÕES NÃO REPETIDAS (OXIGENIO E HIDROGENIOS, POR EX.)
    import PDBTools
    import OrderedCollections
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory.xtc")

    i1 = PDBTools.selindex(atoms(simulation), "index 97 or index 106")

    i2 = PDBTools.selindex(atoms(simulation), "residue 15 and name HB3")

    sum_of_dist = reweight(simulation, (i,j,r) -> r, i1, i2, cutoff = 25.0)
    @test sum_of_dist.energy ≈ [
        1.773896547670759, 1.5923698293115915, 1.716614676290554, 
        1.933003841107648, 1.602329229247863, 1.9639005665480983, 
        3.573986006775934, 2.188798265022823, 2.066180657974777, 
        1.6845109623700647
    ]
end

@testitem "Reweight with small trajectory using all distances and contributions from one group" begin
    import PDBTools
    import OrderedCollections
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory.xtc")

    i1 = PDBTools.selindex(atoms(simulation), "resname TFE and name O")

    i2 = PDBTools.selindex(atoms(simulation), "protein and name O")

    α = 5.e-3

    β = 5.e-3

    cut_off = 12.0

    probs_test = reweight(simulation, (i,j,r) -> gaussian_decay_perturbation(r, α, β), i1, i2; cutoff = cut_off)
    @test probs_test.probability ≈ [
        0.08987791339898044
        0.07326337222373071
        0.0973116226496827
        0.10965810145525891
        0.09829891590498603
        0.0916792371461855
        0.08548699059703141
        0.12480704633057726
        0.09973413264337352
        0.12988266765019355
    ]
end