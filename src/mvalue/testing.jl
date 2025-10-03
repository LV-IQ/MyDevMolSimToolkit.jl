@testitem "mvalue" begin
    import PDBTools

    sasa_1MJC_clean = parse_mvalue_server_sasa(
        """
        ALA 		    5 	 (    95.6)    138.1 [   180.6] 	 | 	 (    -4.6)     35.7 [    75.9] 
        PHE 		    6 	 (   430.8)    491.1 [   557.4] 	 | 	 (    65.0)    119.0 [   173.0] 
        LEU 		    2 	 (   131.2)    147.5 [   163.8] 	 | 	 (   -29.8)    -13.8 [     2.2] 
        ILE 		    4 	 (   247.0)    315.0 [   383.0] 	 | 	 (    24.9)     43.9 [    62.9] 
        VAL 		    5 	 (   403.1)    475.9 [   548.6] 	 | 	 (    76.9)     99.4 [   121.9] 
        PRO 		    2 	 (   134.0)    141.0 [   148.0] 	 | 	 (    27.2)     34.4 [    41.6] 
        MET 		    1 	 (    55.6)     72.7 [    89.8] 	 | 	 (    16.1)     24.6 [    33.2] 
        TRP 		    1 	 (    71.5)     73.3 [    75.2] 	 | 	 (    14.5)     22.9 [    31.4] 
        GLY 		   10 	 (     0.0)      0.0 [     0.0] 	 | 	 (   200.9)    306.4 [   411.9] 
        SER 		    7 	 (   110.3)    157.9 [   205.5] 	 | 	 (   -49.0)     -9.8 [    29.4] 
        THR 		    4 	 (   130.9)    158.7 [   186.5] 	 | 	 (    70.0)     91.8 [   113.6] 
        TYR 		    1 	 (    76.4)     87.0 [    97.7] 	 | 	 (    -2.1)      5.8 [    13.7] 
        GLN 		    2 	 (    89.8)    113.5 [   137.2] 	 | 	 (    18.8)     35.0 [    51.2] 
        ASN 		    3 	 (    54.5)     71.1 [    87.8] 	 | 	 (    43.2)     65.9 [    88.5] 
        ASP 		    6 	 (    70.4)    117.2 [   164.0] 	 | 	 (   -53.0)     -5.6 [    41.8] 
        GLU 		    2 	 (    30.2)     51.3 [    72.4] 	 | 	 (    -4.5)     11.1 [    26.7] 
        HIS 		    1 	 (    42.7)     50.3 [    57.9] 	 | 	 (    13.2)     22.4 [    31.7] 
        LYS 		    7 	 (   260.6)    317.6 [   374.7] 	 | 	 (   -10.0)     44.3 [    98.5] 
        ARG 		    0 	 (     0.0)      0.0 [     0.0] 	 | 	 (     0.0)      0.0 [     0.0] 
        CYS 		    0 	 (     0.0)      0.0 [     0.0] 	 | 	 (     0.0)      0.0 [     0.0] 
        """
    )

    sasa_2RN2_clean = parse_mvalue_server_sasa(
        """
        ALA 		   14 	 (   426.2)    545.2 [   664.2] 	 | 	 (   149.9)    262.6 [   375.3] 
        PHE 		    2 	 (   187.5)    207.6 [   229.7] 	 | 	 (    29.9)     47.9 [    65.9] 
        LEU 		   12 	 (   929.1)   1026.9 [  1124.7] 	 | 	 (    91.9)    187.9 [   283.9] 
        ILE 		    7 	 (   603.9)    722.9 [   841.9] 	 | 	 (    55.7)     89.0 [   122.2] 
        VAL 		    9 	 (   518.5)    649.4 [   780.4] 	 | 	 (    76.4)    116.9 [   157.4] 
        PRO 		    5 	 (   161.5)    179.0 [   196.5] 	 | 	 (    51.9)     69.9 [    87.9] 
        MET 		    4 	 (   187.8)    256.2 [   324.6] 	 | 	 (   -29.5)      4.7 [    38.9] 
        TRP 		    6 	 (   716.2)    727.3 [   738.4] 	 | 	 (    29.7)     80.4 [   131.1] 
        GLY 		   14 	 (     0.0)      0.0 [     0.0] 	 | 	 (   439.1)    586.8 [   734.5] 
        SER 		    4 	 (   197.8)    225.0 [   252.2] 	 | 	 (    54.6)     77.0 [    99.4] 
        THR 		   10 	 (   397.0)    466.5 [   536.0] 	 | 	 (    37.0)     91.5 [   146.0] 
        TYR 		    5 	 (   556.9)    610.1 [   663.4] 	 | 	 (    65.4)    104.9 [   144.4] 
        GLN 		    8 	 (   143.6)    238.4 [   333.2] 	 | 	 (    67.8)    132.6 [   197.4] 
        ASN 		    7 	 (   182.6)    221.5 [   260.3] 	 | 	 (    76.7)    129.6 [   182.4] 
        ASP 		    7 	 (   306.3)    360.9 [   415.5] 	 | 	 (    48.7)    104.0 [   159.3] 
        GLU 		   12 	 (   426.2)    552.8 [   679.4] 	 | 	 (    90.4)    184.0 [   277.6] 
        HIS 		    5 	 (    87.1)    125.1 [   163.1] 	 | 	 (   -20.2)     26.0 [    72.3] 
        LYS 		   11 	 (   317.2)    406.8 [   496.5] 	 | 	 (    48.7)    134.0 [   219.2] 
        ARG 		   10 	 (   546.2)    688.2 [   830.2] 	 | 	 (   110.1)    189.6 [   269.1] 
        CYS 		    3 	 (   183.2)    213.4 [   243.5] 	 | 	 (    32.3)     56.7 [    81.2] 
        """
    )

    # Result from Moeser and Horinek https://doi.org/10.1021/jp409934q
    references = Dict(
        #                    total                  bb                   sc         
        "1MJC" => (-0.8470588235294114, -0.44117647058823506, -0.44117647058823506),
        "2RN2" => (-2.3647058823529408, -1.1999999999999997, -1.164705882352941)
    )

    gmx = Sys.which("gmx")
    if isnothing(gmx)
        @warn "gmx executable not available: some tests won't be run"
    end

    dir = @__DIR__

    #
    # 1MJC
    #
    r_1MJC = mvalue(; pdbname=joinpath(dir, "1MJC_native.pdb"), sasas=sasa_1MJC_clean, type=2)
    @test isapprox(r_1MJC.tot, references["1MJC"][1]; rtol=1e-1)
    @test isapprox(r_1MJC.bb, references["1MJC"][2]; rtol=1e-1)
    @test isapprox(r_1MJC.sc, references["1MJC"][3]; rtol=1e-1)

    # with PDBTools.sasa
    sasa_1MJC_julia = delta_sasa_per_restype(;
        native=PDBTools.read_pdb(joinpath(dir, "1MJC_native.pdb"), "protein"),
        desnat=PDBTools.read_pdb(joinpath(dir, "1MJC_straight.pdb"), "protein"),
    )
    r_1MJC = mvalue(; pdbname=joinpath(dir, "1MJC_native.pdb"), sasas=sasa_1MJC_julia)
    @test isapprox(r_1MJC.tot, -0.937; rtol=1e-1)
    @test isapprox(r_1MJC.bb, -0.398; rtol=1e-1)
    @test isapprox(r_1MJC.sc, -0.539; rtol=1e-1)

    if !isnothing(gmx)
        sasa_1MJC = gmx_delta_sasa_per_restype(;
            native_pdb=joinpath(dir, "1MJC_native.pdb"),
            desnat_pdb=joinpath(dir, "1MJC_straight.pdb"),
        )
        r_1MJC = mvalue(; pdbname=joinpath(dir, "1MJC_native.pdb"), sasas=sasa_1MJC)
        @test isapprox(r_1MJC.tot, -1.24; rtol=1e-1)
        @test isapprox(r_1MJC.bb, -0.777; rtol=1e-1)
        @test isapprox(r_1MJC.sc, -0.463; rtol=1e-1)
    end

    #
    # 2RN2
    #
    r_2RN2 = mvalue(; pdbname=joinpath(dir, "2RN2_native.pdb"), sasas=sasa_2RN2_clean, type=2)
    @test isapprox(r_2RN2.tot, references["2RN2"][1]; rtol=1e-1)
    @test isapprox(r_2RN2.bb, references["2RN2"][2]; rtol=1e-1)
    @test isapprox(r_2RN2.sc, references["2RN2"][3]; rtol=1e-1)

    sasa_2RN2_julia = delta_sasa_per_restype(;
        native=PDBTools.read_pdb(joinpath(dir, "2RN2_native.pdb"), "protein"),
        desnat=PDBTools.read_pdb(joinpath(dir, "2RN2_straight.pdb"), "protein"),
    )
    r_2RN2 = mvalue(; pdbname=joinpath(dir, "2RN2_native.pdb"), sasas=sasa_2RN2_julia)
    @test isapprox(r_2RN2.tot, -2.41; rtol=1e-1)
    @test isapprox(r_2RN2.bb, -1.03; rtol=1e-1)
    @test isapprox(r_2RN2.sc, -1.39; rtol=1e-1)

    if !isnothing(gmx)
        sasa_2RN2 = gmx_delta_sasa_per_restype(;
            native_pdb=joinpath(dir, "2RN2_native.pdb"),
            desnat_pdb=joinpath(dir, "2RN2_straight.pdb"),
        )
        r_2RN2 = mvalue(; pdbname=joinpath(dir, "2RN2_native.pdb"), sasas=sasa_2RN2)
        @test isapprox(r_2RN2.tot, -3.17; rtol=1e-1)
        @test isapprox(r_2RN2.bb, -1.94; rtol=1e-1)
        @test isapprox(r_2RN2.sc, -1.23; rtol=1e-1)
    end

    # Results from running the Auton and Bolen m-value server with 1MJC_clean.pdb
    data_mvalue_server_auton_and_bolen_1MJC = Dict(
        "TMAO" => (224, 1160, 2095),
        "Sarcosine" => (406, 1010, 1614),
        "Betaine" => (-502, 76, 650),
        "Proline" => (-200, 226, 649),
        "Sorbitol" => (378, 780, 1183),
        "Sucrose" => (189, 876, 1559),
        "Urea" => (-293, -711, -1132),
    )

    for (cos, dg) in data_mvalue_server_auton_and_bolen_1MJC
        for ig in 1:3
            m = mvalue(model=AutonBolen, cosolvent=cos, pdbname=joinpath(dir, "1MJC_clean.pdb"), sasas=sasa_1MJC_clean, type=ig)
            @test m.tot ≈ 1e-3 * dg[ig] rtol = 0.1
        end
    end

    data_mvalue_server_auton_and_bolen_2RN2 = Dict(
        "TMAO" => (978, 3215, 5453),
        "Sarcosine" => (1410, 2867, 4323),
        "Betaine" => (-1301, 130, 1560),
        "Proline" => (-568, 464, 1496),
        "Sorbitol" => (830, 1757, 2684),
        "Sucrose" => (952, 2578, 4202),
        "Urea" => (-1217, -2226, -3236)
    )

    for (cos, dg) in data_mvalue_server_auton_and_bolen_2RN2
        for ig in 1:3
            m = mvalue(model=AutonBolen, cosolvent=cos, pdbname=joinpath(dir, "2RN2_clean.pdb"), sasas=sasa_2RN2_clean, type=ig)
            @test m.tot ≈ 1e-3 * dg[ig] rtol = 0.1
        end
    end

    # gmx keyword argument tests
    @test_throws ArgumentError gmx_delta_sasa_per_restype(;
        native_pdb=joinpath(dir, "1MJC_native.pdb"),
        desnat_pdb=joinpath(dir, "1MJC_straight.pdb"),
        gmx="/tmp/gmx_fake",
    )
    if !isnothing(gmx)
        @test length(gmx_delta_sasa_per_restype(;
            native_pdb=joinpath(dir, "2RN2_native.pdb"),
            desnat_pdb=joinpath(dir, "2RN2_straight.pdb"),
            gmx=gmx,
        )) == 20
    end

end