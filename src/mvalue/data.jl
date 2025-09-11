const cosolvent_column = Dict(
    "tmao" => 1,
    "TMAO" => 1,  
    "Tmao" => 1,  
    "Sarcosine"=> 2,    
    "Betaine"=> 3,    
    "Proline" => 4,    
    "Sorbitol" => 5,   
    "Sucrose" => 6,
    "Urea" => 7,
    "urea" => 7,
    "UREA" => 7,
    "UreaWrong" => 8,
    "UreaMH" => 9,
)

#
# Data section
#

#=

Amino acid side-chain and peptide backbone unit transfer free energies (cal/mol) from water to 1M osmolyte
Values for Urea GTFE+ values from Table S1 of https://doi.org/10.1021/jp409934q (https://pubs.acs.org/doi/10.1021/jp409934q#_i14)

Only Urea values are available for this model for now.

=#
const tfe_sc_bb_moeser_and_horinek = Dict(
#                TMAO   Sarcosine     Betaine     Proline    Sorbitol    Sucrose         Urea
    "ALA" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,       1.01),
    "PHE" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,     -68.64),
    "LEU" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,     -40.10),
    "ILE" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,     -23.96),
    "VAL" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,      -7.18),
    "PRO" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,      -3.18),
    "MET" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,     -33.87),
    "TRP" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,    -126.90),
    "GLY" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,       0.00),
    "SER" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,      -6.09),
    "THR" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,      -7.62),
    "TYR" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,     -30.61),
    "GLN" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,     -40.34),
    "ASN" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,     -24.32),
    "ASP" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,      18.02),
    "GLU" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,      15.09),
    "HIS" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,     -36.04),
    "HSD" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,     -36.04),
    "HSE" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,     -36.04),
    "LYS" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,      -8.29),
    "ARG" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,      -6.70),
    "CYS" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,       0.00),
    "BB"  => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,        -39),
)

#=

Amino acid side-chain and peptide backbone unit transfer free energies (cal/mol) from water to 1M osmolyte
Supplementary Table 1 of https://doi.org/10.1073/pnas.0507053102

UreaWrong from GTFE* from Supplementary Table S1 of Moeser and Horinek and originally from https://doi.org/10.1073/pnas.0706251104
UreaMH is the data from Moeser and Horinek

The "Urea" column selection points to "UreaWrong" which is what is output from the server.

=#
const tfe_sc_bb_auton_and_bolen = Dict(
#                TMAO   Sarcosine     Betaine     Proline    Sorbitol    Sucrose      UreaWrong   UreaAPP     UreaMH
    "ALA" => ( -14.64,      10.91,       4.77,      -0.07,      16.57,      22.05,       0.63,       -4.69,     1.01),
    "PHE" => (  -9.32,     -12.64,    -112.93,     -71.26,      26.38,     -96.35,     -42.84,      -83.11,   -68.64),
    "LEU" => (  11.62,      38.33,     -17.73,       4.77,      39.07,      37.11,     -14.30,      -54.57,   -40.10),
    "ILE" => ( -25.43,      39.98,      -1.27,      -2.72,      36.90,      28.12,       1.84,      -38.43,   -23.96),
    "VAL" => (  -1.02,      29.32,     -19.63,       7.96,      24.65,      33.92,      18.62,      -21.65,    -7.18),
    "PRO" => (-137.73,     -34.23,    -125.16,     -63.96,      -4.48,     -73.02,      22.62,      -17.65,    -3.18),
    "MET" => (  -7.65,       8.18,     -14.16,     -35.12,      20.97,      -6.66,      -8.07,      -48.34,   -33.87),
    "TRP" => (-152.87,    -113.03,    -369.93,    -198.37,     -67.23,    -215.27,    -101.19,     -141.46,  -126.90),
    "GLY" => (      0,          0,          0,          0,          0,          0,       0.00,           0,     0.00),
    "SER" => ( -39.04,     -27.98,     -41.85,     -33.49,      -1.58,      -2.79,      19.71,      -20.56,    -6.09),
    "THR" => (   3.57,      -7.54,       0.33,     -18.33,      13.20,      20.82,      18.18,      -22.09,    -7.62),
    "TYR" => (-114.32,     -26.37,    -213.09,    -138.41,     -53.50,     -78.41,      -4.81,      -45.08,   -30.61),
    "GLN" => (  41.41,     -10.19,       7.57,     -32.26,     -23.98,     -40.87,     -14.54,      -54.81,   -40.34),
    "ASN" => (  55.69,     -40.93,      33.17,     -17.71,     -21.21,     -28.28,       1.48,      -38.79,   -24.32),
    "ASP" => ( -66.67,     -14.20,    -116.56,     -90.51,     -83.88,     -37.17,      43.82,        3.55,    18.02),
    "GLU" => ( -83.25,     -12.61,    -112.08,     -89.17,     -70.05,     -41.65,      40.89,        0.62,    15.09),
    "HIS" => (  42.07,     -20.80,     -35.97,     -45.10,     -42.45,    -118.66,     -10.24,      -50.51,   -36.04),
    "HSD" => (  42.07,     -20.80,     -35.97,     -45.10,     -42.45,    -118.66,     -10.24,      -50.51,   -36.04),
    "HSE" => (  42.07,     -20.80,     -35.97,     -45.10,     -42.45,    -118.66,     -10.24,      -50.51,   -36.04),
    "LYS" => (-110.23,     -27.42,    -171.99,     -59.87,     -32.47,     -39.60,      17.51,      -22.76,    -8.29),
    "ARG" => (-109.27,     -32.24,    -109.45,     -60.18,     -24.65,     -79.32,      19.10,      -21.17,    -6.70),
    "CYS" => (      0,          0,          0,          0,          0,          0,       0.00,           0,     0.00), # not reported
    "BB"  => (     90,         52,         67,         48,         35,         62,        -39,         -39,      -39),
)

#= 

Isolated ASA values are from the Supporting Table 2 of https://doi.org/10.1073/pnas.0507053102
(https://www.pnas.org/doi/suppl/10.1073/pnas.0507053102/suppl_file/07053table2.pdf)

=#
const isolated_ASA = Dict{String,Tuple{Float64,Float64}}( 
                # BB      SC   (Å^2)
    "ALA"	=> (46.2,	71.9),
    "PHE"	=> (38.4,	184.4),
    "LEU"	=> (35.3,	157.8),
    "ILE"	=> (30.9,	150.1),
    "VAL"	=> (36.1,	128.4),
    "PRO"	=> (35.6,	111.0),
    "MET"	=> (38.6,	164.8),
    "TRP"	=> (37.4,	228.9),
    "GLY"	=> (88.1,	0),
    "SER"	=> (44.0,	85.8),
    "THR"	=> (37.9,	114.6),
    "TYR"	=> (38.7,	198.1),
    "GLN"	=> (37.8,	155.4),
    "ASN"	=> (40.2,	125.3),
    "ASP"	=> (40.5,	118.2),
    "GLU"	=> (37.8,	148.4),
    "HIS"	=> (40.4,	162.1),
    "LYS"	=> (38.7,	187.1),
    "ARG"	=> (39.1,	216.9),
    "CYS"	=> (42.6,	103.5),
) 

#=

Average SASA values from lower and upper bounds of the denatured state ensemble,
reported in Supplementary Table 2 of https://doi.org/10.1073/pnas.0507053102
These are the average values of Table 1 of https://doi.org/10.1021/bi962819o
Can be used for testing, but **do not** provide the same accuracy as the values
calculated with GROMACS or obtained from the m-value server.
=# 

const sasa_desnat_average = Dict(
    "ALA" => Dict(:bb => 27.9,	:sc => 55.1),
    "PHE" => Dict(:bb => 24.3,	:sc => 128.8),
    "LEU" => Dict(:bb => 22.7,	:sc => 109.6),
    "ILE" => Dict(:bb => 20.0,	:sc => 117.1),
    "VAL" => Dict(:bb => 20.4,	:sc => 96.4),
    "PRO" => Dict(:bb => 22.5,	:sc => 87.0),
    "MET" => Dict(:bb => 25.3,	:sc => 122.4),
    "TRP" => Dict(:bb => 23.6,	:sc => 156.6),
    "GLY" => Dict(:bb => 65.2,	:sc => 0.0),
    "SER" => Dict(:bb => 29.4,	:sc => 66.5),
    "THR" => Dict(:bb => 24.1,	:sc => 84.3),
    "TYR" => Dict(:bb => 25.6,	:sc => 141.7),
    "GLN" => Dict(:bb => 25.3,	:sc => 116.9),
    "ASN" => Dict(:bb => 25.2,	:sc => 90.1),
    "ASP" => Dict(:bb => 26.0,	:sc => 87.0),
    "GLU" => Dict(:bb => 25.7,	:sc => 113.4),
    "HIS" => Dict(:bb => 24.2,	:sc => 111.5),
    "LYS" => Dict(:bb => 26.1,	:sc => 150.7),
    "ARG" => Dict(:bb => 25.1,	:sc => 171.1),
    "CYS" => Dict(:bb => 26.4,	:sc => 73.0),
)



