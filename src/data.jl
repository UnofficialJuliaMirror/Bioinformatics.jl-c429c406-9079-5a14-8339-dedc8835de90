const aa_mass = Dict(
    "monoisotopic" => Dict(
        'A' => 71.03711,
        'R' => 156.10111,
        'N' => 114.04293,
        'D' => 115.02694,
        'C' => 103.00919,
        'E' => 129.04259,
        'Q' => 128.05858,
        'G' => 57.02146,
        'H' => 137.05891,
        'I' => 113.08406,
        'L' => 113.08406,
        'K' => 128.09496,
        'M' => 131.04049,
        'F' => 147.06841,
        'P' => 97.05276,
        'S' => 87.03203,
        'T' => 101.04768,
        'W' => 186.07931,
        'Y' => 163.06333,
        'V' => 99.06841
    ),
    "average" => Dict(
        'A' => 71.0788,
        'R' => 156.1875,
        'N' => 114.1038,
        'D' => 115.0886,
        'C' => 103.1388,
        'E' => 129.1155,
        'Q' => 128.1307,
        'G' => 57.0519,
        'H' => 137.1411,
        'I' => 113.1594,
        'L' => 113.1594,
        'K' => 128.1741,
        'M' => 131.1926,
        'F' => 147.1766,
        'P' => 97.1167,
        'S' => 87.0782,
        'T' => 101.1051,
        'W' => 186.2132,
        'Y' => 163.1760,
        'V' => 99.1326
    )
)

const instability_values = Dict(
    "WW" => 1,
    "WC" => 1,
    "WM" => 24.68,
    "WH" => 24.68,
    "WY" => 1,
    "WF" => 1,
    "WQ" => 1,
    "WN" => 13.34,
    "WI" => 1,
    "WR" => 1,
    "WD" => 1,
    "WP" => 1,
    "WT" => -14.03,
    "WK" => 1,
    "WE" => 1,
    "WV" => -7.49,
    "WS" => 1,
    "WG" => -9.37,
    "WA" => -14.03,
    "WL" => 13.34,
    "CW" => 24.68,
    "CC" => 1,
    "CM" => 33.6,
    "CH" => 33.6,
    "CY" => 1,
    "CF" => 1,
    "CQ" => -6.54,
    "CN" => 1,
    "CI" => 1,
    "CR" => 1,
    "CD" => 20.26,
    "CP" => 20.26,
    "CT" => 33.6,
    "CK" => 1,
    "CE" => 1,
    "CV" => -6.54,
    "CS" => 1,
    "CG" => 1,
    "CA" => 1,
    "CL" => 20.26,
    "MW" => 1,
    "MC" => 1,
    "MM" => -1.88,
    "MH" => 58.28,
    "MY" => 24.68,
    "MF" => 1,
    "MQ" => -6.54,
    "MN" => 1,
    "MI" => 1,
    "MR" => -6.54,
    "MD" => 1,
    "MP" => 44.94,
    "MT" => -1.88,
    "MK" => 1,
    "ME" => 1,
    "MV" => 1,
    "MS" => 44.94,
    "MG" => 1,
    "MA" => 13.34,
    "ML" => 1,
    "HW" => -1.88,
    "HC" => 1,
    "HM" => 1,
    "HH" => 1,
    "HY" => 44.94,
    "HF" => -9.37,
    "HQ" => 1,
    "HN" => 24.68,
    "HI" => 44.94,
    "HR" => 1,
    "HD" => 1,
    "HP" => -1.88,
    "HT" => -6.54,
    "HK" => 24.68,
    "HE" => 1,
    "HV" => 1,
    "HS" => 1,
    "HG" => -9.37,
    "HA" => 1,
    "HL" => 1,
    "YW" => -9.37,
    "YC" => 1,
    "YM" => 44.94,
    "YH" => 13.34,
    "YY" => 13.34,
    "YF" => 1,
    "YQ" => 1,
    "YN" => 1,
    "YI" => 1,
    "YR" => -15.91,
    "YD" => 24.68,
    "YP" => 13.34,
    "YT" => -7.49,
    "YK" => 1,
    "YE" => -6.54,
    "YV" => 1,
    "YS" => 1,
    "YG" => -7.49,
    "YA" => 24.68,
    "YL" => 1,
    "FW" => 1,
    "FC" => 1,
    "FM" => 1,
    "FH" => 1,
    "FY" => 33.6,
    "FF" => 1,
    "FQ" => 1,
    "FN" => 1,
    "FI" => 1,
    "FR" => 1,
    "FD" => 13.34,
    "FP" => 20.26,
    "FT" => 1,
    "FK" => -14.03,
    "FE" => 1,
    "FV" => 1,
    "FS" => 1,
    "FG" => 1,
    "FA" => 1,
    "FL" => 1,
    "QW" => 1,
    "QC" => -6.54,
    "QM" => 1,
    "QH" => 1,
    "QY" => -6.54,
    "QF" => -6.54,
    "QQ" => 20.26,
    "QN" => 1,
    "QI" => 1,
    "QR" => 1,
    "QD" => 20.26,
    "QP" => 20.26,
    "QT" => 1,
    "QK" => 1,
    "QE" => 20.26,
    "QV" => -6.54,
    "QS" => 44.94,
    "QG" => 1,
    "QA" => 1,
    "QL" => 1,
    "NW" => -9.37,
    "NC" => -1.88,
    "NM" => 1,
    "NH" => 1,
    "NY" => 1,
    "NF" => -14.03,
    "NQ" => -6.54,
    "NN" => 1,
    "NI" => 44.94,
    "NR" => 1,
    "ND" => 1,
    "NP" => -1.88,
    "NT" => -7.49,
    "NK" => 24.68,
    "NE" => 1,
    "NV" => 1,
    "NS" => 1,
    "NG" => -14.03,
    "NL" => 1,
    "IW" => 1,
    "IC" => 1,
    "IM" => 1,
    "IH" => 13.34,
    "IY" => 1,
    "IF" => 1,
    "IQ" => 1,
    "IN" => 1,
    "II" => 1,
    "IR" => 1,
    "ID" => 1,
    "IP" => -1.88,
    "IT" => 1,
    "IK" => -7.49,
    "IE" => 44.94,
    "IV" => -7.49,
    "IS" => 1,
    "IG" => 1,
    "IA" => 1,
    "IL" => 20.26,
    "RW" => 58.28,
    "RC" => 1,
    "RM" => 1,
    "RH" => 20.26,
    "RY" => -6.54,
    "RF" => 1,
    "RQ" => 20.26,
    "RN" => 13.34,
    "RI" => 1,
    "RR" => 58.28,
    "RD" => 1,
    "RP" => 20.26,
    "RT" => 1,
    "RK" => 1,
    "RE" => 1,
    "RV" => 1,
    "RS" => 44.94,
    "RG" => -7.49,
    "RA" => 1,
    "RL" => 1,
    "DW" => 1,
    "DC" => 1,
    "DM" => 1,
    "DH" => 1,
    "DY" => 1,
    "DF" => -6.54,
    "DQ" => 1,
    "DN" => 1,
    "DI" => 1,
    "DR" => -6.54,
    "DD" => 1,
    "DP" => 1,
    "DT" => -14.03,
    "DK" => -7.49,
    "DE" => 1,
    "DV" => 1,
    "DS" => 20.26,
    "DG" => 1,
    "DA" => 1,
    "DL" => 1,
    "PW" => -1.88,
    "PC" => -6.54,
    "PM" => -6.54,
    "PH" => 1,
    "PY" => 1,
    "PF" => 20.26,
    "PQ" => 20.26,
    "PN" => 1,
    "PI" => 1,
    "PR" => -6.54,
    "PD" => -6.54,
    "PP" => 20.26,
    "PT" => 1,
    "PK" => 1,
    "PE" => 18.38,
    "PV" => 20.26,
    "PS" => 20.26,
    "PG" => 1,
    "PA" => 20.26,
    "PL" => 1,
    "TW" => -14.03,
    "TC" => 1,
    "TM" => 1,
    "TH" => 1,
    "TY" => 1,
    "TF" => 13.34,
    "TQ" => -6.54,
    "TN" => -14.03,
    "TI" => 1,
    "TR" => 1,
    "TD" => 1,
    "TP" => 1,
    "TT" => 1,
    "TK" => 1,
    "TE" => 20.26,
    "TV" => 1,
    "TS" => 1,
    "TG" => -7.49,
    "TA" => 1,
    "TL" => 1,
    "KW" => 1,
    "KC" => 1,
    "KM" => 33.6,
    "KH" => 1,
    "KY" => 1,
    "KF" => 1,
    "KQ" => 24.68,
    "KN" => 1,
    "KI" => -7.49,
    "KR" => 33.6,
    "KD" => 1,
    "KP" => -6.54,
    "KT" => 1,
    "KK" => 1,
    "KE" => 1,
    "KV" => -7.49,
    "KS" => 1,
    "KG" => -7.49,
    "KA" => 1,
    "KL" => -7.49,
    "EW" => -14.03,
    "EC" => 44.94,
    "EM" => 1,
    "EH" => -6.54,
    "EY" => 1,
    "EF" => 1,
    "EQ" => 20.26,
    "EN" => 1,
    "EI" => 20.26,
    "ER" => 1,
    "ED" => 20.26,
    "EP" => 20.26,
    "ET" => 1,
    "EK" => 1,
    "EE" => 33.6,
    "EV" => 1,
    "ES" => 20.26,
    "EG" => 1,
    "EA" => 1,
    "EL" => 1,
    "VW" => 1,
    "VC" => 1,
    "VM" => 1,
    "VH" => 1,
    "VY" => -6.54,
    "VF" => 1,
    "VQ" => 1,
    "VN" => 1,
    "VI" => 1,
    "VR" => 1,
    "VD" => -14.03,
    "VP" => 20.26,
    "VT" => -7.49,
    "VK" => -1.88,
    "VE" => 1,
    "VV" => 1,
    "VS" => 1,
    "VG" => -7.49,
    "VA" => 1,
    "VL" => 1,
    "SW" => 1,
    "SC" => 33.6,
    "SM" => 1,
    "SH" => 1,
    "SY" => 1,
    "SF" => 1,
    "SQ" => 20.26,
    "SN" => 1,
    "SI" => 1,
    "SR" => 20.26,
    "SD" => 1,
    "SP" => 44.94,
    "ST" => 1,
    "SK" => 1,
    "SE" => 20.26,
    "SV" => 1,
    "SS" => 20.26,
    "SG" => 1,
    "SA" => 1,
    "SL" => 1,
    "GW" => 13.34,
    "GC" => 1,
    "GM" => 1,
    "GH" => 1,
    "GY" => -7.49,
    "GF" => 1,
    "GQ" => 1,
    "GN" => -7.49,
    "GI" => -7.49,
    "GR" => 1,
    "GD" => 1,
    "GP" => 1,
    "GT" => -7.49,
    "GK" => -7.49,
    "GE" => -6.54,
    "GV" => 1,
    "GS" => 1,
    "GG" => 13.34,
    "GA" => -7.49,
    "GL" => 1,
    "AW" => 1,
    "AC" => 44.94,
    "AM" => 1,
    "AH" => -7.49,
    "AY" => 1,
    "AF" => 1,
    "AQ" => 1,
    "AN" => 1,
    "AI" => 1,
    "AR" => 1,
    "AD" => -7.49,
    "AP" => 20.26,
    "AT" => 1,
    "AK" => 1,
    "AE" => 1,
    "AV" => 1,
    "AS" => 1,
    "AG" => 1,
    "AA" => 1,
    "AL" => 1,
    "LW" => 24.68,
    "LC" => 1,
    "LM" => 1,
    "LH" => 1,
    "LY" => 1,
    "LF" => 1,
    "LQ" => 33.6,
    "LN" => 1,
    "LI" => 1,
    "LR" => 20.26,
    "LD" => 1,
    "LP" => 20.26,
    "LT" => 1,
    "LK" => -7.49,
    "LE" => 1,
    "LV" => 1,
    "LS" => 1,
    "LG" => 1,
    "LA" => 1,
    "LL" => 1,
    "NA" => 1
)

const hydropathicity = Dict(
    'A' => 1.800,
    'R' => -4.500,
    'N' => -3.500,
    'D' => -3.500,
    'C' => 2.500,
    'E' => -3.500,
    'Q' => -3.500,
    'G' => -0.400,
    'H' => -3.200,
    'I' => 4.500,
    'L' => 3.800,
    'K' => -3.900,
    'M' => 1.900,
    'F' => 2.800,
    'P' => -1.600,
    'S' => -0.800,
    'T' => -0.700,
    'W' => -0.900,
    'Y' => -1.300,
    'V' => 4.200
)
