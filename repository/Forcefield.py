
IntParams = {	 #Sig (ang), Eps(K), Q (e)   #Here are your atom-atom parameters for your forcefield. It's listed as atomtype_molecule or atomtype_forcefield. Check this before you do anything!
	"Zinc":(2.763, 62.34, 0), #From UFF
	"Carbon":(3.898, 47.81, 0), #From Dreiding
	"Hydrogen":(3.195, 7.642, 0), #Dreiding
	"Oxygen":(3.3, 48.11, 0), #Dreiding
	"Oxygen_clust":(3.66, 0.069, 0), #UFF
	"Oxygen_ligand":(3.3, 48.11, 0), #Dreiding
	"Oxygen_DMF":(3.3, 226, -0.5), #DMF parameters taken from 10.1039/C4CP05961A
	"Me_DMF":(3.80, 69, 0.28),
	"Nitrogen_DMF":(3.20, 144, -0.57),
	"Carbon_DMF":(3.70, 47.3, 0.45),
	"Hydrogen_DMF":(2.20, 7.18, 0.06),
	"Carbon_Kamath":(3.41, 68.94, -0.235), #Chlorfoorm parametrs taken from 10.1021/jp0535238
	"Hydrogen_Kamath":(2.81, 10.06, 0.355),
	"Chlorine_Kamath":(3.45, 138.58, -0.04),
	"Nitrogen_TRAPPE":(3.310, 36.0, -0.482),
	"COMN_TRAPPE":(0, 0, 0.964),
        "Methynyl_OPLS": (3.8, 40.26, 0.42),
        "Chlorine_OPLS": (3.47, 150.98, -0.14),
        "Hydrogen_CDP": (2.81, 10.06, -0.0551),
        "Carbon_CDP": (3.41, 68.94, 0.5609),
        "Chlorine_CDP": (3.45, 138.58, -0.1686),
	"Chlorine_UFF": (3.947, 114.128, 0),
        "Carbon_UFF": (3.851, 52.79, 0),
        "Hydrogen_UFF": (2.886, 22.12, 0),
"C_DMFAA":(3.50, 33.20,-0.110), #from 10.1016/j.molliq.2015.03.004
"H_DMFAA":(2.50, 15.16,0.060),
"N_DMFAA":(3.25, 85.52,0.040),
"O_DMFAA":(2.96, 105.61,-0.680),
"Carb_DMFAA":(3.75, 52.80,0.5),
"AldH_DMFAA":(2.5, 15.16,0.00),
    
"C_DMFYang":(3.50,33.20,-0.110), ##from 10.1016/j.molliq.2015.03.004
"H_DMFYang":(2.5,15.16,0.060),
"N_DMFYang":(3.25,85.52,-0.14),
"O_DMFYang":(2.96,105.61,-0.5),
"Carb_DMFYang":(3.75,52.80,0.5),
"AldH_DMFYang":(2.5,15.16,0),

"C_DMFCaleman":(3.5,33.20,-0.110), #from 10.1016/j.molliq.2015.03.004
"H_DMFCaleman":(2.5,15.16,0.06),
"N_DMFCaleman":(3.25,85.52,-0.14),
"O_DMFCaleman":(2.96,105.61,-0.5),
"Carb_DMFCaleman":(3.75,52.80,0.5),
"AldH_DMFCaleman":(2.5,15.16,0),

"Ar_NB":(0, 0, 0), #from TRAPPE, Q values are 2.42 on central dummy atom, -1.21 on pi cloud dummy atoms
"C_NB":(4.5, 15, 0.14),
"CH_NB":(3.74,48, 0),
"N_NB":(3.31, 40., 0.82),
"O_NB":(2.9, 80., -0.48),

"O_diox":(2.39, 488, -0.38),
"C_diox":(3.91, 52.5, 0.19),

"Carb_Ace":(3.82, 40, 0.424),
"C_Ace":(3.75, 98, 0),
"O_Ace":(3.05, 79, -0.424),

"C_EG":(3.95, 46, 0.265),
"O_EG":(3.02, 93, -0.7),
"H_EG":(0,0,0.435),

"Me_EtOH":(3.75, 98, 0),
"CH2_EtOH":(3.95, 46, 0.265),
"O_EtOH":(3.02, 93, -0.7),
"H_EtOH":(0,0,0.435),

"Me_MeOH":(3.75, 98, 0.265),
"O_MeOH":(3.02, 93, -0.7),
"H_MeOH":(0,0,0.435)
}
