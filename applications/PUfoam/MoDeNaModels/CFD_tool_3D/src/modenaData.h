/** @file modenaCalls.h
	@brief Allocate memory for the surrogate models
*/
#ifndef MODENADATA_H
#define MODENADATA_H

Modena::modenaModel bblgr1("bubbleGrowth1");
Modena::modenaModel bblgr2("bubbleGrowth2");

autoPtr<modenaScalarField> kinetics;

Modena::modenaModel density_reaction_mixturemodel("density_reaction_mixture");

Modena::modenaModel rheologymodel("Rheology_Arrhenius");

Modena::modenaModel strutContentmodel("strutContent");

Modena::modenaModel thermalConductivitymodel("foamConductivity");

size_t kineticTime_Pos;

// strut content
size_t rho_foam_Pos;
size_t strut_content_Pos;

// thermal conductivity
size_t porosity_Pos;
size_t cell_size_Pos;
size_t strut_c_Pos;
size_t temp_Pos;
size_t X_CO2_Pos;
size_t X_Cyp_Pos;
size_t X_O2_Pos;
size_t X_N2_Pos;

// rheology MoDeNa
size_t temp_rheopos;
size_t shear_rheopos;
size_t conv_rheopos;
size_t m0_rheopos;
size_t m1_rheopos;
#endif
