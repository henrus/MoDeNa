/** @file modenaCalls.h
    @brief instantiates the surrogate models
*/
#ifndef MODENACALLS_H
#define MODENACALLS_H
/*
instantiate the surrogate models:
    - bubbleGrowth1,
    - bubbleGrowth2,
    - density_reaction_mixture,
    - rheology,
    - simpleKinetics
*/


for(label pid=0; pid < Pstream::nProcs();++pid) {
  if(Pstream::myProcNo()==pid) {
    Pout << "Allocating modena models for process " << Pstream::myProcNo() << "..." << endl;
    kinetics = autoPtr<modenaScalarField>(new modenaScalarField("RF-1-public"));
    Pout << " finished. " << endl;
  }

  label tmp = Pstream::myProcNo();
  reduce(tmp,sumOp<label>());
}

kineticTime_Pos = kinetics().model().inputs_argPos("'kineticTime'");

kinetics().registerInputField("'Catalyst_1'", Catalyst_1);
kinetics().registerInputField("'CE_A0'", CE_A0);
kinetics().registerInputField("'CE_A1'", CE_A1);
kinetics().registerInputField("'CE_B'", CE_B);
kinetics().registerInputField("'CE_B2'", CE_B2);
kinetics().registerInputField("'CE_I0'", CE_I0);
kinetics().registerInputField("'CE_I1'", CE_I1);
kinetics().registerInputField("'CE_I2'", CE_I2);
kinetics().registerInputField("'CE_PBA'", CE_PBA);
kinetics().registerInputField("'CE_Breac'", CE_Breac);
kinetics().registerInputField("'CE_Areac0'", CE_Areac0);
kinetics().registerInputField("'CE_Areac1'", CE_Areac1);
kinetics().registerInputField("'CE_Ireac0'", CE_Ireac0);
kinetics().registerInputField("'CE_Ireac1'", CE_Ireac1);
kinetics().registerInputField("'CE_Ireac2'", CE_Ireac2);
kinetics().registerInputField("'Bulk'", Bulk);
kinetics().registerInputField("'R_1'", R_1);
kinetics().registerInputField("'R_1_mass'", R_1_mass);
volScalarField R_1_temp_RF1_Celsius
(
    // field has WRONG dimensions
    //R_1_temp_RF1 - dimensionedScalar("offset", dimTemperature, 273.15)
    R_1_temp_RF1 - 273.15
);
kinetics().registerInputField("'R_1_temp'", R_1_temp_RF1_Celsius);
kinetics().registerInputField("'R_1_vol'", R_1_vol_RF1);

kinetics().registerOutputField("source_Catalyst_1", source_Catalyst_1);
kinetics().registerOutputField("source_CE_A0", source_CE_A0);
kinetics().registerOutputField("source_CE_B", source_CE_B);
kinetics().registerOutputField("source_CE_B2", source_CE_B2);
kinetics().registerOutputField("source_CE_I0", source_CE_I0);
kinetics().registerOutputField("source_CE_I1", source_CE_I1);
kinetics().registerOutputField("source_CE_I2", source_CE_I2);
kinetics().registerOutputField("source_CE_PBA", source_CE_PBA);
kinetics().registerOutputField("source_CE_Breac", source_CE_Breac);
kinetics().registerOutputField("source_CE_Areac0", source_CE_Areac0);
kinetics().registerOutputField("source_CE_Areac1", source_CE_Areac1);
kinetics().registerOutputField("source_CE_Ireac0", source_CE_Ireac0);
kinetics().registerOutputField("source_CE_Ireac1", source_CE_Ireac1);
kinetics().registerOutputField("source_CE_Ireac2", source_CE_Ireac2);
kinetics().registerOutputField("source_Bulk", source_Bulk);
kinetics().registerOutputField("source_R_1", source_R_1);
kinetics().registerOutputField("source_R_1_mass", source_R_1_mass);
kinetics().registerOutputField("source_R_1_temp", source_R_1_temp_RF1);
kinetics().registerOutputField("source_R_1_vol", source_R_1_vol_RF1);

kinetics().model().argPos_check();

// strut contents argPos
rho_foam_Pos = strutContentmodel.inputs_argPos("rho");
strut_content_Pos = strutContentmodel.outputs_argPos("fs");
strutContentmodel.argPos_check();

// thermal conductivity argPos
strut_c_Pos = thermalConductivitymodel.inputs_argPos("eps");
cell_size_Pos = thermalConductivitymodel.inputs_argPos("dcell");
strut_c_Pos = thermalConductivitymodel.inputs_argPos("fstrut");
temp_Pos = thermalConductivitymodel.inputs_argPos("T");
X_CO2_Pos = thermalConductivitymodel.inputs_argPos("x[CO2]");
X_O2_Pos = thermalConductivitymodel.inputs_argPos("x[O2]");
X_N2_Pos = thermalConductivitymodel.inputs_argPos("x[N2]");
X_Cyp_Pos = thermalConductivitymodel.inputs_argPos("x[CyP]");
thermalConductivitymodel.argPos_check();

temp_rheopos = rheologymodel.inputs_argPos("T");
shear_rheopos = rheologymodel.inputs_argPos("shear");
conv_rheopos = rheologymodel.inputs_argPos("X");
m0_rheopos = rheologymodel.inputs_argPos("m0");
m1_rheopos = rheologymodel.inputs_argPos("m1");
//rheologymodel.argPos_check();
#endif
