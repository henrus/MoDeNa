dimensionedScalar oneTwo
(
    "oneTwo",
    dimDensity,
    1.0
);
rho_gas = 1.2*oneTwo;

// MoDeNa interface to get the liquid mixture density
volScalarField rhoPolySM
(
    IOobject
    (
        "rhoPolySM",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoPolySM", dimless, 1100.0)
);
    // Calling the model for density reaction mixture
    size_t Td_pos = density_reaction_mixturemodel.inputs_argPos("T");
    size_t XOHd_pos = density_reaction_mixturemodel.inputs_argPos("XOH");
    density_reaction_mixturemodel.argPos_check();

forAll(mesh.C(), celli)
{
    if (liquidMixtureDensitySurrogate)
    {
        // set input vector
        density_reaction_mixturemodel.inputs_set(Td_pos, TS[celli]);
        density_reaction_mixturemodel.inputs_set(XOHd_pos, XOH[celli]);
        density_reaction_mixturemodel.call();
        rhoPolySM[celli] = density_reaction_mixturemodel.outputs_get(0);
    }
    else
    {
        rhoPolySM[celli] = rhoPoly;
    }
}
rhoPolySM.correctBoundaryConditions();

if (PBESwitch)
{
    dimensionedScalar dumVar1 ("dumVar1",dimensionSet(1,-1,-2,-1,0,0,0), 1.0);
    dimensionedScalar dumVar2 ("dumVar2",dimensionSet(0,-2,2,1,0,0,0), 1.0);

    forAll(mesh.C(), celli)
    {
        if (wCO2_g[celli] < 0.0)
        {
            wCO2_g[celli] = ROOTVSMALL;
        }
        if (wBA_g[celli] < 0.0)
        {
            wBA_g[celli]   = ROOTVSMALL;
        }
        if (TS[celli] < 0.0)
        {
            TS[celli]     = 298.0;
        }
    }

    volScalarField rho_bubble
    (
        (p/(RR*TS))*((wCO2_g*M_CO2 + wBA_g*M_B)
        /max((scalar(1000.0)*(wCO2_g + wBA_g)),ROOTVSMALL))
    );

    rho_foam =
        (
            dumVar2*(rho_bubble*(mOne/(scalar(1.0) + mOne))
            + dumVar1*((scalar(1.0) + L0)*rhoPolySM
            *(scalar(1.0) - (mOne/(scalar(1.0) + mOne)))))
        );
}
else
{
    dimensionedScalar dumVar1 ("dumVar1",dimensionSet(0,2,-2,-1,0,0,0), 1.0);
    dimensionedScalar dumVar2 ("dumVar2",dimensionSet(-1,3,0,0,0,0,0), 1.0);

    forAll(mesh.C(), celli)
    {
        scalar temp, henry_coeff;
        temp = TS[celli];
        henry_coeff = henryCoefficient(temp);

        scalar wCO2_g_val, wBA_g_val, pr, bubble_radius, mZero_val, mOne_val;

        mZero_val = mOne_val = 0.0;
        bubble_radius = bubbleRadius(mZero_val, mOne_val);

        wCO2_g_val = wCO2_g[celli];
        wBA_g_val = wBA_g[celli];
        pr = p[celli];

        scalar partialPressure_CO2;
        partialPressure_CO2 =
            partialPressureCO2
            (
                M_CO2, M_B, surfaceTension,
                wCO2_g_val, wBA_g_val, pr, bubble_radius
            );

        scalar wCO2_Max;
        wCO2_Max = wCO2Max (M_CO2, M_liq, partialPressure_CO2, henry_coeff);

        CO2[celli] =
        (
            ((CW_0*XW[celli]*M_CO2)/(scalar(1000.0)*rhoPolySM[celli])) - wCO2_Max
        );

        if (CO2[celli] < 0.0)
        {
            CO2[celli] = 0.0;
        }
    }
    rho_foam =
    (
        (scalar(1.0) + L0)/
        (((CO2*scalar(1000.0)*RR*TS)/(p*M_CO2))*dumVar1
      + ((wBA_g*scalar(1000.0)*RR*TS)/(p*M_B))*dumVar1
      + ((wBA_l)/(rhoBL))*dumVar2
      + (scalar(1.0)/rhoPolySM)*dumVar2)
    );
}
