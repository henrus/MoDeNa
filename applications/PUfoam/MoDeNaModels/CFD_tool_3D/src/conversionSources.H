if (KineticsModel == "generic")
{
    scalar temp_current,
           temp_old,
           Arhns_current,
           Arhns_old,
           XW_val,
           XOH_valOld,
           BA_dilution,
           XOHSource_Old,
           QW,
           rhoPolyS;

    // call surrogate model for density reaction mixture
    size_t T_pos = density_reaction_mixturemodel.inputs_argPos("T");
    size_t XOH_pos = density_reaction_mixturemodel.inputs_argPos("XOH");
    density_reaction_mixturemodel.argPos_check();

    forAll(mesh.C(), celli)
    {
        temp_current = TS[celli];
        temp_old = TS.oldTime()[celli];

        if (liquidMixtureDensitySurrogate)
        {
            // set input vector
            density_reaction_mixturemodel.inputs_set(T_pos, temp_current);
            density_reaction_mixturemodel.inputs_set(XOH_pos, XOH[celli]);
            // call the model
            density_reaction_mixturemodel.call();
            rhoPolyS = density_reaction_mixturemodel.outputs_get(0);
        }
        else
        {
            rhoPolyS = rhoPoly;
        }

        if (XW[celli] < 0.0)
        {
            XW[celli] = ROOTVSMALL;
        }
        if (XW[celli] > 1.0)
        {
            XW[celli] = ROOTVSMALL;
        }
        QW = QKinW (AW, EW, temp_current, L0, rhoPolyS, rhoBL);
        if (QW < 0.0 || CW_0 == 0.0)
        {
            QW = ROOTVSMALL;        }

        if (XOH[celli] < XOH_Gel)
        {
            XWSource[celli] = QW;
        }

        if (XW[celli] < 0.0)
        {
            XW_val = ROOTVSMALL;
        }
        else if (XW[celli] > 1.0)
        {
            XW_val = ROOTVSMALL;
        }
        else
        {
            XW_val = XW[celli];
        }

        if (XOH[celli] < 0.0)
        {
            XOH[celli] = ROOTVSMALL;
        }
        if (XOH[celli] > 1.0)
        {
            XOH[celli] = ROOTVSMALL;
        }

        if (XOH[celli] < XOH_Gel)
        {
            XOH_valOld = XOH[celli];
            Arhns_current = arrhenius (AOH, EOH, temp_current);
            Arhns_old = arrhenius (AOH, EOH, temp_old);
            BA_dilution = BAdilution(L0,rhoPolyS,rhoBL);
            XOHSource_Old =
                (
                    (Arhns_old*COH_0*(1 - XOH_valOld)*
                    (CNCO_0/COH_0 - 2*XW_val*CW_0/COH_0 - XOH_valOld)*BA_dilution)
                );
            if (XOHSource_Old < 0.0)
            {
                ScNonXOH[celli] = ROOTVSMALL;
                Coeff_1stOrd[celli] = ROOTVSMALL;
                SuXOH[celli] = ROOTVSMALL;
                SpXOH[celli] = ROOTVSMALL;
            }
            else
            {
                ScNonXOH[celli] =
                    (
                        BA_dilution*
                        (Arhns_current*CNCO_0 - 2.0*Arhns_current*CW_0*XW_val)
                    );
                Coeff_1stOrd[celli] =
                (
                    BA_dilution*
                    (Arhns_current*COH_0 + Arhns_current*CNCO_0
                    - 2.0*Arhns_current*CW_0*XW_val)
                );

                SuXOH[celli] =
                (
                    -1.0*BA_dilution*Arhns_old*
                    COH_0*Foam::pow(XOH_valOld,2.0)
                );
                SpXOH[celli] =
                (
                    2.0*BA_dilution*Arhns_old*COH_0*XOH_valOld
                );
            }
        }
    }
}
if (KineticsModel == "RF-1")
{
    kinetics().model().inputs_set(kineticTime_Pos, runTime.value());

    // call the model
    kinetics().call();

    forAll(mesh.C(), celli)
    {
        if (CW_0 > scalar(1e-8)) {
            XWSource[celli] = scalar(-1.0)*
                source_CE_B2[celli]/(CW_0/scalar(1000.0));
        } else {
            XWSource[celli] = scalar(0.0);
        }

        XOHSource[celli] = scalar(-1.0)*
            (source_CE_A0[celli] + source_CE_A1[celli])/
            ((COH_0))*scalar(1000.0);
    }
}
