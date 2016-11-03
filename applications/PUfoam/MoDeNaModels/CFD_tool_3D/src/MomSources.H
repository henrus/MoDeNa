if (PBESwitch)
{
    if
    (
        PBEMethod == "QMOM"
     && nNodes == 2
    )
    {
        static int nodes = nNodes;
        int mOrder[2*nodes];
        for (int i = 0; i < 2*nNodes; i++){mOrder[i] = i;}
        double mom[2*nodes], we[nodes], vi[nodes],
        sgBA[2*nodes], sgCO2[2*nodes], sc[2*nodes];
        double sources[2*nodes];

        double radiusGrowthBA[nodes],
               radiusGrowthCO2[nodes],
               volumeGrowthBA[nodes],
               volumeGrowthCO2[nodes],
               nodeRadii[nodes];

        double CO2_l_val, L_l_val, tmp_val, alphafoam;
        double cc1_val = 0.0;
        double LMAX, xBL_val, rhoPolyS;
        // Calling the model for density reaction mixture
        size_t Td_pos = density_reaction_mixturemodel.inputs_argPos("T");
        size_t XOHd_pos = density_reaction_mixturemodel.inputs_argPos("XOH");
        density_reaction_mixturemodel.argPos_check();

        // Calling the bubble surrogate models
        size_t Tbblgr1pos = bblgr1.inputs_argPos("T");
        size_t Rbblgr1pos = bblgr1.inputs_argPos("R");
        // size_t KH1bblgr1pos = bblgr2.inputs_argPos("kH");
        size_t c_1bblgr1pos = bblgr1.inputs_argPos("c");
        size_t p_1bblgr1pos = bblgr1.inputs_argPos("p");
        bblgr1.argPos_check();

        size_t Tbblgr2pos = bblgr2.inputs_argPos("T");
        size_t Rbblgr2pos = bblgr2.inputs_argPos("R");
        // size_t KH2bblgr2pos = bblgr2.inputs_argPos("kH");
        size_t c_2bblgr2pos = bblgr2.inputs_argPos("c");
        size_t p_2bblgr2pos = bblgr2.inputs_argPos("p");
        bblgr2.argPos_check();

        forAll(mesh.C(), celli)
        {
            CO2_l_val = wCO2_l[celli];
            L_l_val = wBA_l[celli];
            tmp_val = TS[celli];
            alphafoam = alpha2[celli];
            double wCO2_g_val, wBA_g_val,prsr;
            wCO2_g_val= wCO2_g[celli];
            wBA_g_val = wBA_g[celli];
            prsr = p[celli];

            if (liquidMixtureDensitySurrogate)
            {
                // set input vector
                density_reaction_mixturemodel.inputs_set(Td_pos, TS[celli]);
                density_reaction_mixturemodel.inputs_set(XOHd_pos, XOH[celli]);
                density_reaction_mixturemodel.call();
                rhoPolyS = density_reaction_mixturemodel.outputs_get(0);
            }
            else
            {
                rhoPolyS = rhoPoly;
            }

            // Using old moments to compute the new weights and nodes
            if (alphafoam > 0.5)
            {
                mom[0] = mZero[celli];
                mom[1] = mOne[celli];
                mom[2] = mTwo[celli];
                mom[3] = mThree[celli];

                // calling the functions, PDA, growthSource, coalescenceSource
                PDA(we, vi, mom, nodes);
                weight0[celli] = we[0];
                weight1[celli] = we[1];
                node0[celli] = vi[0];
                node1[celli] = vi[1];

                // Required variables for wCO2_Max
                double henry_Coefficient, bubble_radius,
                    partialPressure_CO2, wCO2_Max;
                henry_Coefficient = henryCoefficient(tmp_val);
                bubble_radius = bubbleRadius(mom[0], mom[1]);
                partialPressure_CO2 =
                    partialPressureCO2
                    (
                        M_CO2, M_B, surfaceTension, wCO2_g_val,
                        wBA_g_val, prsr, bubble_radius
                    );
                wCO2_Max =
                    wCO2Max
                    (
                        M_CO2, M_liq, partialPressure_CO2, henry_Coefficient
                    );

                if (blowingAgent == "n-pentane")
                {
                    LMAX  = LliqMax (tmp_val);
                }
                if (blowingAgent == "R-11")
                {
                    xBL_val = xBL(tmp_val, dxdT);
                    LMAX = wBL_D(xBL_val, M_B, M_NCO, L0);
                }
                if (blowingAgent == "no")
                {
                    LMAX = 0.0;
                }

                if (bubbleGrowthSurrogateSwitch)
                {
                    // Required variable for bubble growth surrogate models
                    double p_1 =
                    partialPressureBA
                    (
                        M_B, M_CO2, surfaceTension, wBA_g_val, wCO2_g_val,
                        prsr, bubble_radius
                    );
                    double p_2 = partialPressure_CO2;
                    double c_1 = wBA_l[celli]*rhoPolyS*1000.0/M_B;
                    double c_2 = wCO2_l[celli]*rhoPolyS*1000.0/M_CO2;
                    if (bubbleGrowthMode == "twoNodes")
                    {
                        // Info<< "p_1: " << p_1 << endl;
                        // Info<< "p_2: " << p_2 << endl;
                        // Info<< "c_1: " << c_1 << endl;
                        // Info<< "c_2: " << c_2 << endl;
                        // Info<< "bubble_radius: " << bubble_radius << endl;
                        // Info<< "T: " << TS[celli] << endl;
                        for (int i = 0; i < nodes; i++)
                        {
                            nodeRadii[i] = nodeRadius(vi[i]);
                            // set input vector
                            bblgr1.inputs_set(Tbblgr1pos, TS[celli]);
                            bblgr1.inputs_set(Rbblgr1pos, nodeRadii[i]);
                            // bblgr1.inputs_set(KH1bblgr1pos, KH1);
                            bblgr1.inputs_set(c_1bblgr1pos, c_1);
                            bblgr1.inputs_set(p_1bblgr1pos, p_1);
                            // call the bblgr1 model
                            bblgr1.call();
                            radiusGrowthBA[i] = bblgr1.outputs_get(0);

                            // set input vector
                            bblgr2.inputs_set(Tbblgr2pos, TS[celli]);
                            bblgr2.inputs_set(Rbblgr2pos, nodeRadii[i]);
                            // bblgr2.inputs_set(KH2bblgr2pos, KH2);
                            bblgr2.inputs_set(c_2bblgr2pos, c_2);
                            bblgr2.inputs_set( p_2bblgr2pos, p_2);
                            // call the bblgr2 model
                            bblgr2.call();
                            radiusGrowthCO2[i] = bblgr2.outputs_get(0);

                            volumeGrowthBA[i] =
                                (radiusGrowthBA[i]*RR*TS[celli])/(p_1);
                            volumeGrowthCO2[i] =
                                (radiusGrowthCO2[i]*RR*TS[celli])/(p_2);
                            // Info<< "radiusGrowthBA[" << i << "] " << radiusGrowthBA[i] << endl;
                            // Info<< "radiusGrowthCO2[" << i << "] " << radiusGrowthCO2[i] << endl;
                            // Info<< "volumeGrowthBA[" << i << "] " << volumeGrowthBA[i] << endl;
                            // Info<< "volumeGrowthCO2[" << i << "] " << volumeGrowthCO2[i] << endl;
                            if
                            (
                                volumeGrowthBA[i] < 0.0
                             || radiusGrowthBA[i] < 0.0
                             || L0 < 1.0e-8
                             || wBA_l[celli] < 1.0e-8
                            )
                            {
                                volumeGrowthBA[i] = 0.0;
                            }
                            if
                            (
                                volumeGrowthCO2[i] < 0.0
                             || radiusGrowthCO2[i] < 0.0
                             || CW_0 < 1.0e-8
                             || wCO2_l[celli] < 1.0e-8
                            )
                            {
                                volumeGrowthCO2[i] = 0.0;
                            }
                        }
                    }
                    else if (bubbleGrowthMode == "meanRadius")
                    {
                        // Info<< "p_1: " << p_1 << endl;
                        // Info<< "p_2: " << p_2 << endl;
                        // Info<< "c_1: " << c_1 << endl;
                        // Info<< "c_2: " << c_2 << endl;
                        // Info<< "bubble_radius: " << bubble_radius << endl;
                        // Info<< "T: " << TS[celli] << endl;

                        // set input vector
                        bblgr1.inputs_set(Tbblgr1pos, TS[celli]);
                        bblgr1.inputs_set(Rbblgr1pos, bubble_radius);
                        bblgr1.inputs_set(c_1bblgr1pos, c_1);
                        bblgr1.inputs_set(p_1bblgr1pos, p_1);
                        // set input vector
                        bblgr2.inputs_set(Tbblgr2pos, TS[celli]);
                        bblgr2.inputs_set(Rbblgr2pos, bubble_radius);
                        bblgr2.inputs_set(c_2bblgr2pos, c_2);
                        bblgr2.inputs_set(p_2bblgr2pos, p_2);

                        // call the bblgr1 model
                        bblgr1.call();
                        // call the bblgr2 model
                        bblgr2.call();

                        double G1, G2;
                        for (int i = 0; i < nNodes; i++)
                        {
                            volumeGrowthBA[i] = 0.0;
                            volumeGrowthCO2[i] = 0.0;
                        }
                        G1 = bblgr1.outputs_get(0);
                        G2 = bblgr2.outputs_get(0);
                        // Info<< "G1: " << G1 << endl;
                        // Info<< "G2: " << G2 << endl;

                        volumeGrowthBA[0] = (G1*RR*TS[celli])/(p_1);
                        volumeGrowthCO2[0] = (G2*RR*TS[celli])/(p_2);
                        // Info<< "volumeGrowthBA[0]" << volumeGrowthBA[0] << endl;
                        // Info<< "volumeGrowthCO2[0]" << volumeGrowthCO2[0] << endl;
                        if
                        (
                            volumeGrowthBA[0] < 0.0
                         || G1 < 0.0
                         || L0 < 1.0e-8
                         || wBA_l[celli] < 1.0e-8
                        )
                        {
                            volumeGrowthBA[0] = 0.0;
                        }
                        if
                        (
                            volumeGrowthCO2[0] < 0.0
                         || G2 < 0.0
                         || CW_0 < 1.0e-8
                         || wCO2_l[celli] < 1.0e-8
                        )
                        {
                            volumeGrowthCO2[0] = 0.0;
                        }
                        volumeGrowthBA[1] = volumeGrowthBA[0];
                        volumeGrowthCO2[1] = volumeGrowthCO2[0];
                    }
                    else
                    {
                        FatalErrorIn(args.executable())
                        << "Invalid bubble growth mode " << bubbleGrowthMode
                        << "\nValid modes are: "
                        << "\ntwoNodes,"
                        << "\nmeanRadius"
                        << exit(FatalError);
                    }
                }
                else
                {
                    for (int i = 0; i < nodes; i++)
                    {
                        volumeGrowthBA[i] = growthRateBA(L_l_val,LMAX); // growth rate constant due to blowing agent
                        volumeGrowthCO2[i]= growthRateCO2(CO2_l_val, wCO2_Max); // growth rate constant due to CO2
                    }
                }

                growthSource
                    (
                        sgBA, sgCO2, we, vi, nodes, mOrder, CO2_l_val, L_l_val,
                        tmp_val, wCO2_Max, cc1_val, LMAX,
                        volumeGrowthBA, volumeGrowthCO2
                    );

                g1_CO2[celli] = sgCO2[1];
                g1_BA[celli] = sgBA[1];
                cc1[celli] = cc1_val;

                if (XOH[celli] > XOH_Gel)
                {
                    for (int j = 0; j < 2*nodes; j++)
                    {
                        sc[j] = 0.0;
                    }
                }
                else
                {
                    coalescenceSource(sc, we, vi, nodes, mOrder);
                }

                for (int j = 0; j < 2*nodes; j++)
                {
                    sources[j]  = sgBA[j] + sgCO2[j] + 0.5*sc[j];
                }
                mZeroSource[celli] = sources[0];
                mOneSource[celli] = sources[1];
                mTwoSource[celli] = sources[2];
                mThreeSource[celli] = sources[3];
            }
            else
            {
                mZeroSource[celli] = ROOTVSMALL;
                mOneSource[celli] = ROOTVSMALL;
                mTwoSource[celli] = ROOTVSMALL;
                mThreeSource[celli] = ROOTVSMALL;
            }
        } //end of forAll
    }

    if
    (
        PBEMethod == "QMOM"
     && nNodes == 3
    )
    {
        static int nodes = nNodes;
        int mOrder[2*nodes];
        for (int i = 0; i < 2*nNodes; i++){mOrder[i] = i;}
        double mom[2*nodes], we[nodes], vi[nodes], sgBA[2*nodes],
        sgCO2[2*nodes], sc[2*nodes];
        double sources[2*nodes];

        double CO2_l_val, L_l_val, tmp_val, alphafoam;
        double LMAX, xBL_val, rhoPolyS;
        double cc1_val = 0.0;

        double  radiusGrowthBA[nodes],
                radiusGrowthCO2[nodes],
                volumeGrowthBA[nodes],
                volumeGrowthCO2[nodes],
                nodeRadii[nodes];

        // Calling the model for density reaction mixture
        size_t Td_pos = density_reaction_mixturemodel.inputs_argPos("T");
        size_t XOHd_pos = density_reaction_mixturemodel.inputs_argPos("XOH");
        density_reaction_mixturemodel.argPos_check();

        // Calling the bubble surrogate models
        size_t Tbblgr1pos = bblgr1.inputs_argPos("T");
        size_t Rbblgr1pos = bblgr1.inputs_argPos("R");
        // size_t KH1bblgr1pos = bblgr2.inputs_argPos("kH");
        size_t c_1bblgr1pos = bblgr1.inputs_argPos("c");
        size_t p_1bblgr1pos = bblgr1.inputs_argPos("p");
        bblgr1.argPos_check();

        size_t Tbblgr2pos = bblgr2.inputs_argPos("T");
        size_t Rbblgr2pos = bblgr2.inputs_argPos("R");
        // size_t KH2bblgr2pos = bblgr2.inputs_argPos("kH");
        size_t c_2bblgr2pos = bblgr2.inputs_argPos("c");
        size_t p_2bblgr2pos = bblgr2.inputs_argPos("p");
        bblgr2.argPos_check();

        forAll(mesh.C(), celli)
        {
            CO2_l_val = wCO2_l[celli];
            L_l_val = wBA_l[celli];
            tmp_val = TS[celli];
            alphafoam = alpha2[celli];

            // Required variables for wCO2_Max
            double wCO2_g_val, wBA_g_val,prsr;
            wCO2_g_val= wCO2_g[celli];
            wBA_g_val = wBA_g[celli];
            prsr = p[celli];

            if (liquidMixtureDensitySurrogate)
            {
                // set input vector
                density_reaction_mixturemodel.inputs_set(Td_pos, TS[celli]);
                density_reaction_mixturemodel.inputs_set(XOHd_pos, XOH[celli]);
                density_reaction_mixturemodel.call();
                rhoPolyS = density_reaction_mixturemodel.outputs_get(0);
            }
            else
            {
                rhoPolyS = rhoPoly;
            }

            if (alphafoam > 0.5)
            {
                mom[0] = mZero[celli];
                mom[1] = mOne[celli];
                mom[2] = mTwo[celli];
                mom[3] = mThree[celli];
                mom[4] = mFour[celli];
                mom[5] = mFive[celli];

                // calling the functions, PDA, growthSource, coalescenceSource
                PDA(we, vi, mom, nodes);
                weight0[celli] = we[0];
                weight1[celli] = we[1];
                weight2[celli] = we[2];
                node0[celli] = vi[0];
                node1[celli] = vi[1];
                node2[celli] = vi[2];

                // Required variables for wCO2_Max
                double henry_Coefficient, bubble_radius,
                   partialPressure_CO2, wCO2_Max;
                henry_Coefficient = henryCoefficient(tmp_val);
                bubble_radius = bubbleRadius(mom[0], mom[1]);
                partialPressure_CO2 =
                    partialPressureCO2
                    (
                        M_CO2, M_B, surfaceTension, wCO2_g_val, wBA_g_val,
                        prsr, bubble_radius
                    );
                wCO2_Max =
                    wCO2Max
                    (
                        M_CO2, M_liq, partialPressure_CO2, henry_Coefficient
                    );

                if (blowingAgent == "n-pentane")
                {
                    LMAX = LliqMax (tmp_val);
                }
                if (blowingAgent == "R-11")
                {
                    xBL_val = xBL(tmp_val, dxdT);
                    LMAX = wBL_D(xBL_val, M_B, M_NCO, L0);
                }
                if (blowingAgent == "no")
                {
                    LMAX = 0.0;
                }

                if (bubbleGrowthSurrogateSwitch)
                {
                    // Required variable for bubble growth surrogate models
                    double p_1 =
                        partialPressureBA
                        (
                            M_B, M_CO2, surfaceTension, wBA_g_val,
                            wCO2_g_val, prsr, bubble_radius
                        );
                    double p_2 = partialPressure_CO2;
                    double c_1 = wBA_l[celli]*rhoPolyS*1000.0/M_B;
                    double c_2 = wCO2_l[celli]*rhoPolyS*1000.0/M_CO2;

                    for (int i = 0; i < nodes; i++)
                    {
                        nodeRadii[i] = nodeRadius(vi[i]);
                        // set input vector
                        bblgr1.inputs_set(Tbblgr1pos, TS[celli]);
                        bblgr1.inputs_set(Rbblgr1pos, nodeRadii[i]);
                        // bblgr1.inputs_set(KH1bblgr1pos, KH1);
                        bblgr1.inputs_set(c_1bblgr1pos, c_1);
                        bblgr1.inputs_set(p_1bblgr1pos, p_1);
                        // call the bblgr1 model
                        bblgr1.call();
                        radiusGrowthBA[i] = bblgr1.outputs_get(0);

                        // set input vector
                        bblgr2.inputs_set(Tbblgr2pos, TS[celli]);
                        bblgr2.inputs_set(Rbblgr2pos, nodeRadii[i]);
                        // bblgr2.inputs_set(KH2bblgr2pos, KH2);
                        bblgr2.inputs_set(c_2bblgr2pos, c_2);
                        bblgr2.inputs_set( p_2bblgr2pos, p_2);
                        // call the bblgr2 model
                        bblgr2.call();
                        radiusGrowthCO2[i] = bblgr2.outputs_get(0);

                        volumeGrowthBA[i] =
                            (radiusGrowthBA[i]*RR*TS[celli])/(p_1);
                        volumeGrowthCO2[i] =
                            (radiusGrowthCO2[i]*RR*TS[celli])/(p_2);
                        if
                        (
                            volumeGrowthBA[i] < 0.0
                         || radiusGrowthBA[i] < 0.0
                         || L0 < 1.0e-8
                         || wBA_l[celli] < 1.0e-8
                        )
                        {
                            volumeGrowthBA[i] = 0.0;
                        }
                        if
                        (
                            volumeGrowthCO2[i] < 0.0
                         || radiusGrowthCO2[i] < 0.0
                         || CW_0 < 1.0e-8
                         || wCO2_l[celli] < 1.0e-8
                        )
                        {
                            volumeGrowthCO2[i] = 0.0;
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < nodes; i++)
                    {
                        volumeGrowthBA[i] = growthRateBA(L_l_val,LMAX); // growth rate constant due to blowing agent
                        volumeGrowthCO2[i]= growthRateCO2(CO2_l_val, wCO2_Max); // growth rate constant due to CO2
                    }
                }
                growthSource
                (
                    sgBA, sgCO2, we, vi, nodes, mOrder, CO2_l_val, L_l_val,
                    tmp_val, wCO2_Max, cc1_val, LMAX,
                    volumeGrowthBA, volumeGrowthCO2
                );

                g1_CO2[celli] = sgCO2[1];
                g1_BA[celli] = sgBA[1];
                cc1[celli] = cc1_val;

                if (XOH[celli] > XOH_Gel)
                {
                    for (int j = 0; j < 2*nodes; j++)
                    {
                        sc[j] = 0.0;
                    }
                }
                else
                {
                    coalescenceSource(sc, we, vi, nodes, mOrder);
                }

                for (int j = 0; j < 2*nodes; j++)
                {
                    sources[j] = sgBA[j] + sgCO2[j] + 0.5*sc[j];
                }
                mZeroSource[celli] = sources[0];
                mOneSource[celli] = sources[1];
                mTwoSource[celli] = sources[2];
                mThreeSource[celli] = sources[3];
                mFourSource[celli] = sources[4];
                mFiveSource[celli] = sources[5];
            }
            else
            {
                mZeroSource[celli] = ROOTVSMALL;
                mOneSource[celli] = ROOTVSMALL;
                mTwoSource[celli] = ROOTVSMALL;
                mThreeSource[celli] = ROOTVSMALL;
                mFourSource[celli] = ROOTVSMALL;
                mFiveSource[celli] = ROOTVSMALL;
            }
        } //end of forAll
    }
}
