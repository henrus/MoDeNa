!> @file
!! handles input and output
!! @author    Pavel Ferkl
!! @ingroup   bblgr
module in_out
    use foaming_globals_m
    use constants
    use globals
    use ioutils, only:newunit,str
    implicit none
    private
    public set_paths,read_inputs,save_integration_header,&
        save_integration_step,save_integration_close,load_old_results
contains
!********************************BEGINNING*************************************
!> set paths to all files
subroutine set_paths
    fileplacein='../'
    fileplaceout='../results/'
    inputs='unifiedInput.json'
    outputs_1d='outputs_1d.out'
    outputs_GR='outputs_GR.out'
    outputs_c='outputs_c.out'
    outputs_kin='kinetics.out'
    inputs=TRIM(ADJUSTL(fileplacein))//TRIM(ADJUSTL(inputs))
    outputs_1d=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs_1d))
    outputs_GR=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs_GR))
    outputs_c=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs_c))
    outputs_kin=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs_kin))
end subroutine set_paths
!***********************************END****************************************
!> reads input values from a file
!! save them to global variables
subroutine read_inputs
    use fson
    use fson_value_m, only: fson_value_get
    character(len=1024) :: strval
    type(fson_value), pointer :: json_data
    write(*,*) 'loading input file ',TRIM(inputs)
    json_data => fson_parse('../unifiedInput.json')
    call fson_get(json_data, "bubbleGrowth.integrator", integrator)
    call fson_get(json_data, "bubbleGrowth.method", int_meth)
    call fson_get(json_data, "bubbleGrowth.inertialTerm", inertial_term)
    call fson_get(json_data, "bubbleGrowth.solubilityCorrection", solcorr)
    call fson_get(json_data, "bubbleGrowth.meshCoarseningParameter", mshco)
    call fson_get(json_data, "bubbleGrowth.internalNodes", p)
    call fson_get(json_data, "bubbleGrowth.initialTime", tstart)
    if (firstrun) call fson_get(json_data, "bubbleGrowth.finalTime", tend)
    call fson_get(json_data, "bubbleGrowth.outerTimeSteps", its)
    call fson_get(json_data, "bubbleGrowth.maxInnerTimeSteps", maxts)
    call fson_get(json_data, "bubbleGrowth.relativeTolerance", rel_tol)
    call fson_get(json_data, "bubbleGrowth.absoluteTolerance", abs_tol)
    call fson_get(json_data, "bubbleGrowth.numberOfDissolvedGases", ngas)
    allocate(D(ngas),cbl(ngas),xgas(ngas+1),KH(ngas),Mbl(ngas),&
        dHv(ngas),mb(ngas),mb2(ngas),mb3(ngas),avconc(ngas),pressure(ngas),&
        diff_model(ngas),sol_model(ngas),cpblg(ngas),cpbll(ngas),&
        wblpol(ngas),D0(ngas))
    call fson_get(json_data, "bubbleGrowth.carbonDioxidePosition", co2_pos)
    call fson_get(json_data, "physicalProperties.pressure", pamb)
    call fson_get(json_data, "physicalProperties.blowingAgents.PBL.molarMass", Mbl(1))
    call fson_get(json_data, "physicalProperties.blowingAgents.CO2.molarMass", Mbl(2))
    call fson_get(json_data, "physicalProperties.polymer.heatCapacity", cppol)
    call fson_get(json_data, "physicalProperties.blowingAgents.PBL.heatCapacityInLiquidPhase", cpbll(1))
    call fson_get(json_data, "physicalProperties.blowingAgents.CO2.heatCapacityInLiquidPhase", cpbll(2))
    call fson_get(json_data, "physicalProperties.blowingAgents.PBL.heatCapacityInGaseousPhase", cpblg(1))
    call fson_get(json_data, "physicalProperties.blowingAgents.CO2.heatCapacityInGaseousPhase", cpblg(2))
    call fson_get(json_data, "physicalProperties.blowingAgents.PBL.evaporationHeat", dHv(1))
    call fson_get(json_data, "physicalProperties.blowingAgents.CO2.evaporationHeat", dHv(2))
    call fson_get(json_data, "physicalProperties.blowingAgents.PBL.density", rhobl)
    call fson_get(json_data, "initialConditions.temperature", temp0)
    call fson_get(json_data, "initialConditions.bubbleRadius", R0)
    call fson_get(json_data, "initialConditions.numberBubbleDensity", nb0)
    call fson_get(json_data, "initialConditions.concentrations.polyol", OH0)
    call fson_get(json_data, "initialConditions.concentrations.water", W0)
    call fson_get(json_data, "initialConditions.concentrations.isocyanate", NCO0)
    call fson_get(json_data, "initialConditions.concentrations.blowingAgents.PBL", cbl(1))
    call fson_get(json_data, "initialConditions.concentrations.blowingAgents.CO2", cbl(2))
    call fson_get(json_data, "kinetics.kineticModel", strval)
    if (strval=='Baser') then
        kin_model=1
    elseif (strval=='BaserRx') then
        kin_model=3
    elseif (strval=='modena') then
        kin_model=4
    else
        stop 'unknown kinetic model'
    endif
    call fson_get(json_data, "kinetics.useDilution", dilution)
    if (kin_model==1 .or. kin_model==3) then
        call fson_get(json_data, "kinetics.gellingReaction.frequentialFactor", AOH)
        call fson_get(json_data, "kinetics.gellingReaction.activationEnergy", EOH)
        call fson_get(json_data, "kinetics.blowingReaction.frequentialFactor", AW)
        call fson_get(json_data, "kinetics.blowingReaction.activationEnergy", EW)
    endif
    if (kin_model==1 .or. kin_model==3 .or. kin_model==4) then
        call fson_get(json_data, "kinetics.gellingReaction.reactionEnthalpy", dHOH)
        call fson_get(json_data, "kinetics.blowingReaction.reactionEnthalpy", dHW)
    endif
    call fson_get(json_data, "physicalProperties.polymer.polymerDensityModel", rhop_model)
    if (rhop_model==1) then
        call fson_get(json_data, "physicalProperties.polymer.density", rhop)
    elseif (rhop_model==2) then
    else
        stop 'unknown polymer density model'
    endif
    call fson_get(json_data, "physicalProperties.surfaceTensionModel", itens_model)
    if (itens_model==1) then
        call fson_get(json_data, "physicalProperties.surfaceTension", sigma)
    elseif (itens_model==2) then
    else
        stop 'unknown interfacial tension model'
    endif
    call fson_get(json_data, "physicalProperties.blowingAgents.PBL.diffusivityModel", diff_model(1))
    if (diff_model(1)==1) then
        call fson_get(json_data, "physicalProperties.blowingAgents.PBL.diffusivity", D(1))
    elseif (diff_model(1)==2) then
    else
        stop 'unknown diffusivity model'
    endif
    call fson_get(json_data, "physicalProperties.blowingAgents.CO2.diffusivityModel", diff_model(2))
    if (diff_model(2)==1) then
        call fson_get(json_data, "physicalProperties.blowingAgents.CO2.diffusivity", D(2))
    elseif (diff_model(2)==2) then
    else
        stop 'unknown diffusivity model'
    endif
    call fson_get(json_data, "physicalProperties.blowingAgents.PBL.solubilityModel", sol_model(1))
    if (sol_model(1)==1) then
        call fson_get(json_data, "physicalProperties.blowingAgents.PBL.solubility", KH(1))
    elseif (sol_model(1)==2) then
    elseif (sol_model(1)==3) then
    elseif (sol_model(1)==4) then
    elseif (sol_model(1)==5) then
    elseif (sol_model(1)==6) then
    else
        stop 'unknown solubility model'
    endif
    call fson_get(json_data, "physicalProperties.blowingAgents.CO2.solubilityModel", sol_model(2))
    if (sol_model(2)==1) then
        call fson_get(json_data, "physicalProperties.blowingAgents.CO2.solubility", KH(2))
    elseif (sol_model(2)==2) then
    else
        stop 'unknown solubility model'
    endif
    call fson_get(json_data, "physicalProperties.polymer.viscosityModel", visc_model)
    if (visc_model==1) then
        call fson_get(json_data, "physicalProperties.polymer.viscosity", eta)
    elseif (visc_model==2) then
    elseif (visc_model==3) then
    else
        stop 'unknown viscosity model'
    endif
    call fson_get(json_data, "physicalProperties.polymer.maxViscosity", maxeta)
    write(*,*) 'done: inputs loaded'
    write(*,*)
end subroutine read_inputs
!***********************************END****************************************


!********************************BEGINNING*************************************
!> opens output files and writes a header
subroutine save_integration_header
    open (unit=newunit(fi1), file = outputs_1d)
    write(fi1,'(1000A23)') '#time', 'radius','pressure1', 'pressure2',&
        'conversion_of_polyol',&
        'conversion_of_water', 'eq.concentration', 'first_concentration', &
        'viscosity', 'moles_in_polymer', 'moles_in_bubble', 'total_moles', &
        'shell_thickness', 'temperature', 'foam_density', 'weight_fraction1', &
        'weight_fraction2','porosity'
    open (unit=newunit(fi2), file = outputs_GR)
    write(fi2,'(1000A23)') '#GrowthRate1', 'GrowthRate2', 'temperature', &
        'bubbleRadius','c1','c2','p1','p2'
    if (kin_model==4) then
        open (unit=newunit(fi3), file = outputs_kin)
        write(fi3,'(1000A23)') "time","Catalyst_1","CE_A0","CE_A1","CE_B",&
            "CE_B2","CE_I0","CE_I1","CE_I2","CE_PBA","CE_Breac","CE_Areac0",&
            "CE_Areac1","CE_Ireac0","CE_Ireac1","CE_Ireac2","Bulk","R_1",&
            "R_1_mass","R_1_temp","R_1_vol"
    endif
    open (unit=newunit(fi4), file = outputs_c)
end subroutine save_integration_header
!***********************************END****************************************


!********************************BEGINNING*************************************
!> writes an integration step to output file
subroutine save_integration_step(iout)
    integer :: i,iout
    write(fi1,"(1000es23.15)") time,radius,pressure,Y(xOHeq),Y(xWeq),&
        eqconc,Y(fceq),eta,mb(1),mb2(1),mb3(1),st,temp,rhofoam,&
        wblpol,porosity
    write(fi2,"(1000es23.15)") grrate, temp, radius, avconc, pressure
    if (kin_model==4) then
        write(fi3,"(1000es23.15)") time,Y(kineq(1)),Y(kineq(2)),Y(kineq(3)),&
            Y(kineq(4)),Y(kineq(5)),Y(kineq(6)),Y(kineq(7)),Y(kineq(8)),&
            Y(kineq(9)),Y(kineq(10)),Y(kineq(11)),Y(kineq(12)),Y(kineq(13)),&
            Y(kineq(14)),Y(kineq(15)),Y(kineq(16)),Y(kineq(17)),Y(kineq(18)),&
            Y(kineq(19)),Y(kineq(20))
    endif
    write(fi4,"(1000es23.15)") (Y(fceq+i+1),i=0,ngas*p,ngas)
    ! save arrays, which are preserved for future use
    etat(iout,1)=time
    etat(iout,2)=eta
    port(iout,1)=time
    port(iout,2)=porosity
    init_bub_rad(iout,1)=time
    init_bub_rad(iout,2)=radius
end subroutine save_integration_step
!***********************************END****************************************


!********************************BEGINNING*************************************
!> closes output files
subroutine save_integration_close(iout)
    integer :: iout
    real(dp), dimension(:,:), allocatable :: matr
    close(fi1)
    close(fi2)
    close(fi3)
    if (kin_model==4) then
        close(fi4)
    endif
    ! reallocate matrices for eta_rm and bub_vf functions
    ! interpolation doesn't work otherwise
    if (iout /= its) then
        allocate(matr(0:its,2))
        matr=etat
        deallocate(etat)
        allocate(etat(0:iout,2))
        etat=matr(0:iout,:)
        matr=port
        deallocate(port)
        allocate(port(0:iout,2))
        port=matr(0:iout,:)
        matr=init_bub_rad
        deallocate(init_bub_rad)
        allocate(init_bub_rad(0:iout,2))
        init_bub_rad=matr(0:iout,:)
    endif
end subroutine save_integration_close
!***********************************END****************************************


!********************************BEGINNING*************************************
!> loads old results
subroutine load_old_results
    integer :: i,j,ios
    real(dp), dimension(:,:), allocatable :: matrix
    j=0
    open(newunit(fi5),file=outputs_1d)
        do  !find number of points
            read(fi5,*,iostat=ios)
            if (ios/=0) exit
            j=j+1
        enddo
        allocate(matrix(j,18))
        rewind(fi5)
        read(fi5,*)
        do i=2,j
            read(fi5,*) matrix(i,:)
        enddo
    close(fi5)
    allocate(bub_rad(j-1,2))
    bub_rad(:,1)=matrix(2:j,1)
    bub_rad(:,2)=matrix(2:j,2)
    deallocate(matrix)
end subroutine load_old_results
!***********************************END****************************************
end module in_out
