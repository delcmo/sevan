#
#####################################################
# Define some global parameters used in the blocks. #
#####################################################
#
[GlobalParams]
###### Other parameters #######
order = FIRST

###### Stabilization #######
viscosity_name = ENTROPY
diffusion_name = ENTROPY
Cmax = 0.5
Ce = 1.
Cjump_liquid = 5.
Cjump_gas = 5.
Calpha = 5.
isShock = true
isJumpOn = true

###### Initial Conditions #######
pressure_init_left = 1.e6
pressure_init_right = 1.e5
vel_init_left = 0.
vel_init_right = 0.
#temp_init_left = 70600
#temp_init_right = 70600
#temp_init_left_gas = 8.e5
#temp_init_right_gas = 8.e5
rho_init_left_liq = 10.
rho_init_right_liq = 10.
rho_init_left_gas = 1.
rho_init_right_gas = 1.
alpha_init_left = 0.9
alpha_init_right = 0.1
membrane = 0.5
length = 0.
[]

#############################################################################
#                          USER OBJECTS                                     #
#############################################################################
# Define the user object class that store the EOS parameters.               #
#############################################################################

[UserObjects]
  [./eos_gas]
    type = EquationOfState
  	gamma = 1.4
  	Pinf = 0
    q = 0.
  	Cv = 2.5
  	q_prime = 0.
  [../]

  [./eos_liq]
    type = EquationOfState
    gamma = 3.
    Pinf = 0.
    q = 0.
    Cv = 2.5
    q_prime = 0.
  [../]
  
  [./JumpGradPressLiq]
    type = JumpGradientInterface
    variable = pressure_aux_l
    jump_name = jump_grad_press_aux_l
    execute_on = timestep_begin
  [../]

  [./JumpGradPressGas]
    type = JumpGradientInterface
    variable = pressure_aux_g
    jump_name = jump_grad_press_aux_g
    execute_on = timestep_begin
  [../]

  [./JumpGradDensLiq]
    type = JumpGradientInterface
    variable = density_aux_l
    jump_name = jump_grad_dens_aux_l
    execute_on = timestep_begin
  [../]

  [./JumpGradDensGas]
    type = JumpGradientInterface
    variable = density_aux_g
    jump_name = jump_grad_dens_aux_g
    execute_on = timestep_begin
  [../]

  [./JumpGradPressSmoothLiq]
    type = SmoothFunction
    variable = jump_grad_press_aux_l
    var_name = jump_grad_press_smooth_aux_l
    execute_on = timestep_begin
  [../]

  [./JumpGradDensSmoothLiq]
    type = SmoothFunction
    variable = jump_grad_dens_aux_l
    var_name = jump_grad_dens_smooth_aux_l
    execute_on = timestep_begin
  [../]

  [./JumpGradPressSmoothGas]
    type = SmoothFunction
    variable = jump_grad_press_aux_g
    var_name = jump_grad_press_smooth_aux_g
    execute_on = timestep_begin
  [../]

  [./JumpGradDensSmoothGas]
    type = SmoothFunction
    variable = jump_grad_dens_aux_g
    var_name = jump_grad_dens_smooth_aux_g
    execute_on = timestep_begin
  [../]

  [./JumpGradAlpha]
    type = JumpGradientInterface
    variable = alpha_aux_l
    jump_name = jump_grad_alpha_aux_l
    execute_on = timestep_begin
  [../]
[]

###### Mesh #######
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 250
  ny = 1
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  block_id = '0'
[]

#############################################################################
#                             VARIABLES                                     #
#############################################################################
# Define the variables we want to solve for: l=liquid phase,  g=vapor phase.#
#############################################################################
[Variables]
####### LIQUID PHASE ########
  [./alA_l]
    family = LAGRANGE
    scaling = 1e+0
    [./InitialCondition]
        type = ConservativeVariables1DXIC
        area = area
        eos = eos_liq
    [../]
  [../]

  [./alrhoA_l]
    family = LAGRANGE
    scaling = 1e-2
	[./InitialCondition]
        type = ConservativeVariables1DXIC
        area = area
        eos = eos_liq
	[../]
  [../]

  [./alrhouA_l]
    family = LAGRANGE
    scaling = 1e-4
	[./InitialCondition]
        type = ConservativeVariables1DXIC
        area = area
        eos = eos_liq
	[../]
  [../]

  [./alrhoEA_l]
    family = LAGRANGE
    scaling = 1e-6
	[./InitialCondition]
        type = ConservativeVariables1DXIC
        area = area
        eos = eos_liq
	[../]
  [../]

####### VAPOR PHASE ########
  [./alrhoA_g]
    family = LAGRANGE
    scaling = 1e-2
    [./InitialCondition]
        type = ConservativeVariables1DXIC
        area = area
        eos = eos_gas
        isLiquid = false
    [../]
  [../]

  [./alrhouA_g]
    family = LAGRANGE
    scaling = 1e-4
    [./InitialCondition]
        type = ConservativeVariables1DXIC
        area = area
        eos = eos_gas
        isLiquid = false
    [../]
  [../]

  [./alrhoEA_g]
    family = LAGRANGE
    scaling = 1e-6
    [./InitialCondition]
        type = ConservativeVariables1DXIC
        area = area
        eos = eos_gas
        isLiquid = false
    [../]
  [../]
[]

############################################################################################################
#                                            KERNELS                                                       #
############################################################################################################
# Define the kernels for time dependent, convection and viscosity terms. Same index as for variable block. #
############################################################################################################

[Kernels]
######### Liquid phase ##########
  [./VoidFractionTimeLiq]
    type = EelTimeDerivative
    variable = alA_l
  [../]

  [./MassTimeLiq]
    type = EelTimeDerivative
    variable = alrhoA_l
  [../]

  [./MomTimeLiq]
    type = EelTimeDerivative
    variable = alrhouA_l
  [../]

  [./EnerTimeLiq]
    type = EelTimeDerivative
    variable = alrhoEA_l
  [../]

  [./VoidFracConvLiq]
    type = EelVoidFraction
    variable = alA_l
    pressure_liq = pressure_aux_l
    pressure_gas = pressure_aux_g
    vf_liquid = alpha_aux_l
    area = area_aux
  [../]

  [./MassConvLiq]
    type = EelMass
    variable = alrhoA_l
    alrhouA_x = alrhouA_l
    area = area_aux
  [../]

  [./MomConvLiq]
    type = EelMomentum
    variable = alrhouA_l
    vel_x = velocity_x_aux_l
    vel_x_2 = velocity_x_aux_g
    pressure = pressure_aux_l
    area = area_aux
    vf_liquid = alpha_aux_l
  [../]

  [./EnergyConvLiq]
    type = EelEnergy
    variable = alrhoEA_l
    alrhoA = alrhoA_l
    alrhouA_x = alrhouA_l
    vel_x_2 = velocity_x_aux_g
    pressure_liq = pressure_aux_l
    pressure_gas = pressure_aux_g
    area = area_aux
    vf_liquid = alpha_aux_l
    eos = eos_liq
  [../]

  [./VoidFractionViscLiq]
    type = EelArtificialVisc
    variable = alA_l
    equation_name = VOID_FRACTION
    density = density_aux_l
    pressure = pressure_aux_l
    velocity_x = velocity_x_aux_l
    internal_energy = internal_energy_aux_l
    area = area_aux
    vf_liquid = alpha_aux_l
    eos = eos_liq
  [../]

  [./MassViscLiq]
    type = EelArtificialVisc
    variable = alrhoA_l
    equation_name = CONTINUITY
    density = density_aux_l
    pressure = pressure_aux_l
    velocity_x = velocity_x_aux_l
    internal_energy = internal_energy_aux_l
    area = area_aux
    vf_liquid = alpha_aux_l
    eos = eos_liq
  [../]

  [./MomViscLiq]
    type = EelArtificialVisc
    variable = alrhouA_l
    equation_name = XMOMENTUM
    density = density_aux_l
    pressure = pressure_aux_l
    velocity_x = velocity_x_aux_l
    internal_energy = internal_energy_aux_l
    area = area_aux
    vf_liquid = alpha_aux_l
    eos  =eos_liq
  [../]

  [./EnergyViscLiq]
    type = EelArtificialVisc
    variable = alrhoEA_l
    equation_name = ENERGY
    density = density_aux_l
    pressure = pressure_aux_l
    velocity_x = velocity_x_aux_l
    internal_energy = internal_energy_aux_l
    area = area_aux
    vf_liquid = alpha_aux_l
    eos = eos_liq
  [../]

######### Gas phase ##########
  [./MassTimeGas]
    type = EelTimeDerivative
    variable = alrhoA_g
  [../]

  [./MomTimeGas]
    type = EelTimeDerivative
    variable = alrhouA_g
  [../]

  [./EnerTimeGas]
    type = EelTimeDerivative
    variable = alrhoEA_g
  [../]

  [./MassConvGas]
    type = EelMass
    variable = alrhoA_g
    alrhouA_x = alrhouA_g
    area = area_aux
    isLiquid = false
  [../]

  [./MomConvGas]
    type = EelMomentum
    variable = alrhouA_g
    vel_x = velocity_x_aux_g
    vel_x_2 = velocity_x_aux_l
    pressure = pressure_aux_g
    area = area_aux
    vf_liquid = alpha_aux_l
    isLiquid = false
  [../]

  [./EnergyConvGas]
    type = EelEnergy
    variable = alrhoEA_g
    alrhoA = alrhoA_g
    alrhouA_x = alrhouA_g
    vel_x_2 = velocity_x_aux_l
    pressure_liq = pressure_aux_l
    pressure_gas = pressure_aux_g
    area = area_aux
    vf_liquid = alpha_aux_l
    isLiquid = false
    eos = eos_gas
  [../]

  [./MassViscGas]
    type = EelArtificialVisc
    variable = alrhoA_g
    equation_name = CONTINUITY
    density = density_aux_g
    pressure = pressure_aux_g
    velocity_x = velocity_x_aux_g
    internal_energy = internal_energy_aux_g
    area = area_aux
    vf_liquid = alpha_aux_l
    isLiquid = false
    eos  =eos_gas
  [../]

  [./MomViscGas]
    type = EelArtificialVisc
    variable = alrhouA_g
    equation_name = XMOMENTUM
    density = density_aux_g
    pressure = pressure_aux_g
    velocity_x = velocity_x_aux_g
    internal_energy = internal_energy_aux_g
    area = area_aux
    vf_liquid = alpha_aux_l
    isLiquid = false
    eos = eos_gas
  [../]

  [./EnergyViscGas]
    type = EelArtificialVisc
    variable = alrhoEA_g
    equation_name = ENERGY
    density = density_aux_g
    pressure = pressure_aux_g
    velocity_x = velocity_x_aux_g
    internal_energy = internal_energy_aux_g
    area = area_aux
    vf_liquid = alpha_aux_l
    isLiquid = false
    eos = eos_gas
  [../]
[]

##############################################################################################
#                                       AUXILARY VARIABLES                                   #
##############################################################################################
# Define the auxilary variables                                                              #
##############################################################################################

[AuxVariables]

  [./area_aux]
    family = LAGRANGE
  [../]

  [./alpha_aux_l]
    family = LAGRANGE
  [../]

  [./PI_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

#  [./velI_aux]
#    family = MONOMIAL
#    order = CONSTANT
#  [../]

  [./PI_bar_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

#  [./velI_bar_aux]
#    family = MONOMIAL
#    order = CONSTANT
#  [../]

  [./rhoI_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./tempI_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./P_rel_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./vel_rel_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./Omega_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./Aint_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./ht_liq_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./ht_gas_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./EI_liq_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./EI_gas_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./beta_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./beta_max_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./jump_grad_alpha_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]
######### Liquid phase ##########
   [./velocity_x_aux_l]
      family = LAGRANGE
   [../]

   [./density_aux_l]
      family = LAGRANGE
   [../]

   [./total_energy_aux_l]
      family = LAGRANGE
   [../]

   [./internal_energy_aux_l]
      family = LAGRANGE
   [../]

   [./pressure_aux_l]
      family = LAGRANGE
   [../]

  [./mach_aux_l]
    family = LAGRANGE
  [../]

  [./mu_aux_l]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./mu_max_aux_l]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./kappa_aux_l]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./kappa_max_aux_l]
    family = MONOMIAL
    order = CONSTANT
  [../]
  
  [./jump_grad_press_aux_l]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./jump_grad_alpha_aux_l]
    family = MONOMIAL
    order = CONSTANT
  [../]

 [./jump_grad_dens_aux_l]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./jump_grad_press_smooth_aux_l]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./jump_grad_dens_smooth_aux_l]
    family = MONOMIAL
    order = CONSTANT
  [../]
######### Gas phase ##########
  [./velocity_x_aux_g]
    family = LAGRANGE
  [../]

  [./density_aux_g]
    family = LAGRANGE
  [../]

  [./total_energy_aux_g]
    family = LAGRANGE
  [../]

  [./internal_energy_aux_g]
    family = LAGRANGE
  [../]

  [./pressure_aux_g]
    family = LAGRANGE
  [../]

  [./mach_aux_g]
    family = LAGRANGE
  [../]

  [./mu_aux_g]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./kappa_aux_g]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./mu_max_aux_g]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./kappa_max_aux_g]
    family = MONOMIAL
    order = CONSTANT
  [../]
  
  [./jump_grad_press_aux_g]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./jump_grad_dens_aux_g]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./jump_grad_press_smooth_aux_g]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./jump_grad_dens_smooth_aux_g]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

##############################################################################################
#                                       AUXILARY KERNELS                                     #
##############################################################################################
# Define the auxilary kernels for liquid and gas phases. Same index as for variable block.   #
##############################################################################################

[AuxKernels]

  [./AreaAK]
    type = AreaAux
    variable = area_aux
    area = area
  [../]

  [./VoidFractionAKLiq]
    type = VoidFractionAux
    variable = alpha_aux_l
    alA = alA_l
    area = area_aux
  [../]

  [./PIAK]
    type = MaterialRealAux
    variable = PI_aux
    property = interfacial_pressure
  [../]

#  [./VelIAK]
#    type = MaterialRealAux
#    variable = velI_aux
#    property = interfacial_velocity
#  [../]

  [./PIbarAK]
    type = MaterialRealAux
    variable = PI_bar_aux
    property = average_interfacial_pressure
  [../]

#  [./VelIbarAK]
#    type = MaterialRealAux
#    variable = velI_bar_aux
#    property = average_interfacial_velocity
#  [../]

  [./rhoIAK]
    type = MaterialRealAux
    variable = rhoI_aux
    property = interfacial_density
  [../]

  [./tempIAK]
    type = MaterialRealAux
    variable = tempI_aux
    property = interfacial_temperature
  [../]

  [./PrelIAK]
    type = MaterialRealAux
    variable = P_rel_aux
    property = pressure_relaxation
  [../]

  [./VrelIAK]
    type = MaterialRealAux
    variable = vel_rel_aux
    property = velocity_relaxation
  [../]

  [./EIAKLiq]
    type = MaterialRealAux
    variable = EI_liq_aux
    property = liquid_interfacial_energy
  [../]

  [./EIAKGas]
    type = MaterialRealAux
    variable = EI_gas_aux
    property = gas_interfacial_energy
  [../]

  [./htAKLiq]
    type = MaterialRealAux
    variable = ht_liq_aux
    property = liquid_heat_transfer
  [../]

  [./htAKGas]
    type = MaterialRealAux
    variable = ht_gas_aux
    property = gas_heat_transfer
  [../]

  [./MassAKLiq]
    type = MaterialRealAux
    variable = Omega_aux
    property = mass_transfer
  [../]

  [./AintAKGas]
    type = MaterialRealAux
    variable = Aint_aux
    property = interfacial_area
  [../]

  [./BetaAK]
    type = MaterialRealAux
    variable = beta_aux
    property = beta
  [../]

  [./BetaMaxAK]
    type = MaterialRealAux
    variable = beta_max_aux
    property = beta_max
  [../]
####### Liquid phase ##########
  [./VelAKLiq]
    type = VelocityAux
    variable = velocity_x_aux_l
    alrhoA = alrhoA_l
    alrhouA = alrhouA_l
  [../]

  [./DensAKLiq]
    type = DensityAux
    variable = density_aux_l
    alrhoA = alrhoA_l
    vf_liquid = alpha_aux_l
    area = area_aux
  [../]

  [./TotEnerAKLiq]
    type = TotalEnergyAux
    variable = total_energy_aux_l
    alrhoEA = alrhoEA_l
    vf_liquid = alpha_aux_l
    area = area_aux
  [../]

  [./IntEnerAKLiq]
    type = InternalEnergyAux
    variable = internal_energy_aux_l
    alrhoA = alrhoA_l
    alrhouA_x = alrhouA_l
    alrhoEA = alrhoEA_l
    vf_liquid = alpha_aux_l
    area = area_aux
  [../]

  [./PressAKLiq]
    type = PressureAux
    variable = pressure_aux_l
    alrhoA = alrhoA_l
    alrhouA_x = alrhouA_l
    alrhoEA = alrhoEA_l
    vf_liquid = alpha_aux_l
    area = area_aux
    eos = eos_liq
  [../]

  [./MachAKLiq]
    type = MachNumberAux
    variable = mach_aux_l
    alrhoA = alrhoA_l
    alrhouA_x = alrhouA_l
    pressure = pressure_aux_l
    vf_liquid = alpha_aux_l
    area = area_aux
    eos = eos_liq
  [../]

  [./MuAKLiq]
    type = MaterialRealAux
    variable = mu_aux_l
    property = mu_liq
  [../]

  [./KappaAKLiq]
    type = MaterialRealAux
    variable = kappa_aux_l
    property = kappa_liq
  [../]

  [./MuMaxAKLiq]
    type = MaterialRealAux
    variable = mu_max_aux_l
    property = mu_max_liq
  [../]

  [./KappaMaxAKLiq]
    type = MaterialRealAux
    variable = kappa_max_aux_l
    property = kappa_max_liq
  [../]
####### Gas phase ##########
  [./VelAKGas]
    type = VelocityAux
    variable = velocity_x_aux_g
    alrhoA = alrhoA_g
    alrhouA = alrhouA_g
  [../]

  [./DensAKGas]
    type = DensityAux
    variable = density_aux_g
    alrhoA = alrhoA_g
    vf_liquid = alpha_aux_l
    area = area_aux
    isLiquid = false
  [../]

  [./TotEnerAKGas]
    type = TotalEnergyAux
    variable = total_energy_aux_g
    alrhoEA = alrhoEA_g
    vf_liquid = alpha_aux_l
    area = area_aux
    isLiquid = false
  [../]

  [./IntEnerAKGas]
    type = InternalEnergyAux
    variable = internal_energy_aux_g
    alrhoA = alrhoA_g
    alrhouA_x = alrhouA_g
    alrhoEA = alrhoEA_g
    vf_liquid = alpha_aux_l
    area = area_aux
    isLiquid = false
  [../]

  [./PressAKGas]
    type = PressureAux
    variable = pressure_aux_g
    alrhoA = alrhoA_g
    alrhouA_x = alrhouA_g
    alrhoEA = alrhoEA_g
    vf_liquid = alpha_aux_l
    area = area_aux
    eos = eos_gas
    isLiquid = false
  [../]

  [./MachAKGas]
    type = MachNumberAux
    variable = mach_aux_g
    alrhoA = alrhoA_g
    alrhouA_x = alrhouA_g
    vf_liquid = alpha_aux_l
    pressure = pressure_aux_g
    area = area_aux
    eos = eos_gas
    isLiquid = false
  [../]

  [./MuAKGas]
    type = MaterialRealAux
    variable = mu_aux_g
    property = mu_gas
  [../]

  [./KappaAKGas]
    type = MaterialRealAux
    variable = kappa_aux_g
    property = kappa_gas
  [../]

  [./MuMaxAKGas]
    type = MaterialRealAux
    variable = mu_max_aux_g
    property = mu_max_gas
  [../]

  [./KappaMaxAKGas]
    type = MaterialRealAux
    variable = kappa_max_aux_g
    property = kappa_max_gas
  [../]
[]

##############################################################################################
#                                       MATERIALS                                            #
##############################################################################################
# Define functions that are used in the kernels and aux. kernels.                            #
##############################################################################################

[Materials]
  [./ViscCoeffLiq]
    type = ComputeViscCoeff
    block = '0'
    velocity_x = velocity_x_aux_l
    pressure = pressure_aux_l
    density = density_aux_l
    jump_grad_press = jump_grad_press_smooth_aux_l
    jump_grad_dens = jump_grad_dens_smooth_aux_l # jump_grad_dens_aux_l
    jump_grad_alpha = jump_grad_alpha_aux_l # jump_grad_alpha_aux
    vf_liquid = alpha_aux_l
    eos = eos_liq
    rhov2_PPS_name = PpsRhoVel2Liq
    alpha_PPS_name = PpsVolumeFractionLiq
  [../]

  [./ViscCoeffGas]
    type = ComputeViscCoeff
    block = '0'
    velocity_x = velocity_x_aux_g
    pressure = pressure_aux_g
    density = density_aux_g
    jump_grad_press = jump_grad_press_smooth_aux_g # jump_grad_press_aux_g
    jump_grad_dens = jump_grad_dens_smooth_aux_g # jump_grad_dens_aux_g
    vf_liquid = alpha_aux_l
    eos = eos_gas
    isLiquid = false
    rhov2_PPS_name = PpsRhoVel2Gas
    alpha_PPS_name = PpsVolumeFractionLiq
  [../]

  [./InterfacialRelaxationTransfer]
    type = InterfacialRelaxationTransfer
    block = '0'
    Aint = 0.
    velocity_x_liq = velocity_x_aux_l
    pressure_liq = pressure_aux_l
    density_liq = density_aux_l
    vf_liquid = alpha_aux_l
    velocity_x_gas = velocity_x_aux_g
    pressure_gas = pressure_aux_g
    density_gas = density_aux_g
    eos_liq = eos_liq
    eos_gas = eos_gas
  [../]
[]

##############################################################################################
#                                     PPS                                                    #
##############################################################################################
# Define functions that are used in the kernels and aux. kernels.                            #
##############################################################################################
[Postprocessors]
  [./PpsRhoVel2Liq]
    type = ElementAverageMultipleValues
    variable = density_aux_l
    output_type = RHOVEL2
    alA = alA_l
    alrhoA = alrhoA_l
    alrhouA_x = alrhouA_l
    alrhoEA = alrhoEA_l
    area = area_aux
    eos = eos_liq
    #execute_on = timestep_begin
  [../]

  [./PpsRhoVel2Gas]
    type = ElementAverageMultipleValues
    variable = density_aux_g
    output_type = RHOVEL2
    alA = alA_l
    alrhoA = alrhoA_g
    alrhouA_x = alrhouA_g
    alrhoEA = alrhoEA_g
    area = area_aux
    eos = eos_gas
    isLiquid = false
    #execute_on = timestep_begin
  [../]

  [./PpsVolumeFractionLiq]
    type = ElementAverageAbsValue
    variable = alpha_aux_l
    #execute_on = timestep_begin
  [../]
[]

##############################################################################################
#                               BOUNDARY CONDITIONS                                          #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################
[BCs]
######## Void fraction #######
  [./VoidFractionLeftLiq]
    type = EelDirichletBC
    variable = alA_l
    equation_name = VOIDFRACTION
#    value = 0.8 # 0.999
    eos = eos_liq
    boundary = 'left'
  [../]
######## Liquid phase ########
  [./MassLeftLiq]
    type = EelDirichletBC
    variable = alrhoA_l
    equation_name = CONTINUITY
#    value = 999.
    eos = eos_liq
    boundary = 'left'
  [../]

  [./MassRightLiq]
    type = EelDirichletBC
    variable = alrhoA_l
    equation_name = CONTINUITY
    isLeftBC = false
#    value = 1.
    eos = eos_liq
    boundary = 'right'
  [../]

  [./MomLeftLiq]
    type = EelDirichletBC
    variable = alrhouA_l
    equation_name = XMOMENTUM
#    value = 99.9
    eos = eos_liq
    boundary = 'left'
  [../]

  [./MomRightLiq]
    type = EelDirichletBC
    variable = alrhouA_l
    equation_name = XMOMENTUM
    isLeftBC = false
#    value = 100.
    eos = eos_liq
    boundary = 'right'
  [../]

  [./EnergyLeftLiq]
    type = EelDirichletBC
    variable = alrhoEA_l
    equation_name = ENERGY
#    value = 7.80718528e8
    eos = eos_liq
    boundary = 'left'
  [../]

  [./EnergyRightLiq]
    type = EelDirichletBC
    variable = alrhoEA_l
    equation_name = ENERGY
    isLeftBC = false
#    value = 781500.
    eos = eos_liq
    boundary = 'right'
  [../]

######## Gas phase ########
  [./MassLeftGas]
    type = EelDirichletBC
    variable = alrhoA_g
    equation_name = CONTINUITY
    isLiquid = false
    eos = eos_gas
#    value = 0.01
    boundary = 'left'
  [../]

  [./MassRightGas]
    type = EelDirichletBC
    variable = alrhoA_g
    equation_name = CONTINUITY
    isLiquid = false
    isLeftBC = false
#    value = 0.0999
    eos = eos_gas
    boundary = 'right'
  [../]

  [./MomLeftGas]
    type = EelDirichletBC
    variable = alrhouA_g
    equation_name = XMOMENTUM
    isLiquid = false
#    value = 1.
    eos = eos_gas
    boundary = 'left'
  [../]

  [./MomRightGas]
    type = EelDirichletBC
    variable = alrhouA_g
    equation_name = XMOMENTUM
    isLiquid = false
    isLeftBC = false
#    value = 9.99
    eos = eos_gas
    boundary = 'right'
  [../]

  [./EnergyLeftGas]
    type = EelDirichletBC
    variable = alrhoEA_g
    equation_name = ENERGY
    isLiquid = false
#    value = 250.6250
    eos = eos_gas
    boundary = 'left'
  [../]

  [./EnergyRightGas]
    type = EelDirichletBC
    variable = alrhoEA_g
    equation_name = ENERGY
    isLiquid = false
    isLeftBC = false
#    value = 250374.38
    eos = eos_gas
    boundary = 'right'
  [../]
[]

##############################################################################################
#                                       FUNCTIONs                                            #
##############################################################################################
# Define functions that are used in the kernels and aux. kernels.                            #
##############################################################################################

[Functions]

  [./area]
     type = ConstantFunction
     value = 1.
  [../]

[]

##############################################################################################
#                                  PRECONDITIONER                                            #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Preconditioning]
  active = 'FDP_Newton'
  [./FDP_Newton]
    type = FDP
    full = true
    solve_type = 'PJFNK'
    line_search = 'default'
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-mat_fd_coloring_err  -mat_fd_type  -mat_mffd_type'
    petsc_options_value = '1.e-12       ds             ds'
    #petsc_options = '-snes_mf_operator -ksp_converged_reason -ksp_monitor -snes_ksp_ew'
    #petsc_options_iname = '-pc_type'
    #petsc_options_value = 'lu'
  [../]

  [./SMP]
    type=SMP
    full=true
    petsc_options = '-snes_mf_operator'
    petsc_options_iname = '-pc_type'
    petsc_options_value = 'lu'
  [../]
[]

##############################################################################################
#                                     EXECUTIONER                                            #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Executioner]
  type = Transient   # Here we use the Transient Executioner
  scheme = 'bdf2'
  #num_steps = 10
  end_time = 305.e-6
  dt = 1e-7
  dtmin = 1e-9
  l_tol = 1e-8
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-7
  l_max_its = 50
  nl_max_its = 20
  [./TimeStepper]
    type = FunctionDT
    time_t =  '0      1.e-6  7.e-4'
    time_dt = '1.e-7  1.e-7  1.e-6'
  [../]
  [./Quadrature]
    type = GAUSS
    order = SECOND
  [../]
[]

##############################################################################################
#                                        OUTPUT                                              #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Outputs]
  output_initial = true
  interval = 1
  exodus = true
  perf_log = true
[]
