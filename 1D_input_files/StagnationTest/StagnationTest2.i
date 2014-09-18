#
#####################################################
# Define some global parameters used in the blocks. #
#####################################################
#
[GlobalParams]
###### Other parameters #######
order = FIRST
viscosity_name = ENTROPY
diffusion_name = ENTROPY
Ce = 1.
Cjump_liquid = 5.
Cjump_gas = 5.
Calpha = 1.

###### Mass and heat transfer ######
isJumpOn = false
isMassOn = false
isHeatOn = false
isShock = false
Aint = 0.
gravity = '9.8, 0.0, 0.0'

###### Initial Conditions #######
p_bc = 1.e6
T_bc = 300.
gamma_bc = 0.
pressure_init_left = 1.e6
pressure_init_right = 1.e6
vel_init_left = 0.
vel_init_right = 0.
temp_init_left = 300.
temp_init_right = 300.
alpha_init_left = 0.5
alpha_init_right = 0.5
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
  	gamma = 1.43
  	Pinf = 0
  	q = 2030e3
  	Cv = 1040
  	q_prime = -23e3
  [../]

  [./eos_liq]
    type = EquationOfState
    gamma = 2.35
    Pinf = 1.e9
    q = -1167e3
    Cv = 1816
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

  [./SmoothJumpGradDensGas]
    type = SmoothFunction
    variable = jump_grad_dens_aux_g
    var_name = smooth_jump_grad_dens_aux_g
    execute_on = timestep_begin
  [../]

  [./SmoothJumpGradPressGas]
    type = SmoothFunction
    variable = jump_grad_press_aux_g
    var_name = smooth_jump_grad_press_aux_g
    execute_on = timestep_begin
  [../]

  [./SmoothJumpGradDensLiq]
    type = SmoothFunction
    variable = jump_grad_dens_aux_l
    var_name = smooth_jump_grad_dens_aux_l
    execute_on = timestep_begin
  [../]

  [./SmoothJumpGradPressLiq]
    type = SmoothFunction
    variable = jump_grad_press_aux_l
    var_name = smooth_jump_grad_press_aux_l
    execute_on = timestep_begin
  [../]

  [./JumpGradAlphaLiq]
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
  nx = 100
  xmin = 0
  xmax = 1
  block_id = '0'
[]

##############################################################################################
#                                       FUNCTIONs                                            #
##############################################################################################
# Define functions that are used in the kernels and aux. kernels.                            #
##############################################################################################

[Functions]
  [./area]
    type = ParsedFunction
    value = 1.
  [../]
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
    scaling = 1e-2
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
    scaling = 1e-4
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
    scaling = 1e-4
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
    alrhoA = alrhoA_l
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
    velocity_x =   velocity_x_aux_l
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
    alrhoA = alrhoA_g
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

  [./beta_max_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./beta_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./PI_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./PI_bar_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

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

######### Liquid phase ##########
   [./velocity_x_aux_l]
      family = LAGRANGE
   [../]

   [./density_aux_l]
      family = LAGRANGE
   [../]

   [./internal_energy_aux_l]
      family = LAGRANGE
   [../]

   [./pressure_aux_l]
      family = LAGRANGE
   [../]

  [./mach_number_aux_l]
    family = LAGRANGE
  [../]

   [./temperature_aux_l]
    family = LAGRANGE
   [../]

   [./norm_velocity_aux_l]
    family = LAGRANGE
   [../]

  [./mu_max_aux_l]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./kappa_max_aux_l]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./mu_aux_l]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./kappa_aux_l]
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

  [./smooth_jump_grad_press_aux_l]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./smooth_jump_grad_dens_aux_l]
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

  [./internal_energy_aux_g]
    family = LAGRANGE
  [../]

  [./pressure_aux_g]
    family = LAGRANGE
  [../]

  [./mach_number_aux_g]
    family = LAGRANGE
  [../]

  [./temperature_aux_g]
    family = LAGRANGE
  [../]

  [./norm_velocity_aux_g]
    family = LAGRANGE
  [../]

  [./mu_max_aux_g]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./kappa_max_aux_g]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./mu_aux_g]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./kappa_aux_g]
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

  [./smooth_jump_grad_press_aux_g]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./smooth_jump_grad_dens_aux_g]
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

  [./PIbarAK]
    type = MaterialRealAux
    variable = PI_bar_aux
    property = average_interfacial_pressure
  [../]

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
  
  [./MachNumAKLiq]
    type = MachNumberAux
    variable = mach_number_aux_l
    alrhoA = alrhoA_l
    alrhouA_x = alrhouA_l
    vf_liquid = alpha_aux_l
    area = area_aux
    pressure = pressure_aux_l
    eos = eos_liq
  [../]

  [./TempAKLiq]
    type = TemperatureAux
    variable = temperature_aux_l
    pressure = pressure_aux_l
    density = density_aux_l
    eos = eos_liq
  [../]

  [./NormVelAKLiq]
    type = NormVectorAux
    variable = norm_velocity_aux_l
    x_component = velocity_x_aux_l
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

  [./BetaMaxAKLiq]
    type = MaterialRealAux
    variable = beta_max_aux
    property = beta_max
  [../]

  [./BetaAKLiq]
    type = MaterialRealAux
    variable = beta_aux
    property = beta
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

  [./MachNumAKGas]
    type = MachNumberAux
    variable = mach_number_aux_g
    alrhoA = alrhoA_g
    alrhouA_x = alrhouA_g
    vf_liquid = alpha_aux_l
    area = area_aux
    pressure = pressure_aux_g
    eos = eos_gas
    isLiquid = false
  [../]

  [./TempAKGas]
    type = TemperatureAux
    variable = temperature_aux_g
    pressure = pressure_aux_g
    density = density_aux_g
    eos = eos_gas
  [../]

  [./NormVelAKGas]
    type = NormVectorAux
    variable = norm_velocity_aux_g
    x_component = velocity_x_aux_g
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
    jump_grad_press = smooth_jump_grad_press_aux_l
    jump_grad_dens = smooth_jump_grad_dens_aux_l
    jump_grad_alpha = jump_grad_alpha_aux_l
    vf_liquid = alpha_aux_l
    area = area_aux
    eos = eos_liq
    rhov2_PPS_name = AverageRhov2Liq
    alpha_PPS_name = AverageAlphaLiq
  [../]

  [./ViscCoeffGas]
    type = ComputeViscCoeff
    block = '0'
    velocity_x = velocity_x_aux_g
    pressure = pressure_aux_g
    density = density_aux_g
    jump_grad_press = smooth_jump_grad_press_aux_g
    jump_grad_dens = smooth_jump_grad_dens_aux_g
    vf_liquid = alpha_aux_l
    area = area_aux
    eos = eos_gas
    rhov2_PPS_name = AverageRhov2Gas
    alpha_PPS_name = AverageAlphaLiq
    isLiquid = false
  [../]

  [./InterfacialRelaxationTransfer]
    type = InterfacialRelaxationTransfer
    block = '0'
    velocity_x_liq = velocity_x_aux_l
    pressure_liq = pressure_aux_l
    density_liq = density_aux_l
    vf_liquid = alpha_aux_l
    velocity_x_gas = velocity_x_aux_g
    pressure_gas  =pressure_aux_g
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
  [./AverageAlphaLiq]
    type = ElementAverageValue # NodalMaxValue
    variable = alpha_aux_l
    #execute_on = timestep_begin
  [../]

  [./AverageRhov2Liq]
    type = ElementAverageMultipleValues
    variable = norm_velocity_aux_g
    output_type = RHOVEL2
    alA = alA_l
    alrhoA = alrhoA_l
    alrhouA_x = alrhouA_l
    alrhoEA = alrhoEA_l
    eos = eos_liq
    area = area_aux
  [../]

  [./AverageRhov2Gas]
    type = ElementAverageMultipleValues
    variable = norm_velocity_aux_g
    output_type = RHOVEL2
    alA = alA_l
    alrhoA = alrhoA_g
    alrhouA_x = alrhouA_g
    alrhoEA = alrhoEA_g
    eos = eos_gas
    area = area_aux
    isLiquid = false
  [../]
[]

##############################################################################################
#                               BOUNDARY CONDITIONS                                          #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################
[BCs]
    [./VoidFractionRightLiq]
    type = DirichletBC
    variable = alA_l
    value = 0.5
    boundary = 'right'
    [../]
######## Liquid phase ########
  [./MassLeftLiq]
    type = EelWallBC
    variable = alrhoA_l
    equation_name = CONTINUITY
    pressure = pressure_aux_l
    vf_liquid = alpha_aux_l
    area = area_aux
    boundary = 'left'
  [../]

    [./MassRightLiq]
    type = EelStaticPandTBC
    variable = alrhoA_l
    equation_name = CONTINUITY
    vel_x = velocity_x_aux_l
    area = area_aux
    temperature = temperature_aux_l
    vf_liquid = alpha_aux_l
    eos = eos_liq
    boundary = 'right'
    [../]

  [./MomLeftLiq]
    type = EelWallBC
    variable = alrhouA_l
    equation_name = XMOMENTUM
    pressure = pressure_aux_l
    vf_liquid = alpha_aux_l
    area = area_aux
    boundary = 'left'
  [../]

    [./MomRightLiq]
    type = EelStaticPandTBC
    variable = alrhouA_l
    equation_name = XMOMENTUM
    vel_x = velocity_x_aux_l
    area = area_aux
    temperature = temperature_aux_l
    vf_liquid = alpha_aux_l
    eos = eos_liq
    boundary = 'right'
    [../]

  [./EnergyLeftLiq]
    type = EelWallBC
    variable = alrhoEA_l
    equation_name = ENERGY
    pressure = pressure_aux_l
    vf_liquid = alpha_aux_l
    area = area_aux
    boundary = 'left'
  [../]

    [./EnergyRightLiq]
    type = EelStaticPandTBC
    variable = alrhoEA_l
    equation_name = ENERGY
    vel_x = velocity_x_aux_l
    area = area_aux
    temperature = temperature_aux_l
    vf_liquid = alpha_aux_l
    eos = eos_liq
    boundary = 'right'
    [../]

######## Gas phase ########
  [./MassLeftGas]
    type = EelWallBC
    variable = alrhoA_g
    equation_name = CONTINUITY
    pressure = pressure_aux_g
    vf_liquid = alpha_aux_l
    area = area_aux
    isLiquid = false
    boundary = 'left'
  [../]

    [./MassRightGas]
    type = EelStaticPandTBC
    variable = alrhoA_g
    equation_name = CONTINUITY
    vel_x = velocity_x_aux_g
    area = area_aux
    temperature = temperature_aux_g
    vf_liquid = alpha_aux_l
    eos = eos_gas
    isLiquid = false
    boundary = 'right'
    [../]

  [./MomLeftGas]
    type = EelWallBC
    variable = alrhouA_g
    equation_name = XMOMENTUM
    pressure = pressure_aux_g
    vf_liquid = alpha_aux_l
    area = area_aux
    isLiquid = false
    boundary = 'left'
  [../]

    [./MomRightGas]
    type = EelStaticPandTBC
    variable = alrhouA_g
    equation_name = XMOMENTUM
    vel_x = velocity_x_aux_g
    area = area_aux
    temperature = temperature_aux_g
    vf_liquid = alpha_aux_l
    eos = eos_gas
    isLiquid = false
    boundary = 'right'
    [../]

  [./EnergyLeftGas]
    type = EelWallBC
    variable = alrhoEA_g
    equation_name = ENERGY
    pressure = pressure_aux_g
    vf_liquid = alpha_aux_l
    area = area_aux
    isLiquid = false
    boundary = 'left'
  [../]

    [./EnergyRightGas]
    type = EelStaticPandTBC
    variable = alrhoEA_g
    equation_name = ENERGY
    vel_x = velocity_x_aux_g
    area = area_aux
    temperature = temperature_aux_g
    vf_liquid = alpha_aux_l
    eos = eos_gas
    isLiquid = false
    boundary = 'right'
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
    petsc_options_iname = '-mat_fd_coloring_err  -mat_fd_type  -mat_mffd_type'
    petsc_options_value = '1.e-10       ds             ds'
    line_search = 'default'
  [../]

  [./SMP]
    type=SMP
    full=true
    solve_type = 'PJFNK'
    line_search = 'none'
  [../]
[]

##############################################################################################
#                                     EXECUTIONER                                            #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Executioner]
  type = Transient
  scheme = 'bdf2'
  num_steps = 100
  end_time = 1.
  dt = 1.e-2
  dtmin = 1e-9
  l_tol = 1e-8
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-7
  l_max_its = 50
  nl_max_its = 10
[./TimeStepper]
    type = FunctionDT
#    time_t =  '0.      2.e-4    1.e-2   2.e-2   0.56'
#    time_dt = '1.e-5   1.e-4    1.e-4   1.e-3   1.e-3'
#    time_t =  '0.      1.e-2    2.e-2   4.e-2   0.56'
#    time_dt = '1.e-4   1.e-4    1.e-3   1.e-2   1.e-2'
    time_t =  '0.      1.e-4    3.e-2   1.e-1   0.56'
    time_dt = '1.e-4   1.e-2    1.e-2   1.e-2   1.e-2'
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
  console = true
  exodus = true
[]

##############################################################################################
#                                        DEBUG                                               #
##############################################################################################
# Debug                 #
##############################################################################################

#[Debug]
#  show_var_residual_norms = true
#[]
