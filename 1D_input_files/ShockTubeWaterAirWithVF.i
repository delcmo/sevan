#
#####################################################
# Define some global parameters used in the blocks. #
#####################################################
#
[GlobalParams]
###### Other parameters #######
order = FIRST
viscosity_name = FIRST_ORDER
diffusion_name = ENTROPY

###### Boundary Conditions ######
p0_bc = 1.e9
T0_bc = 4.7058824e2
T0_bc_gas = 5.e5
gamma0_bc = 0.
alpha0_bc = 0.9
p_bc = 1.e5
T_bc = 1.765e2
T_bc_gas = 50.
gamma_bc = 0.
alpha_bc = 0.1

###### Initial Conditions #######
pressure_init_left = 1.e9
pressure_init_right = 1.e5
vel_init_left = 0
vel_init_right = 0
temp_init_left = 4.7058824e2
temp_init_right = 1.765e2
temp_init_left_gas = 5.e5
temp_init_right_gas = 50.
alpha_init_left = 0.9
alpha_init_right = 0.1
membrane = 0.5
length = 0.0001
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
  	Cv = 1.e2
  	q_prime = -23e3
  [../]

  [./eos_liq]
    type = EquationOfState
    gamma = 4.4
    Pinf = 6.e8
    q = 0.
    Cv = 1.e3
    q_prime = 0.
  [../]

[]

###### Mesh #######
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 200
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
    scaling = 1e-3
	[./InitialCondition]
        type = ConservativeVariables1DXIC
        area = area
        eos = eos_liq
	[../]
  [../]

  [./alrhouA_l]
    family = LAGRANGE
    scaling = 1e-6
	[./InitialCondition]
        type = ConstantIC
        value = 0.
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
    scaling = 1e-3
    [./InitialCondition]
        type = ConservativeVariables1DXIC
        area = area
        eos = eos_gas
        isLiquid = false
    [../]
  [../]

  [./alrhouA_g]
    family = LAGRANGE
    scaling = 1e-6
    [./InitialCondition]
        type = ConstantIC
        value = 0.
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
#[../]

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
    norm_velocity = norm_velocity_aux_l
    vf_liquid = alpha_aux_l
    eos = eos_liq
  [../]

  [./ViscCoeffGas]
    type = ComputeViscCoeff
    block = '0'
    velocity_x = velocity_x_aux_g
    pressure = pressure_aux_g
    density = density_aux_g
    norm_velocity = norm_velocity_aux_g
    vf_liquid = alpha_aux_l
    eos = eos_gas
    isLiquid = false
  [../]

  [./InterfacialRelaxationTransfer]
    type = InterfacialRelaxationTransfer
    block = '0'
    Aint = 1e6
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
#                               BOUNDARY CONDITIONS                                          #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################
[BCs]
######## Void fraction #######
  [./VoidFractionLeftLiq]
    type = DirichletBC
    variable = alA_l
    value = 0.9
    boundary = 'left'
  [../]
######## Liquid phase ########
  [./MassLeftLiq]
    type = EelStagnationPandTBC
    variable = alrhoA_l
    equation_name = CONTINUITY
    alrhoA = alrhoA_l
    alrhouA_n = alrhouA_l
    area = area_aux
    eos = eos_liq
    #value = 1.5
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
    #value = 0.5
    boundary = 'right'
  [../]

  [./MomLeftLiq]
    type = EelStagnationPandTBC
    variable = alrhouA_l
    equation_name = XMOMENTUM
    alrhoA = alrhoA_l
    alrhouA_n = alrhouA_l
    area = area_aux
    eos = eos_liq
    #value = 0
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
    #value = 0
    boundary = 'right'
  [../]

  [./EnergyLeftLiq]
    type = EelStagnationPandTBC
    variable = alrhoEA_l
    equation_name = ENERGY
    alrhoA = alrhoA_l
    alrhouA_n = alrhouA_l
    area = area_aux
    eos = eos_liq
    #value = 3.75
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
    #value = 1.25
    boundary = 'right'
  [../]

######## Gas phase ########
  [./MassLeftGas]
    type = EelStagnationPandTBC
    variable = alrhoA_g
    equation_name = CONTINUITY
    alrhoA = alrhoA_g
    alrhouA_n = alrhouA_g
    area = area_aux
    eos = eos_gas
    isLiquid = false
    #value = 1.5
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
    #value = 0.5
    boundary = 'right'
  [../]

  [./MomLeftGas]
    type = EelStagnationPandTBC
    variable = alrhouA_g
    equation_name = XMOMENTUM
    alrhoA = alrhoA_g
    alrhouA_n = alrhouA_g
    area = area_aux
    eos = eos_gas
    isLiquid = false
    #value = 0
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
    #value = 0
    boundary = 'right'
  [../]

  [./EnergyLeftGas]
    type = EelStagnationPandTBC
    variable = alrhoEA_g
    equation_name = ENERGY
    alrhoA = alrhoA_g
    alrhouA_n = alrhouA_g
    area = area_aux
    eos = eos_gas
    isLiquid = false
    #value = 0.375
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
    #value = 0.125
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
     type = AreaFunction
     #value = Ao * ( 1 - 0.5*sin((x-left)/l*pi) ) + Bo
     left = 0.0
     length = 1.
     Ao = 0.0
     Bo = 1.0
  [../]

#  [./SaturationTemperature]
#    type = SaturationTemperature
#    eos_liq = eos_liq
#    eos_gas = eos_gas
#[../]

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
    petsc_options = '-snes_mf_operator -snes_ksp_ew'
    petsc_options_iname = '-mat_fd_coloring_err'
    petsc_options_value = '1.e-12'
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
  #scheme = 'explicit-euler'
  string scheme = 'bdf2'
  #petsc_options = '-snes'
  #petsc_options_iname = '-pc_type'
  #petsc_options_value = 'lu'
  #num_steps = 10
  end_time = 7e-5
  dt = 1e-7
  dtmin = 1e-9
  l_tol = 1e-8
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-7
  l_max_its = 50
  nl_max_its = 50
  [./TimeStepper]
    type = FunctionDT
    time_t =  '0     1e-6  1e+5'
    time_dt = '1e-8  1e-7  1e-7'
  [../]
  [./Quadrature]
    type = TRAP
  [../]
[]

##############################################################################################
#                                        OUTPUT                                              #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Output]
  output_initial = true
  interval = 1
  exodus = true
  perf_log = true
[]
