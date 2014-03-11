#
#####################################################
# Define some global parameters used in the blocks. #
#####################################################
#
[GlobalParams]
###### Other parameters #######
order = FIRST
viscosity_name = FIRST_ORDER

###### Initial Conditions #######
pressure_init_left = 3
pressure_init_right = 1
vel_init_left = 0
vel_init_right = 0
temp_init_left = 1
temp_init_right = 1
alpha_init_left = 0.7
alpha_init_right = 0.7
membrane = 0.5
length = 0.
[]

#############################################################################
#                          USER OBJECTS                                     #
#############################################################################
# Define the user object class that store the EOS parameters.               #
#############################################################################

[UserObjects]
  [./eos_liq]
    type = EquationOfState
  	gamma = 1.4
  	Pinf = 0
  	q = 0
  	Cv = 2.5
  	q_prime = 0.
  [../]

  [./eos_gas]
    type = EquationOfState
    gamma = 1.4
    Pinf = 0
    q = 0
    Cv = 2.5
    q_prime = 0.
  [../]

[]

###### Mesh #######
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
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
    scaling = 1e+0
	[./InitialCondition]
        type = ConservativeVariables1DXIC
        area = area
        eos = eos_liq
	[../]
  [../]

  [./alrhouA_l]
    family = LAGRANGE
    scaling = 1e-3
	[./InitialCondition]
        type = ConstantIC
        value = 0.
	[../]
  [../]

  [./alrhovA_l]
    family = LAGRANGE
    scaling = 1e-3
    [./InitialCondition]
    type = ConstantIC
    value = 0.
    [../]
  [../]

  [./alrhoEA_l]
    family = LAGRANGE
    scaling = 1e-3
	[./InitialCondition]
        type = ConservativeVariables1DXIC
        area = area
        eos = eos_liq
	[../]
  [../]

####### VAPOR PHASE ########
  [./alrhoA_g]
    family = LAGRANGE
    scaling = 1e+0
    [./InitialCondition]
        type = ConservativeVariables1DXIC
        area = area
        eos = eos_gas
        isLiquid = false
    [../]
  [../]

  [./alrhouA_g]
    family = LAGRANGE
    scaling = 1e-3
    [./InitialCondition]
        type = ConstantIC
        value = 0.
    [../]
  [../]

  [./alrhovA_g]
    family = LAGRANGE
    scaling = 1e-3
    [./InitialCondition]
    type = ConstantIC
    value = 0.
    [../]
  [../]

  [./alrhoEA_g]
    family = LAGRANGE
    scaling = 1e-3
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

  [./XMomTimeLiq]
    type = EelTimeDerivative
    variable = alrhouA_l
  [../]

  [./YMomTimeLiq]
    type = EelTimeDerivative
    variable = alrhovA_l
  [../]

  [./EnerTimeLiq]
    type = EelTimeDerivative
    variable = alrhoEA_l
  [../]

  [./MassConvLiq]
    type = EelMass
    variable = alrhoA_l
    alrhouA_x = alrhouA_l
    alrhouA_y = alrhovA_l
  [../]

  [./XMomConvLiq]
    type = EelMomentum
    variable = alrhouA_l
    vel_x = velocity_x_aux_l
    vel_y = velocity_y_aux_l
    pressure = pressure_aux_l
    area = area_aux
    vf_liquid = alpha_aux_l
  [../]

  [./YMomConvLiq]
    type = EelMomentum
    variable = alrhovA_l
    vel_x = velocity_x_aux_l
    vel_y = velocity_y_aux_l
    pressure = pressure_aux_l
    area = area_aux
    vf_liquid = alpha_aux_l
    component = 1
  [../]

  [./EnergyConvLiq]
    type = EelEnergy
    variable = alrhoEA_l
    alrhoA = alrhoA_l
    alrhouA_x = alrhouA_l
    alrhouA_y = alrhovA_l
    pressure = pressure_aux_l
    area = area_aux
    vf_liquid = alpha_aux_l
  [../]

  [./MassViscLiq]
    type = EelArtificialVisc
    variable = alrhoA_l
    equation_name = CONTINUITY
    density = density_aux_l
    velocity_x = velocity_x_aux_l
    velocity_y = velocity_y_aux_l
    internal_energy = internal_energy_aux_l
    area = area_aux
    vf_liquid = alpha_aux_l
  [../]

  [./XMomViscLiq]
    type = EelArtificialVisc
    variable = alrhouA_l
    equation_name = XMOMENTUM
    density = density_aux_l
    velocity_x = velocity_x_aux_l
    velocity_y = velocity_y_aux_l
    internal_energy = internal_energy_aux_l
    area = area_aux
    vf_liquid = alpha_aux_l
  [../]

  [./YMomViscLiq]
    type = EelArtificialVisc
    variable = alrhovA_l
    equation_name = YMOMENTUM
    density = density_aux_l
    velocity_x = velocity_x_aux_l
    velocity_y = velocity_y_aux_l
    internal_energy = internal_energy_aux_l
    area = area_aux
    vf_liquid = alpha_aux_l
  [../]

  [./EnergyViscLiq]
    type = EelArtificialVisc
    variable = alrhoEA_l
    equation_name = ENERGY 
    density = density_aux_l
    velocity_x = velocity_x_aux_l
    velocity_y = velocity_y_aux_l
    internal_energy = internal_energy_aux_l
    area = area_aux
    vf_liquid = alpha_aux_l
  [../]

######### Gas phase ##########
  [./MassTimeGas]
    type = EelTimeDerivative
    variable = alrhoA_g
  [../]

  [./XMomTimeGas]
    type = EelTimeDerivative
    variable = alrhouA_g
  [../]

  [./YMomTimeGas]
    type = EelTimeDerivative
    variable = alrhovA_g
  [../]

  [./EnerTimeGas]
    type = EelTimeDerivative
    variable = alrhoEA_g
  [../]

  [./MassConvGas]
    type = EelMass
    variable = alrhoA_g
    alrhouA_x = alrhouA_g
    alrhouA_y = alrhovA_g
  [../]

  [./XMomConvGas]
    type = EelMomentum
    variable = alrhouA_g
    vel_x = velocity_x_aux_g
    vel_y = velocity_y_aux_g
    pressure = pressure_aux_g
    area = area_aux
    vf_liquid = alpha_aux_l
    isLiquid = false
  [../]

  [./YMomConvGas]
    type = EelMomentum
    variable = alrhovA_g
    vel_x = velocity_x_aux_g
    vel_y = velocity_y_aux_g
    pressure = pressure_aux_g
    area = area_aux
    vf_liquid = alpha_aux_l
    isLiquid = false
    component = 1
[../]

  [./EnergyConvGas]
    type = EelEnergy
    variable = alrhoEA_g
    alrhoA = alrhoA_g
    alrhouA_x = alrhouA_g
    alrhouA_y = alrhovA_g
    pressure = pressure_aux_g
    area = area_aux
    vf_liquid = alpha_aux_l
    isLiquid = false
  [../]

  [./MassViscGas]
    type = EelArtificialVisc
    variable = alrhoA_g
    equation_name = CONTINUITY
    density = density_aux_g
    velocity_x = velocity_x_aux_g
    velocity_y = velocity_y_aux_g
    internal_energy = internal_energy_aux_g
    area = area_aux
    vf_liquid = alpha_aux_l
    isLiquid = false
  [../]

  [./XMomViscGas]
    type = EelArtificialVisc
    variable = alrhouA_g
    equation_name = XMOMENTUM
    density = density_aux_g
    velocity_x = velocity_x_aux_g
    velocity_y = velocity_y_aux_g
    internal_energy = internal_energy_aux_g
    area = area_aux
    vf_liquid = alpha_aux_l
    isLiquid = false
  [../]

  [./YMomViscGas]
    type = EelArtificialVisc
    variable = alrhovA_g
    equation_name = YMOMENTUM
    density = density_aux_g
    velocity_x = velocity_x_aux_g
    velocity_y = velocity_y_aux_g
    internal_energy = internal_energy_aux_g
    area = area_aux
    vf_liquid = alpha_aux_l
    isLiquid = false
  [../]

  [./EnergyViscGas]
    type = EelArtificialVisc
    variable = alrhoEA_g
    equation_name = ENERGY
    density = density_aux_g
    velocity_x = velocity_x_aux_g
    velocity_y = velocity_y_aux_g
    internal_energy = internal_energy_aux_g
    area = area_aux
    vf_liquid = alpha_aux_l
    isLiquid = false
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

######### Liquid phase ##########
   [./velocity_x_aux_l]
      family = LAGRANGE
   [../]

   [./velocity_y_aux_l]
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

  [./velocity_y_aux_g]
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

####### Liquid phase ##########
  [./VelXAKLiq]
    type = VelocityAux
    variable = velocity_x_aux_l
    alrhoA = alrhoA_l
    alrhouA = alrhouA_l
  [../]

  [./VelYAKLiq]
    type = VelocityAux
    variable = velocity_y_aux_l
    alrhoA = alrhoA_l
    alrhouA = alrhovA_l
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
    alrhouA_y = alrhovA_l
    alrhoEA = alrhoEA_l
    vf_liquid = alpha_aux_l
    area = area_aux
  [../]

  [./PressAKLiq]
    type = PressureAux
    variable = pressure_aux_l
    alrhoA = alrhoA_l
    alrhouA_x = alrhouA_l
    alrhouA_y = alrhovA_l
    alrhoEA = alrhoEA_l
    vf_liquid = alpha_aux_l
    area = area_aux
    eos = eos_liq
  [../]

  [./NormVelAKLiq]
    type = NormVectorAux
    variable = norm_velocity_aux_l
    x_component = velocity_x_aux_l
    y_component = velocity_y_aux_l
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
  [./VelXAKGas]
    type = VelocityAux
    variable = velocity_x_aux_g
    alrhoA = alrhoA_g
    alrhouA = alrhouA_g
  [../]

  [./VelYAKGas]
    type = VelocityAux
    variable = velocity_y_aux_g
    alrhoA = alrhoA_g
    alrhouA = alrhovA_g
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
    alrhouA_y = alrhovA_g
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
    alrhouA_y = alrhovA_g
    alrhoEA = alrhoEA_g
    vf_liquid = alpha_aux_l
    area = area_aux
    eos = eos_gas
    isLiquid = false
  [../]

  [./NormVelAKGas]
    type = NormVectorAux
    variable = norm_velocity_aux_g
    x_component = velocity_x_aux_g
    y_component = velocity_y_aux_g
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
    velocity_y = velocity_y_aux_l
    pressure = pressure_aux_l
    density = density_aux_l
    norm_velocity = norm_velocity_aux_l
    eos = eos_liq
  [../]

  [./ViscCoeffGas]
    type = ComputeViscCoeff
    block = '0'
    velocity_x = velocity_x_aux_g
    velocity_y = velocity_y_aux_g
    pressure = pressure_aux_g
    density = density_aux_g
    norm_velocity = norm_velocity_aux_g
    eos = eos_gas
    isLiquid = false
  [../]
[]

##############################################################################################
#                               BOUNDARY CONDITIONS                                          #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################
[BCs]
######## Liquid phase ########
  [./MassLeftLiq]
    type = DirichletBC
    variable = alrhoA_l
    value = 2.1
    boundary = 'left'
  [../]

  [./MassRightLiq]
    type = DirichletBC
    variable = alrhoA_l
    value = 0.7
    boundary = 'right'
  [../]

  [./MassTopLiq]
    type = EelWallBC
#    type = DirichletBC
    variable = alrhoA_l
    equation_name = CONTINUITY
    pressure = pressure_aux_l
    vf_liquid = alpha_aux_l
    area = area_aux
#    value = 1.5
    boundary = 'top'
  [../]

  [./MassBottomLiq]
    type = EelWallBC
    variable = alrhoA_l
    equation_name = CONTINUITY
    pressure = pressure_aux_l
    vf_liquid = alpha_aux_l
    area = area_aux
#    type = DirichletBC
#    value = 1.5
    boundary = 'bottom'
  [../]

  [./XMomLeftLiq]
    type = DirichletBC
    variable = alrhouA_l
    value = 0
    boundary = 'left'
  [../]

  [./XMomRightLiq]
    type = DirichletBC
    variable = alrhouA_l
    value = 0
    boundary = 'right'
  [../]

  [./XMomTopLiq]
    type = EelWallBC
    variable = alrhouA_l
    equation_name = XMOMENTUM
    pressure = pressure_aux_l
    vf_liquid = alpha_aux_l
    area = area_aux
#    type = DirichletBC
#    value = 0.
    boundary = 'top'
  [../]

  [./XMomBottomLiq]
    type = EelWallBC
    variable = alrhouA_l
    equation_name = XMOMENTUM
    pressure = pressure_aux_l
    vf_liquid = alpha_aux_l
    area = area_aux
#    type = DirichletBC
#    value = 0.
    boundary = 'bottom'
  [../]

  [./YMomLeftLiq]
    type = DirichletBC
    variable = alrhovA_l
    value = 0
    boundary = 'left'
  [../]

  [./YMomRightLiq]
    type = DirichletBC
    variable = alrhovA_l
    value = 0
    boundary = 'right'
  [../]

  [./YMomTopLiq]
    type = EelWallBC
    variable = alrhovA_l
    equation_name = YMOMENTUM
    pressure = pressure_aux_l
    vf_liquid = alpha_aux_l
    area = area_aux
#    type = DirichletBC
#    value = 0.
    boundary = 'top'
  [../]

  [./YMomBottomLiq]
    type = EelWallBC
    variable = alrhovA_l
    equation_name = YMOMENTUM
    pressure = pressure_aux_l
    vf_liquid = alpha_aux_l
    area = area_aux
#    type = DirichletBC
#    value = 0.
    boundary = 'bottom'
  [../]

  [./EnergyLeftLiq]
    type = DirichletBC
    variable = alrhoEA_l
    value = 5.25
    boundary = 'left'
  [../]

  [./EnergyRightLiq]
    type = DirichletBC
    variable = alrhoEA_l
    value = 1.75
    boundary = 'right'
  [../]

  [./EnergyTopLiq]
    type = EelWallBC
    variable = alrhoEA_l
    equation_name = ENERGY
    pressure = pressure_aux_l
    vf_liquid = alpha_aux_l
    area = area_aux
#    type = DirichletBC
#    value = 3.75.
    boundary = 'top'
  [../]

  [./EnergyBottomLiq]
    type = EelWallBC
    variable = alrhoEA_l
    equation_name = ENERGY
    pressure = pressure_aux_l
    vf_liquid = alpha_aux_l
    area = area_aux
#    type = DirichletBC
#    value = 3.75
    boundary = 'bottom'
  [../]
######## Gas phase ########
  [./MassLeftGas]
    type = DirichletBC
    variable = alrhoA_g
    value = 0.9
    boundary = 'left'
  [../]

  [./MassRightGas]
    type = DirichletBC
    variable = alrhoA_g
    value = 0.3
    boundary = 'right'
  [../]

  [./MassTopGas]
    type = EelWallBC
    variable = alrhoA_g
    equation_name = CONTINUITY
    pressure = pressure_aux_g
    vf_liquid = alpha_aux_l
    area = area_aux
    isLiquid = false
#    type = DirichletBC
#    value = 1.5
    boundary = 'top'
  [../]

  [./MassBottomGas]
    type = EelWallBC
    variable = alrhoA_g
    equation_name = CONTINUITY
    pressure = pressure_aux_g
    vf_liquid = alpha_aux_l
    area = area_aux
    isLiquid = false
#    type = DirichletBC
#    value = 1.5
    boundary = 'bottom'
  [../]

  [./XMomLeftGas]
    type = DirichletBC
    variable = alrhouA_g
    value = 0
    boundary = 'left'
  [../]

  [./XMomRightGas]
    type = DirichletBC
    variable = alrhouA_g
    value = 0
    boundary = 'right'
  [../]

  [./XMomTopGas]
   type = EelWallBC
    variable = alrhouA_g
   equation_name = XMOMENTUM
    pressure = pressure_aux_g
    vf_liquid = alpha_aux_l
    area = area_aux
    isLiquid = false
#    type = DirichletBC
#    value = 0.
    boundary = 'top'
  [../]

  [./XMomBottomGas]
    type = EelWallBC
    variable = alrhouA_g
    equation_name = XMOMENTUM
    pressure = pressure_aux_g
    vf_liquid = alpha_aux_l
    area = area_aux
    isLiquid = false
#    type = DirichletBC
#    value = 0.
    boundary = 'bottom'
  [../]

  [./YMomLeftGas]
    type = DirichletBC
    variable = alrhovA_g
    value = 0
    boundary = 'left'
  [../]

  [./YMomRightGas]
    type = DirichletBC
    variable = alrhovA_g
    value = 0
    boundary = 'right'
  [../]

  [./YMomTopGas]
    type = EelWallBC
    variable = alrhovA_g
    equation_name = YMOMENTUM
    pressure = pressure_aux_g
    vf_liquid = alpha_aux_l
    area = area_aux
    isLiquid = false
#    type = DirichletBC
#    value = 0.
    boundary = 'top'
  [../]

  [./YMomBottomGas]
    type = EelWallBC
    variable = alrhovA_g
    equation_name = YMOMENTUM
    pressure = pressure_aux_g
    vf_liquid = alpha_aux_l
    area = area_aux
    isLiquid = false
#    type = DirichletBC
#    value = 0.
    boundary = 'bottom'
  [../]

  [./EnergyLeftGas]
    type = DirichletBC
    variable = alrhoEA_g
    value = 2.25
    boundary = 'left'
  [../]

  [./EnergyRightGas]
    type = DirichletBC
    variable = alrhoEA_g
    value = 0.75
    boundary = 'right'
  [../]

  [./EnergyTopGas]
    type = EelWallBC
    variable = alrhoEA_g
    equation_name = ENERGY
    pressure = pressure_aux_g
    vf_liquid = alpha_aux_l
    area = area_aux
    isLiquid = false
#    type = DirichletBC
#    value = 0.375
    boundary = 'top'
  [../]

  [./EnergyBottomGas]
    type = EelWallBC
    variable = alrhoEA_g
    equation_name = ENERGY
    pressure = pressure_aux_g
    vf_liquid = alpha_aux_l
    area = area_aux
    isLiquid = false
#    type = DirichletBC
#    value = 0.375
    boundary = 'bottom'
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
  #num_steps = 1
  end_time = 0.1
  dt = 1e-3
  dtmin = 1e-9
  #dtmax = 1e-5
  l_tol = 1e-8
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-6
  l_max_its = 50
  nl_max_its = 30
  #predictor_scale = 0.0
  #e_tol = 0.01
  #e_max = 0.05
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
