[Mesh]
  file = square_cavity_2D_20.exo
[]

[Variables]
  [./temp]
  [../]
[]

[AuxVariables]
  [./heat_flux]
    order = CONSTANT
    family = MONOMIAL
    
  [../]
[]

[AuxKernels]
  [./heat_flux_aux_kernel]
    type = HeatFluxAuxKernel
    variable = heat_flux
    boundary = 'in_left in_bottom in_right in_top'
  [../]
[]

[ICs]
  [./temp_ic]
    variable = temp
    type = ConstantIC
    value = 300
  [../]
[]

[Kernels]
  [./temporal]
    type = HeatConductionTimeDerivative
    variable = temp
  [../]
  [./diff]
    type = HeatConductionKernel
    variable = temp
  [../]
[]

[BCs]
  [./out300]
    type = IsoThermalBC
    variable = temp
    boundary = 'out_left out_top out_bottom '
    value = 300
  [../]
  [./out500]
    type = HeatFluxBC
    variable = temp
    boundary = out_right
    value = 1000000
  [../]
  [./inner]
    type = HeatRadiationBC
    variable = temp
    boundary = 'in_left in_bottom in_right in_top'
  [../]

[]


[UserObjects]
  [./montecarlo_userobject]
    type = ComputeTemperatureBar
    boundary = 'in_left in_bottom in_right in_top'
    max_reflect_count = 10
    particle_count=10000
    absorptivity=0.5
    diffuse_reflectivity=0.5
    mirrors_reflectivity=0.5
    temperature = temp
  [../]
[]

[Materials]
  [./material]
    type = HeatConductionMaterial
    temperature = temp
    block = ANY_BLOCK_ID
    t_list =  '100 200'
    roe_list = '400 400'
    k_list =  '1000 1000'
    cp_list = '100 100'
    sigma = 1
  [../]

  [./material_monte]
    type = MonteCarloRadiationMaterial
    temperature = temp
    monte_carlo = montecarlo_userobject
    boundary = 'in_left in_bottom in_right in_top'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = newton
  dt = 1E-01
  num_steps = 1000

  l_tol = 1e-04
  nl_rel_tol = 1e-05
  l_max_its = 10
  nl_max_its = 10
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
  output_on = 'initial timestep_end'
  [./console]
    type = Console
    perf_log = true
    output_on = 'timestep_end failed nonlinear linear'
  [../]
[]

