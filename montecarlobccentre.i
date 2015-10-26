[Mesh]
  dim = 2
  file = square_cavity_2D_test_20.exo
  uniform_refine = 0
[]
[MeshModifiers]
  [./scale]
    type = Transform
    transform = SCALE
    vector_value = '0.1 0.1 0.1'
  [../]
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
    boundary = 'in_left_1 in_bottom_1 in_right_1 in_top_1 in_left_0 in_bottom_0 in_right_0 in_top_0'
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
    type = HeatFluxBC
    variable = temp
    boundary = 'out_left out_top out_bottom '
    value = 0
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
    boundary = 'in_left_1 in_bottom_1 in_right_1 in_top_1 in_left_0 in_bottom_0 in_right_0 in_top_0'
  [../]

[]


[UserObjects]
  [./montecarlo_userobject]
    type = ComputeTemperatureBar
    boundary = 'in_left_1 in_bottom_1 in_right_1 in_top_1 in_left_0 in_bottom_0 in_right_0 in_top_0'
    temperature = temp
  [../]
[]

[Materials]
  [./material1]
    type = HeatConductionMaterial
    temperature = temp
    block = '1 2 3 4'
    t_list =  '100 200'
    roe_list = '400 400'
    k_list =  '50 50'
    cp_list = '100 100'
    sigma = 1
  [../]
[./material2]
    type = HeatConductionMaterial
    temperature = temp
    block = 5
    t_list =  '100 200'
    roe_list = '400 400'
    k_list =  '0.5 0.5'
    cp_list = '100 100'
    sigma = 1
  [../]

  [./material_monte]
    type = MonteCarloRadiationMaterial
    temperature = temp
    monte_carlo = montecarlo_userobject
    boundary = 'in_left_1 in_bottom_1 in_right_1 in_top_1 in_left_0 in_bottom_0 in_right_0 in_top_0'
    max_reflect_count = 10
    particle_count=10000
    epsilon=1.0
    absorptivity=1.0
    diffuse_reflectivity=0.5
    mirrors_reflectivity=0.5
  [../]
[]

[Executioner]
  type = Transient
  solve_type = newton
  dt = 1E-02
  num_steps = 1000

  l_tol = 1e-04
  nl_rel_tol = 1e-05
  l_max_its = 10
  nl_max_its = 10
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  [./console]
    type = Console
    perf_log = true
    output_on = 'timestep_end failed nonlinear linear'
  [../]
  [./exodus]
    type = Exodus
    perf_log = true
    output_on = 'timestep_end'
  [../]
[]

