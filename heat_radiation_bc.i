[Mesh]
  file = four_blocks.exo
[]

[Variables]
  [./temp]
    family = LAGRANGE
    order = FIRST
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
  [./left]
    type = HeatRadiationBC
    variable = temp
    boundary = left
  [../]
  [./right]
    type = IsoThermalBC
    variable = temp
    boundary = right
    value = 2000
  [../]
  [./top]
    type = IsoThermalBC
    variable = temp
    boundary = top
    value = 2000
  [../]
 [./bottom]
    type = HeatFluxBC
    variable = temp
    boundary = bottom
    value = 0
  [../]
[]

[Materials]
  [./material]
    type = HeatConductionMaterial
    temperature = temp
    block = ANY_BLOCK_ID
    t_list = '0 1'
    roe_list = '2300 2300'
    k_list = '1 1'
    cp_list = '100 100'
    epsilon_real = 1
  [../]
[]

[Executioner]
  type = Transient
  solve_type = newton
  dt = 1E-02
  num_steps = 500

  l_tol = 1e-04
  nl_rel_tol = 1e-05
  l_max_its = 10
  nl_max_its = 10
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  [./exodus]
    type = Exodus
    output_on = 'initial timestep_end'
  [../]
  [./console]
    type = Console
    perf_log = true
    output_on = 'timestep_end failed nonlinear linear'
  [../]
[]

