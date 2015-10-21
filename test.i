[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  #second_order = true
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

[AuxVariables]
  [./sigma]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./aukernel_for_sigma]
    type = GetSigma
    variable = sigma
  []
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

[UserObjects]
  [./montecarlo_material]
    type = ComputeTemperatureBar
    boundary = '1'
    max_reflect_count = 10
    particle_count=10000
    absorptivity=1.0
    diffuse_reflectivity=0.5
    mirrors_reflectivity=0.5
    temperature = temp
  [../]
[]

[BCs]
  [./left]
    type = HeatFluxBC
    variable = temp
    boundary = left
    value = 0
    qc = 0.01
    #data_file = qc.dat
  [../]
  [./right]
    type = HeatRadiationBC
    variable = temp
    boundary = right
    epsilon = 1e-6
    tw0 = 300
    scale = '1 1 1'
    sigma = sigma
    #sigma = 0.1
    flux_data = qc.dat
    fluxcoff = 1
  [../]
  [./top]
    type = DirichletBC
    variable = temp
    boundary = top
    value = 2000
  [../]
 [./bottom]
    type = DirichletBC
    variable = temp
    boundary = bottom
    value = 2000
  [../]
[]

[Materials]
  [./material]
    type = HeatConductionMaterial
    temperature = temp
    block = 0
    t_list = '0 0.5 1'
    roe_list = '20 20 20'
    k_list = '1 1 1'
    cp_list = '1 1 1'
    property_file = 'k(T1).data'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = newton
  dt = 1E-02
  num_steps = 1

  l_tol = 1e-04
  nl_rel_tol = 1e-05
  l_max_its = 10
  nl_max_its = 10
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Postprocessors]
  [./run_time]
    type = RunTime
    time_type = active
  [../]
[]
[Outputs]
  [./exodus]
    type = Exodus
    #refinements = 1
    output_on = 'initial timestep_end'
  [../]
  [./console]
    type = Console
    perf_log = true
    output_on = 'timestep_end failed nonlinear linear'
  [../]
[]

