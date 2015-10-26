/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "OutputTest.h"
#include "FEProblem.h"
#include "Postprocessor.h"
#include "PetscSupport.h"
#include "Executioner.h"
#include "MooseApp.h"
#include "pcrecpp.h"
#include "Moose.h"
#include "Transient.h"
#include <iostream>

template<>
InputParameters validParams<OutputTest>()
{
  InputParameters params = validParams<Console>();
  return params;
}

OutputTest::OutputTest(const InputParameters & parameters) :
    Console(parameters)
{
}

void OutputTest::initialSetup()
{
  Console::initialSetup();

//  mooseError("getOutput test passed");

}

void OutputTest::outputSystemInformation()
{
	mooseError("OutputTest已废弃");
//	 if (!_write_screen && !_write_file)
//	    return;
//
//	  std::stringstream oss;
//
//	  // Framework information
//	  if (_app.getSystemInfo() != NULL)
//	    oss << _app.getSystemInfo()->getInfo();
//
//	  if (_problem_ptr->legacyUoAuxComputation() || _problem_ptr->legacyUoInitialization())
//	  {
//	    oss << "LEGACY MODES ENABLED:\n";
//	    if (_problem_ptr->legacyUoAuxComputation())
//	      oss << "  Computing EXEC_LINEAR AuxKernel types when any UserObject type is executed.\n";
//	    if (_problem_ptr->legacyUoInitialization())
//	      oss << "  Computing all UserObjects during initial setup.\n";
//	  }
//
//	  oss << std::left << '\n'
//	      << "Parallelism:\n"
//	      << std::setw(_field_width) << "  Num Processors: " << static_cast<std::size_t>(n_processors()) << '\n'
//	      << std::setw(_field_width) << "  Num Threads: " << static_cast<std::size_t>(n_threads()) << '\n'
//	      << '\n';
//
//	  MooseMesh & moose_mesh = _problem_ptr->mesh();
//	  MeshBase & mesh = moose_mesh.getMesh();
//	  oss << "Mesh: " << '\n'
//	      << std::setw(_field_width) << "  Distribution: " << (moose_mesh.isParallelMesh() ? "parallel" : "serial")
//	      << (moose_mesh.isDistributionForced() ? " (forced) " : "") << '\n'
//	      << std::setw(_field_width) << "  Mesh Dimension: " << mesh.mesh_dimension() << '\n'
//	      << std::setw(_field_width) << "  Spatial Dimension: " << mesh.spatial_dimension() << '\n'
//	      << std::setw(_field_width) << "  Nodes:" << '\n'
//	      << std::setw(_field_width) << "    Total:" << mesh.n_nodes() << '\n'
//	      << std::setw(_field_width) << "    Local:" << mesh.n_local_nodes() << '\n'
//	      << std::setw(_field_width) << "  Elems:" << '\n'
//	      << std::setw(_field_width) << "    Total:" << mesh.n_elem() << '\n'
//	      << std::setw(_field_width) << "    Local:" << mesh.n_local_elem() << '\n'
//	      << std::setw(_field_width) << "  Num Subdomains: "       << static_cast<std::size_t>(mesh.n_subdomains()) << '\n'
//	      << std::setw(_field_width) << "  Num Partitions: "       << static_cast<std::size_t>(mesh.n_partitions()) << '\n';
//	  if (n_processors() > 1 && moose_mesh.partitionerName() != "")
//	    oss << std::setw(_field_width) << "  Partitioner: "       << moose_mesh.partitionerName()
//	        << (moose_mesh.isPartitionerForced() ? " (forced) " : "")
//	        << '\n';
//	  oss << '\n';
//
//	  EquationSystems & eq = _problem_ptr->es();
//	  unsigned int num_systems = eq.n_systems();
//	  for (unsigned int i=0; i<num_systems; ++i)
//	  {
//	    const System & system = eq.get_system(i);
//	    if (system.system_type() == "TransientNonlinearImplicit")
//	      oss << "Nonlinear System:" << '\n';
//	    else if (system.system_type() == "TransientExplicit")
//	      oss << "Auxiliary System:" << '\n';
//	    else
//	      oss << std::setw(_field_width) << system.system_type() << '\n';
//
//	    if (system.n_dofs())
//	    {
//	      oss << std::setw(_field_width) << "  Num DOFs: " << system.n_dofs() << '\n'
//	          << std::setw(_field_width) << "  Num Local DOFs: " << system.n_local_dofs() << '\n';
//
//	      std::streampos begin_string_pos = oss.tellp();
//	      std::streampos curr_string_pos = begin_string_pos;
//	      oss << std::setw(_field_width) << "  Variables: ";
//	      for (unsigned int vg=0; vg<system.n_variable_groups(); vg++)
//	      {
//	        const VariableGroup &vg_description (system.variable_group(vg));
//
//	        if (vg_description.n_variables() > 1) oss << "{ ";
//	        for (unsigned int vn=0; vn<vg_description.n_variables(); vn++)
//	        {
//	          oss << "\"" << vg_description.name(vn) << "\" ";
//	          curr_string_pos = oss.tellp();
//	          insertNewline(oss, begin_string_pos, curr_string_pos);
//	        }
//
//	        if (vg_description.n_variables() > 1) oss << "} ";
//	      }
//	      oss << '\n';
//
//	      begin_string_pos = oss.tellp();
//	      curr_string_pos = begin_string_pos;
//	      oss << std::setw(_field_width) << "  Finite Element Types: ";
//	#ifndef LIBMESH_ENABLE_INFINITE_ELEMENTS
//	      for (unsigned int vg=0; vg<system.n_variable_groups(); vg++)
//	      {
//	        oss << "\""
//	            << libMesh::Utility::enum_to_string<FEFamily>(system.get_dof_map().variable_group(vg).type().family)
//	            << "\" ";
//	        curr_string_pos = oss.tellp();
//	        insertNewline(oss, begin_string_pos, curr_string_pos);
//	      }
//	      oss << '\n';
//	#else
//	      for (unsigned int vg=0; vg<system.n_variable_groups(); vg++)
//	      {
//	        oss << "\""
//	            << libMesh::Utility::enum_to_string<FEFamily>(system.get_dof_map().variable_group(vg).type().family)
//	            << "\", \""
//	            << libMesh::Utility::enum_to_string<FEFamily>(system.get_dof_map().variable_group(vg).type().radial_family)
//	            << "\" ";
//	        curr_string_pos = oss.tellp();
//	        insertNewline(oss, begin_string_pos, curr_string_pos);
//	      }
//	      oss << '\n';
//
//	      begin_string_pos = oss.tellp();
//	      curr_string_pos = begin_string_pos;
//	      oss << std::setw(_field_width) << "  Infinite Element Mapping: ";
//	      for (unsigned int vg=0; vg<system.n_variable_groups(); vg++)
//	      {
//	        oss << "\""
//	            << libMesh::Utility::enum_to_string<InfMapType>(system.get_dof_map().variable_group(vg).type().inf_map)
//	            << "\" ";
//	        curr_string_pos = oss.tellp();
//	        insertNewline(oss, begin_string_pos, curr_string_pos);
//	      }
//	      oss << '\n';
//	#endif
//
//	      begin_string_pos = oss.tellp();
//	      curr_string_pos = begin_string_pos;
//	      oss << std::setw(_field_width) << "  Approximation Orders: ";
//	      for (unsigned int vg=0; vg<system.n_variable_groups(); vg++)
//	      {
//	#ifndef LIBMESH_ENABLE_INFINITE_ELEMENTS
//	        oss << "\""
//	            << Utility::enum_to_string<Order>(system.get_dof_map().variable_group(vg).type().order)
//	            << "\" ";
//	#else
//	        oss << "\""
//	            << Utility::enum_to_string<Order>(system.get_dof_map().variable_group(vg).type().order)
//	            << "\", \""
//	            << Utility::enum_to_string<Order>(system.get_dof_map().variable_group(vg).type().radial_order)
//	            << "\" ";
//	#endif
//	        curr_string_pos = oss.tellp();
//	        insertNewline(oss, begin_string_pos, curr_string_pos);
//	      }
//	      oss << "\n\n";
//	    }
//	    else
//	      oss << "   *** EMPTY ***\n\n";
//	  }
//
//	  oss << "Execution Information:\n"
//	      << std::setw(_field_width) << "  Executioner: " << demangle(typeid(*_app.getExecutioner()).name()) << '\n';
//
//	  std::string time_stepper = _app.getExecutioner()->getTimeStepperName();
//	  if (time_stepper != "")
//	    oss << std::setw(_field_width) << "  TimeStepper: " << time_stepper << '\n';
//
//	  oss << std::setw(_field_width) << "  Solver Mode: " << Moose::stringify<Moose::SolveType>(_problem_ptr->solverParams()._type) << '\n';
//
//	  const std::string & pc_desc = _problem_ptr->getPreconditionerDescription();
//	  if (!pc_desc.empty())
//	    oss << std::setw(_field_width) << "  Preconditioner: " << pc_desc << '\n';
//	  oss << '\n';
//
//	  // Output information
//	  if (_show_output_on_info)
//	  {
//	    const std::vector<Output *> & outputs = _app.getOutputWarehouse().all();
//	    oss << "Outputs:\n";
//	    for (std::vector<Output *>::const_iterator it = outputs.begin(); it != outputs.end(); ++it)
//	    {
//	      // Display the "output_on" settings
//	      const MultiMooseEnum & output_on = (*it)->outputOn();
//	      oss << "  " << std::setw(_field_width-2) << (*it)->name() <<  "\"" << output_on << "\"\n";
//
//	      // Display the advanced "output_on" settings, only if they are different from "output_on"
//	      if ((*it)->isAdvanced())
//	      {
//	        const OutputOnWarehouse & adv_on = (*it)->advancedOutputOn();
//	        for (std::map<std::string, MultiMooseEnum>::const_iterator adv_it = adv_on.begin(); adv_it != adv_on.end(); ++adv_it)
//	          if (output_on != adv_it->second)
//	            oss << "    " << std::setw(_field_width-4) << adv_it->first + ":" <<  "\"" << adv_it->second << "\"\n";
//	      }
//	    }
//	  }
//
//	oss << "\n\n";
//
////	oss << "--------------------------------------------\n\n";
//	oss << COLOR_CYAN <<"        热传导方程有限元模拟      \n\n";
//	oss << COLOR_DEFAULT;
////	oss << "--------------------------------------------\n\n";
//
//	oss << COLOR_YELLOW << "输入参数: \n";
//	oss << COLOR_DEFAULT;
//
//	write(oss.str());
}
