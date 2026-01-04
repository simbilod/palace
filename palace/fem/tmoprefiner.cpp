// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "tmoprefiner.hpp"

#include "fem/mesh.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"
#include <mfem/fem/tmop_tools.hpp>

namespace palace
{

TMOPRefiner::TMOPRefiner(const IoData &iodata)
  : tmop_data(iodata.model.refinement.tmop_data)
{
}

void TMOPRefiner::Apply(Mesh &mesh, int verbosity) const
{
  // 1. Get the underlying MFEM ParMesh
  auto &par_mesh = mesh.Get();

  // 2. Setup the TargetConstructor
  // TargetType:
  // 1 - Ideal Shape (IDEAL_SHAPE_UNIT_SIZE)
  // 2 - Ideal Shape Given Size (IDEAL_SHAPE_GIVEN_SIZE)
  mfem::TargetConstructor::TargetType target_t;
  if (tmop_data.target_id == 1)
  {
      target_t = mfem::TargetConstructor::IDEAL_SHAPE_UNIT_SIZE; 
  }
  else if (tmop_data.target_id == 2)
  {
      target_t = mfem::TargetConstructor::IDEAL_SHAPE_GIVEN_SIZE;
  }
  else
  {
      MFEM_ABORT("Unknown TMOP Target ID: " << tmop_data.target_id);
      return;
  }
  
  // Create target constructor with MPI communicator for parallel mesh
  auto target_c = std::make_unique<mfem::TargetConstructor>(target_t, par_mesh.GetComm());
  
  // Set initial node positions as target reference
  target_c->SetNodes(*par_mesh.GetNodes());

  // 3. Setup the Quality Metric
  // Available metrics:
  //   2D: 002, 007, 077, 080
  //   3D: 302, 303, 301 (shape+size)
  std::unique_ptr<mfem::TMOP_QualityMetric> metric;
  const int dim = par_mesh.Dimension();
  
  switch (tmop_data.metric_id)
  {
    case 2:  // Default shape metric
    case 302:
      if (dim == 2) { metric = std::make_unique<mfem::TMOP_Metric_002>(); }
      else { metric = std::make_unique<mfem::TMOP_Metric_302>(); }
      break;
    case 7:  // Shape+size balanced
    case 307:
      if (dim == 2) { metric = std::make_unique<mfem::TMOP_Metric_007>(); }
      else { metric = std::make_unique<mfem::TMOP_Metric_302>(); } // 307 not available
      break;
    case 77:  // Size-preserving  
    case 377:
      if (dim == 2) { metric = std::make_unique<mfem::TMOP_Metric_077>(); }
      else { metric = std::make_unique<mfem::TMOP_Metric_302>(); }
      break;
    case 80:  // Shape(gamma) + Size(1-gamma) combo
    case 380:
      if (dim == 2) { metric = std::make_unique<mfem::TMOP_Metric_080>(0.5); }
      else { metric = std::make_unique<mfem::TMOP_Metric_302>(); }
      break;
    case 303:  // Alternative 3D shape
      metric = std::make_unique<mfem::TMOP_Metric_303>();
      break;
    case 301:  // 3D shape+size
      metric = std::make_unique<mfem::TMOP_Metric_301>();
      break;
    default:
      MFEM_ABORT("Unsupported TMOP Metric ID: " << tmop_data.metric_id);
  }

  // 4. Setup the TMOP Integrator
  auto tmop_integ = std::make_unique<mfem::TMOP_Integrator>(metric.get(), target_c.get());
  
  // Set integration rules - use lower order for efficiency
  mfem::IntegrationRules int_rules(0, mfem::Quadrature1D::GaussLobatto);
  int quad_order = 4; // reasonable default for tetrahedral mesh
  tmop_integ->SetIntegrationRules(int_rules, quad_order);

  // 5. Setup Tangential Relaxation for Boundary Preservation
  // Identify all boundary attributes
  //  mfem::Array<int> bdr_attr_marker(par_mesh.bdr_attributes.Max() > 0 ? par_mesh.bdr_attributes.Max() : 0);
  // TODO: Implement boundary preservation correctly. SetBoundaryTangentialRelaxation not found in this MFEM version.
  // Might need EnableSurfaceFitting or manual constraint.
  // bdr_attr_marker = 1; // Mark all boundaries for tangential relaxation
  // tmop_integ->SetBoundaryTangentialRelaxation(bdr_attr_marker);
  
  // Enable limiting to avoid mesh inversion/crossing
  // We use the current nodes as the reference for limiting
  tmop_integ->SetLimitingNodes(*par_mesh.GetNodes());


  // 6. Setup the Nonlinear Form
  // We utilize the mesh's own nodal FiniteElementSpace
  // GetNodalFESpace() returns const, so we go through GetNodes()
  if (!par_mesh.GetNodes())
  {
      par_mesh.EnsureNodes();
  }
  mfem::FiniteElementSpace *fespace = par_mesh.GetNodes()->FESpace();
  
  auto *pfespace = dynamic_cast<mfem::ParFiniteElementSpace*>(fespace);
  MFEM_VERIFY(pfespace, "Nodal FESpace must be a ParFiniteElementSpace!");

  mfem::ParNonlinearForm nlf(pfespace);
  nlf.AddDomainIntegrator(tmop_integ.release()); 
  // Note: Integrator ownership is passed to the NonlinearForm

  // 6b. Boundary Preservation: Fix all boundary DOFs using essential BCs.
  // This prevents boundary nodes from moving at all, preserving CAD geometry.
  // Mark all boundary attributes as essential.
  mfem::Array<int> ess_bdr(par_mesh.bdr_attributes.Max());
  ess_bdr = 1; // Mark all boundaries as essential (fixed)
  nlf.SetEssentialBC(ess_bdr);

  // 7. Setup the Solver using MFEM's TMOPNewtonSolver
  // TMOPNewtonSolver uses LBFGS with built-in line search and mesh inversion detection
  mfem::IntegrationRules int_rules_solver(0, mfem::Quadrature1D::GaussLobatto);
  mfem::TMOPNewtonSolver solver(par_mesh.GetComm(), int_rules_solver.Get(mfem::Geometry::TETRAHEDRON, quad_order), 1);
  // solver_type = 1 means LBFGS
  solver.SetIntegrationRules(int_rules_solver, quad_order);
  
  solver.SetOperator(nlf);
  solver.SetRelTol(tmop_data.tol);
  solver.SetAbsTol(1e-12);
  solver.SetMaxIter(tmop_data.max_it);
  solver.SetPrintLevel(verbosity > 0 ? 1 : 0);

  // 8. Solve
  mfem::Vector b(fespace->GetTrueVSize());
  b = 0.0;
  mfem::Vector x; 
  par_mesh.GetNodes()->GetTrueDofs(x);

  if (verbosity > 0)
  {
      Mpi::Print("\n");
      Mpi::Print("=== TMOP Refinement ===\n");
      Mpi::Print(" Metric ID: {}\n", tmop_data.metric_id);
      Mpi::Print(" Target ID: {}\n", tmop_data.target_id);
      Mpi::Print(" Solver: LBFGS (TMOPNewtonSolver)\n");
      Mpi::Print(" MaxIts: {}\n", tmop_data.max_it);
      Mpi::Print(" RelTol: {:.2e}\n", tmop_data.tol);
      Mpi::Print(" DOFs: {}\n", fespace->GetTrueVSize());
      Mpi::Print(" Essential BCs: {} boundary attributes fixed\n", par_mesh.bdr_attributes.Max());
      Mpi::Print("\n");
  }

  solver.Mult(b, x);

  // 9. Update Mesh
  par_mesh.GetNodes()->SetFromTrueDofs(x);
  
  if (verbosity > 0)
  {
      Mpi::Print("TMOP Refinement Complete.\n\n");
  }
}

}  // namespace palace
