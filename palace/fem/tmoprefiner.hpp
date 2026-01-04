// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_TMOP_REFINER_HPP
#define PALACE_FEM_TMOP_REFINER_HPP

#include <memory>
#include <mfem.hpp>
#include "utils/configfile.hpp"

namespace palace
{

class IoData;
class Mesh;
class Timer;

//
// Class handling TMOP-based r-refinement (mesh node movement).
//
class TMOPRefiner
{
private:
  const config::TMOPData &tmop_data;

public:
  TMOPRefiner(const IoData &iodata);

  // Apply TMOP refinement to the given mesh.
  // This moves the mesh nodes to optimize the metric defined in the configuration.
  void Apply(Mesh &mesh, int verbosity = 0) const;
};

}  // namespace palace

#endif  // PALACE_FEM_TMOP_REFINER_HPP
