// This code conforms with the UFC specification version 2018.1.0
// and was automatically generated by FFC version 2019.1.0.post0.
//
// This code was generated with the option '-l dolfin' and
// contains DOLFIN-specific wrappers that depend on DOLFIN.
//
// This code was generated with the following parameters:
//

//  add_tabulate_tensor_timing:     False
//  convert_exceptions_to_warnings: False
//  cpp_optimize:                   True
//  cpp_optimize_flags:             '-O2'
//  epsilon:                        1e-14
//  error_control:                  False
//  external_include_dirs:          ''
//  external_includes:              ''
//  external_libraries:             ''
//  external_library_dirs:          ''
//  form_postfix:                   True
//  format:                         'dolfin'
//  generate_dummy_tabulate_tensor: False
//  max_signature_length:           0
//  optimize:                       True
//  precision:                      None
//  quadrature_degree:              None
//  quadrature_rule:                None
//  representation:                 'auto'
//  split:                          False

#ifndef __LINEARSPACE2D_H
#define __LINEARSPACE2D_H
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <ufc.h>

class linearspace2d_finite_element_0: public ufc::finite_element
{
public:

  linearspace2d_finite_element_0() : ufc::finite_element()
  {
    // Do nothing
  }

  ~linearspace2d_finite_element_0() override
  {
    // Do nothing
  }

  const char * signature() const final override
  {
    return "FiniteElement('Lagrange', triangle, 1)";
  }

  ufc::shape cell_shape() const final override
  {
    return ufc::shape::triangle;
  }

  std::size_t topological_dimension() const final override
  {
    return 2;
  }

  std::size_t geometric_dimension() const final override
  {
    return 2;
  }

  std::size_t space_dimension() const final override
  {
    return 3;
  }

  std::size_t value_rank() const final override
  {
    return 0;
  }

  std::size_t value_dimension(std::size_t i) const final override
  {
    return 1;
  }

  std::size_t value_size() const final override
  {
    return 1;
  }

  std::size_t reference_value_rank() const final override
  {
    return 0;
  }

  std::size_t reference_value_dimension(std::size_t i) const final override
  {
    return 1;
  }

  std::size_t reference_value_size() const final override
  {
    return 1;
  }

  std::size_t degree() const final override
  {
    return 1;
  }

  const char * family() const final override
  {
    return "Lagrange";
  }

  void evaluate_reference_basis(double * reference_values,
                                std::size_t num_points,
                                const double * X) const final override
  {
    static const double coefficients0[1][3] = { { 0.4714045207910317, -0.2886751345948129, -0.16666666666666666 } };
    static const double coefficients1[1][3] = { { 0.4714045207910317, 0.2886751345948129, -0.16666666666666666 } };
    static const double coefficients2[1][3] = { { 0.4714045207910316, 0.0, 0.3333333333333333 } };
    for (std::size_t k = 0; k < num_points * 3; ++k)
        reference_values[k] = 0.0;
    for (std::size_t ip = 0; ip < num_points; ++ip)
    {
        // Map from UFC reference coordinate X to FIAT reference coordinate Y
        const double Y[2] = { 2.0 * X[ip * 2] - 1.0, 2.0 * X[ip * 2 + 1] - 1.0 };
        // Compute basisvalues for each relevant embedded degree
        double basisvalues1[3] = {};
        basisvalues1[0] = 1.0;
        const double tmp1_1 = (1.0 + 2.0 * Y[0] + Y[1]) / 2.0;
        basisvalues1[1] = tmp1_1;
        basisvalues1[2] = (0.5 + 1.5 * Y[1]) * basisvalues1[0];
        basisvalues1[0] *= std::sqrt(0.5);
        basisvalues1[2] *= std::sqrt(1.0);
        basisvalues1[1] *= std::sqrt(3.0);
        // Accumulate products of coefficients and basisvalues
        for (std::size_t r = 0; r < 3; ++r)
            reference_values[3 * ip] += coefficients0[0][r] * basisvalues1[r];
        for (std::size_t r = 0; r < 3; ++r)
            reference_values[3 * ip + 1] += coefficients1[0][r] * basisvalues1[r];
        for (std::size_t r = 0; r < 3; ++r)
            reference_values[3 * ip + 2] += coefficients2[0][r] * basisvalues1[r];
    }
  }

  void evaluate_reference_basis_derivatives(double * reference_values,
                                            std::size_t order,
                                            std::size_t num_points,
                                            const double * X) const final override
  {
    if (order == 0)
    {
        evaluate_reference_basis(reference_values, num_points, X);
        return;
    }
    const std::size_t num_derivatives = std::pow(2, order);
    std::fill_n(reference_values, num_points * 3 * num_derivatives, 0.0);
    if (order > 1)
        return;
    // Tables of derivatives of the polynomial base (transpose).
    alignas(32) static const double dmats0[2][3][3] =
        { { { 0.0, 0.0, 0.0 },
            { 4.8989794855663495, 0.0, 0.0 },
            { 0.0, 0.0, 0.0 } },
          { { 0.0, 0.0, 0.0 },
            { 2.449489742783182, 0.0, 0.0 },
            { 4.242640687119285, 0.0, 0.0 } } };
    static const double coefficients0[1][3] = { { 0.4714045207910317, -0.2886751345948129, -0.16666666666666666 } };
    static const double coefficients1[1][3] = { { 0.4714045207910317, 0.2886751345948129, -0.16666666666666666 } };
    static const double coefficients2[1][3] = { { 0.4714045207910316, 0.0, 0.3333333333333333 } };
    const std::size_t reference_offset[3] = {};
    const std::size_t num_components[3] = { 1, 1, 1 };
    // Precomputed combinations
    const std::size_t combinations[1][2][1] =
        { { { 0 },
            { 1 } } };
    for (std::size_t ip = 0; ip < num_points; ++ip)
    {
        // Map from UFC reference coordinate X to FIAT reference coordinate Y
        const double Y[2] = { 2.0 * X[ip * 2] - 1.0, 2.0 * X[ip * 2 + 1] - 1.0 };
        // Compute basisvalues for each relevant embedded degree
        double basisvalues1[3] = {};
        basisvalues1[0] = 1.0;
        const double tmp1_1 = (1.0 + 2.0 * Y[0] + Y[1]) / 2.0;
        basisvalues1[1] = tmp1_1;
        basisvalues1[2] = (0.5 + 1.5 * Y[1]) * basisvalues1[0];
        basisvalues1[0] *= std::sqrt(0.5);
        basisvalues1[2] *= std::sqrt(1.0);
        basisvalues1[1] *= std::sqrt(3.0);
        // Loop over all dofs
        for (std::size_t i = 0; i < 3; ++i)
        {
            double derivatives[2] = {};
            switch (i)
            {
            case 0:
                // Compute reference derivatives for dof 0.
                for (std::size_t r = 0; r < num_derivatives; ++r)
                {
                    double aux[3] = {};
                    // Declare derivative matrix (of polynomial basis).
                    double dmats[3][3] = {};
                    // Initialize dmats.
                    std::size_t comb = combinations[order - 1][r][0];
                    std::copy_n(&dmats0[comb][0][0], 9, &dmats[0][0]);
                    // Looping derivative order to generate dmats.
                    for (std::size_t s = 1; s < order; ++s)
                    {
                        // Store previous dmats matrix.
                        double dmats_old[3][3];
                        std::copy_n(&dmats[0][0], 9, &dmats_old[0][0]);
                        // Resetting dmats.
                        std::fill_n(&dmats[0][0], 9, 0.0);
                        // Update dmats using an inner product.
                        comb = combinations[order - 1][r][s];
                        for (std::size_t t = 0; t < 3; ++t)
                            for (std::size_t u = 0; u < 3; ++u)
                                for (std::size_t tu = 0; tu < 3; ++tu)
                                    dmats[t][u] += dmats0[comb][t][tu] * dmats_old[tu][u];
                    }
                    for (std::size_t s = 0; s < 3; ++s)
                        for (std::size_t t = 0; t < 3; ++t)
                            aux[s] += dmats[s][t] * basisvalues1[t];
                    derivatives[r] = 0.0;
                    for (std::size_t s = 0; s < 3; ++s)
                        derivatives[r] += coefficients0[0][s] * aux[s];
                }
                break;
            case 1:
                // Compute reference derivatives for dof 1.
                for (std::size_t r = 0; r < num_derivatives; ++r)
                {
                    double aux[3] = {};
                    // Declare derivative matrix (of polynomial basis).
                    double dmats[3][3] = {};
                    // Initialize dmats.
                    std::size_t comb = combinations[order - 1][r][0];
                    std::copy_n(&dmats0[comb][0][0], 9, &dmats[0][0]);
                    // Looping derivative order to generate dmats.
                    for (std::size_t s = 1; s < order; ++s)
                    {
                        // Store previous dmats matrix.
                        double dmats_old[3][3];
                        std::copy_n(&dmats[0][0], 9, &dmats_old[0][0]);
                        // Resetting dmats.
                        std::fill_n(&dmats[0][0], 9, 0.0);
                        // Update dmats using an inner product.
                        comb = combinations[order - 1][r][s];
                        for (std::size_t t = 0; t < 3; ++t)
                            for (std::size_t u = 0; u < 3; ++u)
                                for (std::size_t tu = 0; tu < 3; ++tu)
                                    dmats[t][u] += dmats0[comb][t][tu] * dmats_old[tu][u];
                    }
                    for (std::size_t s = 0; s < 3; ++s)
                        for (std::size_t t = 0; t < 3; ++t)
                            aux[s] += dmats[s][t] * basisvalues1[t];
                    derivatives[r] = 0.0;
                    for (std::size_t s = 0; s < 3; ++s)
                        derivatives[r] += coefficients1[0][s] * aux[s];
                }
                break;
            case 2:
                // Compute reference derivatives for dof 2.
                for (std::size_t r = 0; r < num_derivatives; ++r)
                {
                    double aux[3] = {};
                    // Declare derivative matrix (of polynomial basis).
                    double dmats[3][3] = {};
                    // Initialize dmats.
                    std::size_t comb = combinations[order - 1][r][0];
                    std::copy_n(&dmats0[comb][0][0], 9, &dmats[0][0]);
                    // Looping derivative order to generate dmats.
                    for (std::size_t s = 1; s < order; ++s)
                    {
                        // Store previous dmats matrix.
                        double dmats_old[3][3];
                        std::copy_n(&dmats[0][0], 9, &dmats_old[0][0]);
                        // Resetting dmats.
                        std::fill_n(&dmats[0][0], 9, 0.0);
                        // Update dmats using an inner product.
                        comb = combinations[order - 1][r][s];
                        for (std::size_t t = 0; t < 3; ++t)
                            for (std::size_t u = 0; u < 3; ++u)
                                for (std::size_t tu = 0; tu < 3; ++tu)
                                    dmats[t][u] += dmats0[comb][t][tu] * dmats_old[tu][u];
                    }
                    for (std::size_t s = 0; s < 3; ++s)
                        for (std::size_t t = 0; t < 3; ++t)
                            aux[s] += dmats[s][t] * basisvalues1[t];
                    derivatives[r] = 0.0;
                    for (std::size_t s = 0; s < 3; ++s)
                        derivatives[r] += coefficients2[0][s] * aux[s];
                }
                break;
            }
            for (std::size_t r = 0; r < num_derivatives; ++r)
                for (std::size_t c = 0; c < num_components[i]; ++c)
                    reference_values[3 * num_derivatives * ip + num_derivatives * i + r + (reference_offset[i] + c)] = derivatives[num_derivatives * c + r];
        }
    }
  }

  void transform_reference_basis_derivatives(double * values,
                                             std::size_t order,
                                             std::size_t num_points,
                                             const double * reference_values,
                                             const double * X,
                                             const double * J,
                                             const double * detJ,
                                             const double * K,
                                             int cell_orientation) const final override
  {
    const std::size_t num_derivatives = std::pow(2, order);
    // Precomputed combinations
    const std::size_t combinations[1][2][1] =
        { { { 0 },
            { 1 } } };
    std::fill_n(values, num_points * 3 * num_derivatives, 0.0);
    const std::size_t reference_offsets[3] = {};
    const std::size_t physical_offsets[3] = {};
    for (std::size_t ip = 0; ip < num_points; ++ip)
    {
        double transform[2][2];
        for (std::size_t r = 0; r < num_derivatives; ++r)
            for (std::size_t s = 0; s < num_derivatives; ++s)
                transform[r][s] = 1.0;
        for (std::size_t r = 0; r < num_derivatives; ++r)
            for (std::size_t s = 0; s < num_derivatives; ++s)
                for (std::size_t k = 0; k < order; ++k)
                    transform[r][s] *= K[2 * 2 * ip + 2 * combinations[order - 1][s][k] + combinations[order - 1][r][k]];
        for (std::size_t d = 0; d < 3; ++d)
        {
            for (std::size_t s = 0; s < num_derivatives; ++s)
            {
                for (std::size_t i = 0; i < 1; ++i)
                {
                    // Using affine transform to map values back to the physical element.
                    const double mapped_value = reference_values[3 * num_derivatives * ip + num_derivatives * d + s + reference_offsets[d]];
                    // Mapping derivatives back to the physical element
                    for (std::size_t r = 0; r < num_derivatives; ++r)
                        values[3 * num_derivatives * ip + num_derivatives * d + r + (physical_offsets[d] + i)] += transform[r][s] * mapped_value;
                }
            }
        }
    }
  }

  void evaluate_basis(std::size_t i,
                      double * values,
                      const double * x,
                      const double * coordinate_dofs,
                      int cell_orientation,
                      const ufc::coordinate_mapping * cm=nullptr
                      ) const final override
  {
    double X[2] = {};
    double J[4];
    double detJ;
    double K[4];
    if (cm)
    {
        cm->compute_reference_geometry(X, J, &detJ, K, 1, x, coordinate_dofs, cell_orientation);
    }
    else
    {
        compute_jacobian_triangle_2d(J, coordinate_dofs);
        compute_jacobian_inverse_triangle_2d(K, detJ, J);
        // Compute constants
        const double C0 = coordinate_dofs[2] + coordinate_dofs[4];
        const double C1 = coordinate_dofs[3] + coordinate_dofs[5];
        // Get coordinates and map to the reference (FIAT) element
        double Y[2] = { (J[1] * (C1 - 2.0 * x[1]) + J[3] * (2.0 * x[0] - C0)) / detJ, (J[0] * (2.0 * x[1] - C1) + J[2] * (C0 - 2.0 * x[0])) / detJ };
        // Map to FFC reference coordinate
        for (std::size_t k = 0; k < 2; ++k)
            X[k] = (Y[k] + 1.0) / 2.0;
    }
    // Evaluate basis on reference element
    double ref_values[3];
    evaluate_reference_basis(ref_values, 1, X);
    // Push forward
    double physical_values[3];
    transform_reference_basis_derivatives(physical_values, 0, 1, ref_values, X, J, &detJ, K, cell_orientation);
    for (std::size_t k = 0; k < 1; ++k)
        values[k] = physical_values[i + k];
  }

  void evaluate_basis_all(double * values,
                          const double * x,
                          const double * coordinate_dofs,
                          int cell_orientation,
                          const ufc::coordinate_mapping * cm=nullptr
                          ) const final override
  {
    // Helper variable to hold value of a single dof.
    double dof_values = 0.0;
    // Loop dofs and call evaluate_basis
    for (std::size_t r = 0; r < 3; ++r)
    {
        evaluate_basis(r, &dof_values, x, coordinate_dofs, cell_orientation);
        values[r] = dof_values;
    }
  }

  void evaluate_basis_derivatives(std::size_t i,
                                  std::size_t n,
                                  double * values,
                                  const double * x,
                                  const double * coordinate_dofs,
                                  int cell_orientation,
                                  const ufc::coordinate_mapping * cm=nullptr
                                  ) const final override
  {
    std::size_t num_derivatives = std::pow(2, n);
    std::fill_n(values, num_derivatives, 0.0);
    // Call evaluate_basis_all if order of derivatives is equal to zero.
    if (n == 0)
    {
        evaluate_basis(i, values, x, coordinate_dofs, cell_orientation);
        return;
    }
    // If order of derivatives is greater than the maximum polynomial degree, return zeros.
    if (n > 1)
        return;
    // Compute Jacobian
    double J[4];
    compute_jacobian_triangle_2d(J, coordinate_dofs);
    // Compute Inverse Jacobian and determinant
    double K[4];
    double detJ;
    compute_jacobian_inverse_triangle_2d(K, detJ, J);
    // Compute constants
    const double C0 = coordinate_dofs[2] + coordinate_dofs[4];
    const double C1 = coordinate_dofs[3] + coordinate_dofs[5];
    // Get coordinates and map to the reference (FIAT) element
    double Y[2] = { (J[1] * (C1 - 2.0 * x[1]) + J[3] * (2.0 * x[0] - C0)) / detJ, (J[0] * (2.0 * x[1] - C1) + J[2] * (C0 - 2.0 * x[0])) / detJ };
    // Precomputed combinations
    const std::size_t combinations[1][2][1] =
        { { { 0 },
            { 1 } } };
    // Declare transformation matrix
    double transform[2][2] =
        { { 1.0, 1.0 },
          { 1.0, 1.0 } };
    // Construct transformation matrix
    for (std::size_t row = 0; row < num_derivatives; ++row)
        for (std::size_t col = 0; col < num_derivatives; ++col)
            for (std::size_t k = 0; k < n; ++k)
                transform[row][col] *= K[2 * combinations[n - 1][col][k] + combinations[n - 1][row][k]];
    switch (i)
    {
    case 0:
        {
            double basisvalues[3] = {};
            basisvalues[0] = 1.0;
            const double tmp1_1 = (1.0 + 2.0 * Y[0] + Y[1]) / 2.0;
            basisvalues[1] = tmp1_1;
            basisvalues[2] = (0.5 + 1.5 * Y[1]) * basisvalues[0];
            basisvalues[0] *= std::sqrt(0.5);
            basisvalues[2] *= std::sqrt(1.0);
            basisvalues[1] *= std::sqrt(3.0);
            // Table(s) of coefficients
            static const double coefficients0[3] = { 0.4714045207910317, -0.2886751345948129, -0.16666666666666666 };
            // Tables of derivatives of the polynomial base (transpose).
            static const double dmats0[3][3] =
                { { 0.0, 0.0, 0.0 },
                  { 4.8989794855663495, 0.0, 0.0 },
                  { 0.0, 0.0, 0.0 } };
            static const double dmats1[3][3] =
                { { 0.0, 0.0, 0.0 },
                  { 2.449489742783182, 0.0, 0.0 },
                  { 4.242640687119285, 0.0, 0.0 } };
            // Compute reference derivatives.
            // Declare array of derivatives on FIAT element.
            double derivatives[2] = {};
            // Declare derivative matrix (of polynomial basis).
            double dmats[3][3] =
                { { 1.0, 0.0, 0.0 },
                  { 0.0, 1.0, 0.0 },
                  { 0.0, 0.0, 1.0 } };
            // Declare (auxiliary) derivative matrix (of polynomial basis).
            double dmats_old[3][3] =
                { { 1.0, 0.0, 0.0 },
                  { 0.0, 1.0, 0.0 },
                  { 0.0, 0.0, 1.0 } };
            // Loop possible derivatives.
            for (std::size_t r = 0; r < num_derivatives; ++r)
            {
                // Reset dmats to identity
                std::fill_n(&dmats[0][0], 9, 0.0);
                for (std::size_t t = 0; t < 3; ++t)
                    dmats[t][t] = 1.0;
                // Looping derivative order to generate dmats.
                for (std::size_t s = 0; s < n; ++s)
                {
                    std::copy_n(&dmats[0][0], 9, &dmats_old[0][0]);
                    std::fill_n(&dmats[0][0], 9, 0.0);
                    // Update dmats using an inner product.
                    // _dmats_product(shape_dmats, comb[r][s], 0)
                    if (combinations[n - 1][r][s] == 0)
                    {
                        for (std::size_t t = 0; t < 3; ++t)
                            for (std::size_t u = 0; u < 3; ++u)
                                for (std::size_t tu = 0; tu < 3; ++tu)
                                    dmats[t][u] += dmats_old[tu][u] * dmats0[t][tu];
                    }
                    // _dmats_product(shape_dmats, comb[r][s], 1)
                    if (combinations[n - 1][r][s] == 1)
                    {
                        for (std::size_t t = 0; t < 3; ++t)
                            for (std::size_t u = 0; u < 3; ++u)
                                for (std::size_t tu = 0; tu < 3; ++tu)
                                    dmats[t][u] += dmats_old[tu][u] * dmats1[t][tu];
                    }
                }
                for (std::size_t s = 0; s < 3; ++s)
                    for (std::size_t t = 0; t < 3; ++t)
                        derivatives[r] += coefficients0[s] * dmats[s][t] * basisvalues[t];
            }
            // Transform derivatives back to physical element
            for (std::size_t r = 0; r < num_derivatives; ++r)
                for (std::size_t s = 0; s < num_derivatives; ++s)
                    values[r] += transform[r][s] * derivatives[s];
        }
        break;
    case 1:
        {
            double basisvalues[3] = {};
            basisvalues[0] = 1.0;
            const double tmp1_1 = (1.0 + 2.0 * Y[0] + Y[1]) / 2.0;
            basisvalues[1] = tmp1_1;
            basisvalues[2] = (0.5 + 1.5 * Y[1]) * basisvalues[0];
            basisvalues[0] *= std::sqrt(0.5);
            basisvalues[2] *= std::sqrt(1.0);
            basisvalues[1] *= std::sqrt(3.0);
            // Table(s) of coefficients
            static const double coefficients0[3] = { 0.4714045207910317, 0.2886751345948129, -0.16666666666666666 };
            // Tables of derivatives of the polynomial base (transpose).
            static const double dmats0[3][3] =
                { { 0.0, 0.0, 0.0 },
                  { 4.8989794855663495, 0.0, 0.0 },
                  { 0.0, 0.0, 0.0 } };
            static const double dmats1[3][3] =
                { { 0.0, 0.0, 0.0 },
                  { 2.449489742783182, 0.0, 0.0 },
                  { 4.242640687119285, 0.0, 0.0 } };
            // Compute reference derivatives.
            // Declare array of derivatives on FIAT element.
            double derivatives[2] = {};
            // Declare derivative matrix (of polynomial basis).
            double dmats[3][3] =
                { { 1.0, 0.0, 0.0 },
                  { 0.0, 1.0, 0.0 },
                  { 0.0, 0.0, 1.0 } };
            // Declare (auxiliary) derivative matrix (of polynomial basis).
            double dmats_old[3][3] =
                { { 1.0, 0.0, 0.0 },
                  { 0.0, 1.0, 0.0 },
                  { 0.0, 0.0, 1.0 } };
            // Loop possible derivatives.
            for (std::size_t r = 0; r < num_derivatives; ++r)
            {
                // Reset dmats to identity
                std::fill_n(&dmats[0][0], 9, 0.0);
                for (std::size_t t = 0; t < 3; ++t)
                    dmats[t][t] = 1.0;
                // Looping derivative order to generate dmats.
                for (std::size_t s = 0; s < n; ++s)
                {
                    std::copy_n(&dmats[0][0], 9, &dmats_old[0][0]);
                    std::fill_n(&dmats[0][0], 9, 0.0);
                    // Update dmats using an inner product.
                    // _dmats_product(shape_dmats, comb[r][s], 0)
                    if (combinations[n - 1][r][s] == 0)
                    {
                        for (std::size_t t = 0; t < 3; ++t)
                            for (std::size_t u = 0; u < 3; ++u)
                                for (std::size_t tu = 0; tu < 3; ++tu)
                                    dmats[t][u] += dmats_old[tu][u] * dmats0[t][tu];
                    }
                    // _dmats_product(shape_dmats, comb[r][s], 1)
                    if (combinations[n - 1][r][s] == 1)
                    {
                        for (std::size_t t = 0; t < 3; ++t)
                            for (std::size_t u = 0; u < 3; ++u)
                                for (std::size_t tu = 0; tu < 3; ++tu)
                                    dmats[t][u] += dmats_old[tu][u] * dmats1[t][tu];
                    }
                }
                for (std::size_t s = 0; s < 3; ++s)
                    for (std::size_t t = 0; t < 3; ++t)
                        derivatives[r] += coefficients0[s] * dmats[s][t] * basisvalues[t];
            }
            // Transform derivatives back to physical element
            for (std::size_t r = 0; r < num_derivatives; ++r)
                for (std::size_t s = 0; s < num_derivatives; ++s)
                    values[r] += transform[r][s] * derivatives[s];
        }
        break;
    case 2:
        {
            double basisvalues[3] = {};
            basisvalues[0] = 1.0;
            const double tmp1_1 = (1.0 + 2.0 * Y[0] + Y[1]) / 2.0;
            basisvalues[1] = tmp1_1;
            basisvalues[2] = (0.5 + 1.5 * Y[1]) * basisvalues[0];
            basisvalues[0] *= std::sqrt(0.5);
            basisvalues[2] *= std::sqrt(1.0);
            basisvalues[1] *= std::sqrt(3.0);
            // Table(s) of coefficients
            static const double coefficients0[3] = { 0.4714045207910316, 0.0, 0.3333333333333333 };
            // Tables of derivatives of the polynomial base (transpose).
            static const double dmats0[3][3] =
                { { 0.0, 0.0, 0.0 },
                  { 4.8989794855663495, 0.0, 0.0 },
                  { 0.0, 0.0, 0.0 } };
            static const double dmats1[3][3] =
                { { 0.0, 0.0, 0.0 },
                  { 2.449489742783182, 0.0, 0.0 },
                  { 4.242640687119285, 0.0, 0.0 } };
            // Compute reference derivatives.
            // Declare array of derivatives on FIAT element.
            double derivatives[2] = {};
            // Declare derivative matrix (of polynomial basis).
            double dmats[3][3] =
                { { 1.0, 0.0, 0.0 },
                  { 0.0, 1.0, 0.0 },
                  { 0.0, 0.0, 1.0 } };
            // Declare (auxiliary) derivative matrix (of polynomial basis).
            double dmats_old[3][3] =
                { { 1.0, 0.0, 0.0 },
                  { 0.0, 1.0, 0.0 },
                  { 0.0, 0.0, 1.0 } };
            // Loop possible derivatives.
            for (std::size_t r = 0; r < num_derivatives; ++r)
            {
                // Reset dmats to identity
                std::fill_n(&dmats[0][0], 9, 0.0);
                for (std::size_t t = 0; t < 3; ++t)
                    dmats[t][t] = 1.0;
                // Looping derivative order to generate dmats.
                for (std::size_t s = 0; s < n; ++s)
                {
                    std::copy_n(&dmats[0][0], 9, &dmats_old[0][0]);
                    std::fill_n(&dmats[0][0], 9, 0.0);
                    // Update dmats using an inner product.
                    // _dmats_product(shape_dmats, comb[r][s], 0)
                    if (combinations[n - 1][r][s] == 0)
                    {
                        for (std::size_t t = 0; t < 3; ++t)
                            for (std::size_t u = 0; u < 3; ++u)
                                for (std::size_t tu = 0; tu < 3; ++tu)
                                    dmats[t][u] += dmats_old[tu][u] * dmats0[t][tu];
                    }
                    // _dmats_product(shape_dmats, comb[r][s], 1)
                    if (combinations[n - 1][r][s] == 1)
                    {
                        for (std::size_t t = 0; t < 3; ++t)
                            for (std::size_t u = 0; u < 3; ++u)
                                for (std::size_t tu = 0; tu < 3; ++tu)
                                    dmats[t][u] += dmats_old[tu][u] * dmats1[t][tu];
                    }
                }
                for (std::size_t s = 0; s < 3; ++s)
                    for (std::size_t t = 0; t < 3; ++t)
                        derivatives[r] += coefficients0[s] * dmats[s][t] * basisvalues[t];
            }
            // Transform derivatives back to physical element
            for (std::size_t r = 0; r < num_derivatives; ++r)
                for (std::size_t s = 0; s < num_derivatives; ++s)
                    values[r] += transform[r][s] * derivatives[s];
        }
        break;
    }
  }

  void evaluate_basis_derivatives_all(std::size_t n,
                                      double * values,
                                      const double * x,
                                      const double * coordinate_dofs,
                                      int cell_orientation,
                                      const ufc::coordinate_mapping * cm=nullptr
                                      ) const final override
  {
    // Call evaluate_basis_all if order of derivatives is equal to zero.
    if (n == 0)
    {
        evaluate_basis_all(values, x, coordinate_dofs, cell_orientation);
        return;
    }
    unsigned int num_derivatives = std::pow(2, n);
    // Set values equal to zero.
    std::fill_n(values, num_derivatives * 3, 0.0);
    // If order of derivatives is greater than the maximum polynomial degree, return zeros.
    if (n > 1)
        return;
    // Helper variable to hold values of a single dof.
    double dof_values[2] = {};
    // Loop dofs and call evaluate_basis_derivatives.
    for (std::size_t r = 0; r < 3; ++r)
    {
        evaluate_basis_derivatives(r, n, dof_values, x, coordinate_dofs, cell_orientation);
        for (std::size_t s = 0; s < num_derivatives; ++s)
            values[num_derivatives * r + s] = dof_values[s];
    }
  }

  double evaluate_dof(std::size_t i,
                      const ufc::function& f,
                      const double * coordinate_dofs,
                      int cell_orientation,
                      const ufc::cell& c,
                      const ufc::coordinate_mapping * cm=nullptr
                      ) const final override
  {
    // Declare variables for result of evaluation
    double vals[1];
    // Declare variable for physical coordinates
    double y[2];
    switch (i)
    {
    case 0:
        {
            y[0] = coordinate_dofs[0];
            y[1] = coordinate_dofs[1];
            f.evaluate(vals, y, c);
            return vals[0];
        }
        break;
    case 1:
        {
            y[0] = coordinate_dofs[2];
            y[1] = coordinate_dofs[3];
            f.evaluate(vals, y, c);
            return vals[0];
        }
        break;
    case 2:
        {
            y[0] = coordinate_dofs[4];
            y[1] = coordinate_dofs[5];
            f.evaluate(vals, y, c);
            return vals[0];
        }
        break;
    }
    return 0.0;
  }

  void evaluate_dofs(double * values,
                             const ufc::function& f,
                             const double * coordinate_dofs,
                             int cell_orientation,
                             const ufc::cell& c,
                             const ufc::coordinate_mapping * cm=nullptr
                             ) const final override
  {
    // Declare variables for result of evaluation
    double vals[1];
    // Declare variable for physical coordinates
    double y[2];
    y[0] = coordinate_dofs[0];
    y[1] = coordinate_dofs[1];
    f.evaluate(vals, y, c);
    values[0] = vals[0];
    y[0] = coordinate_dofs[2];
    y[1] = coordinate_dofs[3];
    f.evaluate(vals, y, c);
    values[1] = vals[0];
    y[0] = coordinate_dofs[4];
    y[1] = coordinate_dofs[5];
    f.evaluate(vals, y, c);
    values[2] = vals[0];
  }

  void interpolate_vertex_values(double * vertex_values,
                                 const double * dof_values,
                                 const double * coordinate_dofs,
                                 int cell_orientation,
                                 const ufc::coordinate_mapping * cm=nullptr
                                 ) const final override
  {
    // Evaluate function and change variables
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[1];
    vertex_values[2] = dof_values[2];
  }

  void tabulate_dof_coordinates(double * dof_coordinates,
                                const double * coordinate_dofs,
                                const ufc::coordinate_mapping * cm=nullptr
                                ) const final override
  {
    dof_coordinates[0] = coordinate_dofs[0];
    dof_coordinates[1] = coordinate_dofs[1];
    dof_coordinates[2] = coordinate_dofs[2];
    dof_coordinates[2 + 1] = coordinate_dofs[3];
    dof_coordinates[2 * 2] = coordinate_dofs[4];
    dof_coordinates[2 * 2 + 1] = coordinate_dofs[5];
  }

  void tabulate_reference_dof_coordinates(double * reference_dof_coordinates) const final override
  {
    static const double dof_X[6] = { 0.0, 0.0, 1.0, 0.0, 0.0, 1.0 };
    std::copy_n(dof_X, 6, reference_dof_coordinates);
  }

  std::size_t num_sub_elements() const final override
  {
    return 0;
  }

  ufc::finite_element * create_sub_element(std::size_t i) const final override
  {
    return nullptr;
  }

  ufc::finite_element * create() const final override
  {
    return new linearspace2d_finite_element_0();
  }

};


class linearspace2d_dofmap_0: public ufc::dofmap
{
public:

  linearspace2d_dofmap_0() : ufc::dofmap()
  {
    // Do nothing
  }

  ~linearspace2d_dofmap_0() override
  {
    // Do nothing
  }

  const char * signature() const final override
  {
    return "FFC dofmap for FiniteElement('Lagrange', triangle, 1)";
  }

  bool needs_mesh_entities(std::size_t d) const final override
  {
    static const bool return_values[3] = { true, false, false };
    if (d >= 3)
        return false;
    return return_values[d];
  }

  std::size_t topological_dimension() const final override
  {
    return 2;
  }

  std::size_t global_dimension(const std::vector<std::size_t>&
                               num_global_entities) const final override
  {
    return num_global_entities[0];
  }

  std::size_t num_global_support_dofs() const final override
  {
    return 0;
  }

  std::size_t num_element_support_dofs() const final override
  {
    return 3;
  }

  std::size_t num_element_dofs() const final override
  {
    return 3;
  }

  std::size_t num_facet_dofs() const final override
  {
    return 2;
  }

  std::size_t num_entity_dofs(std::size_t d) const final override
  {
    static const std::size_t return_values[3] = { 1, 0, 0 };
    if (d >= 3)
        return 0;
    return return_values[d];
  }

  std::size_t num_entity_closure_dofs(std::size_t d) const final override
  {
    static const std::size_t return_values[3] = { 1, 2, 3 };
    if (d >= 3)
        return 0;
    return return_values[d];
  }

  void tabulate_dofs(std::size_t * dofs,
                     const std::vector<std::size_t>& num_global_entities,
                     const std::vector<std::vector<std::size_t>>& entity_indices) const final override
  {
    dofs[0] = entity_indices[0][0];
    dofs[1] = entity_indices[0][1];
    dofs[2] = entity_indices[0][2];
  }

  void tabulate_facet_dofs(std::size_t * dofs,
                           std::size_t facet) const final override
  {
    switch (facet)
    {
    case 0:
        dofs[0] = 1;
        dofs[1] = 2;
        break;
    case 1:
        dofs[0] = 0;
        dofs[1] = 2;
        break;
    case 2:
        dofs[0] = 0;
        dofs[1] = 1;
        break;
    }
  }

  void tabulate_entity_dofs(std::size_t * dofs,
                            std::size_t d, std::size_t i) const final override
  {
    switch (d)
    {
    case 0:
        switch (i)
        {
        case 0:
            dofs[0] = 0;
            break;
        case 1:
            dofs[0] = 1;
            break;
        case 2:
            dofs[0] = 2;
            break;
        }
        break;
    }
  }

  void tabulate_entity_closure_dofs(std::size_t * dofs,
                                    std::size_t d, std::size_t i) const final override
  {
    switch (d)
    {
    case 0:
        switch (i)
        {
        case 0:
            dofs[0] = 0;
            break;
        case 1:
            dofs[0] = 1;
            break;
        case 2:
            dofs[0] = 2;
            break;
        }
        break;
    case 1:
        switch (i)
        {
        case 0:
            dofs[0] = 1;
            dofs[1] = 2;
            break;
        case 1:
            dofs[0] = 0;
            dofs[1] = 2;
            break;
        case 2:
            dofs[0] = 0;
            dofs[1] = 1;
            break;
        }
        break;
    case 2:
        switch (i)
        {
        case 0:
            dofs[0] = 0;
            dofs[1] = 1;
            dofs[2] = 2;
            break;
        }
        break;
    }
  }


  std::size_t num_sub_dofmaps() const final override
  {
    return 0;
  }

  ufc::dofmap * create_sub_dofmap(std::size_t i) const final override
  {
    return nullptr;
  }

  ufc::dofmap * create() const final override
  {
    return new linearspace2d_dofmap_0();
  }

};

// DOLFIN wrappers

// Standard library includes
#include <string>

// DOLFIN includes
#include <dolfin/common/NoDeleter.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MultiMesh.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/Form.h>
#include <dolfin/fem/MultiMeshForm.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/MultiMeshFunctionSpace.h>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/function/CoefficientAssigner.h>
#include <dolfin/function/MultiMeshCoefficientAssigner.h>
#include <dolfin/adaptivity/ErrorControl.h>
#include <dolfin/adaptivity/GoalFunctional.h>
#include <dolfin/la/GenericVector.h>

namespace LinearSpace2D
{

class FunctionSpace: public dolfin::FunctionSpace
{
public:

  // Constructor for standard function space
  FunctionSpace(std::shared_ptr<const dolfin::Mesh> mesh):
    dolfin::FunctionSpace(mesh,
                          std::make_shared<const dolfin::FiniteElement>(std::make_shared<linearspace2d_finite_element_0>()),
                          std::make_shared<const dolfin::DofMap>(std::make_shared<linearspace2d_dofmap_0>(), *mesh))
  {
    // Do nothing
  }

  // Constructor for constrained function space
  FunctionSpace(std::shared_ptr<const dolfin::Mesh> mesh, std::shared_ptr<const dolfin::SubDomain> constrained_domain):
    dolfin::FunctionSpace(mesh,
                          std::make_shared<const dolfin::FiniteElement>(std::make_shared<linearspace2d_finite_element_0>()),
                          std::make_shared<const dolfin::DofMap>(std::make_shared<linearspace2d_dofmap_0>(), *mesh, constrained_domain))
  {
    // Do nothing
  }

};

}

#endif
