// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_APPLY_13_QF_H
#define PALACE_LIBCEED_APPLY_13_QF_H

CEED_QFUNCTION(f_apply_13)(void *, CeedInt Q, const CeedScalar *const *in,
                           CeedScalar *const *out)
{
  const CeedScalar *__restrict__ qd1 = in[0], *__restrict__ qd2 = in[0] + Q,
                                 *__restrict__ u1 = in[1], *__restrict__ u2 = in[2];
  CeedScalar *__restrict__ v1 = out[0], *__restrict__ v2 = out[1];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    v1[i] = qd1[i] * u1[i];

    const CeedScalar u20 = u2[i + Q * 0];
    const CeedScalar u21 = u2[i + Q * 1];
    const CeedScalar u22 = u2[i + Q * 2];
    v2[i + Q * 0] = qd2[i + Q * 0] * u20 + qd2[i + Q * 3] * u21 + qd2[i + Q * 6] * u22;
    v2[i + Q * 1] = qd2[i + Q * 1] * u20 + qd2[i + Q * 4] * u21 + qd2[i + Q * 7] * u22;
    v2[i + Q * 2] = qd2[i + Q * 2] * u20 + qd2[i + Q * 5] * u21 + qd2[i + Q * 8] * u22;
  }
  return 0;
}

#endif  // PALACE_LIBCEED_APPLY_13_QF_H
