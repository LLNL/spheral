#include "gtest/gtest.h"
#include "Kernel/TableKernel.hh"
#include "Kernel/GaussianKernel.hh"
#include "chai/ManagedArray.hpp"
#include "RAJA/RAJA.hpp"

using namespace Spheral;

TEST(TableKernelGPUTest, HostDeviceCopyAndUse) {
    // 1. Create a TableKernel on the host
    using Dim = Dim<3>;
    GaussianKernel<Dim> hostKernel(1.0);
    TableKernel<Dim> hostTableKernel(hostKernel);

    // 2. Create a ManagedArray for TableKernel on the device
    chai::ManagedArray<TableKernel<Dim>> deviceTableKernel(1);

    // 3. Copy the host kernel to the device (chai handles this implicitly)
    deviceTableKernel[0] = hostTableKernel;

    // 4. Create a ManagedArray for the result on the device
    chai::ManagedArray<double> device_result(1);

    // 5. Launch a kernel to use the TableKernel on the device
    RAJA::forall<RAJA::cuda_exec<256>>(RAJA::RangeSegment(0, 1), [=] RAJA_DEVICE (int i) {
        device_result[i] = deviceTableKernel[0].kernelValue(0.5, 1.0);
    });

    // 6. Copy the result back to the host (chai handles this implicitly)
    double host_result = device_result[0];

    // 7. Assert that the result is correct
    double expected_result = hostTableKernel.kernelValue(0.5, 1.0);
    EXPECT_NEAR(host_result, expected_result, 1.0e-12);

    // 8. Clean up device memory (chai handles this automatically)
}