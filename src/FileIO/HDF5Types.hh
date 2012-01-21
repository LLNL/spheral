#include "H5Cpp.h"
using namespace H5;

namespace Spheral {
namespace FileIO {

static CompType H5Vector1d;
static CompType H5Vector2d;
static CompType H5Vector3d;

static CompType H5Tensor1d;
static CompType H5Tensor2d;
static CompType H5Tensor3d;

static CompType H5SymTensor1d;
static CompType H5SymTensor2d;
static CompType H5SymTensor3d;

int initializeSpheralH5Types();

}
}
