#include <boost/python.hpp>
#include "Utilities/iterateIdealH.hh"
#include "Geometry/Dimension.hh"

using namespace boost::python;

namespace Spheral {

//------------------------------------------------------------------------------
// In order to provide our default arguments, we have to write thin wrappers.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
iterateIdealH7(DataBaseSpace::DataBase<Dimension>& dataBase,
               const std::vector<BoundarySpace::Boundary<Dimension>*>& boundaries,
               const KernelSpace::TableKernel<Dimension>& W,
               const int maxIterations,
               const double tolerance,
               const double nPerhForIteration,
               const bool sphericalStart) {
  iterateIdealH(dataBase, boundaries, W, maxIterations, tolerance, nPerhForIteration, sphericalStart);
}

template<typename Dimension>
inline
void
iterateIdealH6(DataBaseSpace::DataBase<Dimension>& dataBase,
               const std::vector<BoundarySpace::Boundary<Dimension>*>& boundaries,
               const KernelSpace::TableKernel<Dimension>& W,
               const int maxIterations,
               const double tolerance,
               const double nPerhForIteration) {
  iterateIdealH(dataBase, boundaries, W, maxIterations, tolerance, nPerhForIteration);
}

template<typename Dimension>
inline
void
iterateIdealH5(DataBaseSpace::DataBase<Dimension>& dataBase,
               const std::vector<BoundarySpace::Boundary<Dimension>*>& boundaries,
               const KernelSpace::TableKernel<Dimension>& W,
               const int maxIterations,
               const double tolerance) {
  iterateIdealH(dataBase, boundaries, W, maxIterations, tolerance);
}

template<typename Dimension>
inline
void
iterateIdealH4(DataBaseSpace::DataBase<Dimension>& dataBase,
               const std::vector<BoundarySpace::Boundary<Dimension>*>& boundaries,
               const KernelSpace::TableKernel<Dimension>& W,
               const int maxIterations) {
  iterateIdealH(dataBase, boundaries, W, maxIterations);
}

template<typename Dimension>
inline
void
iterateIdealH3(DataBaseSpace::DataBase<Dimension>& dataBase,
               const std::vector<BoundarySpace::Boundary<Dimension>*>& boundaries,
               const KernelSpace::TableKernel<Dimension>& W) {
  iterateIdealH(dataBase, boundaries, W);
}

//------------------------------------------------------------------------------
// Expose the ideal H iteration algorithm.
//------------------------------------------------------------------------------
void wrapIterateIdealH() {
  def("iterateIdealH", &iterateIdealH<Dim<1> >, "Apply the ideal H algorithm to the nodes in the DataBase.");
  def("iterateIdealH", &iterateIdealH<Dim<2> >, "Apply the ideal H algorithm to the nodes in the DataBase.");
  def("iterateIdealH", &iterateIdealH<Dim<3> >, "Apply the ideal H algorithm to the nodes in the DataBase.");

  def("iterateIdealH", &iterateIdealH3<Dim<1> >, "Apply the ideal H algorithm to the nodes in the DataBase.");
  def("iterateIdealH", &iterateIdealH3<Dim<2> >, "Apply the ideal H algorithm to the nodes in the DataBase.");
  def("iterateIdealH", &iterateIdealH3<Dim<3> >, "Apply the ideal H algorithm to the nodes in the DataBase.");
                                                          
  def("iterateIdealH", &iterateIdealH4<Dim<1> >, "Apply the ideal H algorithm to the nodes in the DataBase.");
  def("iterateIdealH", &iterateIdealH4<Dim<2> >, "Apply the ideal H algorithm to the nodes in the DataBase.");
  def("iterateIdealH", &iterateIdealH4<Dim<3> >, "Apply the ideal H algorithm to the nodes in the DataBase.");
                                                          
  def("iterateIdealH", &iterateIdealH5<Dim<1> >, "Apply the ideal H algorithm to the nodes in the DataBase.");
  def("iterateIdealH", &iterateIdealH5<Dim<2> >, "Apply the ideal H algorithm to the nodes in the DataBase.");
  def("iterateIdealH", &iterateIdealH5<Dim<3> >, "Apply the ideal H algorithm to the nodes in the DataBase.");
                                                          
  def("iterateIdealH", &iterateIdealH6<Dim<1> >, "Apply the ideal H algorithm to the nodes in the DataBase.");
  def("iterateIdealH", &iterateIdealH6<Dim<2> >, "Apply the ideal H algorithm to the nodes in the DataBase.");
  def("iterateIdealH", &iterateIdealH6<Dim<3> >, "Apply the ideal H algorithm to the nodes in the DataBase.");
                                                          
  def("iterateIdealH", &iterateIdealH7<Dim<1> >, "Apply the ideal H algorithm to the nodes in the DataBase.");
  def("iterateIdealH", &iterateIdealH7<Dim<2> >, "Apply the ideal H algorithm to the nodes in the DataBase.");
  def("iterateIdealH", &iterateIdealH7<Dim<3> >, "Apply the ideal H algorithm to the nodes in the DataBase.");
}

}
