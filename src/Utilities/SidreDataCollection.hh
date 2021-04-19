#ifndef SidreDataCollection_HH
#define SidreDataCollection_HH

#include "axom/sidre.hpp"
#include "Field/Field.hh"
#include "Geometry/Dimension.hh"

namespace Spheral
{

class SidreDataCollection
{
public:
    SidreDataCollection();
    ~SidreDataCollection();

    template<typename Dimension, typename DataType>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, DataType> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, std::string> &field);
    template<typename Dimension, typename DataType>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, std::vector<DataType>> &field);
    template<typename Dimension, typename DataType>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, std::tuple<DataType, DataType, DataType>> &field);
    template<typename Dimension, typename DataType>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, std::tuple<DataType, DataType, DataType, DataType>> &field);
    template<typename Dimension, typename DataType>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, std::tuple<DataType, DataType, DataType, DataType, DataType>> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, Dim<1>::Vector> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, Dim<1>::Vector3d> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, Dim<1>::Tensor> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, Dim<1>::SymTensor> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, Dim<1>::ThirdRankTensor> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, Dim<1>::FourthRankTensor> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, Dim<1>::FifthRankTensor> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, Dim<1>::FacetedVolume> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, Dim<2>::Vector> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, Dim<2>::Tensor> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, Dim<2>::SymTensor> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, Dim<2>::ThirdRankTensor> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, Dim<2>::FourthRankTensor> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, Dim<2>::FifthRankTensor> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, Dim<2>::FacetedVolume> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, Dim<3>::Vector> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, Dim<3>::Tensor> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, Dim<3>::SymTensor> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, Dim<3>::ThirdRankTensor> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, Dim<3>::FourthRankTensor> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, Dim<3>::FifthRankTensor> &field);
    template<typename Dimension>
    axom::sidre::View *alloc_view(const std::string &view_name,
                                  const Spheral::Field<Dimension, Dim<3>::FacetedVolume> &field);
    
    void printDataStore() {m_datastore_ptr->getRoot()->print();};
private:
    axom::sidre::DataStore *m_datastore_ptr;
};

}

#include "Utilities/SidreDataCollectionInline.hh"

#endif
