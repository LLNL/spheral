//---------------------------------Spheral++----------------------------------//
// std::make_shared doesn't play well with initializer_lists, so we provide an
// explicit function to help out
//----------------------------------------------------------------------------//
namespace Spheral {

template<typename T, typename... Args>
inline
std::shared_ptr<T> make_policy(std::initializer_list<std::string> depends, Args&&... args) {
  return std::make_shared<T>(depends, std::forward<Args>(args)...);
}

template<typename T, typename... Args>
inline
std::shared_ptr<T> make_policy(Args&&... args) {
  return std::make_shared<T>(std::forward<Args>(args)...);
}

}
