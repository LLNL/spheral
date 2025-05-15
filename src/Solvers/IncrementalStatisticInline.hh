#include <iostream>
#include <iomanip>
#include <limits>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename DataType>
inline
IncrementalStatistic<DataType>::
IncrementalStatistic(const DataType meanGuess,
                     std::string name,
                     bool print):
  mMeanGuess(meanGuess),
  mName(name),
  mPrint(print),
  mMin(std::numeric_limits<DataType>::max()),
  mMax(std::numeric_limits<DataType>::min()),
  mNumPoints(0),
  mShiftedSum(0),
  mShiftedSum2(0) {
}

//------------------------------------------------------------------------------
// Add a data point
//------------------------------------------------------------------------------
template<typename DataType>
inline
void
IncrementalStatistic<DataType>::
add(const DataType data) {
  if (data < mMin) mMin = data;
  if (data > mMax) mMax = data;
  mNumPoints += 1;
  const DataType delta = data - mMeanGuess;
  mShiftedSum += delta;
  mShiftedSum2 += delta * delta;
  
  if (mPrint) {
    std::cout << mName;
    std::cout << "\tadd: " << data;
    std::cout << std::endl;
  }
}

//------------------------------------------------------------------------------
// Calculate mean
//------------------------------------------------------------------------------
template<typename DataType>
inline
DataType
IncrementalStatistic<DataType>::
mean() const {
  return mShiftedSum / mNumPoints + mMeanGuess;
}

//------------------------------------------------------------------------------
// Calculate variance
//------------------------------------------------------------------------------
template<typename DataType>
inline
DataType
IncrementalStatistic<DataType>::
variance() const {
  return (mShiftedSum2 - mShiftedSum * mShiftedSum / mNumPoints) / mNumPoints;
}

//------------------------------------------------------------------------------
// Calculate total
//------------------------------------------------------------------------------
template<typename DataType>
inline
DataType
IncrementalStatistic<DataType>::
total() const {
  return mShiftedSum + mNumPoints * mMeanGuess;
}

//------------------------------------------------------------------------------
// Get name
//------------------------------------------------------------------------------
template<typename DataType>
inline
std::string
IncrementalStatistic<DataType>::
name() const {
  return mName;
}

//------------------------------------------------------------------------------
// Set name
//------------------------------------------------------------------------------
template<typename DataType>
inline
void
IncrementalStatistic<DataType>::
setName(std::string name) {
  mName = name;
  return;
}

//------------------------------------------------------------------------------
// Get mean guess
//------------------------------------------------------------------------------
template<typename DataType>
inline
DataType
IncrementalStatistic<DataType>::
min() const {
  return mMin;
}

//------------------------------------------------------------------------------
// Get mean guess
//------------------------------------------------------------------------------
template<typename DataType>
inline
DataType
IncrementalStatistic<DataType>::
max() const {
  return mMax;
}

//------------------------------------------------------------------------------
// Get mean guess
//------------------------------------------------------------------------------
template<typename DataType>
inline
DataType
IncrementalStatistic<DataType>::
meanGuess() const {
  return mMeanGuess;
}

//------------------------------------------------------------------------------
// Get numPoints
//------------------------------------------------------------------------------
template<typename DataType>
inline
DataType
IncrementalStatistic<DataType>::
numPoints() const {
  return mNumPoints;
}

//------------------------------------------------------------------------------
// Get shiftedSum
//------------------------------------------------------------------------------
template<typename DataType>
inline
DataType
IncrementalStatistic<DataType>::
shiftedSum() const {
  return mShiftedSum;
}

//------------------------------------------------------------------------------
// Get shiftedSum2
//------------------------------------------------------------------------------
template<typename DataType>
inline
DataType
IncrementalStatistic<DataType>::
shiftedSum2() const {
  return mShiftedSum2;
}

//------------------------------------------------------------------------------
// Print results
//------------------------------------------------------------------------------
template<typename DataType>
inline
void
IncrementalStatistic<DataType>::
print() const {
  typedef std::pair<std::string, DataType> PairType;

  // Get function to make the output consistent
  auto outputPair = [](PairType val) {
    std::cout << std::setw(20) << val.first;
    std::cout << std::setw(14) << std::setprecision(6) << val.second;
  };

  // Print name
  std::cout << mName << std::endl;

  if (numPoints() == 0) {
    // If there is no data, output simpler message
    PairType val = {"points", numPoints()};
    outputPair(val);
    std::cout << "(no data recorded)" << std::endl;
  }
  else {
    // Output in two columns
    std::vector<PairType> vals = {{"points",   numPoints()}, {"total",   total()},
                                  {"mean",     mean()},      {"minimum", mMin},
                                  {"variance", variance()},  {"maximum", mMax}};
    auto index = 0;
    for (auto val : vals) {
      outputPair(val);
      if (index % 2 == 1) {
        std::cout << std::endl;
      }
      index += 1;
    }
  }
}

} // end namespace Spheral
