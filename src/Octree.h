#ifndef UNIBN_OCTREE_H_
#define UNIBN_OCTREE_H_

// Copyright (c) 2015 Jens Behley, University of Bonn
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights  to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

#include <stdint.h>
#include <cmath>
#include <vector>

// needed for gtest access to protected/private members ...
namespace
{
class OctreeTest;
}

namespace unibn
{

/**
 * Some traits to access coordinates regardless of the specific implementation of point
 * inspired by boost.geometry, which needs to be implemented by new points.
 *
 */
namespace traits
{

template <typename PointT, int D>
struct access
{
};

template <class PointT>
struct access<PointT, 0>
{
  static float get(const PointT& p)
  {
    return p.x;
  }
};

template <class PointT>
struct access<PointT, 1>
{
  static float get(const PointT& p)
  {
    return p.y;
  }
};

template <class PointT>
struct access<PointT, 2>
{
  static float get(const PointT& p)
  {
    return p.z;
  }
};
}  // namespace traits

/** convenience function for access of point coordinates **/
template <int D, typename PointT>
inline float get(const PointT& p)
{
  return traits::access<PointT, D>::get(p);
}

/**
 * Some generic distances: Manhattan, (squared) Euclidean, and Maximum distance.
 *
 * A Distance has to implement the methods
 * 1. compute of two points p and q to compute and return the distance between two points, and
 * 2. norm of x,y,z coordinates to compute and return the norm of a point p = (x,y,z)
 * 3. sqr and sqrt of value to compute the correct radius if a comparison is performed using squared norms (see
 *L2Distance)...
 */
template <typename PointT>
struct L1Distance
{
  static inline float compute(const PointT& p, const PointT& q)
  {
    float diff1 = get<0>(p) - get<0>(q);
    float diff2 = get<1>(p) - get<1>(q);
    float diff3 = get<2>(p) - get<2>(q);

    return std::abs(diff1) + std::abs(diff2) + std::abs(diff3);
  }

  static inline float norm(float x, float y, float z)
  {
    return std::abs(x) + std::abs(y) + std::abs(z);
  }

  static inline float sqr(float r)
  {
    return r;
  }

  static inline float sqrt(float r)
  {
    return r;
  }
};

template <typename PointT>
struct L2Distance
{
  static inline float compute(const PointT& p, const PointT& q)
  {
    float diff1 = get<0>(p) - get<0>(q);
    float diff2 = get<1>(p) - get<1>(q);
    float diff3 = get<2>(p) - get<2>(q);

    return std::pow(diff1, 2) + std::pow(diff2, 2) + std::pow(diff3, 2);
  }

  static inline float norm(float x, float y, float z)
  {
    return std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2);
  }

  static inline float sqr(float r)
  {
    return r * r;
  }

  static inline float sqrt(float r)
  {
    return std::sqrt(r);
  }
};

template <typename PointT>
struct MaxDistance
{
  static inline float compute(const PointT& p, const PointT& q)
  {
    float diff1 = std::abs(get<0>(p) - get<0>(q));
    float diff2 = std::abs(get<1>(p) - get<1>(q));
    float diff3 = std::abs(get<2>(p) - get<2>(q));

    float maximum = diff1;
    if (diff2 > maximum) maximum = diff2;
    if (diff3 > maximum) maximum = diff3;

    return maximum;
  }

  static inline float norm(float x, float y, float z)
  {
    float maximum = x;
    if (y > maximum) maximum = y;
    if (z > maximum) maximum = z;
    return maximum;
  }

  static inline float sqr(float r)
  {
    return r;
  }

  static inline float sqrt(float r)
  {
    return r;
  }
};

struct OctreeParams
{
 public:
  OctreeParams(uint32_t bucketSize = 32, bool copyPoints = false, float minExtent = 0.0f)
      : bucketSize(bucketSize), copyPoints(copyPoints), minExtent(minExtent)
  {
  }
  uint32_t bucketSize;
  bool copyPoints;
  float minExtent;
};

/** \brief Index-based Octree implementation offering different queries and insertion/removal of points.
 *
 * The index-based Octree uses a successor relation and a startIndex in each Octant to improve runtime
 * performance for radius queries. The efficient storage of the points by relinking list elements
 * bases on the insight that children of an Octant contain disjoint subsets of points inside the Octant and
 * that we can reorganize the points such that we get an continuous single connect list that we can use to
 * store in each octant the start of this list.
 *
 * Special about the implementation is that it allows to search for neighbors with arbitrary p-norms, which
 * distinguishes it from most other Octree implementations.
 *
 * We decided to implement the Octree using a template for points and containers. The container must have an
 * operator[], which allows to access the points, and a size() member function, which allows to get the size of the
 * container. For the points, we used an access trait to access the coordinates inspired by boost.geometry.
 * The implementation already provides a general access trait, which expects to have public member variables x,y,z.
 *
 * f you use the implementation or ideas from the corresponding paper in your academic work, it would be nice if you
 * cite the corresponding paper:
 *
 *    J. Behley, V. Steinhage, A.B. Cremers. Efficient Radius Neighbor Search in Three-dimensional Point Clouds,
 *    Proc. of the IEEE International Conference on Robotics and Automation (ICRA), 2015.
 *
 * In future, we might add also other neighbor queries and implement the removal and adding of points.
 *
 * \version 0.1-icra
 *
 * \author behley
 */

template <typename PointT, typename ContainerT = std::vector<PointT> >
class Octree
{
 public:
  Octree();
  ~Octree();

  /** \brief initialize octree with all points **/
  void initialize(const ContainerT& pts, const OctreeParams& params = OctreeParams());

  /** \brief initialize octree only from pts that are inside indexes. **/
  void initialize(const ContainerT& pts, const std::vector<uint32_t>& indexes,
                  const OctreeParams& params = OctreeParams());

  /** \brief remove all data inside the octree. **/
  void clear();

  /** \brief radius neighbor queries where radius determines the maximal radius of reported indices of points in
   * resultIndices **/
  template <typename Distance>
  void radiusNeighbors(const PointT& query, float radius, std::vector<uint32_t>& resultIndices) const;

  /** \brief radius neighbor queries with explicit (squared) distance computation. **/
  template <typename Distance>
  void radiusNeighbors(const PointT& query, float radius, std::vector<uint32_t>& resultIndices,
                       std::vector<float>& distances) const;

  /** \brief nearest neighbor queries. Using minDistance >= 0, we explicitly disallow self-matches.
   * @return index of nearest neighbor n with Distance::compute(query, n) > minDistance and otherwise -1.
   **/
  template <typename Distance>
  int32_t findNeighbor(const PointT& query, float minDistance = -1) const;

 protected:
  class Octant
  {
   public:
    Octant();
    ~Octant();

    bool isLeaf;

    // bounding box of the octant needed for overlap and contains tests...
    float x, y, z;  // center
    float extent;   // half of side-length

    uint32_t start, end;  // start and end in succ_
    uint32_t size;        // number of points

    Octant* child[8];
  };

  // not copyable, not assignable ...
  Octree(Octree&);
  Octree& operator=(const Octree& oct);

  /**
   * \brief creation of an octant using the elements starting at startIdx.
   *
   * The method reorders the index such that all points are correctly linked to successors belonging
   * to the same octant.
   *
   * \param x,y,z           center coordinates of octant
   * \param extent          extent of octant
   * \param startIdx        first index of points inside octant
   * \param endIdx          last index of points inside octant
   * \param size            number of points in octant
   *
   * \return  octant with children nodes.
   */
  Octant* createOctant(float x, float y, float z, float extent, uint32_t startIdx, uint32_t endIdx, uint32_t size);

  /** @return true, if search finished, otherwise false. **/
  template <typename Distance>
  bool findNeighbor(const Octant* octant, const PointT& query, float minDistance, float& maxDistance,
                    int32_t& resultIndex) const;

  template <typename Distance>
  void radiusNeighbors(const Octant* octant, const PointT& query, float radius, float sqrRadius,
                       std::vector<uint32_t>& resultIndices) const;

  template <typename Distance>
  void radiusNeighbors(const Octant* octant, const PointT& query, float radius, float sqrRadius,
                       std::vector<uint32_t>& resultIndices, std::vector<float>& distances) const;

  /** \brief test if search ball S(q,r) overlaps with octant
   *
   * @param query   query point
   * @param radius  "squared" radius
   * @param o       pointer to octant
   *
   * @return true, if search ball overlaps with octant, false otherwise.
   */
  template <typename Distance>
  static bool overlaps(const PointT& query, float radius, float sqRadius, const Octant* o);

  /** \brief test if search ball S(q,r) contains octant
   *
   * @param query    query point
   * @param sqRadius "squared" radius
   * @param octant   pointer to octant
   *
   * @return true, if search ball overlaps with octant, false otherwise.
   */
  template <typename Distance>
  static bool contains(const PointT& query, float sqRadius, const Octant* octant);

  /** \brief test if search ball S(q,r) is completely inside octant.
   *
   * @param query   query point
   * @param radius  radius r
   * @param octant  point to octant.
   *
   * @return true, if search ball is completely inside the octant, false otherwise.
   */
  template <typename Distance>
  static bool inside(const PointT& query, float radius, const Octant* octant);

  OctreeParams params_;
  Octant* root_;
  const ContainerT* data_;

  std::vector<uint32_t> successors_;  // single connected list of next point indices...

  friend class ::OctreeTest;
};

}  // namespace unibn

#include "Octree.hpp"
#endif /* OCTREE_H_ */