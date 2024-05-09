#include "Octree.h"

//#include <cassert>
#include <cstring>  // memset.
#include <limits>   // infinity()

namespace unibn
{

/// ----------------------------------------------------------------------------------------------
template <typename PointT, typename ContainerT>
Octree<PointT, ContainerT>::Octant::Octant()
    : isLeaf(true), x(0.0f), y(0.0f), z(0.0f), extent(0.0f), start(0), end(0), size(0)
{
  memset(&child, 0, 8 * sizeof(Octant*));
}

template <typename PointT, typename ContainerT>
Octree<PointT, ContainerT>::Octant::~Octant()
{
  for (uint32_t i = 0; i < 8; ++i) delete child[i];
}


/// ----------------------------------------------------------------------------------------------
template <typename PointT, typename ContainerT>
Octree<PointT, ContainerT>::Octree() : root_(nullptr), data_(nullptr)
{
}

template <typename PointT, typename ContainerT>
Octree<PointT, ContainerT>::~Octree()
{
  // is it enough to only delete root? other leaf nodes?
  // it is OK if the ~Octant() is recursive
  delete root_;
  if (params_.copyPoints) delete data_;
}

template <typename PointT, typename ContainerT>
void Octree<PointT, ContainerT>::initialize(const ContainerT& pts, const OctreeParams& params)
{
  clear();
  params_ = params;

  if (params_.copyPoints)
    data_ = new ContainerT(pts);
  else
    data_ = &pts;

  const uint32_t N = static_cast<uint32_t> (pts.size());
  successors_ = std::vector<uint32_t>(N);

  // determine axis-aligned bounding box.
  float min[3], max[3];
  min[0] = get<0>(pts[0]);
  min[1] = get<1>(pts[0]);
  min[2] = get<2>(pts[0]);
  max[0] = min[0];
  max[1] = min[1];
  max[2] = min[2];

  for (uint32_t i = 0; i < N; ++i)
  {
    // initially each element links simply to the following element.
    successors_[i] = i + 1;

    const PointT& p = pts[i];

    if (get<0>(p) < min[0]) min[0] = get<0>(p);
    if (get<1>(p) < min[1]) min[1] = get<1>(p);
    if (get<2>(p) < min[2]) min[2] = get<2>(p);
    if (get<0>(p) > max[0]) max[0] = get<0>(p);
    if (get<1>(p) > max[1]) max[1] = get<1>(p);
    if (get<2>(p) > max[2]) max[2] = get<2>(p);
  }

  float ctr[3] = {min[0], min[1], min[2]};

  float maxextent = 0.5f * (max[0] - min[0]);
  ctr[0] += maxextent;
  for (uint32_t i = 1; i < 3; ++i)
  {
    float extent = 0.5f * (max[i] - min[i]);
    ctr[i] += extent;
    if (extent > maxextent) maxextent = extent;
  }

  root_ = createOctant(ctr[0], ctr[1], ctr[2], maxextent, 0, N - 1, N);
}

template <typename PointT, typename ContainerT>
void Octree<PointT, ContainerT>::initialize(const ContainerT& pts, const std::vector<uint32_t>& indexes,
                                            const OctreeParams& params)
{
  clear();
  params_ = params;

  if (params_.copyPoints)
    data_ = new ContainerT(pts);
  else
    data_ = &pts;

  const uint32_t N = static_cast<uint32_t> (pts.size());
  successors_ = std::vector<uint32_t>(N);

  if (indexes.size() == 0) return;

  // determine axis-aligned bounding box.
  uint32_t lastIdx = indexes[0];
  float min[3], max[3];
  min[0] = get<0>(pts[lastIdx]);
  min[1] = get<1>(pts[lastIdx]);
  min[2] = get<2>(pts[lastIdx]);
  max[0] = min[0];
  max[1] = min[1];
  max[2] = min[2];

  for (uint32_t i = 1; i < indexes.size(); ++i)
  {
    uint32_t idx = indexes[i];
    // initially each element links simply to the following element.
    successors_[lastIdx] = idx;

    const PointT& p = pts[idx];

    if (get<0>(p) < min[0]) min[0] = get<0>(p);
    if (get<1>(p) < min[1]) min[1] = get<1>(p);
    if (get<2>(p) < min[2]) min[2] = get<2>(p);
    if (get<0>(p) > max[0]) max[0] = get<0>(p);
    if (get<1>(p) > max[1]) max[1] = get<1>(p);
    if (get<2>(p) > max[2]) max[2] = get<2>(p);

    lastIdx = idx;
  }

  float ctr[3] = {min[0], min[1], min[2]};

  float maxextent = 0.5f * (max[0] - min[0]);
  ctr[0] += maxextent;
  for (uint32_t i = 1; i < 3; ++i)
  {
    float extent = 0.5f * (max[i] - min[i]);
    ctr[i] += extent;
    if (extent > maxextent) maxextent = extent;
  }

  root_ = createOctant(ctr[0], ctr[1], ctr[2], maxextent, indexes[0], lastIdx, indexes.size());
}

template <typename PointT, typename ContainerT>
void Octree<PointT, ContainerT>::clear()
{
  delete root_;
  if (params_.copyPoints) delete data_;
  root_ = 0;
  data_ = 0;
  successors_.clear();
}

template <typename PointT, typename ContainerT>
typename Octree<PointT, ContainerT>::Octant* Octree<PointT, ContainerT>::createOctant(float x, float y, float z,
                                                                                      float extent, uint32_t startIdx,
                                                                                      uint32_t endIdx, uint32_t size)
{
  // For a leaf we don't have to change anything; points are already correctly linked or correctly reordered.
  Octant* octant = new Octant;

  octant->isLeaf = true;

  octant->x = x;
  octant->y = y;
  octant->z = z;
  octant->extent = extent;

  octant->start = startIdx;
  octant->end = endIdx;
  octant->size = size;

  static const float factor[] = {-0.5f, 0.5f};

  // subdivide subset of points and re-link points according to Morton codes
  if (size > params_.bucketSize && extent > 2 * params_.minExtent)
  {
    octant->isLeaf = false;

    const ContainerT& points = *data_;
    std::vector<uint32_t> childStarts(8, 0);
    std::vector<uint32_t> childEnds(8, 0);
    std::vector<uint32_t> childSizes(8, 0);

    // re-link disjoint child subsets...
    uint32_t idx = startIdx;

    for (uint32_t i = 0; i < size; ++i)
    {
      const PointT& p = points[idx];

      // determine Morton code for each point...
      uint32_t mortonCode = 0;
      if (get<0>(p) > x) mortonCode |= 1;
      if (get<1>(p) > y) mortonCode |= 2;
      if (get<2>(p) > z) mortonCode |= 4;

      // set child starts and update successors...
      if (childSizes[mortonCode] == 0)
        childStarts[mortonCode] = idx;
      else
        successors_[childEnds[mortonCode]] = idx;
      childSizes[mortonCode] += 1;

      childEnds[mortonCode] = idx;
      idx = successors_[idx];
    }

    // now, we can create the child nodes...
    float childExtent = 0.5f * extent;
    bool firsttime = true;
    uint32_t lastChildIdx = 0;
    for (uint32_t i = 0; i < 8; ++i)
    {
      if (childSizes[i] == 0) continue;

      float childX = x + factor[(i & 1) > 0] * extent;
      float childY = y + factor[(i & 2) > 0] * extent;
      float childZ = z + factor[(i & 4) > 0] * extent;

      octant->child[i] = createOctant(childX, childY, childZ, childExtent, childStarts[i], childEnds[i], childSizes[i]);

      if (firsttime)
        octant->start = octant->child[i]->start;
      else
        successors_[octant->child[lastChildIdx]->end] =
            octant->child[i]->start;  // we have to ensure that also the child ends link to the next child start.

      lastChildIdx = i;
      octant->end = octant->child[i]->end;
      firsttime = false;
    }
  }

  return octant;
}

template <typename PointT, typename ContainerT>
template <typename Distance>
void Octree<PointT, ContainerT>::radiusNeighbors(const Octant* octant, const PointT& query, float radius,
                                                 float sqrRadius, std::vector<uint32_t>& resultIndices) const
{
  const ContainerT& points = *data_;

  // if search ball S(q,r) contains octant, simply add point indexes.
  if (contains<Distance>(query, sqrRadius, octant))
  {
    uint32_t idx = octant->start;
    for (uint32_t i = 0; i < octant->size; ++i)
    {
      resultIndices.push_back(idx);
      idx = successors_[idx];
    }

    return;  // early pruning.
  }

  if (octant->isLeaf)
  {
    uint32_t idx = octant->start;
    for (uint32_t i = 0; i < octant->size; ++i)
    {
      const PointT& p = points[idx];
      float dist = Distance::compute(query, p);
      if (dist < sqrRadius) resultIndices.push_back(idx);
      idx = successors_[idx];
    }

    return;
  }

  // check whether child nodes are in range.
  for (uint32_t c = 0; c < 8; ++c)
  {
    if (octant->child[c] == 0) continue;
    if (!overlaps<Distance>(query, radius, sqrRadius, octant->child[c])) continue;
    radiusNeighbors<Distance>(octant->child[c], query, radius, sqrRadius, resultIndices);
  }
}

template <typename PointT, typename ContainerT>
template <typename Distance>
void Octree<PointT, ContainerT>::radiusNeighbors(const Octant* octant, const PointT& query, float radius,
                                                 float sqrRadius, std::vector<uint32_t>& resultIndices,
                                                 std::vector<float>& distances) const
{
  const ContainerT& points = *data_;

  // if search ball S(q,r) contains octant, simply add point indexes and compute squared distances.
  if (contains<Distance>(query, sqrRadius, octant))
  {
    uint32_t idx = octant->start;
    for (uint32_t i = 0; i < octant->size; ++i)
    {
      resultIndices.push_back(idx);
      distances.push_back(Distance::compute(query, points[idx]));
      idx = successors_[idx];
    }

    return;  // early pruning.
  }

  if (octant->isLeaf)
  {
    uint32_t idx = octant->start;
    for (uint32_t i = 0; i < octant->size; ++i)
    {
      const PointT& p = points[idx];
      float dist = Distance::compute(query, p);
      if (dist < sqrRadius)
      {
        resultIndices.push_back(idx);
        distances.push_back(dist);
      }
      idx = successors_[idx];
    }

    return;
  }

  // check whether child nodes are in range.
  for (uint32_t c = 0; c < 8; ++c)
  {
    if (octant->child[c] == 0) continue;
    if (!overlaps<Distance>(query, radius, sqrRadius, octant->child[c])) continue;
    radiusNeighbors<Distance>(octant->child[c], query, radius, sqrRadius, resultIndices, distances);
  }
}

template <typename PointT, typename ContainerT>
template <typename Distance>
void Octree<PointT, ContainerT>::radiusNeighbors(const PointT& query, float radius,
                                                 std::vector<uint32_t>& resultIndices) const
{
  resultIndices.clear();
  if (root_ == 0) return;

  float sqrRadius = Distance::sqr(radius);  // "squared" radius
  radiusNeighbors<Distance>(root_, query, radius, sqrRadius, resultIndices);
}

template <typename PointT, typename ContainerT>
template <typename Distance>
void Octree<PointT, ContainerT>::radiusNeighbors(const PointT& query, float radius,
                                                 std::vector<uint32_t>& resultIndices,
                                                 std::vector<float>& distances) const
{
  resultIndices.clear();
  distances.clear();
  if (root_ == 0) return;

  float sqrRadius = Distance::sqr(radius);  // "squared" radius
  radiusNeighbors<Distance>(root_, query, radius, sqrRadius, resultIndices, distances);
}

template <typename PointT, typename ContainerT>
template <typename Distance>
bool Octree<PointT, ContainerT>::overlaps(const PointT& query, float radius, float sqRadius, const Octant* o)
{
  // we exploit the symmetry to reduce the test to testing if its inside the Minkowski sum around the positive quadrant.
  float x = get<0>(query) - o->x;
  float y = get<1>(query) - o->y;
  float z = get<2>(query) - o->z;

  x = std::abs(x);
  y = std::abs(y);
  z = std::abs(z);

  float maxdist = radius + o->extent;

  // Completely outside, since q' is outside the relevant area.
  if (x > maxdist || y > maxdist || z > maxdist) return false;

  int32_t num_less_extent = (x < o->extent) + (y < o->extent) + (z < o->extent);

  // Checking different cases:

  // a. inside the surface region of the octant.
  if (num_less_extent > 1) return true;

  // b. checking the corner region && edge region.
  x = std::max(x - o->extent, 0.0f);
  y = std::max(y - o->extent, 0.0f);
  z = std::max(z - o->extent, 0.0f);

  return (Distance::norm(x, y, z) < sqRadius);
}

template <typename PointT, typename ContainerT>
template <typename Distance>
bool Octree<PointT, ContainerT>::contains(const PointT& query, float sqRadius, const Octant* o)
{
  // we exploit the symmetry to reduce the test to test
  // whether the farthest corner is inside the search ball.
  float x = get<0>(query) - o->x;
  float y = get<1>(query) - o->y;
  float z = get<2>(query) - o->z;

  x = std::abs(x);
  y = std::abs(y);
  z = std::abs(z);
  // reminder: (x, y, z) - (-e, -e, -e) = (x, y, z) + (e, e, e)
  x += o->extent;
  y += o->extent;
  z += o->extent;

  return (Distance::norm(x, y, z) < sqRadius);
}

template <typename PointT, typename ContainerT>
template <typename Distance>
int32_t Octree<PointT, ContainerT>::findNeighbor(const PointT& query, float minDistance) const
{
  float maxDistance = std::numeric_limits<float>::infinity();
  int32_t resultIndex = -1;
  if(root_ == 0) return resultIndex;
  
  findNeighbor<Distance>(root_, query, minDistance, maxDistance, resultIndex);

  return resultIndex;
}

template <typename PointT, typename ContainerT>
template <typename Distance>
bool Octree<PointT, ContainerT>::findNeighbor(const Octant* octant, const PointT& query, float minDistance,
                                              float& maxDistance, int32_t& resultIndex) const
{
  const ContainerT& points = *data_;
  // 1. first descend to leaf and check in leafs points.
  if (octant->isLeaf)
  {
    uint32_t idx = octant->start;
    float sqrMaxDistance = Distance::sqr(maxDistance);
    float sqrMinDistance = (minDistance < 0) ? minDistance : Distance::sqr(minDistance);

    for (uint32_t i = 0; i < octant->size; ++i)
    {
      const PointT& p = points[idx];
      float dist = Distance::compute(query, p);
      if (dist > sqrMinDistance && dist < sqrMaxDistance)
      {
        resultIndex = idx;
        sqrMaxDistance = dist;
      }
      idx = successors_[idx];
    }

    maxDistance = Distance::sqrt(sqrMaxDistance);
    return inside<Distance>(query, maxDistance, octant);
  }

  // determine Morton code for each point...
  uint32_t mortonCode = 0;
  if (get<0>(query) > octant->x) mortonCode |= 1;
  if (get<1>(query) > octant->y) mortonCode |= 2;
  if (get<2>(query) > octant->z) mortonCode |= 4;

  if (octant->child[mortonCode] != 0)
  {
    if (findNeighbor<Distance>(octant->child[mortonCode], query, minDistance, maxDistance, resultIndex)) return true;
  }

  // 2. if current best point completely inside, just return.
  float sqrMaxDistance = Distance::sqr(maxDistance);

  // 3. check adjacent octants for overlap and check these if necessary.
  for (uint32_t c = 0; c < 8; ++c)
  {
    if (c == mortonCode) continue;
    if (octant->child[c] == 0) continue;
    if (!overlaps<Distance>(query, maxDistance, sqrMaxDistance, octant->child[c])) continue;
    if (findNeighbor<Distance>(octant->child[c], query, minDistance, maxDistance, resultIndex))
      return true;  // early pruning
  }

  // all children have been checked...check if point is inside the current octant...
  return inside<Distance>(query, maxDistance, octant);
}

template <typename PointT, typename ContainerT>
template <typename Distance>
bool Octree<PointT, ContainerT>::inside(const PointT& query, float radius, const Octant* octant)
{
  // we exploit the symmetry to reduce the test to test
  // whether the farthest corner is inside the search ball.
  float x = get<0>(query) - octant->x;
  float y = get<1>(query) - octant->y;
  float z = get<2>(query) - octant->z;

  x = std::abs(x) + radius;
  y = std::abs(y) + radius;
  z = std::abs(z) + radius;

  if (x > octant->extent) return false;
  if (y > octant->extent) return false;
  if (z > octant->extent) return false;

  return true;
}

}   // end of namespace unibn

