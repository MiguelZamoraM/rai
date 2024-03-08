/*  ------------------------------------------------------------------
    Copyright (c) 2011-2020 Marc Toussaint
    email: toussaint@tu-berlin.de

    This code is distributed under the MIT License.
    Please see <root-path>/LICENSE for details.
    --------------------------------------------------------------  */

#pragma once

// #include "mesh.h"
#include "../Kin/frame.h"

#include <unordered_set>

// namespace fcl {
// class CollisionObject;
// class DynamicAABBTreeCollisionManager;
// class BroadPhaseCollisionManager;
// };

#include <fcl/config.h>
#if FCL_MINOR_VERSION >= 6
#  include "fcl/fcl.h"
typedef fcl::CollisionObject<float> CollObject;
typedef fcl::CollisionGeometry<float> CollGeom;

typedef fcl::Vector3<float> Vec3f;
typedef fcl::Quaternionf Quaternionf;
typedef fcl::BroadPhaseCollisionManager<float> BroadPhaseCollisionManager;
typedef fcl::DynamicAABBTreeCollisionManager<float> DynamicAABBTreeCollisionManager;
typedef fcl::NaiveCollisionManager<float> NaiveCollisionManager;
typedef fcl::CollisionRequest<float> CollisionRequest;
typedef fcl::CollisionResult<float> CollisionResult;
typedef fcl::DistanceRequest<float> DistanceRequest;
typedef fcl::DistanceResult<float> DistanceResult;
typedef fcl::Box<float> Box;
typedef fcl::Sphere<float> Sphere;
typedef fcl::Capsule<float> Capsule;
typedef fcl::Cylinder<float> Cylinder;
#else
#  include <fcl/broadphase/broadphase.h>
#  include <fcl/BVH/BVH_model.h>
#  include <fcl/distance.h>
#  include <fcl/collision.h>
#  include <fcl/collision_data.h>
typedef fcl::CollisionObject CollObject;
typedef fcl::CollisionGeometry CollGeom;
typedef fcl::Vec3f Vec3f;
typedef fcl::Quaternion3f Quaternionf;
typedef fcl::BroadPhaseCollisionManager BroadPhaseCollisionManager;
typedef fcl::DynamicAABBTreeCollisionManager DynamicAABBTreeCollisionManager;
typedef fcl::NaiveCollisionManager NaiveCollisionManager;
typedef fcl::CollisionRequest CollisionRequest;
typedef fcl::CollisionResult CollisionResult;
typedef fcl::DistanceRequest DistanceRequest;
typedef fcl::DistanceResult DistanceResult;

typedef fcl::Box Box;
typedef fcl::Sphere Sphere;
typedef fcl::Capsule Capsule;
typedef fcl::Cylinder Cylinder;

#endif

// class CollObject;
// class BroadPhaseCollisionManager;

namespace rai {

struct FclInterface {
  Array<ptr<struct ConvexGeometryData>> convexGeometryData;

  std::vector<CollObject*> objects;
  shared_ptr<BroadPhaseCollisionManager> manager;

  double cutoff=0.; //0 -> perform fine boolean collision check; >0 -> perform fine distance computations; <0 -> only broadphase
  uintA collisions; //return values!
  arr X_lastQuery;  //memory to check whether an object has moved in consecutive queries

  bool stopEarly{false};

  // hashmap
  std::size_t hash(const uint& l, const uint& r) const {
    //uintmax_t k = std::hash<std::size_t>{}(l);
    //k <<= sizeof(uintmax_t) * 4;
    std::size_t k = l*10000;
    k += r;
    //k ^= std::hash<std::size_t>{}(r);
    return k;
    //return std::hash<uintmax_t>{}(k);
  }
  std::size_t key(const uint& a, const uint& b){
    if (a < b){
      return hash(a, b);
    }

    return hash(b, a);
  }
  std::unordered_set<std::size_t> deactivatedPairs;
  std::unordered_set<std::size_t> temporaryDeactivatedPairs;

  std::unordered_set<std::size_t> relevant_ids;

  // sets of things that do not need to be collision checked against each other.
  // avoids the insertion cost
  std::vector<std::unordered_set<std::size_t>*> *vec = nullptr;

  FclInterface(const Array<Shape*>& shapes, double _cutoff=0.);
  ~FclInterface();

  void step(const arr& X);

  void deactivatePairs(const uintA& collisionExcludePairIDs);
  void temporaryDeactivatePairs(const uintA& collisionExcludePairIDs);
  void temporaryDeactivatePairs(const std::unordered_set<std::size_t>& hashedPairIds);

private: //called by collision callback
  void addCollision(void* userData1, void* userData2);
  static bool BroadphaseCallback(CollObject* o1, CollObject* o2, void* cdata_);
};

struct SplitFclInterface {
  Array<ptr<struct ConvexGeometryData>> convexGeometryData;

  std::vector<CollObject*> robot_objects;
  std::vector<CollObject*> env_objects;
  std::vector<CollObject*> obs_objects;

  std::unordered_set<std::size_t> robot_ids;
  std::unordered_set<std::size_t> env_ids;
  std::unordered_set<std::size_t> obs_ids;

  shared_ptr<BroadPhaseCollisionManager> robot_manager;
  shared_ptr<BroadPhaseCollisionManager> env_manager;
  shared_ptr<BroadPhaseCollisionManager> obs_manager;

  double cutoff=0.; //0 -> perform fine boolean collision check; >0 -> perform fine distance computations; <0 -> only broadphase
  uintA collisions; //return values!
  arr X_lastQuery;  //memory to check whether an object has moved in consecutive queries

  bool stopEarly{false};

  // hashmap
  std::size_t hash(const uint& l, const uint& r) const {
    //uintmax_t k = std::hash<std::size_t>{}(l);
    //k <<= sizeof(uintmax_t) * 4;
    std::size_t k = l*10000;
    k += r;
    //k ^= std::hash<std::size_t>{}(r);
    return k;
    //return std::hash<uintmax_t>{}(k);
  }
  std::size_t key(const uint& a, const uint& b){
    if (a < b){
      return hash(a, b);
    }
    return hash(b, a);
  }
  std::unordered_set<std::size_t> deactivatedPairs;

  bool initialized_frame_state_ {false};

  SplitFclInterface(){};
  ~SplitFclInterface();

  void Init(const Array<Shape*>& geometries, const std::unordered_set<std::size_t>& robot, const std::unordered_set<std::size_t>& obs, const std::unordered_set<std::size_t>& env,  double _cutoff=0.);

  void step(const arr& X, const bool check_robot, const bool check_robot_obs, const bool check_obs_env);

  void deactivatePairs(const uintA& collisionExcludePairIDs);

private: //called by collision callback
  void addCollision(void* userData1, void* userData2);
  static bool BroadphaseCallback(CollObject* o1, CollObject* o2, void* cdata_);
};

}

