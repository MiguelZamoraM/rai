/*  ------------------------------------------------------------------
    Copyright (c) 2011-2020 Marc Toussaint
    email: toussaint@tu-berlin.de

    This code is distributed under the MIT License.
    Please see <root-path>/LICENSE for details.
    --------------------------------------------------------------  */

#pragma once

#include "mesh.h"

#include <unordered_set>

namespace fcl {
class CollisionObject;
class DynamicAABBTreeCollisionManager;
class BroadPhaseCollisionManager;
};

namespace rai {

struct FclInterface {
  Array<ptr<struct ConvexGeometryData>> convexGeometryData;
  std::vector<fcl::CollisionObject*> objects;
  shared_ptr<fcl::BroadPhaseCollisionManager> manager;

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

  FclInterface(const Array<ptr<Mesh>>& geometries, double _cutoff=0.);
  ~FclInterface();

  void step(const arr& X);

  void deactivatePairs(const uintA& collisionExcludePairIDs);
  void temporaryDeactivatePairs(const uintA& collisionExcludePairIDs);
  void temporaryDeactivatePairs(const std::unordered_set<std::size_t>& hashedPairIds);

private: //called by collision callback
  void addCollision(void* userData1, void* userData2);
  static bool BroadphaseCallback(fcl::CollisionObject* o1, fcl::CollisionObject* o2, void* cdata_);
};

struct SplitFclInterface {
  Array<ptr<struct ConvexGeometryData>> convexGeometryData;

  std::vector<fcl::CollisionObject*> robot_objects;
  std::vector<fcl::CollisionObject*> env_objects;
  std::vector<fcl::CollisionObject*> obs_objects;

  std::unordered_set<std::size_t> robot_ids;
  std::unordered_set<std::size_t> env_ids;
  std::unordered_set<std::size_t> obs_ids;

  shared_ptr<fcl::BroadPhaseCollisionManager> robot_manager;
  shared_ptr<fcl::BroadPhaseCollisionManager> env_manager;
  shared_ptr<fcl::BroadPhaseCollisionManager> obs_manager;

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

  void Init(const Array<ptr<Mesh>>& geometries, const std::unordered_set<std::size_t>& robot, const std::unordered_set<std::size_t>& obs, const std::unordered_set<std::size_t>& env,  double _cutoff=0.);

  void step(const arr& X, const bool check_robot, const bool check_robot_obs, const bool check_obs_env);

  void deactivatePairs(const uintA& collisionExcludePairIDs);

private: //called by collision callback
  void addCollision(void* userData1, void* userData2);
  static bool BroadphaseCallback(fcl::CollisionObject* o1, fcl::CollisionObject* o2, void* cdata_);
};

}

