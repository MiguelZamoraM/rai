/*  ------------------------------------------------------------------
    Copyright (c) 2011-2020 Marc Toussaint
    email: toussaint@tu-berlin.de

    This code is distributed under the MIT License.
    Please see <root-path>/LICENSE for details.
    --------------------------------------------------------------  */

#include "fclInterface.h"

#ifdef RAI_FCL

// #include <fcl/broadphase/broadphase.h>
// #include <fcl/BVH/BVH_model.h>
// #include <fcl/distance.h>
// #include <fcl/collision.h>
// #include <fcl/collision_data.h>

namespace rai {
struct ConvexGeometryData {
  arr plane_dis;
  intA polygons;
};
}

rai::FclInterface::FclInterface(const rai::Array<Shape*>& geometries, double _cutoff)
  : cutoff(_cutoff) {
  convexGeometryData.resize(geometries.N);
  for(long int i=0; i<geometries.N; i++) {
    Shape* shape = geometries(i);

    if(shape) {
      // rai::Mesh& mesh = *geometries(i);
      std::shared_ptr<CollGeom> geom;
      if(shape->type()==ST_capsule) {
        geom = make_shared<Capsule>(shape->size(-1), shape->size(-2));
      } else if(shape->type()==ST_cylinder) {
        geom = make_shared<Cylinder>(shape->size(-1), shape->size(-2));
      } else if(shape->type()==ST_sphere) {
        geom = make_shared<Sphere>(shape->size(-1));
      }else if(shape->type()==ST_box){
        geom = make_shared<Box>(shape->size(0), shape->size(1), shape->size(2));
      } else {
        std::cout << shape->type() << std::endl;
        std::cout << shape->frame.name << std::endl;
        rai::Mesh& mesh = shape->mesh();
#if 0
        auto model = make_shared<fcl::BVHModel<fcl::OBBRSS>>();
        model->beginModel();
        for(uint i=0; i<mesh.T.d0; i++)
          model->addTriangle(Vec3f(&mesh.V(mesh.T(i, 0), 0)), Vec3f(&mesh.V(mesh.T(i, 1), 0)), Vec3f(&mesh.V(mesh.T(i, 2), 0)));
        model->endModel();
#elif 1
        mesh.computeNormals();
        std::shared_ptr<ConvexGeometryData> dat = make_shared<ConvexGeometryData>();
        dat->plane_dis = mesh.computeTriDistances();
        copy<int>(dat->polygons, mesh.T);
        dat->polygons.insColumns(0);
        for(uint i=0; i<dat->polygons.d0; i++) {dat->polygons(i, 0) = 3;}

#if FCL_MINOR_VERSION >= 7
        auto verts = make_shared<std::vector<fcl::Vector3<float>>>(mesh.V.d0);
        auto faces = make_shared<std::vector<int>>(mesh.T.N);
        for(uint i=0; i<verts->size(); i++){(*verts)[i] = {(float)mesh.V(i, 0), (float)mesh.V(i, 1), (float)mesh.V(i, 2)};
        std::cout << mesh.V(i, 0) << " " << mesh.V(i, 1) << " " << mesh.V(i, 2) << std::endl;}
        for(uint i=0; i<faces->size(); i++){(*faces)[i] = mesh.T.elem(i) * 1.;
        std::cout << mesh.T.elem(i) << std::endl;}
        geom = make_shared<fcl::Convex<float>>(verts, mesh.T.d0, faces, true);
#else
        geom = make_shared<fcl::Convex>((Vec3f*)mesh.Tn.p, dat->plane_dis.p, mesh.T.d0, (Vec3f*)mesh.V.p, mesh.V.d0, (int*)dat->polygons.p);
#endif

        convexGeometryData(i) = dat;
#else
        const auto model = make_shared<fcl::Sphere>(mesh.getRadius());
#endif
      }
      CollObject* obj = new CollObject(geom, fcl::Transform3f());
      obj->setUserData((void*)(i));
      objects.push_back(obj);
    }
  }

  //manager = make_shared<fcl::IntervalTreeCollisionManager>();
  //manager = make_shared<DynamicAABBTreeCollisionManager>();
  //manager = make_shared<fcl::SaPCollisionManager>();
  manager = make_shared<NaiveCollisionManager>();
  //manager = make_shared<fcl::SpatialHashingCollisionManager()>;
  manager->registerObjects(objects);
  manager->setup();
}

rai::FclInterface::~FclInterface() {
  for(size_t i = 0; i < objects.size(); ++i)
    delete objects[i];
}

void rai::FclInterface::step(const arr& X) {
  CHECK_EQ(X.nd, 2, "");
  CHECK_EQ(X.d0, convexGeometryData.N, "");
  CHECK_EQ(X.d1, 7, "");

  std::vector<CollObject*> updated_objs;

  // #pragma omp parallel for
  for(auto* obj:objects) {
    const uint i = (long int)obj->getUserData();
    if(i<X_lastQuery.d0 && maxDiff(X_lastQuery[i], X[i])<1e-8) {
      continue;
    }
    obj->setTranslation(Vec3f(X(i, 0), X(i, 1), X(i, 2)));
    obj->setQuatRotation(Quaternionf(X(i, 3), X(i, 4), X(i, 5), X(i, 6)));
    obj->computeAABB();

    updated_objs.push_back(obj);
  }

  //manager->update();
  manager->update(updated_objs);
  collisions.clear();
  manager->collide(this, BroadphaseCallback);

  collisions.reshape(collisions.N/2, 2);

  X_lastQuery = X;
}

void rai::FclInterface::deactivatePairs(const uintA& collisionExcludePairIDs){
  for (uint i=0; i<collisionExcludePairIDs.d0; ++i){
    deactivatedPairs.insert(key(collisionExcludePairIDs(i, 0), collisionExcludePairIDs(i, 1)));
  }
}

void rai::FclInterface::temporaryDeactivatePairs(const uintA& collisionExcludePairIDs){
  std::vector<std::size_t> tmp;
  for (uint i=0; i<collisionExcludePairIDs.d0; ++i){
    tmp.push_back(key(collisionExcludePairIDs(i, 0), collisionExcludePairIDs(i, 1)));
  }
  temporaryDeactivatedPairs.insert(tmp.begin(), tmp.end());
}

void rai::FclInterface::temporaryDeactivatePairs(const std::unordered_set<std::size_t>& hashedPairIds){
  temporaryDeactivatedPairs.insert(hashedPairIds.begin(), hashedPairIds.end());
}

void rai::FclInterface::addCollision(void* userData1, void* userData2) {
  const uint a = (long int)userData1;
  const uint b = (long int)userData2;
  collisions.resizeCopy(collisions.N+2);
  collisions.elem(-2) = a;
  collisions.elem(-1) = b;

  // std::cout << a << " " << b << std::endl;
}

// Return value indicates if we can stop early
bool rai::FclInterface::BroadphaseCallback(CollObject* o1, CollObject* o2, void* cdata_) {
  rai::FclInterface* self = static_cast<rai::FclInterface*>(cdata_);

  if (self->relevant_ids.size() > 0) {
    if (self->relevant_ids.count((long int)o1->getUserData()) == 0 &&
        self->relevant_ids.count((long int)o2->getUserData()) == 0){
      return false;
    }
  }

  // If a pair is part of the deactivated pairs, skip it
  const auto pair = self->key((long int)o1->getUserData(), (long int)o2->getUserData());
  if (self->vec != nullptr){
    for (auto tmp: *(self->vec)){
      if (tmp->count(pair) > 0){
        return false;
      }   
    }
  }
  if (self->deactivatedPairs.count(pair) > 0 || self->temporaryDeactivatedPairs.count(pair) > 0){
    return false;
  }

  bool addedCollision = false;

  if(self->cutoff==0.) { //fine boolean collision query
    CollisionRequest request;
    CollisionResult result;
    fcl::collide(o1, o2, request, result);
    if(result.isCollision()) {
      self->addCollision(o1->getUserData(), o2->getUserData());
      addedCollision = true;
    }
  } else if(self->cutoff>0.) { //fine distance query
    DistanceRequest request;
    DistanceResult result;
    fcl::distance(o1, o2, request, result);
    if(result.min_distance<self->cutoff){
      self->addCollision(o1->getUserData(), o2->getUserData());
      addedCollision = true;
    }
  } else { //just broadphase
    self->addCollision(o1->getUserData(), o2->getUserData());
    addedCollision = true;
  }
  
  // We can stop if an collision was added, and the 'stopearly' flag is set
  if (self->stopEarly && addedCollision){
    return true;
  }

  return false;
}

void rai::SplitFclInterface::Init(const Array<Shape*>& geometries, const std::unordered_set<std::size_t>& robot, const std::unordered_set<std::size_t>& obs, const std::unordered_set<std::size_t>& env,  double _cutoff){
  cutoff = _cutoff;
  convexGeometryData.resize(geometries.N);

  robot_ids = robot;
  env_ids = env;
  obs_ids = obs;

  for(long int i=0; i<geometries.N; i++) {
    if(geometries(i)) {
      rai::Mesh& mesh = geometries(i)->mesh();
#if 0
      auto model = make_shared<fcl::BVHModel<fcl::OBBRSS>>();
      model->beginModel();
      for(uint i=0; i<mesh.T.d0; i++)
        model->addTriangle(Vec3f(&mesh.V(mesh.T(i, 0), 0)), Vec3f(&mesh.V(mesh.T(i, 1), 0)), Vec3f(&mesh.V(mesh.T(i, 2), 0)));
      model->endModel();
#elif 1
      mesh.computeNormals();
      std::shared_ptr<ConvexGeometryData> dat = make_shared<ConvexGeometryData>();
      dat->plane_dis = mesh.computeTriDistances();
      copy<int>(dat->polygons, mesh.T);
      dat->polygons.insColumns(0);
      for(uint j=0; j<dat->polygons.d0; j++) {dat->polygons(j, 0) = 3;}

#if FCL_MINOR_VERSION >= 7
      auto verts = make_shared<std::vector<fcl::Vector3<float>>>(mesh.V.d0);
      auto faces = make_shared<std::vector<int>>(mesh.T.N);
      for(uint i=0; i<verts->size(); i++)(*verts)[i] = {(float)mesh.V(i, 0), (float)mesh.V(i, 1), (float)mesh.V(i, 2)};
      for(uint i=0; i<faces->size(); i++)(*faces)[i] = mesh.T.elem(i);
      auto model = make_shared<fcl::Convex<float>>(verts, mesh.T.d0, faces, true);
#else
      const auto model = make_shared<fcl::Convex>((Vec3f*)mesh.Tn.p, dat->plane_dis.p, mesh.T.d0, (Vec3f*)mesh.V.p, mesh.V.d0, (int*)dat->polygons.p);
#endif     

      convexGeometryData(i) = dat;
#else
      const auto model = make_shared<fcl::Sphere>(mesh.getRadius());
#endif
      CollObject* obj = new CollObject(model, fcl::Transform3f());
      obj->setUserData((void*)(i));

      if (robot_ids.count(i) > 0){
        robot_objects.push_back(obj);
        //std::cout << "adding " << i << " to robot" << std::endl;
      }
      if (env_ids.count(i) > 0){
        env_objects.push_back(obj);
        //std::cout << "adding " << i << " to env" << std::endl;
      }
      if (obs_ids.count(i) > 0){
        obs_objects.push_back(obj);
        //std::cout << "adding " << i << " to obs" << std::endl;
      }
    }
  }

  robot_manager = make_shared<NaiveCollisionManager>();
  env_manager = make_shared<NaiveCollisionManager>();
  obs_manager = make_shared<NaiveCollisionManager>();
  
  //robot_manager = make_shared<DynamicAABBTreeCollisionManager>();
  //env_manager = make_shared<DynamicAABBTreeCollisionManager>();
  //obs_manager = make_shared<DynamicAABBTreeCollisionManager>();

  // register objs.
  robot_manager->registerObjects(robot_objects);
  env_manager->registerObjects(env_objects);
  obs_manager->registerObjects(obs_objects);

  robot_manager->setup();
  env_manager->setup();
  obs_manager->setup();

  //manager = make_shared<fcl::IntervalTreeCollisionManager>();
  //manager = make_shared<DynamicAABBTreeCollisionManager>();
  //manager = make_shared<fcl::SaPCollisionManager>();
  //manager = make_shared<fcl::NaiveCollisionManager>();
  //manager = make_shared<fcl::SpatialHashingCollisionManager()>;
  //manager->registerObjects(objects);
  //manager->setup();
}

rai::SplitFclInterface::~SplitFclInterface() {
  for(size_t i = 0; i < robot_objects.size(); ++i){
    delete robot_objects[i];
  }
  for(size_t i = 0; i < env_objects.size(); ++i){
    delete env_objects[i];
  }
  for(size_t i = 0; i < obs_objects.size(); ++i){
    delete obs_objects[i];
  }
}

void rai::SplitFclInterface::step(const arr& X, const bool check_robot, const bool check_robot_obs, const bool check_obs_env) {
  CHECK_EQ(X.nd, 2, "");
  CHECK_EQ(X.d0, convexGeometryData.N, "");
  CHECK_EQ(X.d1, 7, "");

  std::vector<CollObject*> updated_robot_objs;

  arr lhs;
  arr rhs;

  for(auto* obj:robot_objects) {
    const uint i = (long int)obj->getUserData();
    if (X_lastQuery.d0 > 0){
    lhs.referToDim(X_lastQuery, i);
    rhs.referToDim(X, i);
    if(i<X_lastQuery.d0 && maxDiff(lhs, rhs)<1e-8) {
      continue;
    }
    }
    obj->setTranslation(Vec3f(X(i, 0), X(i, 1), X(i, 2)));
    obj->setQuatRotation(Quaternionf(X(i, 3), X(i, 4), X(i, 5), X(i, 6)));
    obj->computeAABB();

    updated_robot_objs.push_back(obj);
    //std::cout << "r updating " << i << std::endl;
  }
  if (updated_robot_objs.size() > 0){
    robot_manager->update(updated_robot_objs);
  }

  // the static environment  only needs to be setup once
  if (!initialized_frame_state_){
    std::vector<CollObject*> updated_env_objs;
    for(auto* obj:env_objects) {
      const uint i = (long int)obj->getUserData();
    if (X_lastQuery.d0 > 0){
      lhs.referToDim(X_lastQuery, i);
      rhs.referToDim(X, i);
      if(i<X_lastQuery.d0 && maxDiff(lhs, rhs)<1e-8) {
        continue;
      }}
      obj->setTranslation(Vec3f(X(i, 0), X(i, 1), X(i, 2)));
      obj->setQuatRotation(Quaternionf(X(i, 3), X(i, 4), X(i, 5), X(i, 6)));
      obj->computeAABB();

      updated_env_objs.push_back(obj);
    }
    if (updated_env_objs.size() > 0){
      env_manager->update(updated_env_objs);
    }
  }

  std::vector<CollObject*> updated_obs_objs;
  for(auto* obj:obs_objects) {
    const uint i = (long int)obj->getUserData();
    if (X_lastQuery.d0 > 0){
    lhs.referToDim(X_lastQuery, i);
    rhs.referToDim(X, i);
    if(i<X_lastQuery.d0 && maxDiff(lhs, rhs)<1e-8) {
      continue;
    }}
    obj->setTranslation(Vec3f(X(i, 0), X(i, 1), X(i, 2)));
    obj->setQuatRotation(Quaternionf(X(i, 3), X(i, 4), X(i, 5), X(i, 6)));
    obj->computeAABB();

    updated_obs_objs.push_back(obj);
    //std::cout << "o updating " << i << std::endl;
  }

  if (updated_obs_objs.size() > 0){
    obs_manager->update(updated_obs_objs);
  }

  collisions.clear();

  //std::cout << "A" << std::endl;
  //std::cout << "checking: " << check_robot << " " << check_robot_obs << " " << check_obs_env << std::endl;

  // robot self-collisions & robot env
  if (check_robot){
    robot_manager->collide(this, BroadphaseCallback);
    robot_manager->collide(env_manager.get(), this, BroadphaseCallback);

    //std::cout << "R: " << collisions.d0 << std::endl;
  }

  // robot-obs
  if (check_robot_obs){
    robot_manager->collide(obs_manager.get(), this, BroadphaseCallback);
    obs_manager->collide(this, BroadphaseCallback);
    //std::cout << "obs: " << collisions.d0 << std::endl;
  }

  // obs-env
  if (check_obs_env){
    obs_manager->collide(env_manager.get(), this, BroadphaseCallback);
    //std::cout << "env: " << collisions.d0 << std::endl;
  }

  collisions.reshape(collisions.N/2, 2);

  X_lastQuery = X;

  initialized_frame_state_ = true;
}

void rai::SplitFclInterface::deactivatePairs(const uintA& collisionExcludePairIDs){
  for (uint i=0; i<collisionExcludePairIDs.d0; ++i){
    deactivatedPairs.insert(key(collisionExcludePairIDs(i, 0), collisionExcludePairIDs(i, 1)));
  }
}

void rai::SplitFclInterface::addCollision(void* userData1, void* userData2) {
  const uint a = (long int)userData1;
  const uint b = (long int)userData2;
  collisions.resizeCopy(collisions.N+2);
  collisions.elem(-2) = a;
  collisions.elem(-1) = b;
}

// Return value indicates if we can stop early
bool rai::SplitFclInterface::BroadphaseCallback(CollObject* o1, CollObject* o2, void* cdata_) {
  rai::SplitFclInterface* self = static_cast<rai::SplitFclInterface*>(cdata_);

  // If a pair is part of the deactivated pairs, skip it
  const auto pair = self->key((long int)o1->getUserData(), (long int)o2->getUserData());
  if (self->deactivatedPairs.count(pair) > 0){
    //std::cout << "Q" << std::endl;
    return false;
  }

  //std::cout << "colliding " << (long int)o1->getUserData() << " " <<(long int)o2->getUserData() << std::endl;

  bool addedCollision = false;

  if(self->cutoff==0.) { //fine boolean collision query
    CollisionRequest request;
    CollisionResult result;
    fcl::collide(o1, o2, request, result);
    if(result.isCollision()) {
      self->addCollision(o1->getUserData(), o2->getUserData());
      addedCollision = true;
    }
  } else if(self->cutoff>0.) { //fine distance query
    DistanceRequest request;
    DistanceResult result;
    fcl::distance(o1, o2, request, result);
    if(result.min_distance<self->cutoff){
      self->addCollision(o1->getUserData(), o2->getUserData());
      addedCollision = true;
    }
  } else { //just broadphase
    self->addCollision(o1->getUserData(), o2->getUserData());
    addedCollision = true;
  }
  
  // We can stop if an collision was added, and the 'stopearly' flag is set
  if (self->stopEarly && addedCollision){
    return true;
  }

  return false;
}

#else //RAI_FCL
rai::FclInterface::FclInterface(const Array<ptr<Mesh>>& _geometries, double _cutoff) { NICO }
rai::FclInterface::~FclInterface() { NICO }
void rai::FclInterface::step(const arr& X) { NICO }
#endif
