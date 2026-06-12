#pragma once

#include "real_type.h"

namespace gpu::hamiltonian {

struct Vec3 {
   real x;
   real y;
   real z;
};

__device__ inline Vec3 make_vec3(real x, real y, real z) {
   return {x, y, z};
}

__device__ inline Vec3 operator+(Vec3 a, Vec3 b) {
   return {a.x + b.x, a.y + b.y, a.z + b.z};
}

__device__ inline Vec3& operator+=(Vec3& a, Vec3 b) {
   a.x += b.x;
   a.y += b.y;
   a.z += b.z;
   return a;
}

__device__ inline real dot(Vec3 a, Vec3 b) {
   return a.x * b.x + a.y * b.y + a.z * b.z;
}

__device__ inline Vec3 dm_field(Vec3 d, Vec3 s) {
   return {
      -d.z * s.y + d.y * s.z,
      -d.x * s.z + d.z * s.x,
      -d.y * s.x + d.x * s.y
   };
}

} // namespace gpu::hamiltonian
